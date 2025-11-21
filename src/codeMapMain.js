const codeMapMain = {
  chip_pre: `
      subroutine pre_process(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT,FLAT_FILE
      integer i
      character*2 chip_name

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

20      format(I2.2)

      do ichip=1,nchip
        call get_chip_id(IMAGE_FILE(ichip),i)
        write(chip_name,20) i
        FLAT_FILE=trim(FLAT_PATH)//'/flat_'
     .//trim(chip_name)//'_weight.fits'

        call chip_pre_process(IMAGE_FILE(ichip),DIR_OUTPUT
     .,FLAT_FILE)
c        call chip_pre_process(IMAGE_FILE(ichip),DIR_OUTPUT)
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chip_pre_process(IMAGE_FILE,DIR_OUTPUT,FLAT_FILE)
c      subroutine chip_pre_process(IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      character*(*) IMAGE_FILE,DIR_OUTPUT
      integer proc_error
      character*(strl) PREFIX,filename,catfile

      integer nx,ny
      real array(npx,npy),normap(npx,npy)
      integer weight(npx,npy)
c      real flat(npx,npy)

      integer i,j,u,v,nxc
      double precision cRPIX(2),cD(2,2),cRVAL(2)

      real sigabc(2,3)
      integer order
      common /sig_pass/ sigabc,order

      real aa,bb,cc

      character*(strl) FLAT_FILE
      integer nxx , nyy
      real flat_weight(npx,npy)


      proc_error=0

c      call readimage(FLAT_FILE,nx,ny,npx,npy,flat)
      call readimage_para(IMAGE_FILE                              
     .,nx,ny,npx,npy,array,cRPIX,cD,cRVAL)

      do i=1,nx
        do j=1,ny
          weight(i,j)=1
c          if (include_FLAT.eq.1) then
c            if (flat(i,j).lt.0.5) then
c              weight(i,j)=0
c            else
c              array(i,j)=array(i,j)*flat(i,j)
c            endif
c          endif
          if (array(i,j).gt.saturation_thresh) weight(i,j)=0
          normap(i,j)=array(i,j)
        enddo
      enddo

c------------------------------------------------------

      nxc=nx/2

      if (ccD_split.eq.2) then

        call set_background(1,nxc,1,ny,npx,npy,normap
     .,blocksize,nct,ncx,proc_error)
        call set_background(nxc+1,nx,1,ny,npx,npy,normap
     .,blocksize,nct,ncx,proc_error)
        call set_sig(1,nxc,1,ny,npx,npy,normap,aa,bb,cc,proc_error)
        sigabc(1,1)=aa
        sigabc(1,2)=bb
        sigabc(1,3)=cc
        call set_sig(nxc+1,nx,1,ny,npx,npy,normap,aa,bb,cc,proc_error)
        sigabc(2,1)=aa
        sigabc(2,2)=bb
        sigabc(2,3)=cc
      else
        call set_background(1,nx,1,ny,npx,npy,normap        
     .,blocksize,nct,ncx,proc_error)
        call set_sig(1,nx,1,ny,npx,npy,normap,aa,bb,cc,proc_error)
        sigabc(1,1)=aa
        sigabc(1,2)=bb
        sigabc(1,3)=cc
      endif


c--------------------------------------------------------------

      call get_PREFIX(IMAGE_FILE,PREFIX)
      filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
      if (ASTROMETRY_trivial.eq.1) then
        call gen_astrometry_data_trivial(cRPIX,cD,cRVAL,filename)
      else
        catfile=ASTROMETRY_cAT
        call generate_gaia_file_name(cRVAL,catfile)
        call gen_astrometry_data(catfile,nx,ny,npx,npy   
     .,normap,weight,cRPIX,cD,cRVAL,filename,proc_error)
      endif

      call locate_defects(nx,ny,npx,npy,array,normap
     .,weight,area_max,area_thresh,proc_error)
      call merge_defects(nx,ny,npx,npy,weight,normap
     .,area_max,source_thresh,area_thresh,proc_error)

c------------------------------------------------------------------

      if (proc_error.eq.0) then
c     Super flat field correction
        if (include_FLAT.eq.1) then
          call readimage(FLAT_FILE,nxx,nyy,npx,npy,flat_weight)
          if (nxx.ne.nx .or. nyy.ne.ny) then
            write(*,*) 'Error: wrong size of flat file!'
            proc_error=1
          endif
          do i=1,nx
            do j=1,ny
              if (flat_weight(i,j).lt.(-900.0)) weight(i,j) = 0
              if (weight(i,j).eq.0) normap(i,j)=-1000. 
            enddo
          enddo    
        else
          do i=1,nx
            do j=1,ny
              if (weight(i,j).eq.0) normap(i,j)=-1000.
            enddo
          enddo
        endif
      else
        do i=1,nx
          do j=1,ny
            normap(i,j)=-1000.
          enddo
        enddo
      endif

      if (proc_error.eq.0) then
        normap(1,1)=-1.
      else
        normap(1,1)=1.
      endif

      do i=1,ccD_split
        do j=1,3
          normap(1+i,j)=sigabc(i,j)
        enddo
      enddo

      filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
     .//'_norm.fits'
      call writeimage_copyhdu(IMAGE_FILE,filename
     .,nx,ny,npx,npy,normap)

c-----------------------------------------------------------
      if (proc_error.eq.0) then
        write(*,*) 'Status of processing ',trim(IMAGE_FILE)
     .,': OK.'
      else
        write(*,*) 'Status of processing ',trim(IMAGE_FILE)
     .,': ERROR!'
      endif

c      pause

      return
      END
`,

  astro: `
      subroutine chip_process_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer proc_error,ichip
      character PREFIX*(strl),filename*(strl)
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)
      double precision cRPIX_b(2),cD_b(2,2),PU_b(2,npd),cRVAL_b(2)
      integer n,i,j,k
      double precision ra,dec,x,y,ra2,dec2

      if (ASTROMETRY_trivial.eq.1) then
        call get_astrometry_trivial(IMAGE_FILE,nchip,DIR_OUTPUT)
      else
        call get_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
c        if (ext_cat.eq.1) call get_astrometry_gal(IMAGE_FILE,nchip
c     .,DIR_OUTPUT)
      endif

      if (ASTROMETRY_trivial.eq.1) return

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_check.dat'

      open(unit=50,file=trim(filename),status='replace')
      rewind 50

      do ichip=1,nchip
        proc_error=0

        call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head'

        call read_astrometry_para(filename,ichip    
     .,cRPIX,cD,cRVAL,PU,npd,proc_error)

        if (proc_error.eq.0) then

          call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
          filename=trim(DIR_OUTPUT)//'/stamps/'
     .//trim(PREFIX)//'_norm.fits'
          call update_para(filename,cRPIX,cD)

          filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'
          open(unit=40,file=trim(filename),status='old')
          rewind 40
          read(40,*)
          read(40,*)
          read(40,*) n,j,k
          do i=1,n
            read(40,*) ra,dec,x,y
            call coordinate_transfer_PU(ra2,dec2,x,y,1
     .,cRPIX,cD,cRVAL,PU,npd)
            write(50,*) ra,dec,ra2,dec2
          enddo
          close(40)
        endif
      enddo

      close(50)


      return
      END
`,

  source: `
      subroutine chip_process_source(IMAGE_FILE,ichip,DIR_OUTPUT
     .,FLAT_FILE)
      implicit none
      include 'para.inc'

      integer ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT,FLAT_FILE
      character*(strl) catfile,PREFIX,PREFIX_head,filename
      character*(strl) sortfile(27)
      integer sortnum
      integer proc_error

      integer nx,ny
      real array(npx,npy),sigmap(npx,npy),normap(npx,npy)
      real flat(npx,npy)

      integer weight(npx,npy)

      integer i,j,u,v,ii
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)

      integer ngal,nstar,nxc

      real sigabc(2,3)


      proc_error=0
      nstar=0

      if (include_FLAT.eq.2)                    
     . call readimage(FLAT_FILE,nx,ny,npx,npy,flat)
      call readimage(IMAGE_FILE(ichip),nx,ny,npx,npy,array)
      call get_PREFIX(IMAGE_FILE(ichip),PREFIX)

      PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
      filename=trim(PREFIX)//'_norm.fits'
      call readimage(filename,nx,ny,npx,npy,normap)
      do i=1,ccD_split
        do j=1,3
          sigabc(i,j)=normap(1+i,j)
        enddo
      enddo
      if (normap(1,1).gt.0) proc_error=1
      do i=1,nx
        do j=1,ny
          weight(i,j)=1
          if (normap(i,j).lt.(-900.)) weight(i,j)=0
          if (include_FLAT.eq.2) then
            if (abs(flat(i,j)).le.flat_thresh) then
              array(i,j)=array(i,j)*(1.-flat(i,j))
            else
              weight(i,j)=0
            endif
          endif
        enddo
      enddo

      nxc=nx/2
      if (ccD_split.eq.2) then
        do i=1,nxc
          ii=i+nxc
          do j=1,ny
            sigmap(i,j)=sigabc(1,1)+sigabc(1,2)*i+sigabc(1,3)*j
            sigmap(i,j)=sqrt(0.5*sigmap(i,j))
            sigmap(ii,j)=sigabc(2,1)+sigabc(2,2)*ii+sigabc(2,3)*j
            sigmap(ii,j)=sqrt(0.5*sigmap(ii,j))
          enddo
        enddo
      else
        do i=1,nx
          do j=1,ny
            sigmap(i,j)=sigabc(1,1)+sigabc(1,2)*i+sigabc(1,3)*j
            sigmap(i,j)=sqrt(0.5*sigmap(i,j))
          enddo
        enddo
      endif

c------------------------------------------------------
      call get_expo_catalog(PREFIX,nx,ny,sigmap,weight,normap   
     .,proc_error)

      if (ext_cat.eq.0) then

        call gen_source_catalog(PREFIX,nx,ny,array,weight
     .,ngal,proc_error)

        if (ext_PSF.eq.0) call gen_star_candidate(PREFIX,nstar
     .,proc_error)

      else

        call get_PREFIX_expo(IMAGE_FILE(1),PREFIX_head)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX_head)//'.head'
        call read_astrometry_para(filename,ichip    
     .,cRPIX,cD,cRVAL,PU,npd,proc_error)

        catfile=SOURcE_cAT

        if (proc_error.eq.0)
     . call generate_gal_cat_file_name(cRVAL,catfile,sortfile,sortnum)

        if (deblending.eq.1) call de_blending(sortfile,sortnum,nx,ny  
     .,weight,cRPIX,cD,cRVAL,PU,proc_error)


        call gen_source_ext_catalog(sortfile,sortnum               
     .,PREFIX,nx,ny,array,weight,sigmap,cRPIX,cD,cRVAL,PU
     .,ngal,proc_error)

        if (ext_PSF.eq.0) call gen_star_candidate_direct(PREFIX    
     .,nx,ny,array,weight,nstar,proc_error)

      endif

      write(*,*) trim(IMAGE_FILE(ichip)),proc_error,nstar,ngal

      return
      END
`,

  fourier: `
      subroutine chip_process_Fourier_T(IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      character*(strl) IMAGE_FILE,DIR_OUTPUT

      integer ierror
      integer i,j,u,v,nn1,nn2
      character*(strl) PREFIX,filename

      integer nsource
      real source_para(ngal_max,npara)
      real source_coll(ngal_max,ns,ns)
      real noise_coll(ngal_max,ns,ns)
      real power_coll(ngal_max,ns,ns)
      real power_ori_coll(ngal_max,ns,ns)

      real source_p(ns,ns),source(ns,ns)
      real noise_p(ns,ns),noise(ns,ns)

      real aa(npara),flux_alt,temp

      real pc
      common /pc_pass/ pc

      call get_PREFIX(IMAGE_FILE,PREFIX)
      PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)

c------------------------------------------------------------------
      nsource=0
      filename=trim(PREFIX)//'_source_info.dat'
      open(unit=10,file=filename,status='old',iostat=ierror)
      rewind 10
      if (ierror.ne.0) then
        write(*,*) filename
        write(*,*) 'catalog file error!!'
        read(*,*)
      endif
      read(10,*)
c      read(10,*) 'ig xc yc sigma peak imax jmax'
c     .,'half_light_flux half_light_area flag'
      do while (ierror.ge.0)
        read(10,*,iostat=ierror) (aa(i),i=1,iflag)
        if (ierror.lt.0) cycle
        nsource=nsource+1
        do i=1,iflag
          source_para(nsource,i)=aa(i)
        enddo
      enddo
      close(10)

      if (nsource.gt.0) then
        nn1=ns*len_g
        nn2=ns*(int(nsource/len_g)+1)
        filename=trim(PREFIX)//'_source.fits'
        call read_stamps(ngal_max,1,nsource,ns,ns
     .,source_coll,nn1,nn2,filename)

        filename=trim(PREFIX)//'_noise.fits'
        call read_stamps(ngal_max,1,nsource,ns,ns
     .,noise_coll,nn1,nn2,filename)

        do i=1,nsource

          do u=1,ns
            do v=1,ns
              source(u,v)=source_coll(i,u,v)
              noise(u,v)=noise_coll(i,u,v)
            enddo
          enddo

          call get_power(ns,ns,source,source_p,2)
          source_para(i,11)=sqrt(max(pc,source_p(ns_2+1,ns_2+1)))
          source_para(i,12)=source_para(i,11)/source_para(i,4)*ns

          call get_power(ns,ns,source,source_p,gal_smooth)
          call get_power(ns,ns,noise,noise_p,gal_smooth)

          call process_powers(ns,source_p,noise_p)
          do u=1,ns
            do v=1,ns
              power_coll(i,u,v)=source_p(u,v)
            enddo
          enddo
        enddo
      endif

      filename=trim(PREFIX)//'_source_info.dat'
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) 'ig xp yp sigma peak imax jmax '
     .,'half_light_flux half_light_area flag flux2 SNR_F'
      if (nsource.gt.0) then
        do i=1,nsource
          write(10,*) (source_para(i,j),j=1,iSNR_F)
        enddo
        close(10)

        nn1=ns*len_g
        nn2=ns*(int(nsource/len_g)+1)
        filename=trim(PREFIX)//'_source_p.fits'
        call write_stamps(ngal_max,1,nsource,ns,ns
     .,power_coll,nn1,nn2,filename)
      else
        close(10)
      endif

c-------------------------------------------------------------------
      if (ext_PSF.eq.1) return
c-------------------------------------------------------------------
      nsource=0
      filename=trim(PREFIX)//'_star_can_info.dat'
      open(unit=10,file=filename,status='old',iostat=ierror)
      rewind 10
      if (ierror.ne.0) then
        write(*,*) filename
        write(*,*) 'catalog file error!!'
        read(*,*)
      endif
      read(10,*)
c      read(10,*) 'ig xp yp SNR'
      do while (ierror.ge.0)
        read(10,*,iostat=ierror) (aa(i),i=1,4)
        if (ierror.lt.0) cycle
        nsource=nsource+1
      enddo
      close(10)

      if (nsource.gt.0) then 
        nn1=ns*len_s
        nn2=ns*(int(nsource/len_s)+1)
        filename=trim(PREFIX)//'_star_can.fits'
        call read_stamps(ngal_max,1,nsource,ns,ns
     .,source_coll,nn1,nn2,filename)

        filename=trim(PREFIX)//'_star_can_noise.fits'
        call read_stamps(ngal_max,1,nsource,ns,ns
     .,noise_coll,nn1,nn2,filename)

        do i=1,nsource
          do u=1,ns
            do v=1,ns
              source(u,v)=source_coll(i,u,v)
              noise(u,v)=noise_coll(i,u,v)
            enddo
          enddo
          call get_power(ns,ns,source,source_p,star_smooth)
          call get_power(ns,ns,noise,noise_p,star_smooth)
          call process_powers(ns,source_p,noise_p)
          call regularize_power(ns,ns,source_p,star_smooth)
          do u=1,ns
            do v=1,ns
              power_coll(i,u,v)=source_p(u,v)
            enddo
          enddo
c        call get_power(ns,ns,source,source_p,0)
c        call get_power(ns,ns,noise,noise_p,0)
c        call process_powers(ns,source_p,noise_p)
c        call regularize_power(ns,ns,source_p,0)
c        do u=1,ns
c          do v=1,ns
c            power_ori_coll(i,u,v)=source_p(u,v)
c          enddo
c        enddo
        enddo

        filename=trim(PREFIX)//'_star_can_power.fits'
        nn1=ns*len_s
        nn2=ns*(int(nsource/len_s)+1)
        call write_stamps(ngal_max,1,nsource,ns,ns
     .,power_coll,nn1,nn2,filename)
      endif

c      filename=trim(PREFIX)//'_star_ori_power.fits'
c      nn1=ns*len_s
c      nn2=ns*(int(nsource/len_s)+1)
c      call write_stamps(ngal_max,1,nsource,ns,ns
c     .,power_ori_coll,nn1,nn2,filename)

c----------------------------------------------------
      write(*,*) trim(PREFIX),nsource



      return
      END
`,

  psf: `
      subroutine proc_PSF(iexpo)
      implicit none
      include 'para.inc'

      integer iexpo,nchip,ichip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer nc
      real p_chip(NMAX_cHIP,4)

      call get_image_list(iexpo,IMAGE_FILE,nchip,DIR_OUTPUT)

      call read_in_candidates(nchip,IMAGE_FILE,DIR_OUTPUT,nc,p_chip)

      call star_selection(nchip)

      call plot_star_expo(nchip,IMAGE_FILE,DIR_OUTPUT)

      call plot_stars(nchip,IMAGE_FILE,DIR_OUTPUT
     .,nc,p_chip)

      if (PSF_type.eq.1) then
        call make_PSF_local_fit(nchip,IMAGE_FILE,DIR_OUTPUT)
      elseif (PSF_type.eq.2) then
        call make_PSF_hybrid(nchip,IMAGE_FILE,DIR_OUTPUT)
      else
        pause 'Invalid PSF fitting method!'
      endif


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_in_candidates(nchip,IMAGE_FILE,DIR_OUTPUT
     .,nc,p_chip)
      implicit none
      include 'para.inc'

      integer nchip,nc
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename,headname

      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      common /star_info_pass/ star_para,nstar

      integer ierror,ntot,i,nn1,nn2,u,v,k,j
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)
      double precision x,y,xx,yy,step
      real p_chip(NMAX_cHIP,4),aa(npara),source_p(ns,ns)
      real star(nstar_max,ns,ns),ee(2),size,temp,FWHM
      real map1(ns,ns),map2(ns,ns)

      real chi_d(NMAX_cHIP,nstar_max,nstar_max)
      common /chi_d_pass/ chi_d

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      headname=trim(DIR_OUTPUT)//'/astrometry/'//trim(PREFIX)//'.head'

      nc=0
      do k=1,nchip

        nstar(k)=0

        ierror=0
        call read_astrometry_para(headname,k,cRPIX,cD,cRVAL,PU
     .,npd,ierror)

        if (ierror.eq.1) cycle
        nc=nc+1
        x=1d0
        y=1d0
        call xy_to_xxyy(x,y,xx,yy,cRPIX,cD)
        p_chip(nc,1)=xx
        p_chip(nc,2)=yy
        x=2046d0
        y=4094d0
        call xy_to_xxyy(x,y,xx,yy,cRPIX,cD)
        p_chip(nc,3)=xx
        p_chip(nc,4)=yy

        call get_PREFIX(IMAGE_FILE(k),PREFIX)
        PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
        filename=trim(PREFIX)//'_star_can_info.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) filename
          stop 'catalog file error!!'
        endif
        read(10,*)
c        read(10,*) 'ig xp yp SNR'
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) (aa(i),i=1,4)
          if (ierror.lt.0) cycle
          nstar(k)=nstar(k)+1
          do i=1,4
            star_para(k,nstar(k),i)=aa(i)
          enddo
        enddo
        close(10)

        if (nstar(k).gt.0) then
          nn1=ns*len_s
          nn2=ns*(int(nstar(k)/len_s)+1)
          filename=trim(PREFIX)//'_star_can_power.fits'
          call read_stamps(nstar_max,1,nstar(k),ns,ns,star
     .,nn1,nn2,filename)

          do i=1,nstar(k)
            star_para(k,i,5)=1.
            do u=1,ns
              do v=1,ns
                source_p(u,v)=star(i,u,v)
              enddo
            enddo

            x=star_para(k,i,2)
            y=star_para(k,i,3)
            call xy_to_xxyy(x,y,xx,yy,cRPIX,cD)
            star_para(k,i,6)=xx
            star_para(k,i,7)=yy
            call get_power_all(ns,ns,source_p,ee,size,0.02)
            star_para(k,i,8)=size
            star_para(k,i,9)=ee(1)
            star_para(k,i,10)=ee(2)
            call get_PSF_FWHM(source_p,FWHM)
            star_para(k,i,11)=FWHM

            temp=0.
            do u=1,ns
              do v=1,ns
                temp=temp+source_p(u,v)
              enddo
            enddo
            star_para(k,i,12)=1./temp
          enddo

          do i=1,nstar_max
            do j=1,nstar_max
              chi_d(k,i,j)=0.
            enddo
          enddo

          do i=1,nstar(k)-1
            do j=i+1,nstar(k)
              do u=1,ns
                do v=1,ns
                  map1(u,v)=star(i,u,v)*star_para(k,i,12)
                  map2(u,v)=star(j,u,v)*star_para(k,j,12)
                enddo
              enddo
              call ana_chi2(ns,map1,map2,temp)
              temp=sqrt(temp)
              chi_d(k,i,j)=temp
              chi_d(k,j,i)=temp
            enddo
          enddo
        endif

c        write(*,*) 'reading star cans:',ntot
      enddo

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine star_selection(nchip)
      implicit none
      include 'para.inc'

      integer nchip
      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      common /star_info_pass/ star_para,nstar

      integer npsam
      parameter (npsam=NMAX_cHIP*nstar_max*nstar_max)

      integer i,j,k,u,v,ntot,w,n,u2,v2
      real tmp(npsam),temp
      real thresh,peak,sig,chimin(NMAX_cHIP,nstar_max)
      real chi_d(NMAX_cHIP,nstar_max,nstar_max)
      common /chi_d_pass/ chi_d

      integer id(nstar_max),group_size(nstar_max)
      integer group_id(nstar_max),max_group_id
      integer max_gsize,max_gid,max2_gsize,max2_gid

      ntot=0
      do k=1,nchip
        ntot=ntot+nstar(k)
        do i=1,nstar_max
          chimin(k,i)=1000.
        enddo
      enddo
      if (ntot.lt.nstar_min*2) then
        do k=1,nchip
          do i=1,nstar(k)
            star_para(k,i,5)=-1
          enddo
        enddo
        return
      endif


      ntot=0
      do k=1,nchip
        do i=1,nstar(k)
          ntot=ntot+1
          tmp(ntot)=star_para(k,i,8)
        enddo
      enddo
      call sort(ntot,npsam,tmp)
      thresh=tmp((ntot*2)/3)

!----Determine the threshold for chi^2 ------------------------------

      ntot=0
      do k=1,nchip
        do i=1,nstar(k)-1
          do j=i+1,nstar(k)
            temp=chi_d(k,i,j)
            chimin(k,i)=min(temp,chimin(k,i))
            chimin(k,j)=min(temp,chimin(k,j))
            if (star_para(k,i,8).lt.thresh) cycle
            if (star_para(k,j,8).lt.thresh) cycle
            ntot=ntot+1
            tmp(ntot)=temp
          enddo
        enddo
      enddo

      call get_peak_width_low_side(npsam,ntot,tmp,peak,sig)

      thresh=peak+4.*sig

      do k=1,nchip
        ntot=0
        do i=1,nstar(k)
          if (chimin(k,i).gt.thresh) then
            star_para(k,i,5)=-1
            cycle
          endif
          ntot=ntot+1
        enddo
        if (ntot.lt.nstar_min_local) then
          do i=1,nstar(k)
            star_para(k,i,5)=-1
          enddo
          cycle
        endif

        ntot=0
        do i=1,nstar(k)
          if (star_para(k,i,5).lt.0) cycle
          ntot=ntot+1
          id(ntot)=i
        enddo
        do i=1,ntot
          group_id(i)=0
          group_size(i)=0
        enddo
        max_group_id=0
        do i=1,ntot
          if (group_id(i).eq.0) then
            max_group_id=max_group_id+1
            group_id(i)=max_group_id
          endif
          do j=i+1,ntot
            if (group_id(j).eq.group_id(i)) cycle
            if (chi_d(k,id(i),id(j)).gt.thresh) cycle
            if (group_id(j).eq.0) then
              group_id(j)=group_id(i)
            else
              u=group_id(j)
              v=group_id(i)
              do w=1,ntot
                if (group_id(w).eq.u) group_id(w)=v
              enddo
            endif
          enddo
        enddo
!-----------------Find the largest group ---------------------------
        do i=1,ntot
          u=group_id(i)
          group_size(u)=group_size(u)+1
        enddo

        if (nstar(k).gt.0) then
          max_gsize=group_size(1)
          max_gid=1
          max2_gsize=0
          max2_gid=0

          do i=2,max_group_id
            if (group_size(i).gt.max_gsize) then
              max2_gsize=max_gsize
              max2_gid=max_gid
              max_gsize=group_size(i)
              max_gid=i
            elseif (group_size(i).gt.max2_gsize) then
              max2_gsize=group_size(i)
              max2_gid=i
            endif
          enddo
        endif
c        if (max2_gsize.lt.nstar_min_local/2) then
          do i=1,ntot
            j=id(i)
            if (group_id(i).ne.max_gid) star_para(k,j,5)=-1
          enddo
c        else
c          do i=1,ntot
c            j=id(i)
c            if (group_id(i).ne.max_gid.and.group_id(i).ne.max2_gid)
c    .star_para(k,j,5)=-1
c          enddo
c        endif

      enddo

      do k=1,nchip
        n=0
        do i=1,nstar(k)
          if (star_para(k,i,5).lt.0) cycle
          n=n+1
        enddo
        if (n.lt.nstar_min_local) then
          do i=1,nstar(k)
            star_para(k,i,5)=-1
          enddo
        endif
      enddo


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plot_star_expo(nchip,IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename

      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      real star_test(NMAX_cHIP*nstar_max,ns,ns)
      common /star_info_pass/ star_para,nstar

      integer nmax_stamp,opt(NMAX_cHIP*nstar_max)
      parameter (nmax_stamp=5000)

      integer ntot ,ichip, w, nn1, nn2, start
      integer k, i, u, v
      real star(nstar_max,ns,ns)

      real chi_d(NMAX_cHIP,nstar_max,nstar_max)
      common /chi_d_pass/ chi_d

      do i=1,nstar_max*NMAX_cHIP
        opt(i)=0
      enddo

      ntot = 0
      w = 0
      start = 0
      do ichip=1,nchip
        if (nstar(ichip).eq.0) cycle
        nn1=ns*len_s
        nn2=ns*(int(nstar(ichip)/len_s)+1)

        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
        filename=trim(PREFIX)//'_star_can_power.fits'
        call read_stamps(nstar_max,1,nstar(ichip),ns,ns
     .,star,nn1,nn2,filename)
        do i=1,nstar(ichip)
          do u=1,ns
            do v=1,ns
              star_test(start+i,u,v)=star(i,u,v)
            enddo
          enddo
          w = start + i
          if (star_para(ichip,i,5).le.0) cycle
          ntot = ntot + 1
          if (ntot.lt.nmax_stamp) opt(w)=1
        enddo
        start = start + nstar(ichip)
      enddo

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)

      if (ntot.gt.0) then 
        nn1=ns*len_sam
        nn2=ns*(int(min(ntot,nmax_stamp)/len_sam)+1)
        filename=trim(PREFIX)//'_star_power_expo.fits'
        call write_stamps_2(nstar_max*NMAX_cHIP,w,ns,ns,star_test
     .,opt,1,nn1,nn2,filename)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plot_stars(nchip,IMAGE_FILE
     .,DIR_OUTPUT,nc,p_chip)
      implicit none
      include 'para.inc'

      integer nchip,nc
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename
      real p_chip(NMAX_cHIP,4)

      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      common /star_info_pass/ star_para,nstar

      integer nm
      parameter (nm=1000)
      real PSFmap(nm,nm),sk(NMAX_cHIP*nstar_max,5),source_p(ns,ns)
      integer nmax_stamp,opt(NMAX_cHIP*nstar_max)
      parameter (nmax_stamp=5000)

      integer i,j,k,u,v,ntot,w,nums,nn1,nn2
      real FWHM,FWHM_ave,e1_ave,e2_ave,chi_d_ave
      real chi_d(NMAX_cHIP,nstar_max,nstar_max)
      common /chi_d_pass/ chi_d


      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      PREFIX=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)

      filename=trim(PREFIX)//'_star_info_expo.dat'
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) '# ichip nstar FWHM e1 e2 chi_d'

      ntot=0

      do k=1,nchip

        FWHM_ave=0.
        e1_ave=0.
        e2_ave=0.
        chi_d_ave=0.
        nums=0

        do i=1,nstar(k)
          if (star_para(k,i,5).le.0) cycle
          nums=nums+1
          ntot=ntot+1
          sk(ntot,1)=star_para(k,i,6)
          sk(ntot,2)=star_para(k,i,7)
          sk(ntot,3)=star_para(k,i,8)
          sk(ntot,4)=star_para(k,i,9)
          sk(ntot,5)=star_para(k,i,10)
          FWHM=star_para(k,i,11)

          FWHM_ave=FWHM_ave+FWHM
          e1_ave=e1_ave+star_para(k,i,9)
          e2_ave=e2_ave+star_para(k,i,10)
          if (nums.ge.2) chi_d_ave=chi_d_ave+chi_d(k,i,j)
          j=i
        enddo

        if (nums.ge.nstar_min_local) then
          FWHM_ave=FWHM_ave/nums
          e1_ave=e1_ave/nums
          e2_ave=e2_ave/nums
          chi_d_ave=chi_d_ave/(nums-1.)
          write(10,*) k,nums,FWHM_ave,e1_ave,e2_ave,chi_d_ave
        else
          write(10,*) k,0,-99.,-99.,-99.,-99.
        endif

      enddo

      write(*,*) trim(PREFIX),' total no. of stars:',ntot
      close(10)


      call draw_shear_expo(nm,PSFmap,nchip,NMAX_cHIP,nc,p_chip
     .,NMAX_cHIP*nstar_max,ntot,sk,200.,1.)
      filename=trim(PREFIX)//'_PSF_source.fits'
      call writeimage(filename,nm,nm,nm,nm,PSFmap)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine make_PSF_local_fit(nchip,IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      character*(strl) PREFIX,filename

      integer nstar(NMAX_cHIP)
      double precision star_para(NMAX_cHIP,nstar_max,npara)
      common /star_info_pass/ star_para,nstar

      integer ntot,nc,k,i,j,w,u,v,nums
      double precision posi(nstar_max,2)
      double precision xx,yy,PSF_coe_l(ns,ns,npl)
      real model(ns,ns),px,py,sshape(nstar_max,3),msshape(nstar_max,3)
      real ee(2),size,star(nstar_max,ns,ns)

      integer nn1,nn2

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      filename=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)//
     .'_star_comp_expo.dat'
      open(unit=90,file=filename,status='replace')
      rewind 90

      ntot=0
      do k=1,nchip
        nums=0

        if (nstar(k).gt.0) then
          call get_PREFIX(IMAGE_FILE(k),PREFIX)
          nn1=ns*len_s
          nn2=ns*(int(nstar(k)/len_s)+1)
          filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
     .//'_star_can_power.fits'
          call read_stamps(nstar_max,1,nstar(k),ns,ns,star
     .,nn1,nn2,filename)
        endif
        
        filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)//
     .'_PSF_coe_local.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10

        do i=1,nstar(k)
          if (star_para(k,i,5).lt.0) cycle
          ntot=ntot+1
          nums=nums+1
          posi(nums,1)=star_para(k,i,2)
          posi(nums,2)=star_para(k,i,3)
          sshape(nums,1)=star_para(k,i,8)
          sshape(nums,2)=star_para(k,i,9)
          sshape(nums,3)=star_para(k,i,10)
          do u=1,ns
            do v=1,ns
              star(nums,u,v)=star(i,u,v)
            enddo
          enddo
        enddo

        if (nums.ge.nstar_min_local) then
          call interpolate_PSF(nums,nstar_max,star,posi,ns
     .,npl,nplx,PSF_coe_l)
          write(10,*) nums,1
          do i=1,ns
            do j=1,ns
              write(10,*) (PSF_coe_l(i,j,u),u=1,npl)
            enddo
          enddo
          write(90,*) k,nums,1
          do i=1,nums
            xx=posi(i,1)
            yy=posi(i,2)
            call get_PSF_model(ns,npl,nplx,PSF_coe_l,xx,yy,model)
            call get_power_all(ns,ns,model,ee,size,0.02)
            msshape(i,1)=size
            msshape(i,2)=ee(1)
            msshape(i,3)=ee(2)
            px=posi(i,1)
            py=posi(i,2)
            write(90,*) px,py,(sshape(i,u),u=1,3),(msshape(i,v),v=1,3)
          enddo
        else
          write(10,*) nums,-1
          write(90,*) k,nums,-1
        endif
        close(10)
      enddo

      close(90)

      return
      end
`,

  shear: `
      subroutine expo_shear(nchip,IMAGE_FILE,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer ichip,nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      integer proc_error,nstar,status,nstar_tot
      character*(strl) PREFIX,headname,filename,PREFIX1,PREFIX2
      character*(strl) PREFIX_head
      integer ngal,ierror,nn1,nn2,i,j,u,v,nx,ny,k
      real g1,g2,de,h1,h2
      real gal_p_coll(ngal_max,ns,ns),gal_para(ngal_max,npara)
      double precision cRPIX(2),cD(2,2),cRVAL(2)
      double precision PU(2,npd)
      real gal_p(ns,ns),psf_model(ns,ns),aa(npara)
      double precision ra,dec,gf1,gf2,temp,cos2,sin2,cos4,sin4
      double precision x,y,xx,yy
      integer parity
      double precision local_coe(ns,ns,npl),PSF_coe(ns,ns,npo)
      real ePSF(ns,ns),ePSF_p(ns,ns),psf_FWHM,psfmap(npx,npy)
      real px,py,dstar

      proc_error=0

      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      headname=trim(DIR_OUTPUT)//'/astrometry/'//trim(PREFIX)//'.head'

      if (ext_PSF.eq.1) then
        filename=trim(PSF_PATH)//'/PSF.fits'
        call readimage(filename,nx,ny,ns,ns,ePSF)
        call get_power(ns,ns,ePSF,ePSF_p,0)
        do i=1,ns
          do j=1,ns
            local_coe(i,j,1)=ePSF_p(i,j)
          enddo
        enddo
      endif


      do ichip=1,nchip
        proc_error=0
        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        PREFIX1=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)
        PREFIX2=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)
        if (ext_PSF.ne.1 .and. PSF_type.eq.1) then
          filename=trim(PREFIX1)//'_PSF_coe_local.dat'
          open(unit=10,file=filename,status='old',iostat=ierror)
          rewind 10
          read(10,*) nstar,status
          if (status.eq.-1) then
            proc_error=1
          else
            do i=1,ns
              do j=1,ns
                read(10,*) (local_coe(i,j,k),k=1,npl)
              enddo
            enddo
          endif
          close(10)
        elseif (ext_PSF.ne.1 .and. PSF_type.eq.2) then
          filename=trim(PREFIX1)//'_PSF_local.fits'
          call readimage(filename,nx,ny,npx,npy,psfmap)
          nstar=int(psfmap(step_psf-1,step_psf-1)+0.5)
          if (psfmap(step_psf,step_psf).lt.-1.) proc_error=1
        endif
        if (proc_error.eq.1) then
          ngal=0
          goto 50
        endif
        call read_astrometry_para(headname,ichip,cRPIX,cD,cRVAL,PU,npd
     .,proc_error)
        if (proc_error.eq.1) then
          ngal=0
          goto 50
        endif

        ngal=0
        filename=trim(PREFIX1)//'_source_info.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) filename
          stop 'catalog file error in shear_proc!!'
        endif
        read(10,*)
c        read(10,*) 'ig xc yc sigma peak imax jmax half_light_flux half_light_area flag flux2 SNR_F'
        do while (ierror.ge.0)
          read(10,*,iostat=ierror) (aa(i),i=1,iSNR_F)
          if (ierror.lt.0) cycle
          ngal=ngal+1
          do i=1,iSNR_F
            gal_para(ngal,i)=aa(i)
          enddo
        enddo
        close(10)

        if (ngal.eq.0) then
          goto 50
        endif

        nn1=ns*len_g
        nn2=ns*(int(ngal/len_g)+1)
        filename=trim(PREFIX1)//'_source_p.fits'
        call read_stamps(ngal_max,1,ngal,ns,ns,gal_p_coll
     .,nn1,nn2,filename)

50        filename=trim(PREFIX2)//'_shear.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) 'ig xc yc sigma nstar imax jmax ' 
     .,'half_light_flux half_light_area flag psf_FWHM SNR_F '  
     .,'ra dec gf1 gf2 g1 g2 de h1 h2 cos2 sin2 parity'

        do i=1,ngal
          do u=1,ns
            do v=1,ns
              gal_p(u,v)=gal_p_coll(i,u,v)
            enddo
          enddo
          x=gal_para(i,2)
          y=gal_para(i,3)
          if (ext_PSF.eq.1) then
            call get_PSF_model(ns,1,1,local_coe,x,y,psf_model)
          else
            if (PSF_type.eq.1) then
              call get_PSF_model(ns,npl,nplx,local_coe,x,y,psf_model)
            elseif (PSF_type.eq.2) then
              px=x
              py=y
              call get_PSF_model_very_local(psfmap,px,py
     .,psf_model,dstar)
            endif
          endif
          call get_PSF_area(psf_model,psf_FWHM)

          gal_para(i,iPSF)=psf_FWHM
          gal_para(i,istar)=nstar
          call coordinate_transfer_PU(ra,dec,x,y,1,cRPIX,cD
     .,cRVAL,PU,npd)
          gal_para(i,ira)=ra
          gal_para(i,idec)=dec
          call field_distortion_PU(x,y,npd,PU,cD,cRPIX,gf1,gf2
     .,cos2,sin2,parity)
          gal_para(i,igf1)=gf1
          gal_para(i,igf2)=gf2
          call get_shear(ns,gal_p,psf_model,g1,g2,de,h1,h2)
          gal_para(i,ig1)=g1*cos2+g2*sin2
          gal_para(i,ig2)=g2*cos2-g1*sin2
          gal_para(i,ide)=de
          cos4=cos2*cos2-sin2*sin2
          sin4=2d0*sin2*cos2
          gal_para(i,ih1)=h1*cos4+h2*sin4
          gal_para(i,ih2)=h2*cos4-h1*sin4
          if (parity.eq.-1) then
            gal_para(i,ig2)=-gal_para(i,ig2)
            gal_para(i,ih2)=-gal_para(i,ih2)
          endif
          gal_para(i,icos2)=cos2
          gal_para(i,isin2)=sin2
          gal_para(i,iparity)=parity

          write(10,*) (gal_para(i,j),j=1,iparity)

        enddo
        close(10)
      enddo

      return
      END
`,

  combine: `
      subroutine combine_expo_catalog(nchip,IMAGE_FILE
     .,DIR_OUTPUT,chi2)
      implicit none
      include 'para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT

      real chi2
      character*(strl) PREFIX,filename
      character*1000 cat_content,cat_list1,cat_list2
      integer ichip,ierror,u,i,m,n,chip_index
      real cat(npara),g1c,g2c


      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      filename=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)//'_all.cat'
      open(unit=20,file=filename,status='replace',iostat=ierror)
      rewind 20

      if (ext_cat.eq.1) then
        do ichip=1,nchip
          call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
          filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)//
     .'_orig.cat'
          open(unit=15,file=filename,status='old',iostat=ierror)
          rewind 15
          if (ierror.ne.0) then
            write(*,*) trim(filename),' is missing!'
            stop
          endif
          read(15,'(A)') cat_list2
          read(15,'(A)',iostat=ierror)  cat_content
          if (ierror.lt.0) then
            close(15)
            cycle
          else
            close(15)
            exit
          endif
        enddo
      endif

      n=0
      m=0
      do ichip=1,nchip
        call get_chip_id(IMAGE_FILE(ichip),chip_index)
        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        filename=trim(DIR_OUTPUT)//'/result/'//trim(PREFIX)//
     .'_shear.dat'
        open(unit=10,file=filename,status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) trim(filename),' is missing!'
          stop
        endif
        read(10,'(A)') cat_list1

        if (ext_cat.eq.1) then
          filename=trim(DIR_OUTPUT)//'/stamps/'//trim(PREFIX)//
     .'_orig.cat'
          open(unit=15,file=filename,status='old',iostat=ierror)
          rewind 15
          if (ierror.ne.0) then
            write(*,*) trim(filename),' is missing!'
            stop
          endif
          read(15,*)
          if (ichip.eq.1) then
            write(20,*) trim(cat_list2),' ccD_NUM ',trim(cat_list1)
     .,'Chi2'
            if (chi2.gt.chi2_thresh) then
              close(10)
              close(15)
              close(20)
              write(*,*) trim(PREFIX)//' contains no valid sources!'
              return
            endif
          endif
          do while (ierror.ge.0)
            read(10,*,iostat=ierror) (cat(u),u=1,iparity)
            if (ierror.lt.0) cycle
            read(15,'(A)') cat_content
            if (cat(i_imax).ge.ns.or.cat(i_jmax).ge.ns) then
              m=m+1
              cycle
            endif
            n=n+1
            ! g1c=cat(igf1)+g1_c
            ! g2c=cat(igf2)+g2_c
            g1c = 0.
            g2c = 0.
            cat(ig1)=cat(ig1)-g1c*cat(ide)+g1c*cat(ih1)+g2c*cat(ih2)
            cat(ig2)=cat(ig2)-g2c*cat(ide)+g1c*cat(ih2)-g2c*cat(ih1)
            write(20,*) trim(cat_content),chip_index
     .,(cat(u),u=1,iparity),chi2
          enddo
          close(15)
        else
          if (ichip.eq.1) then
            write(20,*) ' ccD_NUM ',trim(cat_list1)
            if (chi2.gt.chi2_thresh) then
              close(10)
              close(20)
              write(*,*) trim(PREFIX)//' contains no valid sources!'
              return
            endif
          endif
          do while (ierror.ge.0)
            read(10,*,iostat=ierror) (cat(u),u=1,iparity)
            if (ierror.lt.0) cycle
            if (cat(i_imax).ge.ns.or.cat(i_jmax).ge.ns) then
              m=m+1
              cycle
            endif
            n=n+1
            g1c=cat(igf1)+g1_c
            g2c=cat(igf2)+g2_c
            cat(ig1)=cat(ig1)-g1c*cat(ide)+g1c*cat(ih1)+g2c*cat(ih2)
            cat(ig2)=cat(ig2)-g2c*cat(ide)+g1c*cat(ih2)-g2c*cat(ih1)
            write(20,*) chip_index,(cat(u),u=1,iparity)
          enddo
        endif
        close(10)
      enddo
      write(*,*) trim(PREFIX),n,m
      close(20)

      return
      end
`,
};
export default codeMapMain;