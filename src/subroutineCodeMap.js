const functionCodeMap = {
  initialize: `
      subroutine initialize(EXPO_LIST)
      implicit none
      include 'para.inc'

      character*(strl) EXPO_LIST
      character*(strl) EXPO_FILE(NMAX_EXPO)
      integer N_EXPO
      common /filename_pass/ EXPO_FILE,N_EXPO

      character*(strl) expo_name
      integer ierror,nchip

      N_EXPO=0

      open(unit=10,file=EXPO_LIST,status='old',iostat=ierror)
      rewind 10
      if (ierror.ne.0) then
        write(*,*) 'EXPO_LIST reading error!!'
        stop
      endif
      do while (ierror.ge.0)
        read(10,*,iostat=ierror) expo_name , nchip
        if (ierror.lt.0) cycle
        N_EXPO=N_EXPO+1
        EXPO_FILE(N_EXPO)=trim(expo_name)
      enddo
      close(10)
      write(*,*) 'Total number of EXPOSURE: ',N_EXPO


      return
      end
  `,
  mpi_distribute:`
      subroutine mpi_distribute(num_job,job,message)
      use mpi
      implicit none

      integer my_id,num_procs
      common /MPIpar/ my_id,num_procs

      external job

      character*(*) message

      integer num_job,i,j,k,complete
      integer source,tag,ierr,status(mpi_status_size)


      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) ! synchronize all nodes

      i=0
      if (my_id.eq.0) then
        i=1
        j=num_job
      endif
      complete=0
      do while (complete.eq.0)
        if (my_id.ne.0) then
          call MPI_SEND(i,1,mpi_int,0,0,MPI_cOMM_WORLD,ierr)
          call MPI_REcV(i,1,mpi_int,0,0,MPI_cOMM_WORLD,status,ierr)
          if (i.eq.0) then
            complete=1
          else
            call job(i)
          endif
        else
          call MPI_REcV(k,1,mpi_int,MPI_ANY_SOURcE
     .,MPI_ANY_TAG,MPI_cOMM_WORLD,status,ierr)
          tag=status(MPI_TAG)
          source=status(MPI_SOURcE)
          call MPI_SEND(i,1,mpi_int,source,tag,MPI_cOMM_WORLD,ierr)
          if (k.gt.0) then
            j=j-1
          endif
          write(*,*) source,i,j,trim(message)
          if (i.ne.0) i=i+1
          if (i.gt.num_job) i=0
          if (j.eq.0) complete=1
        endif
      enddo

      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) ! synchronize all nodes


      return
      end
  `,
  readimage: `
      SUBROUTINE readimage(filename,nx,ny,npx,npy,array)
! read a 2D fits images from a fits file.
	    IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER group,firstpix,nbuffer,npixels,i,status2
      INTEGER naxes(2)
      INTEGER nx,ny,npx,npy
      REAL array(npx,npy)
      REAL nullval,anyf
      LOGICAL anynull
      CHARACTER filename*(*)
      
!  The STATUS parameter must always be initialized.
      status=0
      nullval=0.
      
!  Get an unused Logical Unit Number to use to open the FITS file.
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)

!  Determine the size of the image.
      CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

!  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the NAXISn keywords of:'
     .,filename	 
          return
      end if

      group=1
      nx=naxes(1)
      ny=naxes(2)

      status2=0
	  
      CALL FTG2DE(unit,group,nullval,npx,nx,ny,array,anyf,status)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      return
      END 
  `,
  writeimage: `
      SUBROUTINE writeimage(filename,nx,ny,npx,npy,array)

!  Create a FITS primary array containing a 2-D image
      IMPLICIT NONE
      
      INTEGER status,unit,blocksize,bitpix,naxis
      INTEGER i,j,group
      INTEGER nx,ny,npx,npy
      INTEGER naxes(2)
      REAL array(npx,npy)
      CHARACTER filename*(*)      
      LOGICAL simple,extend
      
!  The STATUS parameter must be initialized before using FITSIO.  A
!  positive value of STATUS is returned whenever a serious error occurs.
!  FITSIO uses an "inherited status" convention, which means that if a
!  subroutine is called with a positive input value of STATUS, then the
!  subroutine will exit immediately, preserving the status value. For 
!  simplicity, this program only checks the status value at the end of 
!  the program, but it is usually better practice to check the status 
!  value more frequently.

      status=0

!  Delete the file if it already exists, so we can then recreate it.
!  The deletefile subroutine is listed at the end of this file.
      CALL deletefile(filename,status)

!  Get an unused Logical Unit Number to use to open the FITS file.
!  This routine is not required;  programmers can choose any unused
!  unit number to open the file.
      CALL ftgiou(unit,status)

!  Create the new empty FITS file.  The blocksize parameter is a
!  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      CALL ftinit(unit,filename,blocksize,status)

!  Initialize parameters about the FITS image.
!  The size of the image is given by the NAXES values. 
!  The EXTEND = TRUE parameter indicates that the FITS file
!  may contain extensions following the primary array.
      simple=.true.
      bitpix=-32
      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      extend=.false.

!  Write the required header keywords to the file
      CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!  Write the array to the FITS file.
!  The last letter of the subroutine name defines the datatype of the
!  array argument; in this case the 'J' indicates that the array has an
!  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
!  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
!  total number of pixels.  GROUP is seldom used parameter that should
!  almost always be set = 1.
      group=1
      CALL FTP2DE(unit,group,npx,nx,ny,array,status)

!  Write another optional keyword to the header
!  The keyword record will look like this in the FITS file:
!
!  EXPOSURE=                 1500 / Total Exposure Time
!
!     CALL ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

!  Check for any errors, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      return
      END
  `,
  readimage_para: ` 
! This is a simple module contains subroutine that can read and write 2D fits images
! with no extensions to the primary array. It is stored in single precision.
!IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE readimage_para(filename,nx,ny,npx,npy,array
     .,CRPIX,CD,CRVAL)
! read a 2D fits images from a fits file.
	  IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER group,firstpix,nbuffer,npixels,i,status2,j
      INTEGER naxes(2)
      INTEGER nx,ny,npx,npy
      double precision CRPIX(2),CD(2,2),CRVAL(2),temp(2)
      REAL array(npx,npy)
      REAL nullval,anyf
      LOGICAL anynull
      CHARACTER filename*(*),comment*100
      
!  The STATUS parameter must always be initialized.
      status=0
      nullval=0.
      
!  Get an unused Logical Unit Number to use to open the FITS file.
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)


!  Determine the coordinate parameters of the image.
      CALL ftgknd(unit,'CRPIX',1,2,temp,nfound,status)
!  Check that it found both CRPIX1 and CRPIX2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the CRPIXn keywords of: '
     .,filename
          return
      end if
      CRPIX(1)=temp(1)
      CRPIX(2)=temp(2)

      CALL ftgknd(unit,'CRVAL',1,2,temp,nfound,status)
!  Check that it found both CRVAL1 and CRVAL2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the CRVALn keywords of: '
     .,filename
          return
      end if
      CRVAL(1)=temp(1)
      CRVAL(2)=temp(2)

!      CALL ftgknd(unit,'CDELT',1,2,temp,nfound,status)
!  Check that it found both CDELT1 and CDELT2 keywords.
!      if (nfound .ne. 2)then
!          print *,'READIMAGE failed to read the CDELTn keywords.'
!          return
!      end if
!      CDELT(1)=temp(1)
!      CDELT(2)=temp(2)

      CALL ftgknd(unit,'CD1_',1,2,temp,nfound,status)
!  Check that it found both CD1_1 and CD1_2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the CD1_n keywords of: '
     .,filename
          return
      end if
      CD(1,1)=temp(1)
      CD(1,2)=temp(2)

      CALL ftgknd(unit,'CD2_',1,2,temp,nfound,status)
!  Check that it found both CD2_1 and CD2_2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the CD2_n keywords of: '
     .,filename
          return
      end if
      CD(2,1)=temp(1)
      CD(2,2)=temp(2)

c      CALL ftgkyd(unit,'ZP',ZP,comment,status)

c      CALL ftgkyd(unit,'EXPTIME',EXPTIME,comment,status)
      

!  Determine the size of the image.
      CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

!  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 2)then
          print *,'READIMAGE failed to read the NAXISn keywords of: '
     .,filename
          return
      end if

      group=1
      nx=naxes(1)
      ny=naxes(2)

      status2=0
	  
      CALL FTG2DE(unit,group,nullval,npx,nx,ny,array,anyf,status)     

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      return
      END 
  `,
  set_background: `
      subroutine set_background(nx1,nx2,ny1,ny2,npx,npy,image
     .,blocksize,nct,ncx,ierror)
      implicit none

c The purpose of this code is to achieve a fine adjustment of the background.

c input & output:
      integer nx1,nx2,ny1,ny2,npx,npy,ierror,blocksize,nct,ncx
      real image(npx,npy)

c parameters:
      integer npp,nfit
      parameter (npp=1000)
      parameter (nfit=5000)

c local variables:
      real pix(npp),arr(npp,3),c(nct)
      real mean_sam(nfit,3),sam(nfit),tempo(nfit),tmp_fit(nfit,3)
      integer indx(npp),nsam,nsam1,changed
      integer ix,iy,i,j,k,nbx,nby,xmin,xmax,ymin,ymax,nct_min
      real i2,j2,ij,mean,sig,aa,bb,cc,med,x,y
      real bound1,bound2
c uses:
      real ran1,func_val,ratio

      if (ierror.eq.1) return

c a rough flattening of the field first:
      call flatten_chip(nx1,nx2,ny1,ny2,npx,npy,image,4,2,ierror)

c      return

      if (ierror.eq.1) return

      ratio=1./max((nx2-nx1),(ny2-ny1))

c divide the chip into nbx*nby blocks, each with sidelength of "blocksize":
      nbx=max((nx2-nx1)/blocksize,1)
      nby=max((ny2-ny1)/blocksize,1)

c get the median of each block:
      nsam=0
      do i=1,nbx
        xmin=(i-1)*blocksize+nx1
        xmax=min(xmin+blocksize,nx2)
        do j=1,nby
          ymin=(j-1)*blocksize+ny1
          ymax=min(ymin+blocksize,ny2)

          do k=1,npp
            ix=int(ran1()*(xmax-xmin)+xmin)
            iy=int(ran1()*(ymax-ymin)+ymin)
            arr(k,1)=ix*ratio
            arr(k,2)=iy*ratio
            arr(k,3)=image(ix,iy)
            pix(k)=arr(k,3)
          enddo

          call indexx(npp,npp,pix,indx)

          k=indx(npp/2)
          nsam=nsam+1

          mean_sam(nsam,1)=arr(k,1)
          mean_sam(nsam,2)=arr(k,2)
          mean_sam(nsam,3)=arr(k,3)
        enddo
      enddo

c remove the blocks that are outliers:

      nct_min=nct*3/2

      nsam1=nsam
      changed=1
      do while (changed.eq.1 .and. nsam1.ge.nct_min)
        call find_slope_2D(nfit,nsam1,mean_sam,aa,bb,cc)

        do i=1,nsam1
          sam(i)=mean_sam(i,3)-aa-bb*mean_sam(i,1)-cc*mean_sam(i,2)
          tempo(i)=sam(i)
        enddo
        call sort(nsam1,nfit,tempo)
        mean=tempo(nsam1/2)
        sig=0.5*(tempo(nsam1*5/6)-tempo(nsam1/6))
        nsam=0
        changed=0
        do i=1,nsam1
          if (abs(sam(i)-mean).lt.3.*sig) then
            nsam=nsam+1
            mean_sam(nsam,1)=mean_sam(i,1)
            mean_sam(nsam,2)=mean_sam(i,2)
            mean_sam(nsam,3)=mean_sam(i,3)
          else
            changed=1
          endif
        enddo
        nsam1=nsam
      enddo

      if (nsam1.lt.nct_min) then
        write(*,*) 'Background not stable enough!',nsam1
        ierror=1
        return
      endif

      call fit_2D(nfit,nsam1,mean_sam,nct,ncx,c)

      do i=nx1,nx2
        do j=ny1,ny2

          x=i*ratio
          y=j*ratio

          image(i,j)=image(i,j)-func_val(x,y,nct,ncx,c)

        enddo
      enddo


      return
      end
  `,
  set_sig: ` 
      subroutine set_sig(nx1,nx2,ny1,ny2,npx,npy,image
     .,aa,bb,cc,ierror)
      implicit none

c The purpose of this code is to get the rms of noise as a function of position.

c input & output:
      integer nx1,nx2,ny1,ny2,npx,npy,ierror
      real image(npx,npy)

c parameters:
      integer npp
      parameter (npp=2000)

c local variables:
      real pix(npp),arr(npp,3),arr2(npp,3),c(6)
      integer ix,iy,i,j,nr
      real aa,bb,cc,arr_min,arr_max,temp
c uses:
      real ran1

      if (ierror.eq.1) return

c get a sample of points:

      do i=1,npp
        ix=int(ran1()*(nx2-nx1-1)+nx1)
        iy=int(ran1()*(ny2-ny1-1)+ny1)
        pix(i)=0.5*((image(ix,iy)-image(ix+1,iy))**2
     .+(image(ix,iy)-image(ix,iy+1))**2)
        arr(i,1)=ix
        arr(i,2)=iy
        arr(i,3)=pix(i)
      enddo
      call sort(npp,npp,pix)

      arr_min=pix(npp/3)
      arr_max=pix(2*npp/3)

c remove the points of extreme values:
      nr=0
      do i=1,npp
        if (arr(i,3).ge.arr_min.and.arr(i,3).le.arr_max) then
          nr=nr+1
          arr2(nr,1)=arr(i,1)
          arr2(nr,2)=arr(i,2)
          arr2(nr,3)=arr(i,3)
        endif
      enddo

c fit a 2D linear function with the remaining points:

      if (nr.lt.npp/10) then
        ierror=1
        return
      endif

      call find_slope_2D(npp,nr,arr2,aa,bb,cc)

      do i=nx1,nx2
        do j=ny1,ny2
          temp=aa+bb*i+cc*j
          image(i,j)=image(i,j)/sqrt(0.5*temp)
        enddo
      enddo


      return
      end
  `,
  generate_gaia_file_name: ` 
      subroutine generate_gaia_file_name(cRVAL,filename)
      implicit none

      double precision cRVAL(2)
      character filename*(*),c_ra*2,c_dec
      integer ra,dec

10    format(I1.1)
20    format(I2.2)

      dec=min(int(abs(cRVAL(2))/10d0)+1,9)
      write(c_dec,10) dec

      ra=min(int(cRVAL(1)/10d0),35)
      write(c_ra,20) ra

      if (cRVAL(2).ge.0) then
        filename=trim(filename)//'/gaia_p'
      else
        filename=trim(filename)//'/gaia_m'
      endif

      filename=trim(filename)//c_dec
      if (dec.eq.9) then
        filename=trim(filename)//'.cat'
      else
        filename=trim(filename)//'_'//c_ra//'.cat'
      endif

      return
      end
  `,
  gen_astrometry_data: `
      subroutine gen_astrometry_data(cat_standard,nx,ny,npx,npy
     .,map,weight,cRPIX,cD,cRVAL,filename,proc_error)
      implicit none

      real astrometry_shift_ratio
      parameter (astrometry_shift_ratio=0.2)

      character cat_standard*(*),filename*(*)
      integer nx,ny,npx,npy,ierror,n_user,n_ref,proc_error
      real map(npx,npy),flux_min
      integer weight(npx,npy)
      integer i,j,u,v,k
      double precision ra,dra,diffra,dec(2),a,d,flux,mag
      integer nss_max,nss,n_user_max
      parameter (nss_max=10000)
      parameter (n_user_max=200)

      double precision xs(nss_max),ys(nss_max)
      double precision xr(nss_max),yr(nss_max)
      double precision ra_r(nss_max),dec_r(nss_max)
      integer box(nss_max),astrometry_shift_range
      double precision ra2(nss_max),dec2(nss_max)
      double precision x2(nss_max),y2(nss_max)

      double precision cRPIX(2),cD(2,2),cRVAL(2)


      if (proc_error.eq.1) then
        nss=0
        n_user=0
        n_ref=0
        goto 40
      endif


c Get the ranges of ra & dec:
      call get_ra_dec_range(nx,ny,ra,dec,dra,cRPIX,cD,cRVAL
     .,astrometry_shift_ratio)

c------------------------------------------------------------------
c Get the positions of the stars in the reference catalog:

      n_ref=0
      open(unit=20,file=trim(cat_standard)
     .,status='old',iostat=ierror)
      rewind 20
      if (ierror.ne.0) then
        write(*,'(A)') cat_standard
        write(*,*) ierror
        pause 'catalog file error!!'
      endif
      read(20,*)
      do while (ierror.ge.0)
        read(20,*,iostat=ierror) a,d
        if (ierror.lt.0) cycle
        if (abs(diffra(a,ra)).gt.dra*0.5) cycle
        if (d.lt.dec(1).or.d.gt.dec(2)) cycle
        n_ref=n_ref+1
        if (n_ref.gt.nss_max) then
          write(*,*) 'n_ref is too large!!',filename
          close(20)
          nss=0
          n_user=0
          n_ref=0
          goto 40
        endif
        ra_r(n_ref)=a
        dec_r(n_ref)=d
        call coordinate_transfer_simple(a,d,xr(n_ref),yr(n_ref),-1
     .,cRPIX,cD,cRVAL)
      enddo
      close(20)


      call get_astrometry_catalog(nx,ny,npx,npy,map
     .,weight,n_user_max,nss_max,nss,xs,ys)


      n_user=nss

c Match the point sources in the "user" and "ref" catalogs:
      astrometry_shift_range=int(max(nx,ny)*astrometry_shift_ratio)

      call pattern_matching(nss_max,n_ref,xr,yr,nss_max,n_user,xs,ys
     .,astrometry_shift_range,box)

c Get the Astrometric calibration parameters:


      nss=0
      do i=1,n_ref
        if (box(i).eq.0) cycle
        nss=nss+1
        ra2(nss)=ra_r(i)
        dec2(nss)=dec_r(i)
        j=box(i)
        x2(nss)=xs(j)
        y2(nss)=ys(j)
      enddo

      write(*,*) nss,n_ref,n_user,trim(filename)

40    open(unit=10,file=trim(filename),status='replace')
      rewind 10
      write(10,*) cRPIX(1),cRPIX(2),cRVAL(1),cRVAL(2)
      write(10,*) cD(1,1),cD(1,2),cD(2,1),cD(2,2)
      write(10,*) nss,n_user,n_ref
      do i=1,nss
        write(10,*) ra2(i),dec2(i),x2(i),y2(i)
      enddo
      close(10)


      return
      end
  `,
  locate_defects: `
      subroutine locate_defects(nx,ny,npx,npy,array,normap,weight
     .,area_max,area_thresh,ierror)
      implicit none

      integer nx,ny,npx,npy,ierror,area_max,area_thresh
      real normap(npx,npy),map(nx,ny),array(npx,npy)
      integer weight(npx,npy)
      real diffx(nx,ny),diffy(nx,ny)
      integer i,j,ix,iy
      character filename*100

      integer margin
      parameter (margin=10)
      real defect_halo_thresh
      parameter (defect_halo_thresh=1.)
      integer y_smooth,x_smooth
      parameter (y_smooth=200)
      parameter (x_smooth=100)

      real sig,med,sigx,medx,sigy,medy

      real loga,iden
      external loga,iden

      do j=1,ny
        do i=nx/2-margin,nx/2+margin
          weight(i,j)=0
        enddo
        do i=1,margin
          weight(i,j)=0
          weight(nx+1-i,j)=0
        enddo
      enddo
      do i=1,nx
        do j=1,margin
          weight(i,j)=0
          weight(i,ny+1-j)=0
        enddo
      enddo

      if (ierror.eq.1) return

      do i=1,nx
        do j=1,ny
c          map(i,j)=loga(normap(i,j),1)
c          map(i,j)=log(normap(i,j)
          map(i,j)=loga(array(i,j),1)
        enddo
      enddo


      call remove_continuous(nx,ny,nx,ny,map,iden,4)
c      filename='map1.fits'
c      call writeimage(filename,nx,ny,nx,ny,map)
c      call get_sig_med(nx,ny,map,sig,med)
c      do i=1,nx
c        do j=1,ny
c          if (abs(map(i,j)).gt.5.*sig) weight(i,j)=0
c        enddo
c      enddo
      do i=1,nx
        ix=mod(i,nx)+1
        do j=1,ny
          iy=mod(j,ny)+1
          diffx(i,j)=map(i,j)-map(ix,j)
          diffy(i,j)=map(i,j)-map(i,iy)
        enddo
      enddo
      call get_sig_med(nx,ny,diffx,sigx,medx)
      do i=1,nx
        do j=1,ny
          if (abs(diffx(i,j)).gt.8.*sigx) weight(i,j)=0
        enddo
      enddo
      call get_sig_med(nx,ny,diffy,sigy,medy)
      do i=1,nx
        do j=1,ny
          if (abs(diffy(i,j)).gt.8.*sigy) weight(i,j)=0
        enddo
      enddo

      call mask_source_regions(nx,ny,npx,npy
     .,weight,normap,area_max,defect_halo_thresh*2.,area_thresh)

      call detect_stripes(nx,ny,npx,npy,normap,weight
     .,x_smooth,y_smooth)

      call detect_artificial_stripes(nx,ny,npx,npy,weight
     .,diffx,diffy,sigx,sigy,medx,medy)

      do i=1,nx
        do j=1,ny
          if (weight(i,j).gt.1) weight(i,j)=1
        enddo
      enddo

      call detect_stellar_halo(nx,ny,npx,npy,normap,weight
     .,area_max,defect_halo_thresh)

      call detect_dent(nx,ny,npx,npy,normap,weight
     .,area_max,defect_halo_thresh)



      return
      end
`,
  merge_defects: ` 
      subroutine merge_defects(nx,ny,npx,npy
     .,weight,normap,area_max,source_thresh,area_thresh
     .,ierror)
      implicit none

      integer nx,ny,npx,npy,area_max,area_thresh,ierror
      real normap(npx,npy),source_thresh
      integer weight(npx,npy),mark(nx,ny)

      integer nb,nbb,toobig,buffer(area_max,2)
      integer i,j,ix,iy,jx,jy,u,v,k1,k2,k


      if (ierror.eq.1) return

      do i=1,nx
        do j=1,ny
          if (normap(i,j).ge.source_thresh.and.weight(i,j).eq.1) then
            mark(i,j)=1
          else
            mark(i,j)=0
          endif
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          if (mark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            mark(i,j)=-1
            toobig=0
            do while (nb.gt.nbb)
              k1=nbb+1
              k2=nb
              nbb=nb
              do k=k1,k2
                ix=buffer(k,1)
                iy=buffer(k,2)
                do u=max(ix-1,1),min(ix+1,nx)
                  do v=max(iy-1,1),min(iy+1,ny)
                    if (mark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      mark(u,v)=-1
                      if (nb.eq.area_max) then
                        toobig=1
                        goto 20
                      endif
                    elseif (mark(u,v).gt.1
     . .or. (mark(u,v).eq.0 .and. weight(u,v).eq.0)) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                mark(buffer(k,1),buffer(k,2))=area_max
                weight(buffer(k,1),buffer(k,2))=0
              enddo
            else
              do k=1,nb
                mark(buffer(k,1),buffer(k,2))=nb
              enddo
            endif
          endif
        enddo
      enddo


      return
      end
  `,
  writeimage_copyhdu: ` 
      SUBROUTINE writeimage_copyhdu(file1,filename,nx,ny,npx,npy,array)

!  Create a FITS primary array containing a 2-D image
      IMPLICIT NONE
      
      INTEGER status,unit,blocksize,bitpix,naxis,readwrite
      INTEGER i,j,group,unit1
      INTEGER nx,ny,npx,npy
      INTEGER naxes(2)
      REAL array(npx,npy)
      CHARACTER filename*(*),file1*(*)      
      LOGICAL simple,extend

      status=0      

!  Delete the file if it already exists, so we can then recreate it.
!  The deletefile subroutine is listed at the end of this file.
      CALL deletefile(filename,status)

!  Get an unused Logical Unit Number to use to open the FITS file.
!  This routine is not required;  programmers can choose any unused
!  unit number to open the file.
      CALL ftgiou(unit,status)

!  Create the new empty FITS file.  The blocksize parameter is a
!  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      CALL ftinit(unit,filename,blocksize,status)

      status=0
      CALL ftgiou(unit1,status)
      readwrite=0
      CALL ftopen(unit1,file1,readwrite,blocksize,status)
      CALL ftcphd(unit1,unit,status)

      CALL ftmkyj(unit,'bitpix',-32,'new',status)	
c      CALL ftmkyd(unit,'bscale',1d0,'new',status)
c      CALL ftmkyd(unit,'bzero',0d0,'new',status)	

!  Initialize parameters about the FITS image.
!  The size of the image is given by the NAXES values. 
!  The EXTEND = TRUE parameter indicates that the FITS file
!  may contain extensions following the primary array.

c	simple=.true.
c      bitpix=-32
c      naxis=2
c      naxes(1)=nx
c      naxes(2)=ny
c      extend=.false.

!  Write the required header keywords to the file
c      CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!  Write the array to the FITS file.
!  The last letter of the subroutine name defines the datatype of the
!  array argument; in this case the 'J' indicates that the array has an
!  integer*4 datatype. ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
!  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
!  total number of pixels.  GROUP is seldom used parameter that should
!  almost always be set = 1.
      group=1
      CALL FTP2DE(unit,group,npx,nx,ny,array,status)

!  Write another optional keyword to the header
!  The keyword record will look like this in the FITS file:
!
!  EXPOSURE=                 1500 / Total Exposure Time
!
!     CALL ftpkyj(unit,'EXPOSURE',1500,'Total Exposure Time',status)


!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)
      CALL ftclos(unit1, status)
      CALL ftfiou(unit1, status)

!  Check for any errors, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      IF (status .gt. 0) CALL printerror(status)
      return
      END
  `,
// ==============================================================================
  get_astrometry: `
      subroutine get_astrometry(IMAGE_FILE,nchip,DIR_OUTPUT)
      implicit none
      include 'para.inc'

      integer nchip
      character*(strl) IMAGE_FILE(NMAX_cHIP),DIR_OUTPUT
      integer ichip

      character*(strl) PREFIX,filename
      integer nss(nchip),n_user,n_ref,i,k
      integer nss_max
      parameter (nss_max=10000)
      double precision ra2(nchip,nss_max),dec2(nchip,nss_max)
      double precision x2(nchip,nss_max),y2(nchip,nss_max)

      double precision cRPIX2(nchip,2),cD2(nchip,2,2)
      double precision PU(2,npd),cRVAL2(2)
      integer valid(nchip),tot_valid,tot_source

      tot_valid=0
      tot_source=0

      do ichip=1,nchip

        call get_PREFIX(IMAGE_FILE(ichip),PREFIX)
        filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'_astro.dat'

        open(unit=10,file=filename,status='old')
        rewind 10
        read(10,*) cRPIX2(ichip,1),cRPIX2(ichip,2)
     .,cRVAL2(1),cRVAL2(2)
        read(10,*) cD2(ichip,1,1),cD2(ichip,1,2)
     .,cD2(ichip,2,1),cD2(ichip,2,2)
        read(10,*) nss(ichip),n_user,n_ref
        do i=1,nss(ichip)
          read(10,*) ra2(ichip,i),dec2(ichip,i)
     .,x2(ichip,i),y2(ichip,i)
        enddo
        close(10)
        if (nss(ichip).ge.10) then
          valid(ichip)=1
          tot_valid=tot_valid+1
          tot_source=tot_source+nss(ichip)
        else
          valid(ichip)=0
        endif
      enddo


      if (tot_source.ge.(npd+tot_valid*3)*3) then
        call measure_astrometry_global(nss_max,nss,nchip,ra2
     .,dec2,x2,y2,cRPIX2,cD2,cRVAL2,PU,npd,valid)
        call check_astrometry_global(nss_max,nss,nchip,ra2
     .,dec2,x2,y2,cRPIX2,cD2,cRVAL2,PU,npd,valid)
      else
        do i=1,nchip
          valid(i)=0
        enddo
      endif


      call get_PREFIX_expo(IMAGE_FILE(1),PREFIX)
      filename=trim(DIR_OUTPUT)//'/astrometry/'
     .//trim(PREFIX)//'.head'
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) cRVAL2(1),cRVAL2(2)
      do i=1,npd
        write(10,*) PU(1,i),PU(2,i)
      enddo
      do k=1,nchip
        write(10,*) k,valid(k),cRPIX2(k,1),cRPIX2(k,2)
     .,cD2(k,1,1),cD2(k,1,2),cD2(k,2,1),cD2(k,2,2)
      enddo
      close(10)


      return
      end
`,
  get_PREFIX_expo: `
      subroutine get_PREFIX_expo(imagefile,PREFIX)
      implicit none

      character*(*) imagefile,PREFIX
      integer i,p_dot,p_slash,n

      p_dot=0
      p_slash=0

      n=len(trim(imagefile))
      do i=n,1,-1
        if (imagefile(i:i).eq.'_' .and. p_dot.eq.0) p_dot=i
        if (imagefile(i:i).eq.'/' .and. p_slash.eq.0) p_slash=i
      enddo
      if (p_dot.eq.0 .or. p_slash.eq.0. .or. p_dot.le.p_slash+1) 
     .stop 'Image_file name is NOT normal !'

      PREFIX=imagefile(p_slash+1:p_dot-1)

      return
      end
  `,
  get_PREFIX: `
      subroutine get_PREFIX(imagefile,PREFIX)
      implicit none

      character*(*) imagefile,PREFIX
      integer i,p_dot,p_slash,n

      p_dot=0
      p_slash=0

      n=len(trim(imagefile))
      do i=n,1,-1
        if (imagefile(i:i).eq.'.' .and. p_dot.eq.0) p_dot=i
        if (imagefile(i:i).eq.'/' .and. p_slash.eq.0) p_slash=i
      enddo
      if (p_dot.eq.0 .or. p_slash.eq.0. .or. p_dot.le.p_slash+1) then
        write(*,*) 'Image_file name is NOT normal !'
        read(*,*)
      endif

      PREFIX=imagefile(p_slash+1:p_dot-1)

      return
      end
`,

  read_astrometry_para: ` 
      subroutine read_astrometry_para(filename,ichip
     .,cRPIX,cD,cRVAL,PU,npd,proc_error)
      implicit none

      integer ichip,npd,proc_error
      character filename*(*)
      double precision cRPIX(2),cD(2,2),PU(2,npd),cRVAL(2)
      integer valid,j,k,i

      if (proc_error.eq.1) return

      open(unit=11,file=filename,status='old')
      rewind 11
      read(11,*) cRVAL(1),cRVAL(2)
      do i=1,npd
        read(11,*) PU(1,i),PU(2,i)
      enddo
      do k=1,ichip-1
        read(11,*)
      enddo
      read(11,*) j,valid,cRPIX(1),cRPIX(2)
     .,cD(1,1),cD(1,2),cD(2,1),cD(2,2)

      close(11)

      if (valid.eq.0) proc_error=1

      return
      end
  `,
  update_para: `
      SUBROUTINE update_para(filename,CRPIX,CD)
      IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize
      double precision CRPIX(2),CD(2,2)
      CHARACTER filename*(*),comment*100
      
      status=0
      comment='replaced'
      CALL ftgiou(unit,status)
      readwrite=1
      CALL ftopen(unit,filename,readwrite,blocksize,status)
  
      call ftukyd(unit,'CRPIX1',CRPIX(1),9,comment,status)
      call ftukyd(unit,'CRPIX2',CRPIX(2),9,comment,status)

      call ftukyd(unit,'CD1_1',CD(1,1),9,comment,status)
      call ftukyd(unit,'CD1_2',CD(1,2),9,comment,status)
      call ftukyd(unit,'CD2_1',CD(2,1),9,comment,status)
      call ftukyd(unit,'CD2_2',CD(2,2),9,comment,status)
	

      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

      IF (status .gt. 0) CALL printerror(status)
      return
      END 
`,
  coordinate_transfer_PU: `
      subroutine coordinate_transfer_PU(a,d,x,y,direc
     .,cRPIX,cD,cRVAL,PU,npd)
      implicit none

      integer direc,npd
      double precision a,d,x,y,xx,yy,xxx,yyy,rr,dxx,dyy
      double precision cRPIX(2),cD(2,2),cD_1(2,2),cRVAL(2)
      double precision PU(2,npd)

      double precision da,dd,const1,const2,ds,xi,eta,diffra
      double precision cosda,tandc,tandd,tanda,temp,sumra
      double precision pi
      parameter (pi=3.1415926d0)

c      integer tmp_sig
c      common /temp_pass/ tmp_sig

      const1=pi/180d0
      tandc=tan(cRVAL(2)*const1)

      if (direc.eq.1) then                            ! direc 是什么？
        xx=cD(1,1)*(x-cRPIX(1))+cD(1,2)*(y-cRPIX(2))  ! 转换到天球坐标
        yy=cD(2,1)*(x-cRPIX(1))+cD(2,2)*(y-cRPIX(2))

        call mapping_PU(xx,yy,xi,eta,npd,PU,1)

        xxx=xi*const1
        yyy=eta*const1

        da=xxx/(cos(cRVAL(2)*const1)*(1.-yyy*tandc))
        da=da-da*da*da*0.33333333333
        a=sumra(da/const1,cRVAL(1))
        cosda=1.-da*da*0.5+da**4/24.

        dd=(yyy*(cosda+tandc**2)+tandc*(cosda-1.))
     ./(yyy*tandc*(cosda-1.)+1.+cosda*tandc**2)
        dd=dd-dd*dd*dd*0.3333333333
        d=dd/const1+cRVAL(2)

      else

        da=diffra(a,cRVAL(1))*const1
        dd=(d-cRVAL(2))*const1
        tandd=tan(dd)
        cosda=cos(da)
        tanda=tan(da)

        yy=tandc*(cosda-1.)-(1.+cosda*tandc**2)*tandd
        yy=yy/(tandd*tandc*(cosda-1.)-(cosda+tandc**2))
        xx=tanda*(cos(cRVAL(2)*const1)*(1.-yy*tandc))

        xi=xx/const1
        eta=yy/const1

c        if (tmp_sig.eq.1.and.direc.eq.-1) then
c          write(*,*) da/const1,dd/const1,xi,eta
c        endif

        call mapping_PU(xx,yy,xi,eta,npd,PU,2)

        temp=cD(1,1)*cD(2,2)-cD(1,2)*cD(2,1)
        temp=1./temp

        cD_1(1,1)=cD(2,2)*temp
        cD_1(2,2)=cD(1,1)*temp
        cD_1(1,2)=-cD(1,2)*temp
        cD_1(2,1)=-cD(2,1)*temp

        x=xx*cD_1(1,1)+yy*cD_1(1,2)+cRPIX(1)
        y=xx*cD_1(2,1)+yy*cD_1(2,2)+cRPIX(2)

      endif

      return
      end
`,
 get_expo_catalog: `
      subroutine get_expo_catalog(PREFIX,nx,ny,sigmap
     .,weight,normap,ierror)
      implicit none
      include 'para.inc'

      character*(strl) catname,PREFIX
      integer nx,ny,ierror
      real sigmap(npx,npy),normap(npx,npy)
      integer weight(npx,npy),mark(npx,npy)

      integer nb,nbb,toobig,buffer(area_max,2)
      integer i,j,ix,iy,jx,jy,u,v,k1,k2,k
      real xp,yp,r2,rmax,sig,temp,peak,thresh,xc,yc,SNR
      real total_flux,half_light_flux
      integer half_light_area,total_area


      if (ierror.eq.1) return

      do i=1,nx
        do j=1,ny
          if (normap(i,j).ge.source_thresh.and.weight(i,j).ge.1) then
            mark(i,j)=1
          else
            mark(i,j)=0
          endif
        enddo
      enddo

      catname=trim(PREFIX)//'.cat'

      open(unit=10,file=catname,status='replace')
      rewind 10

      write(10,*) ' xp ',' yp '
     .,' total_area ',' half_light_area ',' sig '
     .,' sig_normed_total_flux ',' sig_normed_half_light_flux '
     .,' sig_normed_peak ',' rmax '

      do i=1,nx
        do j=1,ny
          if (mark(i,j).eq.1) then
            nbb=0
            nb=1
            buffer(nb,1)=i
            buffer(nb,2)=j
            mark(i,j)=-1
            toobig=0
            do while (nb.gt.nbb)
              k1=nbb+1
              k2=nb
              nbb=nb
              do k=k1,k2
                ix=buffer(k,1)
                iy=buffer(k,2)
                do u=max(ix-1,1),min(ix+1,nx)
                  do v=max(iy-1,1),min(iy+1,ny)
                    if (mark(u,v).eq.1) then
                      nb=nb+1
                      buffer(nb,1)=u
                      buffer(nb,2)=v
                      mark(u,v)=-1
                      if (nb.eq.area_max) then
                        toobig=1
                        goto 20
                      endif
                    elseif (mark(u,v).gt.1) then
                      toobig=1
                      goto 20
                    endif
                  enddo
                enddo
              enddo
            enddo
20          if (toobig.eq.1) then
              do k=1,nb
                mark(buffer(k,1),buffer(k,2))=area_max
                weight(buffer(k,1),buffer(k,2))=2
              enddo
            else
              if (nb.ge.area_thresh) then
                xc=0.
                yc=0.
                total_flux=0.
                sig=0.
                peak=-100000.
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  mark(ix,iy)=nb
                  weight(ix,iy)=2

                  xc=xc+ix*normap(ix,iy)
                  yc=yc+iy*normap(ix,iy)
                  total_flux=total_flux+normap(ix,iy)
                  sig=sig+sigmap(ix,iy)

                  if (normap(ix,iy).gt.peak) then
                    peak=normap(ix,iy)
                    xp=ix
                    yp=iy
                  endif
                enddo
                xc=xc/total_flux
                yc=yc/total_flux
                total_area=nb
                sig=sig/total_area

                thresh=peak*0.5
                rmax=0
                half_light_area=0
                half_light_flux=0.
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  r2=(ix-xc)**2+(iy-yc)**2
                  rmax=max(rmax,r2)
                  if (normap(ix,iy).ge.thresh) then
                    half_light_area=half_light_area+1
                    half_light_flux=half_light_flux+normap(ix,iy)
                  endif
                enddo
                rmax=sqrt(rmax)
                if (peak.ge.core_thresh) then
                  write(10,*) xp,yp,total_area,half_light_area
     .,sig,total_flux,half_light_flux,peak,rmax
                endif
              else
                do k=1,nb
                  ix=buffer(k,1)
                  iy=buffer(k,2)
                  mark(ix,iy)=nb
                enddo
              endif
            endif
          endif
        enddo
      enddo


      close(10)

      return
      end
`,

generate_gal_cat_file_name: `
      subroutine generate_gal_cat_file_name(cRVAL,filename
     .,sortfile,sortnum)
      implicit none
      include 'para.inc'

      double precision cRVAL(2) , m_ra , m_dec
      character filename*(*)
      character*(strl) sortfile(27)
      integer dec1,dec2
      character c_dec*2
      character c_ra*3
      integer ra1,ra2
      integer ra,dec,sortnum,mra,mdec

10    format(I2.2)
20    format(I3.3)

      m_dec=1.0d0
      ! first judge the dec
      if (abs(cRVAL(2)).lt.30.) then
        m_ra=1.4d0
      elseif (abs(cRVAL(2)).lt.40.) then
        m_ra=1.6d0
      elseif (abs(cRVAL(2)).lt.50.) then
        m_ra=1.8d0
      elseif (abs(cRVAL(2)).lt.60.) then
        m_ra=2.5d0
      else
        m_ra=3.0d0
      endif
      ! --------------dec-------------------------
      dec1=int(floor(cRVAL(2)-m_dec))
      dec2=int(floor(cRVAL(2)+m_dec))
      !--------------- ra -----------------------
      ra1=int(floor(cRVAL(1)-m_ra))
      ra2=int(floor(cRVAL(1)+m_ra))
      !------------------------------------------
      ! write in the file name
      sortnum=0
      do dec=dec1,dec2
        do ra=ra1,ra2
          sortnum=sortnum+1
          mra=ra
          mdec=dec
          if (ra.lt.0) then
            mra=mra+360
          elseif (ra.ge.360) then
            mra=mra-360
          endif 
          write(c_ra,20) mra

          if (dec.ge.0) then
            sortfile(sortnum)=trim(filename)//'/gal_p'
            write(c_dec,10) mdec
          elseif (dec.lt.0) then
            sortfile(sortnum)=trim(filename)//'/gal_m'
            mdec = -mdec - 1
            write(c_dec,10) mdec
          endif
          sortfile(sortnum)=trim(sortfile(sortnum))//c_dec
     .//'_'//c_ra//'.cat'
        enddo
      enddo
      return
      end
`,

  de_blending: `
      subroutine de_blending(sortfile,sortnum,nx,ny,weight
     .,cRPIX,cD,cRVAL,PU,proc_error)
      IMPLICIT NONE
      include 'para.inc'

      character*(strl) sortfile(27)
      integer n,sortnum
      double precision cRPIX(2),cD(2,2),cRVAL(2)
      double precision PU(2,npd),z1
      double precision xx,yy,ra,dec,z,dra,ra_c,dec_bound(2),diffra
      real astrometry_shift_ratio
      parameter (astrometry_shift_ratio=0.2)

      integer nx,ny,proc_error,ierr
      integer weight(npx,npy)
      integer i,j,ix,iy,iz,old_w,new_w
      double precision mag_g,mag_r,mag_i,mag_z,mag_y
      double precision magerr_g,magerr_r,magerr_i,magerr_z,magerr_y
      double precision zp,zperr

      if (proc_error.eq.1) return

      call get_ra_dec_range_fine(nx,ny,ra_c,dec_bound,dra  
     .,cRPIX,cD,cRVAL,PU,npd,astrometry_shift_ratio)

      do n=1,sortnum
        open(unit=10,file=trim(sortfile(n)),status='old',iostat=ierr)
        rewind 10
        if (ierr.ne.0) then
          write(*,*) trim(sortfile(n))
     .,'  catalog file error in deblending!!'
          cycle
        endif
        read(10,*)

        do while (ierr.ge.0)
          read(10,*,iostat=ierr) ra,dec,mag_g,magerr_g,mag_r,magerr_r,
     .mag_i,magerr_i,mag_z,magerr_z,mag_y,magerr_y,zp,zperr
          if (ierr.lt.0) cycle
          if (zp.lt.0 .or. zp.gt.5.) cycle

          z=zp
          if (abs(diffra(ra,ra_c)).gt.dra*0.5) cycle
          if (dec.lt.dec_bound(1).or.dec.gt.dec_bound(2)) cycle
          call coordinate_transfer_PU(ra,dec,xx,yy,-1
     .,cRPIX,cD,cRVAL,PU,npd)
          ix=int(xx+0.5)
          iy=int(yy+0.5)
          if (ix.lt.1 .or. ix.gt.nx .or. iy.lt.1 .or. iy.gt.ny) cycle
          iz=int(z*1000.+10)

          if (weight(ix,iy).lt.2) then
            cycle
          elseif (weight(ix,iy).gt.2) then
            z1=0.001*(weight(ix,iy)-10.)
            if (abs(z-z1).gt.dz_thresh) then
              new_w=0
              old_w=weight(ix,iy)
              call fill_patch(nx,ny,weight,ix,iy,old_w,new_w)
            endif
          else
            old_w=2
            new_w=iz
            call fill_patch(nx,ny,weight,ix,iy,old_w,new_w)
          endif
        enddo
        close(10)
      enddo

      do i=1,nx
        do j=1,ny
          if (weight(i,j).gt.2) weight(i,j)=2
        enddo
      enddo

      return
      end
`,

  gen_source_ext_catalog: `
      subroutine gen_source_ext_catalog(sortfile,sortnum,PREFIX
     .,nx,ny,array,weight,sigmap,cRPIX,cD,cRVAL,PU,ngal,proc_error)
      implicit none
      include 'para.inc'

c      The purpose of this subroutine to extract sources from an existing external catalog.

      character*(strl) sortfile(27),filename,PREFIX
      integer n,sortnum
      integer ierror,flag,ngal,proc_error

      character*500 cat_content,cat_list

      integer nx,ny
      real array(npx,npy),sigmap(npx,npy)
      integer weight(npx,npy)

      integer i,j,u,v,nn1,nn2,ig,igal

      integer sid(ngal_max)

      real source_para(ngal_max,npara)
      real source_collect(ngal_max,ns,ns)
      real noise_collect(ngal_max,ns,ns)
      real source(ns,ns),noise(ns,ns)

      real xp,yp,sig,total_flux,half_light_flux,peak,rf
      integer total_area,half_light_area
      common /stamp_pass/ xp,yp,total_area,half_light_area,sig
     .,total_flux,half_light_flux,peak,rf

      real SNR,temp

      double precision cRPIX(2),cD(2,2),cRVAL(2)
      double precision PU(2,npd)
      double precision xx,yy,ra,dec,dra,ra_c,dec_bound(2),diffra
      real astrometry_shift_ratio
      parameter (astrometry_shift_ratio=0.2)

      integer imax,jmax
      common /defect_pass/ imax,jmax

      integer orighead

      ngal=0
      if (proc_error.eq.1) goto 40

      call get_ra_dec_range_fine(nx,ny,ra_c,dec_bound,dra
     .,cRPIX,cD,cRVAL,PU,npd,astrometry_shift_ratio)

      ig=0
      do n=1,sortnum
        open(unit=10,file=trim(sortfile(n)),status='old',iostat=ierror)
        rewind 10
        if (ierror.ne.0) then
          write(*,*) trim(sortfile(n))
     .,'  catalog file error in source_ext!!'
          cycle
        endif
        read(10,*)

        do while (ierror.ge.0)
          read(10,*,iostat=ierror) ra,dec
          if (ierror.lt.0) cycle
          ig=ig+1

          if (abs(diffra(ra,ra_c)).gt.dra*0.5) cycle
          if (dec.lt.dec_bound(1).or.dec.gt.dec_bound(2)) cycle

          call coordinate_transfer_PU(ra,dec,xx,yy,-1
     .,cRPIX,cD,cRVAL,PU,npd)

          xp=xx
          yp=yy

          if (xp-nl_2.lt.chip_margin.or.xp+nl_2.gt.nx-chip_margin
     . .or. yp-nl_2.lt.chip_margin.or.yp+nl_2.gt.ny-chip_margin) cycle

          sig=sigmap(int(xp+0.5),int(yp+0.5))

          call find_noise(flag,noise,nx,ny,array,weight)
          if (flag.lt.0) cycle
          call check_source(flag,source,nx,ny,array,weight)
          if (flag.lt.0) cycle

          ngal=ngal+1
          do u=1,ns
            do v=1,ns
              source_collect(ngal,u,v)=source(u,v)
              noise_collect(ngal,u,v)=noise(u,v)
            enddo
          enddo
          source_para(ngal,1)=ig
          source_para(ngal,2)=xp
          source_para(ngal,3)=yp
          source_para(ngal,4)=sig
          source_para(ngal,5)=peak
          source_para(ngal,6)=imax
          source_para(ngal,7)=jmax
          source_para(ngal,8)=half_light_flux
          source_para(ngal,9)=half_light_area
          source_para(ngal,10)=flag

          sid(ngal)=ig

          if (ngal.ge.ngal_max) then
            close(10)
            goto 40
          endif
        enddo
        close(10)
      enddo

40    if (ngal.gt.0) then
        nn1=ns*len_g
        nn2=ns*(int(ngal/len_g)+1)
        filename=trim(PREFIX)//'_source.fits'
        call write_stamps(ngal_max,1,ngal,ns,ns
     .,source_collect,nn1,nn2,filename)

        filename=trim(PREFIX)//'_noise.fits'
        call write_stamps(ngal_max,1,ngal,ns,ns
     .,noise_collect,nn1,nn2,filename)
      endif

      filename=trim(PREFIX)//'_source_info.dat'
      open(unit=10,file=filename,status='replace')
      rewind 10
      write(10,*) 'ig xp yp sigma peak imax '
     .,'jmax half_light_flux half_light_area flag'
      do i=1,ngal
        write(10,*) (source_para(i,j),j=1,iflag)
      enddo
      close(10)

      orighead = 0
      filename=trim(PREFIX)//'_orig.cat'
      open(unit=15,file=filename,status='replace')
      rewind 15
      if (proc_error.eq.1 .or. ngal.eq.0) then
        write(15,*) 'No sources!!'
      else
        i=0
        igal=1
        ig=sid(igal)
        do n=1,sortnum
        open(unit=10,file=trim(sortfile(n)),status='old',iostat=ierror)
          rewind 10
          if (ierror.ne.0) then
            close(10)
            cycle
          endif
          read(10,'(A)',iostat=ierror) cat_list
          if (orighead.eq.0) then
            write(15,'(A)') trim(cat_list)
            orighead = 1
          endif
          do while (ierror.ge.0 .and. igal.le.ngal)
            read(10,'(A)',iostat=ierror) cat_content
            if (ierror.lt.0) cycle
            i=i+1
            if (i.eq.ig) then
              write(15,'(A)') trim(cat_content)
              igal=igal+1
              if (igal.le.ngal) ig=sid(igal)
            endif
          enddo
          close(10)
        enddo
      endif

      close(15)

      return
      END
`,

  gen_star_candidate_direct: `
      subroutine gen_star_candidate_direct(PREFIX,nx,ny,array,weight
     .,nstar,proc_error)
      implicit none
      include 'para.inc'

      character*(strl) catname,filename,PREFIX
      integer ierror,flag,proc_error,nstar

      integer nx,ny
      real array(npx,npy)
      integer weight(npx,npy)

      integer i,j,u,v,nn1,nn2,ig

      real star_para(nstar_max,npara),aa(npara)
      real source_coll(ngal_max,ns,ns)
      real noise_coll(ngal_max,ns,ns)

      real source(ns,ns),noise(ns,ns)

      real SNR,temp

      real xp,yp,sig,total_flux,half_light_flux,peak,rf
      integer total_area,half_light_area
      common /stamp_pass/ xp,yp,total_area,half_light_area,sig
     .,total_flux,half_light_flux,peak,rf

      nstar=0

      if (proc_error.eq.1) goto 40

      catname=trim(PREFIX)//'.cat'

      open(unit=10,file=catname,status='old',iostat=ierror)
      rewind 10

      if (ierror.ne.0) then
        write(*,*) catname
        pause 'catalog file error!!'
      endif
      read(10,*)

      do while (ierror.ge.0)

        read(10,*,iostat=ierror) xp,yp,total_area
     .,half_light_area,sig,total_flux,half_light_flux,peak,rf

        if (ierror.lt.0) cycle

        temp=half_light_area
        SNR=half_light_flux/sqrt(temp)
        if (SNR.lt.SNR_PSF*0.5) cycle

        call find_noise(flag,noise,nx,ny,array,weight)
        if (flag.lt.0) cycle

        call check_source(flag,source,nx,ny,array,weight)
        if (flag.lt.0) cycle

        temp=half_light_area
        SNR=half_light_flux/sqrt(temp)

        if (SNR.lt.SNR_PSF) cycle

        nstar=nstar+1
        do u=1,ns
          do v=1,ns
            source_coll(nstar,u,v)=source(u,v)
            noise_coll(nstar,u,v)=noise(u,v)
          enddo
        enddo
        star_para(nstar,1)=nstar
        star_para(nstar,2)=xp
        star_para(nstar,3)=yp
        star_para(nstar,4)=SNR

        if (nstar.ge.nstar_max) exit

      enddo
      close(10)

40    filename=trim(PREFIX)//'_star_can_info.dat'
      open(unit=20,file=filename,status='replace')
      rewind 20
      write(20,*) 'ig xp yp SNR'
      do i=1,nstar
        write(20,*) (star_para(i,j),j=1,4)
      enddo
      close(20)

      if (nstar.gt.0) then
        filename=trim(PREFIX)//'_star_can.fits'
        nn1=ns*len_s
        nn2=ns*(int(nstar/len_s)+1)
        call write_stamps(ngal_max,1,nstar,ns,ns
     .,source_coll,nn1,nn2,filename)

        filename=trim(PREFIX)//'_star_can_noise.fits'
        nn1=ns*len_s
        nn2=ns*(int(nstar/len_s)+1)
        call write_stamps(ngal_max,1,nstar,ns,ns
     .,noise_coll,nn1,nn2,filename)
      endif


      return
      END
`,
  read_stamps: `
	subroutine read_stamps(np,nstart,n,nx,ny,stamp,n1,n2
     .,filename)
	implicit none

	integer np,nstart,n,nx,ny,n1,n2
	integer NMAX
	parameter (NMAX=7000)
	real stamp(np,nx,ny),large_stamp(NMAX,NMAX)
	character filename*(*)
	integer i,j,offx,offy,k

      if ((n2*ny).ge.NMAX) then
        write(*,*) 'large_stamp is too small!!'
        stop
      endif
      call readimage(filename,n1,n2,NMAX,NMAX,large_stamp)

	offx=0	
	offy=0
	
	do k=nstart,n
	  do i=1,nx
	    do j=1,ny
	      stamp(k,i,j)=large_stamp(i+offx,j+offy)
	    enddo
	  enddo
	  offx=offx+nx
	  if (offx+nx.gt.n1) then
	    offx=0
	    offy=offy+ny
	  endif
	enddo


	return
	end
  `,
  get_power: `
      subroutine get_power(n1,n2,map,power,smooth)
      implicit none

c      n1 and n2 must be even numbers!!!!

      integer n1,n2,smooth
      real map(n1,n2),power(n1,n2)
      complex arr(n1,n2)
      integer i,j,ii,jj,n1_2,n2_2
      real pc
      common /pc_pass/ pc


      n1_2=n1/2
      n2_2=n2/2

      do i=1,n1
        do j=1,n2
          arr(i,j)=complex(real(map(i,j)),0.)
        enddo
      enddo

      call FFT2D(n1,n2,arr,1)

      do i=1,n1
        ii=mod(i+n1_2-1,n1)+1
        do j=1,n2
          jj=mod(j+n2_2-1,n2)+1
          power(ii,jj)=abs(arr(i,j))**2
        enddo
      enddo

      pc=power(n1_2+1,n2_2+1)

      if (smooth.eq.1) then
        call smooth_image55_hole(n1,n2,power)
      elseif (smooth.eq.2) then
        call smooth_image55_hole_ln(n1,n2,power)
      endif

      return
      end
`,
  process_powers: `
      subroutine process_powers(n,sourcep,noisep)
      implicit none

      integer n,offc,cc
      real sourcep(n,n),noisep(n,n)
      integer i,j
      real temp

      do i=1,n
        do j=1,n
          sourcep(i,j)=sourcep(i,j)-noisep(i,j)
        enddo
      enddo

      temp=0.
      do i=2,n-1
        temp=temp+sourcep(i,1)+sourcep(i,n)+sourcep(1,i)+sourcep(n,i)
      enddo

      temp=temp/(4.*(n-2.))

      do i=1,n
        do j=1,n
          sourcep(i,j)=sourcep(i,j)-temp
        enddo
      enddo

      return
      end
`,
  write_stamps: `
      subroutine write_stamps(np,nstart,n,nx,ny,stamp,n1,n2
     .,filename)
      implicit none

      integer np,nstart,n,nx,ny,n1,n2
      integer NMAX
      parameter (NMAX=7000)
      real stamp(np,nx,ny),large_stamp(NMAX,NMAX)
      character filename*(*)
      integer i,j,offx,offy,k

      do i=1,n1
        do j=1,n2
          large_stamp(i,j)=0.
        enddo
      enddo

      offx=0	
      offy=0

      do k=nstart,n
              if (offy+ny.gt.n2) then
                write(*,*) 'large_stamp is too small!!'
                stop
              endif
        do i=1,nx
          do j=1,ny
            large_stamp(i+offx,j+offy)=stamp(k,i,j)
          enddo
        enddo
        offx=offx+nx
        if (offx+nx.gt.n1) then
          offx=0
          offy=offy+ny
        endif
      enddo

      call writeimage(filename,n1,n2,NMAX,NMAX,large_stamp)

      return
      end  
  `,
  regularize_power: `
      subroutine regularize_power(nx,ny,power,star_smooth)
      implicit none

      integer nx,ny,star_smooth
      real power(nx,ny),temp
      integer i,j,cx,cy

      cx=nx/2+1
      cy=ny/2+1

      if (star_smooth.ge.1) then

        temp=1./power(cx,cy)

        do i=1,nx
          do j=1,ny
            power(i,j)=power(i,j)*temp
          enddo
        enddo

      else

        temp=4./(power(cx+1,cy)+power(cx-1,cy)
     .+power(cx,cy+1)+power(cx,cy-1))

        do i=1,nx
          do j=1,ny
            power(i,j)=power(i,j)*temp
          enddo
        enddo
        power(cx,cy)=1.

      endif

      return
      end
  `,
  xy_to_xxyy: `
      subroutine xy_to_xxyy(x,y,xx,yy,cRPIX,cD)
      implicit none

      double precision x,y,xx,yy
      double precision cRPIX(2),cD(2,2)

      xx=cD(1,1)*(x-cRPIX(1))+cD(1,2)*(y-cRPIX(2))
      yy=cD(2,1)*(x-cRPIX(1))+cD(2,2)*(y-cRPIX(2))

      return
      end
`,
get_power_all: `
      subroutine get_power_all(nx,ny,power,e,size,thresh_ratio)
      implicit none

      integer nx,ny
      real e(2),thresh_ratio,power(nx,ny),size
      integer area

      call get_power_area(nx,ny,power,area,thresh_ratio)
      size=area
      call get_power_e(nx,ny,power,e,thresh_ratio)

      return
      end
`,

get_PSF_FWHM: `
      subroutine get_PSF_FWHM(power,FWHM)
      implicit none
      include 'para.inc'

      real power(ns,ns),thresh
      integer i,j
      real area,FWHM,beta

      thresh=power(ns/2+1,ns/2+1)*exp(-1.)

      area=0.
      do i=1,ns
        do j=1,ns
          if (power(i,j).ge.thresh) area=area+1.
        enddo
      enddo

      beta=ns/(2.*pi)/sqrt(area/pi)
      FWHM=beta*2.*sqrt(2.*log(2.))*pixel_size

      return
      end
`,

ana_chi2: `
      subroutine ana_chi2(n,map1,map2,p)
      implicit none

      integer n
      real map1(n,n),map2(n,n),p,flux
      integer i,j,n1,n2

      n1=n/4
      n2=(n/4)*3
      p=0.
      flux=0.
      do i=n1,n2
        do j=n1,n2
          flux=flux+(map1(i,j)+map2(i,j))*0.5
          p=p+(map1(i,j)-map2(i,j))**2
        enddo
      enddo
      p=p/flux

      return
      end
`,

sort: `
      SUBROUTINE sort(n,np,arr)
      INTEGER n,np,M,NSTACK
      REAL arr(np)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
          write(*,*) 'NSTACK too small in sort'
          read(*,*)
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
`,

get_peak_width_low_side: `
      subroutine get_peak_width_low_side(np,n,arr,p,sig)
      implicit none

      integer np,n
      real arr(np),p,sig,den(n),den2(n),thresh
      integer i,j,ip

      call sort(n,np,arr)

      do i=2,n-1
        den(i)=((arr(i+1)-arr(i-1))/2.)**2
      enddo

      den(1)=(arr(2)-arr(1))**2
      den(n)=(arr(n)-arr(n-1))**2

      den2=0.
      do i=3,n-2
        do j=i-2,i+2
          den2(i)=den2(i)+den(j)
        enddo
        den2(i)=1./sqrt(den2(i)/5.)
      enddo

      den2(1)=1./sqrt((den(1)+den(2)+den(3))/3.)
      den2(2)=1./sqrt((den(1)+den(2)+den(3)+den(4))/4.)
      den2(n)=1./sqrt((den(n)+den(n-1)+den(n-2))/3.)
      den2(n-1)=1./sqrt((den(n)+den(n-1)+den(n-2)+den(n-3))/4.)

      ip=1
      p=den2(1)

      do i=2,n
        if (den2(i).gt.p) then
          p=den2(i)
          ip=i
        endif
      enddo

      thresh=p*0.5
      do i=ip+1,n
        if (den2(i).lt.thresh) exit
      enddo

      p=arr(ip)
      sig=arr(i)-arr(ip)

      return
      end
`,

write_stamps_2: `
      subroutine write_stamps_2(np,n,nx,ny,stamp,opt,val,n1,n2
          .,filename)
      implicit none

      integer np,n,nx,ny,n1,n2
      integer NMAX,opt(np),val
      parameter (NMAX=7000)
      real stamp(np,nx,ny),large_stamp(NMAX,NMAX)
      character filename*(*)
      integer i,j,offx,offy,k

      do i=1,n1
        do j=1,n2
          large_stamp(i,j)=0.
        enddo
      enddo

      offx=0	
      offy=0

      do k=1,n
              if (offy+ny.gt.n2) then
                write(*,*) 'large_stamp is too small!!'
                read(*,*)
              endif
              if (opt(k).ne.val) cycle
        do i=1,nx
          do j=1,ny
            large_stamp(i+offx,j+offy)=stamp(k,i,j)
          enddo
        enddo
        offx=offx+nx
        if (offx+nx.gt.n1) then
          offx=0
          offy=offy+ny
        endif
      enddo

            call writeimage(filename,n1,n2,NMAX,NMAX,large_stamp)

      return
      end
`,

draw_shear_expo: `
      subroutine draw_shear_expo(n,map,ncp,ncpp,nc,pc,nsp,ns
     .,sk,intensity,thickness)
      implicit none

      integer n,ncp,nc,ns,nsp,ncpp
      real map(n,n),pc(ncpp,4),sk(nsp,5),e(ns),tt(ns)
      real intensity,thickness,xmin,xmax,ymin,ymax
      real ratiox,ratioy,dx,dy,margin,dd,tmp
      real kmin,kmax,emax,ratio,cos2t,cost,sint
      real x1,x2,y1,y2
      parameter (margin=0.05)
      integer i,j

      do i=1,n
        do j=1,n
          map(i,j)=0.
        enddo
      enddo
      xmin=1e10
      xmax=-xmin
      ymin=1e10
      ymax=-ymin
      do i=1,nc
        xmin=min(pc(i,1),xmin)
        xmin=min(pc(i,3),xmin)
        xmax=max(pc(i,1),xmax)
        xmax=max(pc(i,3),xmax)
        ymin=min(pc(i,2),ymin)
        ymin=min(pc(i,4),ymin)
        ymax=max(pc(i,2),ymax)
        ymax=max(pc(i,4),ymax)
      enddo
      dx=(xmax-xmin)*margin
      dy=(ymax-ymin)*margin
      xmin=xmin-dx
      xmax=xmax+dx
      ymin=ymin-dy
      ymax=ymax+dy
      ratiox=(n-1.)/(xmax-xmin)
      ratioy=(n-1.)/(ymax-ymin)
      do i=1,nc
        x1=(pc(i,1)-xmin)*ratiox+1.
        x2=(pc(i,3)-xmin)*ratiox+1.
        y1=(pc(i,2)-ymin)*ratioy+1.
        y2=(pc(i,4)-ymin)*ratioy+1.
        call draw_rectangle(n,n,map,x1,y1,x2,y2,intensity,thickness)
      enddo
      do i=1,ns
        sk(i,1)=(sk(i,1)-xmin)*ratiox+1.
        sk(i,2)=(sk(i,2)-ymin)*ratioy+1.
      enddo
      tmp=ns
      dd=n/(sqrt(tmp))*0.25

      kmin=10000.
      kmax=-kmin
      emax=0.
      do i=1,ns
        kmin=min(kmin,sk(i,3))
        kmax=max(kmax,sk(i,3))
        e(i)=sqrt(sk(i,4)**2+sk(i,5)**2)
        tt(i)=e(i)
      enddo
      call sort(ns,ns,tt)
      tmp=max(abs(kmax),abs(kmin))*0.01
      ratio=intensity*0.7/(kmax-kmin+tmp)

      do i=1,ns
        tmp=(sk(i,3)-kmin)*ratio+1.
        call draw_box_fill(n,n,map,sk(i,1)-dd,sk(i,2)-dd
     .,sk(i,1)+dd,sk(i,2)+dd,tmp)
      enddo

      ratio=dd/tt(ns*3/4)
      do i=1,ns
        if (e(i).le.0.) cycle
        cos2t=sk(i,4)/e(i)
        cost=sqrt((1.+cos2t)*0.5)
        sint=sqrt((1.-cos2t)*0.5)
        if (sk(i,5).gt.0) sint=-sint
        x1=sk(i,1)-e(i)*ratio*cost
        x2=sk(i,1)+e(i)*ratio*cost
        y1=sk(i,2)-e(i)*ratio*sint
        y2=sk(i,2)+e(i)*ratio*sint
        call draw_line(n,n,map,x1,y1,x2,y2,intensity,0.)
      enddo

      call reverse_color(n,n,map)

      return
      end
`,

interpolate_PSF: `
      subroutine interpolate_PSF(nsam,npsam,image,posi,ns
     .,npp,nppx,PSF_coe)
      implicit none

      integer ns,npp,nppx,nsam,npsam
      double precision posi(npsam,2),PSF_coe(ns,ns,npp)
      real image(npsam,ns,ns)
      integer i,j,k
      real arr(nsam,3),coe(npp)

      do i=1,ns
        do j=1,ns
          do k=1,nsam
            arr(k,1)=posi(k,1)
            arr(k,2)=posi(k,2)
            arr(k,3)=image(k,i,j)
          enddo
c          call fit_2D(nsam,nsam,arr,npp,nppx,coe)
          call fit_2D_2(nsam,nsam,arr,npp,coe)
          do k=1,npp
            PSF_coe(i,j,k)=coe(k)
          enddo
        enddo
      enddo

      return
      END
`,

get_PSF_model: `
      subroutine get_PSF_model(ns,npp,nppx,PSF_coe,xx,yy,model)
      implicit none

      integer ns,npp,nppx
      double precision xx,yy,PSF_coe(ns,ns,npp)
      real model(ns,ns),x,y,coe(npp),func_val,func_val_2
      integer i,j,k

      x=xx
      y=yy

      do i=1,ns
        do j=1,ns
          do k=1,npp
            coe(k)=PSF_coe(i,j,k)
          enddo
c          model(i,j)=func_val(x,y,npp,nppx,coe)
          model(i,j)=func_val_2(x,y,npp,coe)
        enddo
      enddo

      return
      end
`,
  get_PSF_area: `
      subroutine get_PSF_area(model,FWHM)
      implicit none
      include 'para.inc'

      real model(ns,ns),area,thresh,FWHM,r,beta
      integer i,j

      thresh=exp(-1.) * model(ns/2+1,ns/2+1)

      area=-1e-5
      do i=1,ns
        do j=1,ns
          if (model(i,j).ge.thresh) area=area+1.
        enddo
      enddo
      if (area.le.0.) then
        FWHM=-1.
        return
      endif
      r=sqrt(area/pi)
      beta=ns/(2.*pi)/r
      FWHM=beta*2.*sqrt(2.*log(2.))*0.2628

      return
      end
`,
  field_distortion_PU: `
      subroutine field_distortion_PU(x,y,npd
     .,PU,cD,cRPIX,g1,g2,cos2,sin2,parity)
      implicit none

      integer npd
      double precision xi,eta,x,y,xx,yy,g1,g2,cos2,sin2,aa,bb
      double precision cRPIX(2),cD(2,2),PU(2,npd),cos1,sin1
      integer parity

      double precision mat(2,2),dm(2,2),temp,det,sqrt_det
      double precision dx_dxx,dx_dyy,dy_dxx,dy_dyy
      double precision dxx_dxi,dxx_deta,dyy_dxi,dyy_deta
      double precision r,dr_dxi,dr_deta
      integer px,py,n,order

      ! error will occur if both x==cRPIX(1) and y==cRPIX(2) !
      xx=cD(1,1)*(x-cRPIX(1))+cD(1,2)*(y-cRPIX(2))
      yy=cD(2,1)*(x-cRPIX(1))+cD(2,2)*(y-cRPIX(2))

      call mapping_PU(xx,yy,xi,eta,npd,PU,1)
 
      temp=1d0/(cD(1,1)*cD(2,2)-cD(1,2)*cD(2,1))

      dx_dxx=cD(2,2)*temp
      dx_dyy=-cD(1,2)*temp
      dy_dxx=-cD(2,1)*temp
      dy_dyy=cD(1,1)*temp

      dxx_dxi=1d0
      dyy_deta=1d0
      dxx_deta=0d0
      dyy_dxi=0d0

c      r=sqrt(xi**2+eta**2)
c      if (r.gt.0d0) then
c        dr_dxi=xi/r
c        dr_deta=eta/r
c      else
c        dr_dxi=0d0
c        dr_deta=0d0
c      endif

      px=0
      py=1
      order=1
      n=0
      do while (n.lt.npd)
        if (py.eq.order) then
c          if (mod(order,2).eq.1) then
c            n=n+1
c            dxx_dxi=dxx_dxi-PU(1,n)*order*r**(order-1)*dr_dxi
c            dxx_deta=dxx_deta-PU(1,n)*order*r**(order-1)*dr_deta
c            dyy_deta=dyy_deta-PU(2,n)*order*r**(order-1)*dr_deta
c            dyy_dxi=dyy_dxi-PU(2,n)*order*r**(order-1)*dr_dxi
c            if (n.eq.npd) cycle
c          endif
          order=order+1
          px=order
          py=0
        else
          px=px-1
          py=py+1
        endif
        n=n+1
        dxx_dxi=dxx_dxi-PU(1,n)*px*xi**(px-1)*eta**py
        dxx_deta=dxx_deta-PU(1,n)*xi**px*py*eta**(py-1)
        dyy_deta=dyy_deta-PU(2,n)*px*eta**(px-1)*xi**py
        dyy_dxi=dyy_dxi-PU(2,n)*eta**px*py*xi**(py-1)

      enddo


      mat(1,1)=dx_dxx*dxx_dxi+dx_dyy*dyy_dxi
      mat(1,2)=dx_dxx*dxx_deta+dx_dyy*dyy_deta
      mat(2,1)=dy_dxx*dxx_dxi+dy_dyy*dyy_dxi
      mat(2,2)=dy_dxx*dxx_deta+dy_dyy*dyy_deta

      det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
      temp=1d0/det

      dm(1,1)=mat(2,2)*temp
      dm(1,2)=-mat(1,2)*temp
      dm(2,1)=-mat(2,1)*temp
      dm(2,2)=mat(1,1)*temp

      parity=1
      if (det.lt.0) then
        dm(1,1)=-dm(1,1)
        dm(1,2)=-dm(1,2)
        parity=-1
      endif

      sqrt_det=sqrt(dm(1,1)*dm(2,2)-dm(1,2)*dm(2,1))

      dm(1,1)=dm(1,1)/sqrt_det
      dm(1,2)=dm(1,2)/sqrt_det
      dm(2,1)=dm(2,1)/sqrt_det
      dm(2,2)=dm(2,2)/sqrt_det

      cos1=0.5d0*(dm(1,1)+dm(2,2))
      sin1=0.5d0*(dm(1,2)-dm(2,1))

      aa=-0.5d0*(dm(1,2)+dm(2,1))
      bb=0.5d0*(dm(2,2)-dm(1,1))

      g1=aa*sin1+bb*cos1
      g2=aa*cos1-bb*sin1

      cos2=cos1*cos1-sin1*sin1
      sin2=2d0*sin1*cos1

      if (parity.eq.-1) g2=-g2

      return
      end
`,
  get_shear: `
      SUBROUTINE get_shear(n,gal,psf,g1,g2,de,h1,h2)
      implicit none

      integer n
      real g1,g2,de,h1,h2,gal(n,n),psf(n,n)
      integer i,j,n_2,cc
      real ks,peak,thresh,area,ks_2,kx,kx2,ky,ky2,k2,k,temp,temp1
      real filter_deriv,filter,ff,norm

      real pi,PSFr_ratio
      parameter (pi=3.1415926)
      parameter (PSFr_ratio=0.75)

      peak=psf(1,1)
      do i=1,n
        do j=1,n
          if (psf(i,j).gt.peak) peak=psf(i,j)
        enddo
      enddo

      thresh=exp(-1.)*peak

      area=0.
      do i=1,n
        do j=1,n
          if (psf(i,j).ge.thresh) area=area+1.
        enddo
      enddo

      ks=sqrt(area/pi)
      ks_2=(ks*PSFr_ratio)**(-2)

      thresh=peak*1e-5
      n_2=n/2
      cc=1+n_2

      g1=0.
      g2=0.
      de=0.
      h1=0.
      h2=0.

      do i=1,n
        kx=i-cc
        kx2=kx*kx
        do j=1,n
          ky=j-cc
          ky2=ky*ky
          k2=kx2+ky2
          k=sqrt(k2)
          if (psf(i,j).gt.thresh) then
            ff=k2*ks_2
            temp=exp(-ff)/psf(i,j)              ! 相当于T
            temp1=temp*gal(i,j)                 ! 相当于T*M

            g1=g1-temp1*(kx2-ky2)
            g2=g2-temp1*2.*kx*ky
            de=de+temp1*k2*(2.-ff)              ! N
            h1=h1+temp1*ks_2*(k2*k2-8.*kx2*ky2) ! U
            h2=h2+temp1*ks_2*4.*kx*ky*(kx2-ky2)

          endif
        enddo
      enddo


      return
      END
  `,
  get_chip_id: `
      subroutine get_chip_id(imagefile,id)
      implicit none

      character*(*) imagefile
      character*2 id_letter
      integer i,p1,p2,n,id

      call read_ccDNUM(imagefile,id)

      return
      end
  `,
  read_ccDNUM: `
      SUBROUTINE read_CCDNUM(filename,ccd_num)
! read the size of the 2D fits image from a fits file.
	  IMPLICIT NONE

      INTEGER status,unit,readwrite,blocksize,nfound
      INTEGER ccd_num
      CHARACTER filename*(*)
      character comment*100
      
      status=0
      
      CALL ftgiou(unit,status)

      readwrite=0
      CALL ftopen(unit,filename,readwrite,blocksize,status)
      CALL ftgkyj(unit,'CCDNUM',ccd_num,comment,status)

      CALL ftclos(unit, status)
      CALL ftfiou(unit, status)

      IF (status .gt. 0) CALL printerror(status)
      return
      END   
  `,
  
};
export default functionCodeMap;