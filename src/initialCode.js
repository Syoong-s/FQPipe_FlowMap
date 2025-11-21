const initialCode = `
      call MPI_Init(ierr)
      call MPI_comm_rank(MPI_cOMM_WORLD, my_id, ierr )
      call MPI_comm_size(MPI_cOMM_WORLD, num_procs, ierr )
      ! my_id = 0,1,...,num_procs-1
      call MPI_BARRIER(MPI_cOMM_WORLD, ierr ) ! synchronize all nodes

      call getarg(1,EXPO_LIST)
      ! ----------------------------------------------------------------------------
      if (my_id.eq.0) call initialize(EXPO_LIST)
      call MPI_Bcast(N_EXPO,1,mpi_int,0,MPI_cOMM_WORLD,ierr)
      call MPI_Bcast(EXPO_FILE,NMAX_EXPO*strl,mpi_character,0
     .,MPI_cOMM_WORLD,ierr)

      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)
      ! ----------------------------------------------------------------------------
      if (mod(PROcESS_stage,2).eq.0)
     .call mpi_distribute(N_EXPO,pre_process,'Pre-process...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,3).eq.0)
     .call mpi_distribute(N_EXPO,proc_astrometry,'Astrometry...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,5).eq.0)
     .call mpi_distribute(N_EXPO,proc_source,'Sources ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,7).eq.0)
     .call mpi_distribute(N_EXPO,proc_FourierT,'FourierT ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,11).eq.0)
     .call mpi_distribute(N_EXPO,proc_PSF,'PSF ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      if (mod(PROcESS_stage,13).eq.0)
     .call mpi_distribute(N_EXPO,proc_shear,'Shear ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

      do iexpo=1,NMAX_EXPO
        do i=1,6
          expo_para(i,iexpo)=0.
        enddo
      enddo

      if (mod(PROcESS_stage,17).eq.0)
     .call mpi_distribute(N_EXPO,proc_info,'Info ...')

      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)
      call MPI_AllReduce(expo_para,expo_para_t,NMAX_EXPO*6
     .,mpi_REAL,MPI_SUM,MPI_cOMM_WORLD,ierr)

      do iexpo=1,NMAX_EXPO
        do i=1,6
          expo_para(i,iexpo)=expo_para_t(i,iexpo)
        enddo
      enddo

      if (my_id.eq.0) then
        call get_dir(EXPO_LIST,root_dir,1)
        filename=trim(root_dir)//'/expo_info.dat'
        open(unit=10,file=filename,status='replace')
        rewind 10
        write(10,*) 'N-valid-chip PSF-FWHM(arcsec) chi_d-stars'    
     .,' nstar-per-chip cRVAL1 cRVAL2 expo_name '
        do i=1,N_EXPO
          write(10,*) (expo_para(j,i),j=1,6),trim(EXPO_FILE(i))
        enddo
        close(10)
      endif

      if (mod(PROcESS_stage,19).eq.0)
     .call mpi_distribute(N_EXPO,proc_comb,'combine ...')
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)

c      !----- MPI finalization -----------------------------------------------
      call MPI_BARRIER(MPI_cOMM_WORLD,ierr)
      call MPI_Finalize (ierr)
`;

export default initialCode;