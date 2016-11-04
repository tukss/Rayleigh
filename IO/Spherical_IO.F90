Module Spherical_IO
    Use Spherical_Buffer
	Use Parallel_Framework
	Use SendReceive
    Use ISendReceive
    Use General_MPI
    Use Fourier_Transform
    Use Legendre_Transforms, Only : Legendre_Transform
    Use BufferedOutput
	Implicit None
	! This module contains routines for outputing spherical data as:
	! 1. Slices of sphere
	! 2. Phi-Averages over sphere (f(r,theta))
	! 3. Phi/theta Averages over sphere ((f(r))
	! 4. Volume averages over sphere
    ! 5. 3-D Output on Spherical grid
    ! 6. Spectra on slices of the sphere
    ! 7. PDFs taken on slices of the sphere

	!////////////////////////////////////////////
    Integer, Parameter :: nqmax=800, nshellmax=100
    Integer, Parameter :: endian_tag = 314      ! first 4 bits of each diagnostic file - used for assessing endianness on read-in

    !Each diagnostic type has an associated version number that is written to the file
    !following the endian tag.  Hopefully, reading routines can be backward compatible if
    !significant changes are made to the output structure (and reflected in the version number)
  
    Integer, Parameter :: shellslice_version = 3
    Integer, Parameter :: azavg_version = 3
    Integer, Parameter :: shellavg_version = 4
    Integer, Parameter :: globalavg_version = 3
    Integer, Parameter :: shellspectra_version = 3
    Integer, Parameter :: full3d_version = 3    !currently unused
    Type, Public :: DiagnosticInfo
        ! Need to see if we can make these allocatable, but for now..
        ! Each instance of this class has two static arrays used for reading in namelist input
        ! These arrays are large enough to hold nqmax and nshellmax values, but 
        ! typically only a small fraction of that amount will be specified at input.
        Integer :: values(1:nqmax) = -1 ! The list of values specified in an input namelist
        Integer :: levels(1:nshellmax) = -1 ! The radial indices output (shell slices, spectra, and histograms only)
        Integer :: compute(1:nqmax)= -1 ! compute(i) = 1 if i was specified in values.  -1 otherwise
        
        Integer :: nq, nlevels ! Number of nonzero elements of values and levels
        Integer :: my_nlevels  ! Number of nonzero elements of levels that are in process

        Integer :: file_unit = 15
        Character*120 :: file_prefix = 'None'

        Integer, Allocatable :: oqvals(:)   ! Array of size nq used by I/O process to record output ordering of diagnostics

        Integer :: frequency = 90000000 ! How often we write this diagnostic type
        Integer :: rec_per_file =1     ! How many of these records we write to a file
        Integer :: current_rec = 1       ! Which record we are on within a file
        Integer :: file_header_size =0 ! Size of file header in bytes
        Integer :: file_record_size = 0 ! Size of each record in bytes
        Integer :: file_position = 1   ! Keep track of where we are (byte-wise) within the file
        Integer :: avg_level = 0       ! global_averages = 3, shell_averages = 1, az_averages = 1, all others = 0

        Integer :: ind = 1              ! index of current diagnostic being stored (advances throughout an output calculation)
        Logical :: begin_output = .false.
        Integer :: mpi_tag = 1          ! For use when communicating before writing the file

        Integer :: output_version = -1  ! See note at top

        !This flag changes throughout the diagnostic evaluation
        !If the qty being output at the current iteration is supposed to be
        !output as this diagnostic type, the grab_this_q will be set to true.
        !Otherwise grab_this_q remains false.
        Logical :: grab_this_q = .false.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Additional Organizational Arrays are required for managing shell-slice like outputs
        Integer, Allocatable :: my_shell_levs(:), have_shell(:), my_shell_ind(:)
        Integer, Allocatable :: shell_r_ids(:), nshells_at_rid(:)
        Integer :: nshell_r_ids

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !Communicatory Info for parallel writing (if used)
        Integer :: ocomm, orank, onp
        Logical :: master = .false.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Methods
        Contains
        Procedure :: Init => Initialize_Diagnostic_Info
        Procedure :: AdvanceInd
        Procedure :: reset => diagnostic_output_reset
        Procedure :: Shell_Balance
        Procedure :: OpenFile
        Procedure :: OpenFile_Par
        Procedure :: Set_File_Info
        Procedure :: CloseFile
        Procedure :: CloseFile_Par
        Procedure :: update_position
        Procedure :: getq_now
        Procedure :: init_ocomm
        Procedure :: CleanUp
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! 


    End Type DiagnosticInfo

    Type(DiagnosticInfo) :: Shell_Averages, Shell_Slices, Global_Averages, AZ_Averages, Full_3D, Shell_Spectra

    Integer :: current_averaging_level = 0
    Integer :: current_qval = 0

    Integer :: averaging_level(1:nqmax) = 0, compute_q(1:nqmax) = 0
    Integer :: shellavg_values(1:nqmax)=-1, globalavg_values(1:nqmax)=-1
    Integer :: shellslice_values(1:nqmax) =-1, shellslice_levels(1:nshellmax)=-1, azavg_values(1:nqmax)=-1
    Integer :: full3d_values(1:nqmax) = -1, shellspectra_values(1:nqmax)=-1, shellspectra_levels(1:nshellmax)=-1
    Integer :: histo_values(1:nqmax) = -1, histo_levels(1:nshellmax)=-1

    Integer :: globalavg_nrec = 1, shellavg_nrec = 1, azavg_nrec = 1, shellslice_nrec =1, shellspectra_nrec =1

    Integer :: globalavg_frequency = 90000000, shellavg_frequency = 90000000
    Integer :: azavg_frequency = 90000000, shellslice_frequency = 90000000
    Integer :: shellspectra_frequency=90000000

    Integer :: full3d_frequency= 90000000
    Character*120 :: local_file_path=''
    Logical :: mem_friendly = .false.

    Namelist /output_namelist/shellavg_values, globalavg_values, &
        & shellslice_values, shellslice_levels, azavg_values, &
        & full3d_values, &
        & full3d_frequency, globalavg_nrec, shellavg_nrec, azavg_nrec, shellslice_nrec, &
        & globalavg_frequency, shellavg_frequency, azavg_frequency, shellslice_frequency, &
        & shellspectra_nrec, shellspectra_frequency, shellspectra_levels, shellspectra_values, &
        & mem_friendly 




    Integer :: integer_zero = 0
    Real*8, Private, Allocatable :: circumference(:,:), qty(:,:,:), f_of_r_theta(:,:)
    Real*8, Private, Allocatable :: azav_outputs(:,:,:), f_of_r(:), rdtheta_total(:)
    Real*8, Private, Allocatable :: shellav_outputs(:,:,:), globav_outputs(:), shell_slice_outputs(:,:,:,:)
    type(SphericalBuffer) :: spectra_buffer
    Real*8, Private :: da_total, int_vol, int_dphi, int_rsquared_dr, int_sintheta_dtheta
    Real*8, Private, Allocatable :: sintheta_dtheta(:), rsquared_dr(:)
    Character*6, Public :: i_ofmt = '(i8.8)', i_pfmt = '(i5.5)'
    integer :: output_ireq(1), output_status(1)

    !////////////////////////////////////////////////////////////////////
    Logical :: start_az_average, start_shell_average  !, start_shell_slice
    Logical :: start_pdf, start_global_average, start_spectra

    Integer :: io_node = 0

    Integer, Private :: current_iteration
    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ! We store some redundant information just to keep the IO level as independent of the Physics level as possible
    ! These variables are all private.
    Real*8, Allocatable,Private :: theta_integration_weights(:), r_integration_weights(:)
    Integer, Private :: my_theta_min, my_theta_max, my_rmin, my_rmax, my_mp_max, my_mp_min
    Integer, Private :: nr, nphi, ntheta, my_ntheta
    Integer, Private :: nproc1, nproc2, myid, nproc, my_row_rank, my_column_rank, my_nr
    Real*8, Allocatable, Private :: radius(:),sintheta(:),costheta(:) 
    Real*8 :: over_nphi_double

    !////////////////////////////////////////////////////////////
    ! And here are some variables used for storing averages so that moments may be taken
    ! These two arrays contain the ell0 and m0 representation of all quantities
    !   computed in shell- and (eventually) azimuthal- averaging.
    ! These two averages are used to compute moments
    Real*8, Private, Allocatable :: IOell0_values(:,:), IOm0_values(:,:,:)
    Integer, Private :: num_avg_store ! number of averages to store in the two arrays above
    Integer, Private :: IOavg_flag = -1
Contains

    Subroutine Begin_Outputting(iter)
		  Implicit None
		  Integer, Intent(In) :: iter
		  current_iteration = iter

        Call Global_Averages%reset()
        Call Shell_Averages%reset()
        Call AZ_Averages%reset()
        Call Shell_Slices%reset()
        Call Shell_Spectra%reset()

        Allocate(f_of_r_theta(my_rmin:my_rmax,my_theta_min:my_theta_max))
        Allocate(f_of_r(my_rmin:my_rmax))

        num_avg_store = shell_averages%nq
        Allocate(IOm0_values(my_rmin:my_rmax,my_theta_min:my_theta_max,1:num_avg_store))
        Allocate(IOell0_values(my_rmin:my_rmax,1:num_avg_store))
        IOm0_values(:,:,:) = 0.0d0
        IOell0_values(:,:) = 0.0d0
    End Subroutine Begin_Outputting



	Subroutine Initialize_Spherical_IO(rad_in,sintheta_in, rw_in, tw_in, costheta_in,file_path)
		Implicit None
		Integer :: k, fcount(3,2), ntot, fcnt, master_rank
		Real*8, Intent(In) :: rad_in(:), sintheta_in(:), rw_in(:), tw_in(:), costheta_in(:)
        Character*120 :: fdir
        Character*120, Intent(In) :: file_path
        local_file_path = file_path
		! Handles output bookkeeping

		my_theta_min = pfi%my_2p%min
		my_theta_max = pfi%my_2p%max
		my_rmin = pfi%my_1p%min
		my_rmax = pfi%my_1p%max
        
        my_mp_min = pfi%my_3s%min
        my_mp_max = pfi%my_3s%max

		my_nr = pfi%my_1p%delta
        my_ntheta = pfi%my_2p%delta

        nproc  = pfi%gcomm%np
		nproc1 = pfi%ccomm%np		! processor's per column (radius split among these)
		nproc2 = pfi%rcomm%np      ! processor's per row  (theta split among these)
		myid   = pfi%gcomm%rank
		my_row_rank = pfi%rcomm%rank
		my_column_rank = pfi%ccomm%rank
		nphi = pfi%n3p
		ntheta = pfi%n2p
		nr = pfi%n1p

        over_nphi_double = 1.0d0/nphi

		Allocate(sintheta(1:ntheta))
		Allocate(costheta(1:ntheta))
		sintheta(:) = sintheta_in(:)
		costheta(:) = costheta_in(:)
		Allocate(radius(1:nr))
		radius(:) = rad_in(:)

        Allocate(r_integration_weights(1:nr))
        Allocate(theta_integration_weights(1:ntheta))

        r_integration_weights(:) = rw_in(:)
        theta_integration_weights(:) = tw_in(:)		


        ! Map the various quantity lists etc. into their associated diagnostic structures

        !Numbers here are the mpi_tags used in communication for each output
        !In theory they can be the same, but it's probably a good idea to keep them unique
        Call            Full_3D%Init(averaging_level,compute_q,myid, &
            & 54,values = full3d_values)

        Call    Global_Averages%Init(averaging_level,compute_q,myid, &
            & 55,avg_level = 3,values = globalavg_values)

		Call        AZ_Averages%Init(averaging_level,compute_q,myid, &
            & 56,avg_level = 1,values = azavg_values)

        Call     Shell_Averages%Init(averaging_level,compute_q,myid, &
            & 57,avg_level = 2,values = shellavg_values)

        Call       Shell_Slices%Init(averaging_level,compute_q,myid, &
            & 58,values = shellslice_values, levels = shellslice_levels)

        Call       Shell_Spectra%Init(averaging_level,compute_q,myid, &
            & 59,values = shellspectra_values, levels = shellspectra_levels)

        !Outputs involve saving and communicating partial shell slices (e.g. Shell_Slices or spectra)
        !require an additional initialization step to load-balance the shells
        Call Shell_Slices%Shell_Balance()
        Call Shell_Spectra%Shell_Balance()
        if (my_row_rank .eq. 0) Then

            If (Shell_Slices%nshell_r_ids .gt. 0) Then
                master_rank = shell_slices%shell_r_ids(1)
                Call Shell_Slices%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,master_rank) ! For parallel IO
            Endif
            Call AZ_Averages%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,0) ! 0 handles file headers etc. for AZ Average output

            If (Shell_Spectra%nshell_r_ids .gt. 0) Then
                master_rank = shell_spectra%shell_r_ids(1)
                Call Shell_Spectra%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,master_rank) 
            Endif
        Endif
        
        If (Shell_Spectra%nlevels .gt. 0) Then
            !Shell spectra require an additional step
            !Initialize the buffer object that we use for transposing spectra
            !Some row ranks might have no buffer initialized.
            ntot = Shell_Spectra%nq*Shell_Spectra%my_nlevels
            fcnt = ntot/my_nr
            k = Mod(ntot,my_nr)
            if (k .gt. 0) fcnt = fcnt+1
            If (fcnt .gt. 0) Then
                fcount(:,:) = fcnt
           		Call spectra_buffer%init(field_count = fcount, config = 'p3b')		
            Endif	
        Endif 
        !Next, provide the file directory, frequency info, output versions etc.

        ! Global Averages
        fdir = 'G_Avgs/'
        Call Global_Averages%set_file_info(globalavg_version,globalavg_nrec,globalavg_frequency,fdir)    

        ! Shell Averages
        fdir = 'Shell_Avgs/'
        Call Shell_Averages%set_file_info(shellavg_version,shellavg_nrec,shellavg_frequency,fdir)    

        ! Azimuthal Averages
        fdir = 'AZ_Avgs/'
        Call AZ_Averages%set_file_info(azavg_version,azavg_nrec,azavg_frequency,fdir)    

        ! Shell Slices
        fdir = 'Shell_Slices/'
        Call Shell_Slices%set_file_info(shellslice_version,shellslice_nrec,shellslice_frequency,fdir)   

        ! Shell Spectra
        fdir = 'Shell_Spectra/'
        Call Shell_Spectra%set_file_info(shellspectra_version,shellspectra_nrec,shellspectra_frequency,fdir) 


        ! Full 3D (special because it only cares about the frequency, not nrec)
        fdir = 'Spherical_3D/'
        Call Full_3D%set_file_info(full3d_version,shellslice_nrec,full3d_frequency,fdir)    

        !Write(6,*)'Shell F: ', Shell_Slices%frequency
   End Subroutine Initialize_Spherical_IO




	Subroutine Get_Shell_Slice(qty)
		Implicit None
		Integer :: j, ilocal, shell_lev_ind, shell_ind
		Real*8, Intent(In) :: qty(:,:,my_theta_min:)
        If (Shell_Slices%nlevels .gt. 0) Then
            shell_ind = Shell_Slices%ind
            !If (myid .eq. 0) THen
                !NOTE:  Same as other remark when allocating.  Really should use master node here
                Shell_Slices%oqvals(shell_ind) = current_qval
            !Endif
        

		    If (Shell_Slices%my_nlevels .gt. 0) Then

		      If (Shell_Slices%begin_output) Then
			    Allocate(shell_slice_outputs(1:nphi,my_theta_min:my_theta_max,1:Shell_Slices%my_nlevels,1:Shell_Slices%nq))
		      Endif

            
		      shell_lev_ind =1
		      Do j = 1, Shell_Slices%nlevels
		        ilocal = Shell_Slices%levels(j)-my_rmin+1
		        If (Shell_Slices%have_shell(j) .eq. 1) Then ! my processor has this radius
				    shell_slice_outputs(:,my_theta_min:my_theta_max,shell_lev_ind,shell_ind) = &
                        &  qty(:,ilocal,my_theta_min:my_theta_max)
				    shell_lev_ind = shell_lev_ind +1
		        Endif
		      Enddo


		    Endif
            ! advance counter for next quantity to store (if any)
            Call Shell_Slices%AdvanceInd()
        Endif
	End Subroutine Get_Shell_Slice

	Subroutine Get_Shell_Spectra(qty)
		Implicit None
		Integer :: j, ilocal, shell_ind, field_ind, rind, counter
        Integer :: k, jj
		Real*8, Intent(In) :: qty(1:,1:,my_theta_min:)

        If (Shell_Spectra%nlevels .gt. 0) Then
            shell_ind = Shell_Spectra%ind
            !If (myid .eq. 0) Then
                Shell_Spectra%oqvals(shell_ind) = current_qval
            !Endif
        

		    If (Shell_Spectra%my_nlevels .gt. 0) Then

		        If (Shell_Spectra%begin_output) Then
			        Call spectra_buffer%construct('p3b')
                    spectra_buffer%p3b(:,:,:,:) = 0.0d0
		        Endif

                            

		        Do j = 1, Shell_Spectra%my_nlevels



		            ilocal = Shell_Spectra%my_shell_levs(j)-my_rmin+1

                    counter = (shell_ind-1)*Shell_Spectra%my_nlevels+ j-1 

                    field_ind = counter/my_nr+1
                    rind = MOD(counter,my_nr)+my_rmin

                    Do k = 1, nphi
                    Do jj = my_theta_min, my_theta_max
				        spectra_buffer%p3b(k,rind,jj,field_ind) = &
                        & qty(k, ilocal, jj)
                    Enddo
                    Enddo

		        Enddo
            Endif
            Call Shell_Spectra%AdvanceInd()
		Endif

	End Subroutine Get_Shell_Spectra

	Subroutine Write_Shell_Spectra(this_iter,simtime)
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:,:), all_spectra(:,:,:,:,:)
        Real*8, Allocatable :: sendbuffer(:,:,:,:,:), out_radii(:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, m, mp, lmax,rind,field_ind,f,r
        Integer :: rone,  p,  counter, nf
		Integer :: n, nn, this_nshell, nq_shell, shell_spectra_tag, nmodes
        Integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp
        Integer(kind=MPI_OFFSET_KIND)  :: qsize, qdisp, rec_size
		Integer :: your_mp_min, your_mp_max, your_nm, your_id
		Integer :: nelem, m_ind, m_val, current_rec
        Integer :: funit, error, sirq, inds(5), dims(3)
        Integer :: my_nlevels, nlevels
        Integer :: lp1, nrirqs, ind5
        Integer :: ierr, rcount, buffsize
        Integer, Allocatable :: rirqs(:)
        Integer :: mstatus(MPI_STATUS_SIZE)       

        nlevels = Shell_Spectra%nlevels             ! The total number of spectra levels that needs to be output
        my_nlevels = Shell_Spectra%my_nlevels       ! The number of radial levels that this rank needs to write out
        nq_shell = Shell_Spectra%nq                 ! The number of quantities 
        shell_spectra_tag = Shell_Spectra%mpi_tag
        funit = Shell_Spectra%file_unit
        lmax = maxval(pfi%inds_3s)
        lp1 = lmax+1
        nmodes = lp1*lp1
		responsible = 0
		If ( (my_row_rank .eq. 0) .and. (my_nlevels .gt. 0) )  responsible = 1


        
        !/////////////
        If (my_nlevels .gt. 0) Then
            !//////////////////////
            ! First thing we do is FFT/reform the buffer/Legendre Transform
            !
            Call FFT_To_Spectral(spectra_buffer%p3b, rsc = .true.)
            spectra_buffer%config ='p3b'
            Call spectra_buffer%reform()
            Call spectra_buffer%construct('s2b')
            Call Legendre_Transform(spectra_buffer%p2b,spectra_buffer%s2b)
            Call spectra_buffer%deconstruct('p2b')

            Allocate(sendbuffer(0:lmax,my_nlevels,nq_shell,2, my_mp_min:my_mp_max ))
            sendbuffer = 0.0d0 


            
            nf = spectra_buffer%nf2b
            Do p = 1, 2  ! Real and imaginary parts
            Do mp = my_mp_min,my_mp_max
                m = pfi%inds_3s(mp)
                    counter = 0
                    Do f = 1, nq_shell

                        field_ind = counter/my_nr+1
                        Do r = 1, shell_spectra%my_nlevels   
                                
                            rind = MOD(counter,my_nr)+my_rmin
                            sendbuffer(m:lmax,r,f,p,mp) = &
                                & spectra_buffer%s2b(mp)%data(m:lmax,rind,p,field_ind)
                            counter = counter+1
                        Enddo
                    Enddo
                Enddo

            Enddo
            call spectra_buffer%deconstruct('s2b')


        Endif

        If (responsible .eq. 1) Then
            ! Rank 0 in reach row receives  all pieces of the shell spectra from the other nodes

            Allocate(all_spectra(0:lmax,0:lmax, my_nlevels,nq_shell, 1:2))
            Allocate(buff(0:lmax,my_nlevels,nq_shell,1:2,1:lp1))  !note - indexing starts at 1 not zero for mp_min etc.
            all_spectra(:,:,:,:,:) = 0.0d0
            buff(:,:,:,:,:) = 0.0d0

            nrirqs = nproc2-1
            Allocate(rirqs(1:nrirqs))
            rirqs(:) = 0
            ind5 = pfi%all_3s(0)%delta+1
            Do nn = 1, nrirqs
                !Write(6,*)'Ind5: ', ind5
                your_id = nn

                your_nm     = pfi%all_3s(nn)%delta
                your_mp_min = pfi%all_3s(nn)%min
                your_mp_max = pfi%all_3s(nn)%max


                nelem = your_nm*my_nlevels*2*lp1*nq_shell

                inds(:) = 1
                inds(5) = ind5  !This is the mp_index here.

                Call Ireceive(buff, rirqs(nn), n_elements = nelem,source= your_id, &
                    &  tag=shell_spectra_tag,grp = pfi%rcomm, indstart = inds)
                ind5 = ind5+your_nm
            Enddo

            ! Stipe my own data into the receive buffer


            Do mp = my_mp_min,  my_mp_max
                m = pfi%inds_3s(mp)
                Do p = 1,2
                    Do f = 1, nq_shell
                        Do r = 1, my_nlevels   

                            buff(m:lmax,r,f,p,mp) = sendbuffer(m:lmax,r,f,p,mp) 

                        Enddo
                    Enddo
                Enddo

            Enddo
            !DeAllocate(sendbuffer)

            Call IWaitAll(nrirqs,rirqs)

            !Stripe the receiver buffer into the spectra buffer
           
            Do mp = 1,lp1
                m = pfi%inds_3s(mp)
                Do p = 1, 2  ! Real and imaginary parts
                    Do f = 1, nq_shell
                        Do r = 1, my_nlevels   
                            all_spectra(m:lmax,m,r,f,p) = buff(m:lmax,r,f,p,mp)  
                        Enddo
                    Enddo
                Enddo

            Enddo

            DeAllocate(sendbuffer)
            DeAllocate(buff)
            DeAllocate(rirqs)
        Else
			!  Non responsible nodes send their info
			If (my_nlevels .gt. 0) Then
                inds(:) = 1

				Call Isend(sendbuffer,sirq, dest = 0,tag=shell_spectra_tag, grp = pfi%rcomm, indstart = inds)
                Call IWait(sirq)
                DeAllocate(sendbuffer)
			Endif
		Endif


        If (my_row_rank .eq. 0) Call Shell_Spectra%OpenFile_Par(this_iter, error)



        If (responsible .eq. 1) Then   
            !Write(6,*)'I am responsible: ', my_column_rank
            funit = shell_spectra%file_unit
            current_rec = Shell_Spectra%current_rec  ! Note that we have to do this after the file is opened
            If  ( (current_rec .eq. 1) .and. (shell_spectra%master) ) Then                
                !Write(6,*)'I am master: ', my_column_rank
                dims(1) =  lmax
                dims(2) =  nlevels
                dims(3) =  nq_shell
                buffsize = 3
                call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_shell
                call MPI_FILE_WRITE(funit,Shell_Spectra%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                allocate(out_radii(1:nlevels))
                Do i = 1, nlevels
                    out_radii(i) = radius(Shell_Spectra%levels(i))
                Enddo
                buffsize = nlevels
	            call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
                DeAllocate(out_radii)
                
	            call MPI_FILE_WRITE(funit, Shell_Spectra%levels, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

            Endif

            ! I might really (really) want to look into file views later.
            ! Depending on the offset size mpi type, disp will crap out past 2GB
            hdisp = 24 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_shell*4 ! nq
            hdisp = hdisp+nlevels*12  ! level indices and level values
            


            ! The file is striped with time step slowest, followed by q


            rcount = 0
            Do p = 1, Shell_Spectra%nshell_r_ids
                if (Shell_Spectra%shell_r_ids(p) .lt. my_column_rank) Then
                    rcount = rcount+ Shell_Spectra%nshells_at_rid(p)
                Endif
            Enddo
            my_rdisp = rcount*nmodes*8

                

            ! This is the LOCAL number ELEMENTS in the real or imaginary component of
            ! of a single quantity  (This is not in bytes)
            buffsize = my_nlevels*nmodes 

            !This is the half-size (bytes) of a single quantity's information
            !Each quantity has real/imaginary components, and
            ! so the full size is twice this value.  THIS IS GLOBAL
            qsize = nlevels*nmodes*8

            !This is the size (bytes) of a single iteration's record
            rec_size = qsize*2*nq_shell+12  ! 12 is for the simtime+iteration at the end

            disp = hdisp+rec_size*(current_rec-1)

            !new_disp = disp+my_rdisp
            Do p = 1, 2
                new_disp = disp+my_rdisp +(p-1)*qsize*nq_shell
                !write(6,*)'new_disp: ', new_disp, my_column_rank
                Do i = 1, nq_shell
                         
          
                    Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                    
                    Call MPI_FILE_WRITE(funit, all_spectra(0,0,1,i,p), buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)

                    new_disp = new_disp+qsize
                Enddo
            Enddo
            disp = hdisp+rec_size*current_rec
            disp = disp-12
            Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)


            If (shell_spectra%master) Then

                buffsize = 1
                Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
                Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                       MPI_INTEGER, mstatus, ierr)
            Endif


			DeAllocate(all_spectra)
        Endif  ! Responsible

        If (my_row_rank .eq. 0) Call shell_spectra%closefile_par()

	End Subroutine Write_Shell_Spectra


	Subroutine Write_Shell_Spectra_MEM(this_iter,simtime)
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:,:), all_spectra(:,:,:,:,:)
        Real*8, Allocatable :: sendbuffer(:,:,:,:,:), out_radii(:)
        Real*8, Allocatable :: bsendbuffer(:,:,:,:,:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, m, mp, lmax,rind,field_ind,f,r
        Integer :: rone,  p,  counter, nf
		Integer :: n, nn, this_nshell, nq_shell, shell_spectra_tag, nmodes
        Integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp
        Integer(kind=MPI_OFFSET_KIND)  :: qsize, qdisp, rec_size
		Integer :: your_mp_min, your_mp_max, your_nm, your_id
		Integer :: nelem, m_ind, m_val, current_rec
        Integer :: funit, error, sirq, inds(5), dims(3)
        Integer :: my_nlevels, nlevels, qindex
        Integer :: lp1, nrirqs, ind5
        Integer :: ierr, rcount, buffsize
        Integer, Allocatable :: rirqs(:)
        Integer :: mstatus(MPI_STATUS_SIZE)       

        nlevels = Shell_Spectra%nlevels             ! The total number of spectra levels that needs to be output
        my_nlevels = Shell_Spectra%my_nlevels       ! The number of radial levels that this rank needs to write out
        nq_shell = Shell_Spectra%nq                 ! The number of quantities 
        shell_spectra_tag = Shell_Spectra%mpi_tag
        funit = Shell_Spectra%file_unit
        lmax = maxval(pfi%inds_3s)
        lp1 = lmax+1
        nmodes = lp1*lp1
		responsible = 0
		If ( (my_row_rank .eq. 0) .and. (my_nlevels .gt. 0) )  Then
            responsible = 1
            Allocate(all_spectra(0:lmax,0:lmax, my_nlevels,1, 1:2))
            Allocate(buff(0:lmax,my_nlevels,1,1:2,1:lp1))  !note - indexing starts at 1 not zero for mp_min etc.
            nrirqs = nproc2-1
            Allocate(rirqs(1:nrirqs))
        Endif


        !Before we start the main communication, all processes that contribute to the
        ! spectral output must get their buffers in the correct form
        If (my_nlevels .gt. 0) Then
            !//////////////////////
            ! First thing we do is FFT/reform the buffer/Legendre Transform
            !
            Call FFT_To_Spectral(spectra_buffer%p3b, rsc = .true.)
            spectra_buffer%config ='p3b'
            Call spectra_buffer%reform()
            Call spectra_buffer%construct('s2b')
            Call Legendre_Transform(spectra_buffer%p2b,spectra_buffer%s2b)
            Call spectra_buffer%deconstruct('p2b')

            Allocate(bsendbuffer(0:lmax,my_nlevels,nq_shell,2, my_mp_min:my_mp_max ))
            Allocate(sendbuffer(0:lmax,my_nlevels,1,2, my_mp_min:my_mp_max ))
            bsendbuffer = 0.0d0 
            sendbuffer = 0.0d0
            nf = spectra_buffer%nf2b
            Do p = 1, 2  ! Real and imaginary parts
            Do mp = my_mp_min,my_mp_max
                m = pfi%inds_3s(mp)
                    counter = 0
                    Do f = 1, nq_shell

                        field_ind = counter/my_nr+1
                        Do r = 1, shell_spectra%my_nlevels   
                                
                            rind = MOD(counter,my_nr)+my_rmin
                            bsendbuffer(m:lmax,r,f,p,mp) = &
                                & spectra_buffer%s2b(mp)%data(m:lmax,rind,p,field_ind)
                            counter = counter+1
                        Enddo
                    Enddo
                Enddo

            Enddo
            call spectra_buffer%deconstruct('s2b')

        Endif


        If (my_row_rank .eq. 0) Call Shell_Spectra%OpenFile_Par(this_iter, error)

        If (responsible .eq. 1) Then
            ! Processes that take part in the write have some extra work to do
            funit = shell_spectra%file_unit
            current_rec = Shell_Spectra%current_rec  ! Note that we have to do this after the file is opened
            If  ( (current_rec .eq. 1) .and. (shell_spectra%master) ) Then                
                !Write(6,*)'I am master: ', my_column_rank
                dims(1) =  lmax
                dims(2) =  nlevels
                dims(3) =  nq_shell
                buffsize = 3
                call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_shell
                call MPI_FILE_WRITE(funit,Shell_Spectra%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                allocate(out_radii(1:nlevels))
                Do i = 1, nlevels
                    out_radii(i) = radius(Shell_Spectra%levels(i))
                Enddo
                buffsize = nlevels
	            call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
                DeAllocate(out_radii)
                
	            call MPI_FILE_WRITE(funit, Shell_Spectra%levels, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

            Endif

            hdisp = 24 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_shell*4 ! nq
            hdisp = hdisp+nlevels*12  ! level indices and level values
            
            rcount = 0
            Do p = 1, Shell_Spectra%nshell_r_ids
                if (Shell_Spectra%shell_r_ids(p) .lt. my_column_rank) Then
                    rcount = rcount+ Shell_Spectra%nshells_at_rid(p)
                Endif
            Enddo
            my_rdisp = rcount*nmodes*8

                

            ! This is the LOCAL number ELEMENTS in the real or imaginary component of
            ! of a single quantity  (This is not in bytes)
            buffsize = my_nlevels*nmodes 

            !This is the half-size (bytes) of a single quantity's information
            !Each quantity has real/imaginary components, and
            ! so the full size is twice this value.  THIS IS GLOBAL
            qsize = nlevels*nmodes*8

            !This is the size (bytes) of a single iteration's record
            rec_size = qsize*2*nq_shell+12  ! 12 is for the simtime+iteration at the end

            disp = hdisp+rec_size*(current_rec-1)

        Endif


        Do qindex = 1, nq_shell  ! Q LOOP starts here!

            !Load the current quantity into the sendbuffer
            If (my_nlevels .gt. 0) Then
                sendbuffer(:,:,1,:,:) = & 
                    & bsendbuffer(:,:,qindex,:,:)
            Endif



            If (responsible .eq. 1) Then
                ! Rank 0 in reach row receives  all pieces of the shell spectra from the other nodes

                all_spectra(:,:,:,:,:) = 0.0d0
                buff(:,:,:,:,:) = 0.0d0


                rirqs(:) = 0
                ind5 = pfi%all_3s(0)%delta+1
                Do nn = 1, nrirqs
                    !Write(6,*)'Ind5: ', ind5
                    your_id = nn

                    your_nm     = pfi%all_3s(nn)%delta
                    your_mp_min = pfi%all_3s(nn)%min
                    your_mp_max = pfi%all_3s(nn)%max


                    nelem = your_nm*my_nlevels*2*lp1

                    inds(:) = 1
                    inds(5) = ind5  !This is the mp_index here.

                    Call Ireceive(buff, rirqs(nn), n_elements = nelem,source= your_id, &
                        &  tag=shell_spectra_tag,grp = pfi%rcomm, indstart = inds)
                    ind5 = ind5+your_nm
                Enddo

                ! Stripe my own data into the receive buffer

                Do mp = my_mp_min,  my_mp_max
                    m = pfi%inds_3s(mp)
                    Do p = 1,2
                        Do r = 1, my_nlevels   
                            buff(m:lmax,r,1,p,mp) = sendbuffer(m:lmax,r,1,p,mp) 
                        Enddo
                    Enddo
                Enddo

                Call IWaitAll(nrirqs,rirqs)

                !Stripe the receiver buffer into the spectra buffer
               
                Do mp = 1,lp1
                    m = pfi%inds_3s(mp)
                    Do p = 1, 2  ! Real and imaginary parts
                        Do r = 1, my_nlevels   
                            all_spectra(m:lmax,m,r,1,p) = buff(m:lmax,r,1,p,mp)  
                        Enddo
                    Enddo
                Enddo

                !Write the slice we just received
                Do p = 1, 2
                    new_disp = disp+my_rdisp +(p-1)*qsize*nq_shell +(qindex-1)*qsize        
                    Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                    
                    Call MPI_FILE_WRITE(funit, all_spectra(0,0,1,1,p), buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                Enddo

            Else
			    !  Non-responsible nodes send their info
			    If (my_nlevels .gt. 0) Then
                    inds(:) = 1
				    Call Isend(sendbuffer,sirq, dest = 0,tag=shell_spectra_tag, grp = pfi%rcomm, indstart = inds)
                    Call IWait(sirq)
			    Endif
		    Endif

        Enddo  ! Q-LOOP

        If (responsible .eq. 1) Then
            disp = hdisp+rec_size*current_rec
            disp = disp-12
            Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)

            If (shell_spectra%master) Then

                buffsize = 1
                Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
                Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                       MPI_INTEGER, mstatus, ierr)
            Endif

            DeAllocate(all_spectra)
            DeAllocate(buff)
            DeAllocate(rirqs)
        Endif


        If (my_row_rank .eq. 0) Call shell_spectra%closefile_par()
        If (my_nlevels .gt. 0) Then 
            DeAllocate(sendbuffer, bsendbuffer)
        Endif



	End Subroutine Write_Shell_Spectra_MEM




	Subroutine Write_Shell_Slices(this_iter,simtime)
        USE MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:), all_shell_slices(:,:,:,:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, p, sizecheck, t
		Integer :: n, nn, this_nshell, nq_shell, shell_slice_tag
		Integer :: your_theta_min, your_theta_max, your_ntheta, your_id
		Integer :: nelem, buffsize, sirq, nrirqs, inds(1:4)
        Integer :: file_pos, funit, error, dims(1:3), first_shell_rank
        Real*8, Allocatable :: out_radii(:)
        Integer, Allocatable :: level_inds(:), rirqs(:)
        
        integer :: ierr, rcount
		integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp, qdisp, full_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
        sizecheck = sizeof(disp)
        if (sizecheck .lt. 8) Then
            if (myid .eq. 0) Then
            Write(6,*)"Warning, MPI_OFFSET_KIND is less than 8 bytes on your system."
            Write(6,*)"Your size (in bytes) is: ", sizecheck
            Write(6,*)"A size of 4 bytes means that shell slices files are effectively limited to 2 GB in size."
            Endif
        endif 

        nq_shell = Shell_Slices%nq
        shell_slice_tag = Shell_Slices%mpi_tag
        funit = Shell_Slices%file_unit

		responsible = 0
		If (my_row_rank .eq. 0) Then
            If (Shell_Slices%my_nlevels .gt. 0) Then
                responsible = 1
            Endif
        Endif




        this_nshell = Shell_Slices%my_nlevels
		If (responsible .eq. 1) Then
			! Responsible node receives  all the pieces of the shell slices from the other nodes
			Allocate(all_shell_slices(1:nphi,1:ntheta,1:Shell_Slices%my_nlevels,nq_shell))

            Allocate(buff(1:nphi,1:this_nshell,1:nq_shell, 1:ntheta))            
			all_shell_slices(:,:,:,:) = 0.0d0
            buff(:,:,:,:) = 0.0d0

            ! Post Ireceives
            nrirqs = nproc2-1
            Allocate(rirqs(1:nrirqs))          

            Do nn = 1, nproc2-1
				your_id = nn

				your_ntheta    = pfi%all_2p(nn)%delta
				your_theta_min = pfi%all_2p(nn)%min

                inds(1) = 1
                inds(2) = 1
                inds(3) = 1
                inds(4) = your_theta_min

				nelem = nphi*your_ntheta*this_nshell*nq_shell

         		Call IReceive(buff, rirqs(nn),n_elements = nelem, source= your_id,tag=shell_slice_tag,grp = pfi%rcomm, indstart = inds)

            Enddo

            ! Stripe my own data into buff

            Do k = 1, nq_shell
                Do j = 1, this_nshell
                    Do t = my_theta_min, my_theta_max
                        Do i = 1, nphi
                            buff(i,j,k,t) = shell_slice_outputs(i,t,j,k)
                        Enddo
                    Enddo
                Enddo
            Enddo
            

            Call IWaitAll(nrirqs,rirqs)


            Do k = 1, nq_shell
                Do j = 1, this_nshell
                    Do t = 1, ntheta
                        Do i = 1, nphi
                            all_shell_slices(i,t,j,k) = buff(i,j,k,t)
                        Enddo
                    Enddo
                Enddo
            Enddo

            DeAllocate(buff)
            DeAllocate(rirqs)
            DeAllocate(shell_slice_outputs)
		Else
			!  Non responsible nodes send their info
			If (Shell_Slices%my_nlevels .gt. 0) Then
                !Everyone needs to restripe their data before sending it down the row
                !Stripe so that theta is slowest 
                Allocate(buff(1:nphi,1:this_nshell,1:nq_shell, my_theta_min:my_theta_max))
                Do t = my_theta_min, my_theta_max
                    Do k = 1, nq_shell
                        Do j = 1, this_nshell
                            Do i = 1, nphi
                                buff(i,j,k,t) = shell_slice_outputs(i,t,j,k)
                            Enddo
                        Enddo
                    Enddo
                Enddo
                nelem = nphi*my_ntheta*this_nshell*nq_shell
                inds(:) = 1
				Call Isend(buff,sirq,n_elements = nelem,dest = 0,tag=shell_slice_tag, grp = pfi%rcomm, indstart = inds)
            
                Call IWait(sirq)
                DeAllocate(shell_slice_outputs)
                DeAllocate(buff)
			Endif
		Endif


        ! Communication is complete.  Now we open the file using MPI-IO
        

        ! For the moment, every process in column 0 participates in the mpi operation
        ! The plan is to tune this later so that 
        if (my_row_rank .eq. 0) Call Shell_Slices%OpenFile_Par(this_iter, error)

        If (responsible .eq. 1) Then   
           funit = shell_slices%file_unit
        If (Shell_Slices%current_rec .eq. 1) Then                
            
            If (shell_slices%master) Then            
                ! The master rank (whoever owns the first output shell level) writes the header
                dims(1) = ntheta
                dims(2) = Shell_Slices%nlevels
                dims(3) =  nq_shell
                buffsize = 3
                call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_shell
                call MPI_FILE_WRITE(funit,Shell_Slices%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                allocate(out_radii(1:Shell_Slices%nlevels))
                Do i = 1, Shell_Slices%nlevels
                    out_radii(i) = radius(Shell_Slices%levels(i))
                Enddo
                buffsize = Shell_Slices%nlevels
	            call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
                DeAllocate(out_radii)
                

                allocate(level_inds(1:Shell_Slices%nlevels))
                Do i = 1, Shell_Slices%nlevels
                    level_inds(i) = Shell_Slices%levels(i)
                Enddo

	            call MPI_FILE_WRITE(funit, Shell_Slices%levels, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 
                DeAllocate(level_inds)
                buffsize = ntheta
	            call MPI_FILE_WRITE(funit, costheta, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 

            Endif
        Endif

            ! I might really (really) want to look into file views later.
            ! Depending on the offset size mpi type, disp will crap out past 2GB
            hdisp = 24 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_shell*4 ! nq
            hdisp = hdisp+Shell_Slices%nlevels*12  ! level indices and level values
            hdisp = hdisp+ ntheta*8  ! costheta

            qdisp = ntheta*Shell_Slices%nlevels*nphi*8
            full_disp = qdisp*nq_shell+12  ! 12 is for the simtime+iteration at the end
            disp = hdisp+full_disp*(Shell_Slices%current_rec-1)
            
            buffsize = Shell_Slices%my_nlevels*ntheta*nphi
            ! The file is striped with time step slowest, followed by q


            rcount = 0
            Do p = 1, Shell_Slices%nshell_r_ids
                if (Shell_Slices%shell_r_ids(p) .lt. my_column_rank) Then
                    rcount = rcount+ Shell_Slices%nshells_at_rid(p)
                Endif
            Enddo
            my_rdisp = rcount*ntheta*nphi*8
            Do i = 1, nq_shell
                new_disp = disp+qdisp*(i-1)+my_rdisp                
                Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                
                Call MPI_FILE_WRITE(funit, all_shell_slices(1,1,1,i), buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
            Enddo
            disp = hdisp+full_disp*Shell_Slices%current_rec
            disp = disp-12
            Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)


            If (shell_slices%master) Then
                buffsize = 1
                Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
                Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                       MPI_INTEGER, mstatus, ierr)
            Endif


			DeAllocate(all_shell_slices)
        Endif  ! Responsible

        If (my_row_rank .eq. 0) Call shell_slices%closefile_par()


	End Subroutine Write_Shell_Slices


	Subroutine Write_Shell_Slices_MEM(this_iter,simtime)
        ! A "more" memory friendly version of write_shell_Slices.  
        ! Writes one quantity at a time
        USE MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:), all_shell_slices(:,:,:,:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, p, sizecheck, t
		Integer :: n, nn, this_nshell, nq_shell, shell_slice_tag
		Integer :: your_theta_min, your_theta_max, your_ntheta, your_id
		Integer :: nelem, buffsize, sirq, nrirqs, inds(1:4),qbuffsize
        Integer :: file_pos, funit, error, dims(1:3), first_shell_rank
        Real*8, Allocatable :: out_radii(:)
        Integer, Allocatable :: level_inds(:), rirqs(:)
        
        integer :: ierr, rcount, qindex
		integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp, qdisp, full_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
        sizecheck = sizeof(disp)
        if (sizecheck .lt. 8) Then
            if (myid .eq. 0) Then
            Write(6,*)"Warning, MPI_OFFSET_KIND is less than 8 bytes on your system."
            Write(6,*)"Your size (in bytes) is: ", sizecheck
            Write(6,*)"A size of 4 bytes means that shell slices files are effectively limited to 2 GB in size."
            Endif
        endif 
        nq_shell = Shell_Slices%nq
        shell_slice_tag = Shell_Slices%mpi_tag
        

		responsible = 0
        this_nshell = Shell_Slices%my_nlevels
		If (my_row_rank .eq. 0) Then
            Call Shell_Slices%OpenFile_Par(this_iter, error)
            If (Shell_Slices%my_nlevels .gt. 0) Then
                responsible = 1
			    Allocate(all_shell_slices(1:nphi,1:ntheta,1:Shell_Slices%my_nlevels,1))
                Allocate(buff(1:nphi,1:this_nshell,1:1, 1:ntheta))  
                nrirqs = nproc2-1
                Allocate(rirqs(1:nrirqs))  

                ! Some displacements for accessing the file
                hdisp = 24 ! dimensions+endian+version+record count
                hdisp = hdisp+nq_shell*4 ! nq
                hdisp = hdisp+Shell_Slices%nlevels*12  ! level indices and level values
                hdisp = hdisp+ ntheta*8  ! costheta

                qdisp = ntheta*Shell_Slices%nlevels*nphi*8
                full_disp = qdisp*nq_shell+12  ! 12 is for the simtime+iteration at the end
                disp = hdisp+full_disp*(Shell_Slices%current_rec-1)
                ! The file is striped with time step slowest, followed by q
                rcount = 0
                Do p = 1, Shell_Slices%nshell_r_ids
                    if (Shell_Slices%shell_r_ids(p) .lt. my_column_rank) Then
                        rcount = rcount+ Shell_Slices%nshells_at_rid(p)
                    Endif
                Enddo
                my_rdisp = rcount*ntheta*nphi*8

                qbuffsize = Shell_Slices%my_nlevels*ntheta*nphi ! Number of elements in one q's worth of shells

            Endif
        Endif
        If (responsible .eq. 0) Then
            If (Shell_Slices%my_nlevels .gt. 0) Then
                Allocate(buff(1:nphi,1:this_nshell,1:1, my_theta_min:my_theta_max))
            Endif
        Endif
        funit = Shell_Slices%file_unit
        !////////////////////////////
        !Write a header
        If (Shell_Slices%current_rec .eq. 1) Then                
            
            If (shell_slices%master .and. (responsible .eq. 1)) Then            
                ! The master rank (whoever owns the first output shell level) writes the header
                dims(1) = ntheta
                dims(2) = Shell_Slices%nlevels
                dims(3) =  nq_shell
                buffsize = 3
                call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_shell
                call MPI_FILE_WRITE(funit,Shell_Slices%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                allocate(out_radii(1:Shell_Slices%nlevels))
                Do i = 1, Shell_Slices%nlevels
                    out_radii(i) = radius(Shell_Slices%levels(i))
                Enddo
                buffsize = Shell_Slices%nlevels
	            call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
                DeAllocate(out_radii)
                

                allocate(level_inds(1:Shell_Slices%nlevels))
                Do i = 1, Shell_Slices%nlevels
                    level_inds(i) = Shell_Slices%levels(i)
                Enddo

	            call MPI_FILE_WRITE(funit, Shell_Slices%levels, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 
                DeAllocate(level_inds)
                buffsize = ntheta
	            call MPI_FILE_WRITE(funit, costheta, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 

            Endif
        Endif


        !///////////////////////////////////////////////////////////////////////////
        ! In this revised version, we send one quantity's worth of shells at a time
        Do qindex = 1, nq_shell

		    If (responsible .eq. 1) Then
			    ! Responsible node receives  all the pieces of the shell slices from the other nodes
          
			    all_shell_slices(:,:,:,:) = 0.0d0
                buff(:,:,:,:) = 0.0d0
            
                Do nn = 1, nproc2-1
				    your_id = nn

				    your_ntheta    = pfi%all_2p(nn)%delta
				    your_theta_min = pfi%all_2p(nn)%min

                    inds(1) = 1
                    inds(2) = 1
                    inds(3) = 1
                    inds(4) = your_theta_min

				    nelem = nphi*your_ntheta*this_nshell

             		Call IReceive(buff, rirqs(nn),n_elements = nelem, source= your_id, &
                        & tag=shell_slice_tag,grp = pfi%rcomm, indstart = inds)

                Enddo

                ! Stripe my own data into the receive buffer
                Do j = 1, this_nshell
                    Do t = my_theta_min, my_theta_max
                        Do i = 1, nphi
                            buff(i,j,1,t) = shell_slice_outputs(i,t,j,qindex)
                        Enddo
                    Enddo
                Enddo
                
                Call IWaitAll(nrirqs,rirqs)

                !Re-organize the buffer
                Do j = 1, this_nshell
                    Do t = 1, ntheta
                        Do i = 1, nphi
                            all_shell_slices(i,t,j,1) = buff(i,j,1,t)
                        Enddo
                    Enddo
                Enddo


		    Else
			    !  Non responsible nodes send their info
			    If (Shell_Slices%my_nlevels .gt. 0) Then
                    !Everyone needs to restripe their data before sending it down the row
                    !Stripe so that theta is slowest 

                    Do t = my_theta_min, my_theta_max
                        Do j = 1, this_nshell
                            Do i = 1, nphi
                                buff(i,j,1,t) = shell_slice_outputs(i,t,j,qindex)
                            Enddo
                        Enddo
                    Enddo
                    nelem = nphi*my_ntheta*this_nshell
                    inds(:) = 1
				    Call Isend(buff,sirq,n_elements = nelem,dest = 0,tag=shell_slice_tag, &
                        & grp = pfi%rcomm, indstart = inds)
                
                    Call IWait(sirq)

			    Endif
		    Endif


            ! Communication is complete.  Write this q-value using MPI-IO

            If (responsible .eq. 1) Then   

                new_disp = disp+qdisp*(qindex-1)+my_rdisp                
                Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                
                Call MPI_FILE_WRITE(funit, all_shell_slices(1,1,1,1), qbuffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)

            Endif  ! Responsible
        Enddo

        If (responsible .eq. 1) Then
			DeAllocate(all_shell_slices)
            DeAllocate(rirqs)
            disp = hdisp+full_disp*Shell_Slices%current_rec
            disp = disp-12
            Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)


            If (shell_slices%master) Then
                buffsize = 1
                Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
                Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                       MPI_INTEGER, mstatus, ierr)
            Endif
        Endif
        If (my_row_rank .eq. 0) Call shell_slices%closefile_par()
        If (Shell_Slices%my_nlevels .gt. 0) Then
            DeAllocate(shell_slice_outputs)
            DeAllocate(buff)
        Endif
	End Subroutine Write_Shell_Slices_MEM



    Function Compute_Quantity(qval) result(yesno)
        integer, intent(in) :: qval 
        Logical             :: yesno 
        current_qval = qval
        yesno = .false.

        If (IOAvg_Flag .eq. 1) Then
            ! We check only the shell_averages so we can compute averages for moments
            If (Shell_Averages%compute(qval) .eq. 1) Then
                Call Shell_Averages%getq_now(yesno)
            Endif
        Else
            ! Otherwise, check everything - normal output
            If (compute_q(qval) .eq. 1) Then

                ! We check all diagnostic types and their respective frequencies
                ! While doing so, we modify the averaging level
                Call Full_3D%getq_now(yesno)
                Call Shell_Slices%getq_now(yesno)
                Call Shell_Spectra%getq_now(yesno)
                Call AZ_Averages%getq_now(yesno)
                Call Shell_Averages%getq_now(yesno)
                Call Global_Averages%getq_now(yesno)

            Endif
        Endif
    End function Compute_quantity

    Function Time_To_Output(iter) result(yesno)
        integer, intent(in) :: iter 
        Logical             :: yesno 
        yesno = .false.


         If (Mod(iter,Global_Averages%Frequency) .eq. 0) yesno = .true.
         If (Mod(iter,Shell_Averages%Frequency) .eq. 0) yesno = .true.        
         If (Mod(iter,Shell_Spectra%Frequency) .eq. 0) yesno = .true.            
         If (Mod(iter,AZ_Averages%Frequency) .eq. 0) yesno = .true.
         If (Mod(iter,Shell_Slices%Frequency) .eq. 0) yesno = .true. 
         If (Mod(iter,Full_3D%Frequency) .eq. 0) yesno = .true.
    End function Time_To_Output

    Subroutine Finalize_Averages()
        ! Use the azimuthal averages to compute the ell=0 mean
        Call IOComputeEll0(IOm0_values,IOell0_values)
        Call Shell_Averages%reset() !need to reset the index counter
    End Subroutine Finalize_Averages
    Subroutine Set_Avg_Flag(flag_val)
        Integer, Intent(In) :: flag_val
        IOavg_flag = flag_val
    End Subroutine Set_Avg_Flag
	Subroutine Add_Quantity(qty)
		Implicit None
		Real*8, Intent(In) :: qty(:,:,:)

        If (IOavg_flag .eq. 1) Then
            !Compute and store the azimuthal average of qty
            If (shell_averages%grab_this_q) Then
                Call IOComputeM0(qty)
            Endif
        Else

		    If (Shell_Slices%grab_this_q) Call get_shell_slice(qty)
		    If (Shell_Spectra%grab_this_q) Call get_shell_spectra(qty)
            Call Get_Averages(qty)

		    If (full_3d%grab_this_q) Call write_full_3d(qty)
        Endif
		



	End Subroutine	Add_Quantity

	Subroutine Complete_Output(iter, sim_time)
		Integer, Intent(In) :: iter
		Real*8, Intent(In) :: sim_time

	    If (Mod(iter,Shell_Slices%frequency) .eq. 0 ) Then
            If (mem_friendly) Then
                Call Write_Shell_Slices_MEM(iter,sim_time)
            Else
                Call Write_Shell_Slices(iter,sim_time)
            Endif
        Endif
	    If (Mod(iter,Shell_Spectra%frequency) .eq. 0 ) Then
            If (mem_friendly) Then
                Call Write_Shell_Spectra_MEM(iter,sim_time)
            else
                Call Write_Shell_Spectra(iter,sim_time)
            Endif
        Endif
	    If (Mod(iter,AZ_Averages%frequency) .eq. 0 ) Call Write_Azimuthal_Average(iter,sim_time)
	    If (Mod(iter,Shell_Averages%frequency) .eq. 0 ) Call Write_Shell_Average(iter,sim_time)
	    If (Mod(iter,Global_Averages%frequency) .eq. 0 ) Call Write_Global_Average(iter,sim_time)

        DeAllocate(f_of_r_theta)
        DeAllocate(f_of_r)
        DeAllocate(IOm0_values)
        DeAllocate(IOell0_values)
	End Subroutine Complete_Output

	Subroutine Get_Averages(qty)
        ! Takes azimuthal, shellular, and partial global averages of 3D array qty
		Implicit None
		Integer :: azav_ind,shellav_ind, globav_ind
        Integer :: i, r,t,p,m
        Real*8 :: this_average, wght
		Real*8, Intent(In) :: qty(1:,my_rmin:,my_theta_min:)

        !//////////////////////////////
        !First the azimuthal average
        If (current_averaging_level .ge. 1) Then

		    f_of_r_theta(:,:) = 0.0D0
		
		    Do t = my_theta_min, my_theta_max
			    Do r = my_rmin, my_rmax
				    f_of_r_theta(r,t) = sum(qty(:,r,t))
			    Enddo
		    Enddo

		    f_of_r_theta = f_of_r_theta*over_nphi_double     ! average in phi

		    If (AZ_Averages%grab_this_q) Then
                If (.not. Allocated(azav_outputs)) Then
			        Allocate(azav_outputs(my_rmin:my_rmax,my_theta_min:my_theta_max,1:AZ_Averages%nq))
                Endif
    		    azav_ind = AZ_Averages%ind
			    azav_outputs(:,:,azav_ind) = f_of_r_theta
                If (myid .eq. 0) AZ_Averages%oqvals(azav_ind) = current_qval
			    Call AZ_Averages%AdvanceInd()
		    Endif
        Endif

        !/////////////////////
        !Now Average over the partial sphere (shell_averages)
		If (current_averaging_level .ge. 2) Then
		    f_of_r(:) = 0.0D0

            Do t = my_theta_min, my_theta_max
                f_of_r(:) = f_of_r(:) + f_of_r_theta(:,t)*theta_integration_weights(t)
            Enddo


		    If (Shell_Averages%grab_this_q) Then

                If (.not. Allocated(shellav_outputs)) Then
			        Allocate(shellav_outputs(my_rmin:my_rmax,1:4,1:Shell_Averages%nq))  ! four moments
                    shellav_outputs(:,:,:) = 0.0d0
                Endif

                shellav_ind = Shell_Averages%ind

                !First, add the partially integrated spherically symmetric mean ( f_of_r )
                Do r = my_rmin, my_rmax
                    shellav_outputs(r,1,shellav_ind) = f_of_r(r)
                Enddo

                ! Now, the moments - rms, skewness, kurtosis
                Do m = 2, 4
                    Do t = my_theta_min, my_theta_max
                        wght = theta_integration_weights(t)*over_nphi_double
                        Do r = my_rmin, my_rmax
                            Do p = 1, nphi
                            shellav_outputs(r,m,shellav_ind) = shellav_outputs(r,m,shellav_ind) + &
                                & wght*(qty(p,r,t)-IOell0_values(r,shellav_ind))**m
                            Enddo
                            
                        Enddo
                    Enddo
                Enddo

                If (myid .eq. 0) Shell_Averages%oqvals(shellav_ind) = current_qval
                Call Shell_Averages%AdvanceInd()
		    Endif
        Endif

        !////////////////////////
        ! Finally, integrate in radius to get the partial global average
        If (current_averaging_level .ge. 3) Then

            this_average =0.0d0
            do i = my_rmin, my_rmax
                this_average = this_average+f_of_r(i)*r_integration_weights(i)
            enddo

		    If (Global_Averages%grab_this_q) Then
    		    If (.not. Allocated(globav_outputs)) Then			
    			    Allocate(globav_outputs(1:Global_Averages%nq))			
    		    Endif
                globav_ind = Global_Averages%ind
			    globav_outputs(globav_ind) = this_average
                If (myid .eq. 0) Global_Averages%oqvals(globav_ind) = current_qval
                Call Global_Averages%AdvanceInd()
		    Endif

        Endif

	END Subroutine Get_Averages





	Subroutine Write_Azimuthal_Average(this_iter,simtime)
        USE MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:), all_azavgs(:,:,:)
		Integer :: responsible, current_rec, s_start, s_end, this_rid
		Integer :: i, j, k,qq, p, sizecheck
		Integer :: n, nn, this_nshell, nq_azav, az_avg_tag
		Integer :: your_theta_min, your_theta_max, your_ntheta
		Integer :: nelem, buffsize
        Integer :: file_pos, funit, error, dims(1:3)
        Integer :: inds(3), nirq,sirq
        Integer, Allocatable :: rirqs(:)
        
        integer :: ierr, rcount
		integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp, qdisp, full_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
        sizecheck = sizeof(disp)
        if (sizecheck .lt. 8) Then
            if (myid .eq. 0) Then
            Write(6,*)"Warning, MPI_OFFSET_KIND is less than 8 bytes on your system."
            Write(6,*)"Your size (in bytes) is: ", sizecheck
            Write(6,*)"A size of 4 bytes means that shell slices files are effectively limited to 2 GB in size."
            Endif
        endif 


		responsible = 0
        nq_azav     = AZ_Averages%nq
        az_avg_tag  = AZ_Averages%mpi_tag
        funit       = AZ_Averages%file_unit
        

		If (my_row_rank .eq. 0) Then
                responsible = 1
        Endif

        !Everyone needs to restripe their data for the sends so that theta is slowest
        Allocate(buff(my_rmin:my_rmax,1:nq_azav, my_theta_min:my_theta_max))
        Do k = my_theta_min, my_theta_max
            Do j = 1, nq_azav
                Do i = my_rmin, my_rmax
                    buff(i,j,k) = azav_outputs(i,k,j)  ! restripe
                Enddo
            Enddo
        Enddo
        DeAllocate(azav_outputs)


		If (responsible .eq. 1) Then
			! Rank 0 in reach row receives from all other row members
		
			Allocate(all_azavgs(my_rmin:my_rmax,1:nq_azav, 1:ntheta))
			all_azavgs(:,:,:) = 0.0d0

            nirq = nproc2-1
            Allocate(rirqs(1:nirq))

            Do nn = 1, nproc2-1

				your_ntheta    = pfi%all_2p(nn)%delta
				your_theta_min = pfi%all_2p(nn)%min
				your_theta_max = pfi%all_2p(nn)%max
                inds(1) = 1
                inds(2) = 1
                inds(3) = your_theta_min
                nelem = your_ntheta*my_nr*nq_azav

                Call IReceive(all_azavgs, rirqs(nn),n_elements = nelem, &
                            &  source= nn,tag = az_avg_tag, grp = pfi%rcomm,indstart = inds)				
			Enddo

            all_azavgs(my_rmin:my_rmax,1:nq_azav,my_theta_min:my_theta_max) = &
                & buff(my_rmin:my_rmax,1:nq_azav,my_theta_min:my_theta_max)

            Call IWaitAll(nirq, rirqs)
            DeAllocate(rirqs)
		Else
			!  Rest of the row sends to process 0 within the row
            inds(1) = 1 !my_rmin
            inds(2) = 1
            inds(3) = 1 ! my_theta_min
            nelem = my_nr*my_ntheta*nq_azav
            Call ISend(buff, sirq,n_elements = nelem, dest = 0, tag = az_avg_tag, & 
                grp = pfi%rcomm, indstart = inds)
            Call IWait(sirq)
		Endif
        DeAllocate(buff)

        ! Communication is complete.  Now we open the file using MPI-IO
        

      

        If (responsible .eq. 1) Then   
            Call AZ_Averages%OpenFile_Par(this_iter, error)
            current_rec = AZ_Averages%current_rec
            funit = AZ_Averages%file_unit
            !before we do anything else, we need to restripe the data yet again (might be able to work around this later)

            Allocate(buff( 1:ntheta,my_rmin:my_rmax, 1:nq_azav))
            Do k = 1, ntheta
                Do j = 1, nq_azav
                    Do i = my_rmin, my_rmax    
                        buff(k,i,j) = all_azavgs(i,j,k)
                    Enddo
                Enddo
            Enddo        
            DeAllocate(all_azavgs)
            

            If ((my_column_rank .eq. 0) .and. (current_rec .eq. 1) ) Then            
                ! Rank 0 in column and row writes the header
                dims(1) =  nr
                dims(2) =  ntheta
                dims(3) =  nq_azav
                buffsize = 3
                call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_azav
                call MPI_FILE_WRITE(funit,AZ_Averages%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nr
                call MPI_FILE_WRITE(funit, radius, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 

                buffsize = ntheta
                call MPI_FILE_WRITE(funit, costheta, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
            Endif


            hdisp = 24 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_azav*4 ! nq
            hdisp = hdisp+nr*8  ! The radius array
            hdisp = hdisp+ ntheta*8  ! costheta

            qdisp = ntheta*nr*8
            full_disp = qdisp*nq_azav+12  ! 12 is for the simtime+iteration at the end
            disp = hdisp+full_disp*(current_rec-1)
            
            buffsize = my_nr*ntheta
            ! The file is striped with time step slowest, followed by q

            my_rdisp = (my_rmin-1)*ntheta*8

            Do i = 1, nq_azav
                new_disp = disp+qdisp*(i-1)+my_rdisp                
                Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                
                Call MPI_FILE_WRITE(funit, buff(1,my_rmin,i), buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
            Enddo

            disp = hdisp+full_disp*current_rec
            disp = disp-12
            Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)


            If (my_column_rank .eq. 0) Then
                buffsize = 1
                Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
                Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                       MPI_INTEGER, mstatus, ierr)
            Endif

			DeAllocate(buff)
            Call AZ_Averages%closefile_par()
        Endif  ! Responsible

	End Subroutine Write_Azimuthal_Average


    Subroutine Write_Global_Average(this_iter,simtime)
        Implicit None
        Real*8, Allocatable :: buff(:), full_avg(:)

        Integer :: responsible
        Integer :: i,n, nq_globav, global_avg_tag, error, file_pos
        Integer :: funit
        Integer, Intent(In) :: this_iter
        Real*8, Intent(In) :: simtime

        nq_globav = Global_Averages%nq
        global_avg_tag = Global_Averages%mpi_tag
        funit = Global_Averages%file_unit


        !///////////////////////////////
        ! Sum across rows, and then across the first column

        Allocate(full_avg(nq_globav))
        Allocate(buff(nq_globav))

        full_avg(:) = 0.0d0
        buff(:) = 0.0d0
        Call dsum1d(globav_outputs,buff,pfi%rcomm)
        If (my_row_rank .eq. 0) Then
            Call dsum1d(buff,full_avg,pfi%ccomm)
        Endif

        !//////////////////////////////
       

        If (myid .eq. io_node) Then

            Call Global_Averages%OpenFile(this_iter, error)
            If (error .eq. 0) Then
                If (Global_Averages%current_rec .eq. 1) Then
                    Write(funit)nq_globav
                    Write(funit)(Global_Averages%oqvals(i),i=1,nq_globav)
                    Call Global_Averages%update_position
                Endif
                file_pos = Global_Averages%file_position

                Write(funit)(full_avg(i),i=1,nq_globav)
                Write(funit)simtime
                Write(funit)this_iter
                Call Global_Averages%CloseFile 
            Endif

        Endif

        DeAllocate(globav_outputs)
        DeAllocate(buff)
        DeAllocate(full_avg)

    End Subroutine Write_Global_Average


	Subroutine Write_Shell_Average(this_iter, simtime)
		Implicit None
        Integer, Intent(In) :: this_iter
		Real*8, Intent(In) :: simtime
		Integer :: i,j, k, n, nn, nq_shellav, shell_avg_tag, m 
        Real*8, Allocatable :: full_shellavg(:,:,:), buff(:,:,:), buff2(:,:,:)		
        Integer :: your_r_min, your_r_max, your_nr, your_id
        Integer :: funit, error, inds(3), ncount, sirq,nirq
        Integer, Allocatable :: rirqs(:)

        shell_avg_tag = Shell_Averages%mpi_tag
        nq_shellav = Shell_Averages%nq
        funit = Shell_Averages%file_unit

        !Sum across the row to complete integration in theta
        Allocate(buff(my_rmin:my_rmax,1:4,nq_shellav))
        buff(:,:,:) = 0.0d0
        Call dsum3d(shellav_outputs,buff,pfi%rcomm)

        If (my_row_rank .eq. 0) Then
            !now set up a series of isends/ireceives along the column
            If (my_column_rank .eq. 0) Then
                !Post ireceives before anything else is done
                nirq = nproc1-1
                Allocate(rirqs(1:nirq))
                Allocate(full_shellavg(1:nq_shellav,1:4,1:nr))
                Do n = 1, nproc1-1
                    your_r_min = pfi%all_1p(n)%min
                    your_nr = pfi%all_1p(n)%delta
                    ncount = your_nr*nq_shellav*4
                    inds(1) = 1
                    inds(2) = 1
                    inds(3) = your_r_min
                    Call IReceive(full_shellavg, rirqs(n),n_elements = ncount, &
                            &  source= n,tag = shell_avg_tag, grp = pfi%ccomm,indstart = inds)
                Enddo
                ! Load my buff into the full_shellavg array
                Do j = my_rmin,my_rmax
                    Do m = 1, 4
                        Do i = 1, nq_shellav
                            full_shellavg(i,m,j) = buff(j,m,i)
                        Enddo
                    Enddo
                Enddo
                DeAllocate(buff)
                Call IWaitAll(nirq, rirqs)
                DeAllocate(rirqs)
            Else
                Allocate(buff2(nq_shellav, 1:4,my_rmin:my_rmax)) ! Transpose for easier send logic
                Do j = my_rmin,my_rmax
                    Do m = 1, 4
                        Do i = 1, nq_shellav
                            buff2(i,m,j) = buff(j,m,i)
                        Enddo
                    Enddo
                Enddo
                ncount = my_nr*nq_shellav*4
                Call ISend(buff2, sirq,ncount, dest = 0, tag = shell_avg_tag, grp = pfi%ccomm)

                Call IWait(sirq)
                DeAllocate(buff2)
                DeAllocate(buff)
            Endif
        Endif

        If (myid .eq. io_node) Then ! Rank Zero writes the file
            Call Shell_Averages%OpenFile(this_iter, error)
            If (error .eq. 0) Then
                If (Shell_Averages%current_rec .eq. 1) Then
                    Write(funit)nr,nq_shellav
                    Write(funit)(Shell_Averages%oqvals(i),i=1,nq_shellav)
                    Write(funit)(radius(i),i=1,nr)
                    Call Shell_Averages%update_position()
                Endif
                    Write(funit)(((full_shellavg(k,m,i),i=1,nr),m=1,4),k=1,nq_shellav)
                Write(funit) simtime

                Write(funit)this_iter


                Call Shell_Averages%CloseFile
                DeAllocate(full_shellavg)
            Endif
		Endif

        DeAllocate(shellav_outputs)

	End Subroutine Write_Shell_Average

	Subroutine Write_Full_3D(qty)
		Use MPI_BASE ! Doing this here for now.  No other routine above sees MPI_Base, and I may want to keep it that way.
		Implicit None		
		Real*8, Intent(In) :: qty(:,my_rmin:,my_theta_min:)
		Real*8, Allocatable :: my_shells(:,:,:), buff(:,:,:)
		Integer :: i, j
		Character*2 :: qstring
		Character*8 :: iterstring
		Character*120 :: cfile

		Integer :: your_theta_min, your_theta_max, your_ntheta
		Integer :: np, buffsize, p
		Integer(kind=MPI_OFFSET_KIND) :: my_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
		Integer :: funit, ierr, full_3d_tag

		! qty is dimensioned 1:n_phi, my_rmin:my_rmax, my_theta_min:my_theta_max

        full_3d_tag = Full_3D%mpi_tag
		If (my_row_rank .eq. 0) Then
			! Everyone in the row communicates to row-rank zero.
			! Each row owns a specific rank of radii
			Allocate(my_shells(1:nphi, 1:ntheta, my_rmin:my_rmax))
			my_shells(:,:,:) = 0.0d0					

			! First each rank stripes its own data into the new array
			! Not that striping is not in the "native" order
			Do j = my_theta_min, my_theta_max
				Do i = my_rmin, my_rmax
					my_shells(:,j,i) = qty(:,i,j)
				Enddo
			Enddo
			np = pfi%rcomm%np
			Do p = 1, np-1	
				your_theta_min = pfi%all_2p(p)%min
				your_theta_max = pfi%all_2p(p)%max
				your_ntheta    = pfi%all_2p(p)%delta
				Allocate(buff(1:nphi,my_rmin:my_rmax, your_theta_min:your_theta_max))
				Call receive(buff, source= p,tag=full_3d_tag,grp = pfi%rcomm)
				Do j = your_theta_min, your_theta_max
					Do i = my_rmin, my_rmax
						my_shells(:,j,i) = buff(:,i,j)
					Enddo
				Enddo
				DeAllocate(buff)
			Enddo
			! Now do the MPI write
			np = pfi%ccomm%np
			my_disp = 0
			Do p = 1, my_column_rank
				my_disp = my_disp+pfi%all_1p(p-1)%delta
			Enddo
			my_disp = my_disp*ntheta*nphi*8	! Displacment of this rank in the MPI File
			buffsize = my_nr*nphi*ntheta		! Number of elements to write


			write(iterstring,i_ofmt) current_iteration
			write(qstring,'(i2.2)') current_qval
         cfile = trim(local_file_path)//'Spherical_3D/'//trim(iterstring)//'_'//qstring
			!Write(6,*)cfile
			call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, & 
                   MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                   MPI_INFO_NULL, funit, ierr) 

			call MPI_FILE_SET_VIEW(funit, my_disp, MPI_DOUBLE_PRECISION, & 
                   MPI_DOUBLE_PRECISION, 'native', & 
                   MPI_INFO_NULL, ierr) 
			call MPI_FILE_WRITE(funit, my_shells, buffsize, MPI_DOUBLE_PRECISION, & 
                   mstatus, ierr) 

			call MPI_FILE_CLOSE(funit, ierr) 
			If (my_column_rank .eq. 0) Then
				! row/column 0 writes out a file with the grid, etc.
				! This file should contain everything that needs to be known for processing later
	         write(iterstring,'(i8.8)') current_iteration
            cfile = trim(local_file_path)//'Spherical_3D/'//trim(iterstring)//'_'//'grid'
	         open(unit=15,file=cfile,form='unformatted', status='replace')
	         Write(15)nr
				Write(15)ntheta
				Write(15)nphi
	         Write(15)(radius(i),i=1,nr)
				Write(15)(acos(costheta(i)),i = 1, ntheta)
	         Close(15)
			Endif


		Else
			! Send an array that's indexed starting at 1.  Shouldn't be necessary, but just in case.
			!Allocate(buff(1:nphi,1:my_ntheta,1:my_nr))
			
			!Do j = my_theta_min, my_theta_max
			!	jind = j-my_theta_min+1
			!	Do i = my_rmin, my_rmax
			!		buff(:,jind,i) = qty(:,i,j)
			!	Enddo
			!Enddo
			Call send(qty, dest= 0,tag=full_3d_tag,grp = pfi%rcomm)
			!DeAllocate(buff)

		Endif	


	End Subroutine Write_Full_3D


    !////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !       Diagnostic Class Methods

    Subroutine Initialize_Diagnostic_Info(self,avg_levels,computes,pid,mpi_tag,avg_level,values, levels)
        Implicit None
        Integer :: i,ind
        Integer, Intent(In) :: pid, mpi_tag
        Integer, Optional, Intent(In) :: avg_level
        Integer, Optional, Intent(In) :: values(1:)
        Integer, Optional, Intent(In) :: levels(1:)
        Integer, Intent(InOut) :: computes(1:), avg_levels(1:)
        Class(DiagnosticInfo) :: self 
        If (present(avg_level)) Then
            if (avg_level .gt. 0) self%avg_level = avg_level
        Endif
        self%mpi_tag = mpi_tag
        self%nq = 0
        If (present(values)) Then
            self%values(:) = values(:)  ! This is clunky - will look into getting the object attributes directly into a namelist later
            
            Do i = 1, nqmax
                if(self%values(i) .gt. 0) Then 
                    self%nq = self%nq+1
                    ind = self%values(i)
                    self%compute(ind) = 1
                    computes(ind) = 1
                    If (avg_levels(ind) .lt. self%avg_level) Then
                        avg_levels(ind) = self%avg_level
                    Endif
                endif 
            Enddo
        Endif
        self%my_nlevels = 0
        If (present(levels)) Then
            self%levels(:) = levels(:)
            Do i = 1, nshellmax
                if( (self%levels(i) .gt. 0) .and. (self%levels(i) .le. nr) ) Then
                    self%nlevels = self%nlevels+1
                Endif
            Enddo
        Endif

        !if (pid .eq. 0) Then
            !NOTE:  Later, we may want to do this only on the master node later, but this isn't a huge memory issue
            Allocate(self%oqvals(1:self%nq))
            self%oqvals(:) = nqmax+100
        !Endif
    End Subroutine Initialize_Diagnostic_Info
    Subroutine AdvanceInd(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        self%ind = self%ind+1
        self%begin_output = .false.
    End Subroutine AdvanceInd
    Subroutine Diagnostic_Output_Reset(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        self%ind = 1
        self%begin_output = .true.
    End Subroutine Diagnostic_Output_Reset

    Subroutine Shell_Balance(self)
        ! I'm being a little sloppy here.  This method of the diagnostic class still uses
        ! a number of module-wide variables.
        ! Eventually, the code might be cleaner if this class were moved to its own module.
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer :: j, ilocal, i, pcount, your_rmax, your_rmin
        Integer, Allocatable :: ptemp(:), ptemp2(:), ptemp3(:)
		self%my_nlevels = 0
		Allocate(self%my_shell_levs(1:self%nlevels))
		Allocate(self%have_shell(1:self%nlevels))
		Allocate(self%my_shell_ind(1:self%nlevels))
		self%have_shell(:) = 0		
		Do j = 1, self%nlevels
			ilocal = self%levels(j)
			If ((ilocal .ge. my_rmin) .and. (ilocal .le. my_rmax)) Then ! my processor has this radius
			   self%have_shell(j) = 1
			   self%my_nlevels = self%my_nlevels+1
			   self%my_shell_levs(self%my_nlevels) = self%levels(j)
			   self%my_shell_ind(self%my_nlevels) = j
			Endif
		Enddo

        !/// ID 0 has a little more work to do
		If (my_row_rank .eq. 0) Then
			!Use process zero to figure out which radial processors will
			! actually output shells.  This will make it easier to read in the input later.
			Allocate(ptemp(1:nproc1))
			Allocate(ptemp2(1:nproc1))
			Allocate(ptemp3(1:nproc1))
			ptemp(:) = 0
			ptemp2(:) = 0
			do i = 0, nproc1-1
				your_rmin = pfi%all_1p(i)%min
				your_rmax = pfi%all_1p(i)%max
				Do j = 1, self%nlevels
					ilocal = self%levels(j)
					If ( (ilocal .ge. your_rmin) .and. (ilocal .le. your_rmax) ) Then
						ptemp(i+1) = ptemp(i+1)+1	!  radial id "i" has another shell that we want to output
					Endif	
				Enddo
			enddo
			pcount =0
			do i = 0, nproc1-1
				if (ptemp(i+1) .ge. 1) then 
					pcount = pcount+1					! pcount is the number of unique radial id's that have shells
					ptemp2(pcount) = i				! ptemp2 is for storing the radial id's of that have shells to output
					ptemp3(pcount) = ptemp(i+1)	! ptemp3 is the number of shells this radial id has
				endif
			enddo
            
			!   Resize the temporary arrays
			Allocate(self%shell_r_ids(1:pcount))
			Allocate(self%nshells_at_rid(1:pcount))
			self%shell_r_ids(:) = ptemp2(1:pcount)	! These are the radial ids of the processors that have shells
			self%nshells_at_rid(:) = ptemp3(1:pcount) ! How many shells this rid has
			self%nshell_r_ids = pcount
			DeAllocate(ptemp3)
			DeAllocate(ptemp2)
			DeAllocate(ptemp)
		Endif

    End Subroutine Shell_Balance

    Subroutine Init_OComm(self,pcomm,pnp,prank,mrank)
        Implicit None
        Integer, Intent(In) :: pcomm, pnp, prank, mrank
        Class(DiagnosticInfo) :: self
        self%ocomm = pcomm
        self%orank = prank
        self%onp = pnp
        If (prank .eq. mrank) then
            self%master = .true.    ! This process handles file headers in parallel IO
        Endif
    End Subroutine Init_OComm

    Subroutine set_file_info(self,oversion,rcount, freq, fpref,funit)
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer :: oversion, rcount, freq
        Integer, Optional, Intent(In) :: funit
        Character*120 :: fpref
        If (present(funit)) Then
            self%file_unit = funit
        Endif
        self%file_prefix = fpref
        self%output_version = oversion
        self%rec_per_file = rcount
        self%frequency = freq
    End Subroutine set_file_info

    Subroutine OpenFile(self,iter,errcheck)
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer, Intent(In) :: iter
        Integer, Intent(InOut) :: errcheck
        Character*8 :: iterstring,istr
        Character*120 :: filename, omsg
        Integer :: modcheck, imod, file_iter, next_iter, ibelong
        ! Note - we should do something to make sure that the file has been started before we start writing to it...
        ! possibly look at self%rec_count*self%frequency
        ! otherwise, if someone output with a strange cadence relative to
        modcheck = self%frequency*self%rec_per_file
        imod = Mod(iter,modcheck) 

        if (imod .eq. 0) then
            ibelong = iter
        else
            ibelong = iter-imod+modcheck    ! This iteration belongs in a file with number= ibelong
        endif
        write(iterstring,i_ofmt) ibelong
        filename = trim(local_file_path)//trim(self%file_prefix)//trim(iterstring)

        If ( (imod .eq. self%frequency) .or. (self%rec_per_file .eq. 1) ) Then   ! time to begin a new file 

            !omsg = ' Creating Filename: '//trim(filename)
            !Write(6,*)'Creating Filename: ', filename
            !Call stdout%print(omsg)
            Call stdout%print(' Creating Filename: '//trim(filename))
            Open(unit=self%file_unit,file=filename,form='unformatted', status='replace',access='stream',iostat = errcheck)
            Write(self%file_unit)endian_tag
            Write(self%file_unit)self%output_version
            Write(self%file_unit)integer_zero ! We write zero initially - only update nrec after the data is actually written
            self%current_rec = 1            
            If (errcheck .ne. 0) Then
                next_iter =file_iter+modcheck
                call stdout%print(' Unable to create file!!: '//trim(filename))
            Endif
        Else
            Open(unit=self%file_unit,file=filename,form='unformatted', status='old',access='stream', &
                & iostat = errcheck, POSITION = 'APPEND')    

            ! This looks redundant, but it allows partial files to be continued following restart
            !Read(self%file_unit,POS = 9)self%current_rec  
            self%current_rec = self%current_rec+1
            If (errcheck .ne. 0) Then
                next_iter =file_iter+modcheck
                Call stdout%print(' --Failed to find needed file: '//trim(filename))
                Call stdout%print(' --Partial diagnostic files are not currently supported.')
                Write(istr,'(i8.8)')ibelong+self%frequency
               !Write(6,*)'No data will be written until a new file is created at iteration: ', ibelong+self%frequency
                Call stdout%print(' --No data will be written until a new file is created at iteration: '//trim(istr))
            Endif
        Endif

    End Subroutine OpenFile


    Subroutine OpenFile_Par(self,iter,ierr)
        !Performs the same tasks as OpenFile, but uses MPI-IO
        ! Opens file, advances record count, writes header etc.
        Use MPI_BASE
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer, Intent(In) :: iter
        Integer, Intent(InOut) :: ierr
        Character*8 :: iterstring
        Character*120 :: filename
        Integer :: modcheck, imod, file_iter, next_iter, ibelong
        Integer :: buffsize, funit
        Integer :: mstatus(MPI_STATUS_SIZE)
        integer(kind=MPI_OFFSET_KIND) :: disp


        modcheck = self%frequency*self%rec_per_file
        imod = Mod(iter,modcheck) 

        if (imod .eq. 0) then
            ibelong = iter
        else
            ibelong = iter-imod+modcheck    ! This iteration belongs in a file with number= ibelong
        endif

        write(iterstring,i_ofmt) ibelong
        filename = trim(local_file_path)//trim(self%file_prefix)//trim(iterstring)

        If ( (imod .eq. self%frequency) .or. (self%rec_per_file .eq. 1) ) Then   ! time to begin a new file 

            

    		call MPI_FILE_OPEN(self%ocomm, filename, & 
                 MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                 MPI_INFO_NULL, funit, ierr) 
            self%file_unit = funit
            If (self%master) Then     
                Write(6,*)'Creating Filename: ', filename      
                buffsize = 1
                call MPI_FILE_WRITE(self%file_unit, endian_tag, buffsize, MPI_INTEGER, & 
                    mstatus, ierr)

                buffsize = 1
                call MPI_FILE_WRITE(self%file_unit, self%Output_Version, buffsize, MPI_INTEGER, & 
                    mstatus, ierr)
                buffsize = 1
                call MPI_FILE_WRITE(self%file_unit, integer_zero, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) ! We write zero initially - only update nrec after the data is actually written
            Endif

            
            self%current_rec = 1            
            If (ierr .ne. 0) Then
                next_iter =file_iter+modcheck
                Write(6,*)'Unable to create file!!: ',filename
            Endif
        Else
    		call MPI_FILE_OPEN(self%ocomm, filename, & 
                 MPI_MODE_RDWR, & 
                 MPI_INFO_NULL, funit, ierr) 
            self%file_unit = funit

            !disp = 8
            !Call MPI_File_Seek(self%file_unit,disp,MPI_SEEK_SET,ierr)
            !Read the current record

            !call MPI_FILE_READ(self%file_unit, self%current_rec, 1, MPI_INTEGER, & 
            !mstatus, ierr)

            self%current_rec = self%current_rec+1
            If (ierr .ne. 0) Then
                next_iter =file_iter+modcheck
                Write(6,*)'Failed to find needed file: ', filename
                Write(6,*)'Partial diagnostic files are not currently supported.'
                Write(6,*)'No data will be written until a new file is created at iteration: ', ibelong+self%frequency
            Endif
        Endif

    End Subroutine OpenFile_Par




    Subroutine CloseFile(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        Call self%update_position()
        Write(self%file_unit,POS = 9)self%current_rec
        Close(self%file_unit)
    End Subroutine CloseFile

    Subroutine CloseFile_Par(self)
        USE MPI_BASE
        Implicit None
        integer :: ierr, buffsize
        Integer :: mstatus(MPI_STATUS_SIZE)
        integer(kind=MPI_OFFSET_KIND) :: disp
        !Parallel File Close
        !Peforms the same task as closefile, but using MPI-IO
        Class(DiagnosticInfo) :: self
        disp = 8
        Call MPI_File_Seek(self%file_unit,disp,MPI_SEEK_SET,ierr)
            If (ierr .ne. 0) Then
                Write(6,*)'Error rewinding to header.  Error code: ', ierr, myid
            Endif
        If (self%master) Then  

            buffsize = 1

            call MPI_FILE_WRITE(self%file_unit,self%current_rec , buffsize, MPI_INTEGER, & 
                mstatus, ierr) 
            If (ierr .ne. 0) Write(6,*)'Error writing to header.  Error code: ', ierr
        Endif
        Call MPI_FILE_CLOSE(self%file_unit, ierr)
        If (ierr .ne. 0) Write(6,*)'Error closing file.  Error code: ',ierr
    End Subroutine CloseFile_Par



    Subroutine Update_Position(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        INQUIRE(UNIT=self%file_unit, POS=self%file_position)
    End Subroutine Update_Position
    Subroutine getq_now(self,yesno)
        Implicit None
        Logical, Intent(InOut) :: yesno
        Integer :: modcheck
        Class(DiagnosticInfo) :: self
        self%grab_this_q = .false.
        If(self%compute(current_qval) .eq. 1) Then
            modcheck = Mod(current_iteration,self%frequency)
            If (modcheck .eq. 0) Then
                yesno = .true.
                current_averaging_level = Max(current_averaging_level,self%avg_level)
                self%grab_this_q = .true.
            Endif
        Endif
    End Subroutine getq_now

    !/////////////////////////////////////////
    ! (Presently) Random Utility routines
    Subroutine ComputeEll0(inbuff,outbuff)
        Real*8, Intent(In) :: inbuff(1:,my_rmin:,my_theta_min:,1:)
        Real*8, Intent(InOut) :: outbuff(my_rmin:,1:)
        Real*8, Allocatable :: tmp_buffer(:,:)
        Integer :: bdims(1:4)
        Integer :: q,nq,r,t,p
        !Averages over theta and phi to get the spherically symmetric mean of all
        ! fields in inbuff at each radii
        ! inbuff is expected to be dimensioned as (1:nphi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,1:nfields)
        ! ** Note that this routine should be used sparingly because it requires 
        ! **   a collective operation (allreduce) across process rows
        ! ** One way to do this is to aggregate several fields into the inbuff when calling this routine


        bdims = shape(inbuff)
        nq = bdims(4)

        Allocate(tmp_buffer(my_rmin:my_rmax,1:nq))
        tmp_buffer(:,:) = 0.0d0
        outbuff(:,:) = 0.0d0

        ! Perform phi-integration and partial averaging in theta
        Do q = 1, nq
            Do t = my_theta_min, my_theta_max
                Do r = my_rmin, my_rmax
                    Do p = 1, nphi
                        tmp_buffer(r,q) = tmp_buffer(r,q)+inbuff(p,r,t,q) &
                            & *theta_integration_weights(t)
                    Enddo
                Enddo
            Enddo
        Enddo

        ! Turn phi-integration into an average
        tmp_buffer(:,:) = tmp_buffer(:,:)*over_nphi_double

        ! Complete the averaging process in theta
        Call DALLSUM2D(tmp_buffer, outbuff, pfi%rcomm)

        DeAllocate(tmp_buffer)

    End Subroutine ComputeEll0

    Subroutine IOComputeEll0(inbuff,outbuff)
        !Works exactly like computeEll0, but inbuff has already been averaged in phi
        Real*8, Intent(In) :: inbuff(my_rmin:,my_theta_min:,1:)
        Real*8, Intent(InOut) :: outbuff(my_rmin:,1:)
        Real*8, Allocatable :: tmp_buffer(:,:)
        Integer :: bdims(1:3)
        Integer :: q,nq,r,t,p
        !Averages over theta and phi to get the spherically symmetric mean of all
        ! fields in inbuff at each radii
        ! inbuff is expected to be dimensioned as (my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,1:nfields)
        ! ** Note that this routine should be used sparingly because it requires 
        ! **   a collective operation (allreduce) across process rows
        ! ** One way to do this is to aggregate several fields into the inbuff when calling this routine


        bdims = shape(inbuff)
        nq = bdims(3)

        Allocate(tmp_buffer(my_rmin:my_rmax,1:nq))
        tmp_buffer(:,:) = 0.0d0
        outbuff(:,:) = 0.0d0

        ! Perform partial averaging in theta
        Do q = 1, nq
            Do t = my_theta_min, my_theta_max
                Do r = my_rmin, my_rmax
                    tmp_buffer(r,q) = tmp_buffer(r,q)+inbuff(r,t,q) &
                        & *theta_integration_weights(t)
                Enddo
            Enddo
        Enddo

        !This next line was erroneous.  inbuff has already been averaged in phi
        ! Turn phi-integration into an average
        !tmp_buffer(:,:) = tmp_buffer(:,:)*over_nphi_double

        ! Complete the averaging process in theta
        Call DALLSUM2D(tmp_buffer, outbuff, pfi%rcomm)

        DeAllocate(tmp_buffer)

    End Subroutine IOComputeEll0


    Subroutine Compute_Radial_Average(inbuff,outbuff)
        Real*8, Intent(In) :: inbuff(my_rmin:,1:)
        Real*8, Intent(InOut) :: outbuff(1:)
        Real*8, Allocatable :: tmp_buffer(:)
        Integer :: bdims(1:2)
        Integer :: q,nq,r,t,p
        !Averages over radius for all fields contained in inbuff
        ! fields in inbuff at each radii
        ! inbuff is expected to be dimensioned as (1:nphi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,1:nfields)
        ! ** Note that this routine should be used sparingly because it requires 
        ! **   a collective operation (allreduce) across process rows
        ! ** One way to do this is to aggregate several fields into the inbuff when calling this routine


        bdims = shape(inbuff)
        nq = bdims(2)

        Allocate(tmp_buffer(1:nq))
        tmp_buffer(:) = 0.0d0
        outbuff(:) = 0.0d0

        ! Perform partial averaging in r
        Do q = 1, nq
            Do r = my_rmin, my_rmax
                tmp_buffer(q) = tmp_buffer(q)+inbuff(r,q) &
                    & *r_integration_weights(r)
            Enddo
        Enddo

        ! Complete the averaging process in theta
        Call DALLSUM1D(tmp_buffer, outbuff, pfi%ccomm)

        DeAllocate(tmp_buffer)

    End Subroutine Compute_Radial_Average

    Subroutine IOComputeM0(qty)
        Real*8, Intent(In) :: qty(1:,my_rmin:,my_theta_min:)


        Integer :: q,nq,r,t,p, ind
        !Averages over phi to get the azimuthally symmetric mean of all
        ! fields in inbuff at each radii and theta value.
        ! inbuff is expected to be dimensioned as (1:nphi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,my_theta_min:my_theta_max,1:nfields)



        ind = Shell_Averages%ind

        ! Perform phi-integration

        Do t = my_theta_min, my_theta_max
            Do r = my_rmin, my_rmax
                Do p = 1, nphi
                    IOm0_values(r,t,ind) = IOm0_values(r,t,ind)+qty(p,r,t)
                Enddo
                IOm0_values(r,t,ind) = IOm0_values(r,t,ind)*over_nphi_double ! Turn integration into an average
            Enddo
        Enddo

        Call Shell_Averages%AdvanceInd()

    End Subroutine IOComputeM0

    Subroutine ComputeM0(inbuff,outbuff)
        Real*8, Intent(In) :: inbuff(1:,my_rmin:,my_theta_min:,1:)
        Real*8, Intent(InOut) :: outbuff(my_rmin:,my_theta_min:,1:)

        Integer :: bdims(1:4)
        Integer :: q,nq,r,t,p
        !Averages over phi to get the azimuthally symmetric mean of all
        ! fields in inbuff at each radii and theta value.
        ! inbuff is expected to be dimensioned as (1:nphi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,my_theta_min:my_theta_max,1:nfields)

        bdims = shape(inbuff)
        nq = bdims(4)

        outbuff(:,:,:) = 0.0d0

        ! Perform phi-integration and partial averaging in theta
        Do q = 1, nq
            Do t = my_theta_min, my_theta_max
                Do r = my_rmin, my_rmax
                    Do p = 1, nphi
                        outbuff(r,t,q) = outbuff(r,t,q)+inbuff(p,r,t,q)
                    Enddo
                Enddo
            Enddo
        Enddo

        ! Turn phi-integration into an average
        outbuff = outbuff*over_nphi_double

    End Subroutine ComputeM0

    !/////////////////////////////////////////
    ! Cleanup Routines
    Subroutine CleanUp(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        ! DeAllocates all allocated components of the info structure.
        ! Reinitializes all static components to their initial values (if they had one).
        self%values(1:nqmax) = -1 
        self%levels(1:nshellmax) = -1 
        self%compute(1:nqmax)= -1
        self%nq = 0
        self%nlevels = 0
        self%my_nlevels = 0
        self%file_unit = 15
        self%file_prefix = 'None'
        If (Allocated(self%oqvals))  DeAllocate(self%oqvals)

        self%frequency = 90000000 
        self%rec_per_file =1     
        self%current_rec = 1       
        self%file_header_size =0 
        self%file_record_size = 0 
        self%file_position = 1   
        self%avg_level = 0       

        self%ind = 1              
        self%begin_output = .false.
        self%mpi_tag = 1          

        self%output_version = -1 
        self%grab_this_q = .false.       

        If (Allocated(self%my_shell_levs)) DeAllocate(self%my_shell_levs)
        If (Allocated(self%have_shell)) DeAllocate(self%have_shell)
        If (Allocated(self%my_shell_ind)) DeAllocate(self%my_shell_ind)
        If (Allocated(self%shell_r_ids)) DeAllocate(self%shell_r_ids)
        If (Allocated(self%nshells_at_rid)) DeAllocate(self%nshells_at_rid)
        self%master = .false.
    End Subroutine CleanUP

    Subroutine CleanUP_Spherical_IO
        Implicit None
        !DeAllocates all allocated module arrays
        !Re-initializes module variables to their default values (if any)
        current_averaging_level = 0
        current_qval = 0

        averaging_level(1:nqmax) = 0
        compute_q(1:nqmax) = 0
        shellavg_values(1:nqmax)=-1 
        globalavg_values(1:nqmax)=-1
        shellslice_values(1:nqmax) =-1
        shellslice_levels(1:nshellmax)=-1
        azavg_values(1:nqmax)=-1
        full3d_values(1:nqmax) = -1 
        shellspectra_values(1:nqmax)=-1
        shellspectra_levels(1:nshellmax)=-1
        histo_values(1:nqmax) = -1
        histo_levels(1:nshellmax)=-1

        globalavg_nrec = 1
        shellavg_nrec = 1
        azavg_nrec = 1
        shellslice_nrec =1
        shellspectra_nrec =1

        globalavg_frequency = 90000000
        shellavg_frequency = 90000000
        azavg_frequency = 90000000
        shellslice_frequency = 90000000
        shellspectra_frequency=90000000

        full3d_frequency= 90000000
        local_file_path=''

        integer_zero = 0
        If (Allocated(circumference)) DeAllocate(circumference)
        If (Allocated(qty)) DeAllocate(qty)
        If (Allocated(f_of_r_theta)) DeAllocate(f_of_r_theta)
        If (Allocated(azav_outputs)) DeAllocate(azav_outputs)
        If (Allocated(f_of_r)) DeAllocate(f_of_r)
        If (Allocated(rdtheta_total)) DeAllocate(rdtheta_total)
        If (Allocated(shellav_outputs)) DeAllocate(shellav_outputs)
        If (Allocated(globav_outputs)) DeAllocate(globav_outputs)
        If (Allocated(shell_slice_outputs)) DeAllocate(shell_slice_outputs)
    
        
        If (Allocated(sintheta_dtheta)) DeAllocate(sintheta_dtheta)
        If (Allocated(rsquared_dr)) DeAllocate(rsquared_dr)
        i_ofmt = '(i8.8)'  ! These should never change during a run, but just in case...
        i_pfmt = '(i5.5)'
        io_node = 0
        If (Allocated(theta_integration_weights)) DeAllocate(theta_integration_weights)
        If (Allocated(r_integration_weights))   DeAllocate(r_integration_weights)
        If (Allocated(radius)) DeAllocate(radius)
        If (Allocated(sintheta)) DeAllocate(sintheta)
        If (Allocated(costheta)) DeAllocate(costheta)


        Call Full_3D%cleanup
        Call Global_Averages%cleanup
        Call AZ_Averages%cleanup
        Call Shell_Averages%cleanup
        Call Shell_Slices%cleanup
        Call Shell_Spectra%cleanup
        !Call spectra_buffer%cleanup

    End Subroutine CleanUP_Spherical_IO




End Module Spherical_IO
