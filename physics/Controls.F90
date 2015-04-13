Module Controls
	! Things that control how the simulation runs
    Use BufferedOutput
	Implicit None
    
    !////////////////////////////////////////////////////////////////////////////////
    ! Multiple run controls,  These are not set in a namelist, but are used through command line options.
    Character*120, Allocatable :: rundirs(:)
    Logical :: multi_run_mode = .false.
    Integer :: nruns = 0 ! only set to non-zero value if multi_run_mode is True
    Integer, Allocatable :: run_cpus(:)
    Character*120 :: my_path = ''

    !////////////////////////////////////////////////////////////////////////////////
    ! Numerical Controls
    ! Flats that control details of the parallelization/data layout and 
    ! how the equations are solved (not what equations are solved).
	Logical :: chebyshev = .false.          ! Set to true to use chebyshev polynomials in radius (default is finite-difference)
	Logical :: bandsolve = .false.          ! Set to true to use band solves with the finite-differences
	Logical :: static_transpose = .false.   ! When true, transpose buffers for sending/receiving are never de-allocated for spherical buffer objects
	Logical :: static_config = .false.      ! When true, configuration buffers (p3a, s1b, etc.) ar enever de-allocated for spherical buffer objects
	Logical :: use_parity = .true.          ! Possibly defunct - should always be true
	Logical :: deriv_cluge = .true.         ! Use modified 2nd derivative in radius for finite-differences (leave true for stability...for now)
	Logical :: pad_alltoall = .false.       ! Normally all-to-allv is used.  Standard alltoall with zero padded buffers can be used when this flag is on.



    Namelist /Numerical_Controls_Namelist/ chebyshev, bandsolve, static_transpose, static_config, &
            & use_parity, deriv_cluge, pad_alltoall

    !////////////////////////////////////////////////////////////////////////////////
    ! Physical Controls
    ! Flags that control various fundamental aspects of the physics employed
	Logical :: magnetism = .false.          ! Turn magnetism on or off
	Logical :: nonlinear = .true.           ! Nonlinear terms can be turned off (calculated but zeroed out - for debugging)
	Logical :: Rotation = .false.           ! Rotate or not
	Logical :: lorentz_forces = .true.     ! Turn Lorentz forces on or off
    Logical :: viscous_heating = .true.     ! Turns viscous heating on/off
    Logical :: ohmic_heating = .true.

    Namelist /Physical_Controls_Namelist/ magnetism, nonlinear, rotation, lorentz_forces, &
                & viscous_heating, ohmic_heating

    !///////////////////////////////////////////////////////////////////////////
    !   Temporal Controls
    !   Flags that control details of the time-stepping (some relate to the numerics, but we keep the time-related things together).
	Real*8  :: alpha_implicit = 0.51d0            ! Crank Nicolson Implict/Explicit weighting factor (1.0 is fully implicit)
	Integer :: max_iterations = 1000000         ! The maximum number of iterations to be run in a given session
    Integer :: check_frequency = 20000          ! Number of iterations between checkpoint dumps
    Real*8  :: cflmax = 0.4d0, cflmin = 0.6d0  ! Limits for the cfl condition
	Real*8  :: max_time_step = 1.0d0            ! Maximum timestep to take, whatever CFL says (should always specify this in main_input file)
    Real*8  :: min_time_step = 1.0d-13
    Integer :: chk_type = 1                     ! Set to 2 for memory friendly IO.  In development
    Integer :: diagnostic_reboot_interval = -1
    Namelist /Temporal_Controls_Namelist/ alpha_implicit, max_iterations, check_frequency, &
                & cflmax, cflmin, max_time_step,chk_type, diagnostic_reboot_interval, min_time_step



    !///////////////////////////////////////////////////////////////////////////
    ! I/O Controls
    ! What is normally sent to standard out can, if desired, be sent to a file instead
    Integer :: stdout_flush_interval = 50  ! Lines stored before stdout buffer is flushed to stdout_unit
    Character*120 :: stdout_file = 'nofile'
    Type(OutputBuffer) :: stdout
    Namelist /IO_Controls_Namelist/ stdout_flush_interval,stdout_file


Contains
    Subroutine Initialize_Controls()
        Implicit None
        character*120 :: ofilename
        !Set default for diagnostic_reboot_interval (if necessary)
        If (diagnostic_reboot_interval .le. 0) Then
            diagnostic_reboot_interval = check_frequency
        Endif

        !Initialize the stdout buffer -- by default, write to unit 6 with frequency of 1

		Select Case(stdout_file)
			Case('stdout')	! Standard out, but flush with user-defined frequency
                Call stdout%init(6,line_count = stdout_flush_interval)
			Case('nofile')
                Call stdout%init(6) ! Standard out, with effectively no buffering (line_count = 1)
		    Case Default
			    ! All stdout written to file, flushed at user-defined flush interval
                ofilename = Trim(my_path)//Trim(stdout_file)
                Call stdout%init(116,line_count = stdout_flush_interval,filename=ofilename)
		End Select

    End Subroutine Initialize_Controls

  

End Module Controls
