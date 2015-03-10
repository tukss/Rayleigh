Module Controls
	! Things that control how the simulation runs
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

    !////////////////////////////////////////////////////
    !   Temporal Controls
    !   Flags that control details of the time-stepping (some relate to the numerics, but we keep the time-related things together).
	Real*8  :: alpha_implicit = 0.51d0            ! Crank Nicolson Implict/Explicit weighting factor (1.0 is fully implicit)
	Integer :: max_iterations = 1000000         ! The maximum number of iterations to be run in a given session
    Integer :: check_frequency = 20000          ! Number of iterations between checkpoint dumps
    Real*8  :: cflmax = 0.4d0, cflmin = 0.6d0  ! Limits for the cfl condition
	Real*8  :: max_time_step = 1.0d0            ! Maximum timestep to take, whatever CFL says (should always specify this in main_input file)
    Integer :: chk_type = 1                     ! Set to 2 for memory friendly IO.  In development
    Integer :: diagnostic_reboot_interval = -1
    Namelist /Temporal_Controls_Namelist/ alpha_implicit, max_iterations, check_frequency, &
                & cflmax, cflmin, max_time_step,chk_type, diagnostic_reboot_interval

Contains
    Subroutine Initialize_Controls()
        Implicit None
        !Set default for diagnostic_reboot_interval (if necessary)
        If (diagnostic_reboot_interval .le. 0) Then
            diagnostic_reboot_interval = check_frequency
        Endif
    End Subroutine Initialize_Controls

End Module Controls
