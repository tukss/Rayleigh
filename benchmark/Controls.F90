Module Controls
	! Things that control how the simulation runs
	Implicit None
	!Real*8 :: Ra, Ek, Pr, Pm
	Real*8 :: alpha_implicit = 0.5d0
	Logical :: chebyshev = .false.
	Logical :: bandsolve = .false.
	Logical :: magnetism = .false.
	Logical :: nonlinear = .true.
	Logical :: Rotation = .false.
	Integer :: max_iterations = 1000000
   Integer :: check_frequency = 2
	Logical :: static_transpose = .false.
	Logical :: static_config = .false.
	Logical :: use_parity = .true.
	Logical :: test_reduce = .false.
	Logical :: lorentz_forces = .false.
	Logical :: deriv_cluge = .false.
	Integer :: gpower = 1
    Logical :: Conserve_L = .false.
	Logical :: pad_alltoall = .false.
	Logical :: read_argv = .false.
    Logical :: no_slip_boundaries = .false.
    Real*8  :: cflmax = 0.1d0, cflmin = 0.05d0
	Real*8 :: max_time_step = 5.0d-4
	Namelist /Controls_Namelist/ max_iterations, chebyshev, nonlinear, &
			alpha_implicit, check_frequency, rotation, static_transpose, static_config, &
			use_parity, test_reduce,magnetism, gpower, lorentz_forces, deriv_cluge, &
            cflmax, cflmin, Conserve_L, pad_alltoall, bandsolve, max_time_step, read_argv, &
            no_slip_boundaries
End Module Controls
