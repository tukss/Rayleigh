Module Controls
	! Things that control how the simulation runs
	Implicit None
	Real*8 :: Ra, Ek, Pr, Pm
	Real*8 :: alpha_implicit = 0.5d0
	Logical :: chebyshev = .false.
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
	Namelist /Controls_Namelist/ Ra, Ek, Pr, Pm,max_iterations, chebyshev, nonlinear, &
			alpha_implicit, check_frequency, rotation, static_transpose, static_config, &
			use_parity, test_reduce,magnetism, gpower, lorentz_forces, deriv_cluge
End Module Controls
