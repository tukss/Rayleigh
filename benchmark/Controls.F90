Module Controls
	! Things that control how the simulation runs
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
	Namelist /Controls_Namelist/ Ra, Ek, Pr, max_iterations, chebyshev, nonlinear, &
			alpha_implicit, check_frequency, rotation, static_transpose, static_config, &
			use_parity
End Module Controls
