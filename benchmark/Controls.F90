Module Controls
	! Things that control how the simulation runs
	Real*8 :: Ra, Ek, Pr, Pm
	Real*8 :: alpha_implicit = 0.5d0
	Logical :: chebyshev = .false.
	Logical :: magnetism = .false.
	Logical :: nonlinear = .true.
	Integer :: max_iterations = 1000000
    Integer :: checkpoint_frequency = 2
	Namelist /Controls_Namelist/ Ra, Ek, Pr, max_iterations, chebyshev, nonlinear, alpha_implicit, checkpoint_frequency
End Module Controls
