Module BoundaryConditions

	Implicit None

	Logical :: Fix_T_Top    = .True.
	Logical :: Fix_T_Bottom = .True.
	Real*8  :: T_Bottom     = 0.0d0
	Real*8  :: T_Top        = 0.0d0

	Namelist /Boundary_Conditions_Namelist/ Fix_T_Top, Fix_T_Bottom, T_Bottom, T_Top

End Module BoundaryConditions
