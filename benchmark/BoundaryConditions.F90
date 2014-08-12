Module BoundaryConditions

	Implicit None

	Logical :: Fix_Tvar_Top    = .True.
	Logical :: Fix_Tvar_Bottom = .True.
	Logical :: Fix_dTdr_Top    = .False.
	Logical :: Fix_dTdr_Bottom = .False.
	Real*8  :: T_Bottom     = 1.0d0
	Real*8  :: T_Top        = 0.0d0
    Real*8  :: dTdr_Top     = 0.0d0
    Real*8  :: dTdr_Bottom  = 0.0d0

	Namelist /Boundary_Conditions_Namelist/ Fix_Tvar_Top, Fix_Tvar_Bottom, T_Bottom, T_Top, dTdr_top, dTdr_bottom, &
		fix_dtdr_bottom, fix_dtdr_top

End Module BoundaryConditions
