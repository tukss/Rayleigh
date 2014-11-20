Module BoundaryConditions

	Implicit None

	Logical :: Fix_Tvar_Top    = .True.
	Logical :: Fix_Tvar_Bottom = .True.
	Logical :: Fix_dTdr_Top    = .False.
	Logical :: Fix_dTdr_Bottom = .False.
    Logical :: Fix_divrt_top = .False.
    Logical :: Fix_divt_top = .False.
    Logical :: Fix_divrfc_top = .False.
    Logical :: Fix_divfc_top = .False.
	Real*8  :: T_Bottom     = 1.0d0
	Real*8  :: T_Top        = 0.0d0
    Real*8  :: dTdr_Top     = 0.0d0
    Real*8  :: dTdr_Bottom  = 0.0d0

    Logical :: Strict_L_Conservation = .false.         ! (In-Progress) Turn on to enforce angular momentum conservation abous x,y, and z-axes
    Logical :: no_slip_boundaries = .false. ! Set to true to use no-slip boundaries.  Stree-free boundaries are the default.

	Namelist /Boundary_Conditions_Namelist/ Fix_Tvar_Top, Fix_Tvar_Bottom, T_Bottom, T_Top, dTdr_top, dTdr_bottom, &
		fix_dtdr_bottom, fix_dtdr_top, fix_divrt_top, fix_divt_top, fix_divrfc_top, fix_divfc_top, &
        no_slip_boundaries, strict_L_Conservation

End Module BoundaryConditions
