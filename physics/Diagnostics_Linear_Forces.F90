#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

Module Diagnostics_Linear_Forces
    Use Diagnostics_Base
    Implicit None

Contains

    Subroutine Compute_Linear_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Call Compute_Buoyancy_Force(buffer)
        Call Compute_Coriolis_Force(buffer)
        Call Compute_Pressure_Force(buffer)
        Call Compute_Viscous_Force(buffer)
    End Subroutine Compute_Linear_Forces

    Subroutine Compute_Buoyancy_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        ! -- full buoyancy 
        If (compute_quantity(buoyancy_force)) Then
            DO_PSI
                qty(PSI) = ref%gravity_term_s(r)*(buffer(PSI,tvar)-&
                           & ell0_values(r,tvar))
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! -- fluctuating buoyancy (ell = 0, m =0 already subtracted)
        If (compute_quantity(buoyancy_pforce)) Then
            DO_PSI
                qty(PSI) = ref%gravity_term_s(r)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! -- mean buoyancy
        If (compute_quantity(buoyancy_pforce)) Then
            DO_PSI
                qty(PSI) = ref%gravity_term_s(r)*(m0_values(PSI2,tvar)-&
                           & ell0_values(r,tvar))
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Buoyancy_Force

    Subroutine Compute_Coriolis_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)

    End Subroutine Compute_Coriolis_Force

    Subroutine Compute_Pressure_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)

    End Subroutine Compute_Pressure_Force

    Subroutine Compute_Viscous_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)

    End Subroutine Compute_Viscous_Force

End Module Diagnostics_Linear_Forces
