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
        Call Compute_Viscous_Force(buffer)
    End Subroutine Compute_Linear_Forces

    Subroutine Compute_Buoyancy_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        ! Buoyancy forces as they contribute to the ell .ne. 0 components of the momentum equation

        ! -- full buoyancy 
        If (compute_quantity(buoyancy_force)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*(buffer(PSI,tvar)-&
                           & ell0_values(r,tvar))
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! -- fluctuating buoyancy (ell = 0, m =0 already subtracted)
        If (compute_quantity(buoyancy_pforce)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! -- mean buoyancy
        If (compute_quantity(buoyancy_mforce)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*(m0_values(PSI2,tvar)-&
                           & ell0_values(r,tvar))
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Buoyancy_Force

    Subroutine Compute_Coriolis_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: coriolis_term
        coriolis_term = ref%Coriolis_Coeff

        If(compute_quantity(Coriolis_Force_r)) Then
            DO_PSI
                qty(PSI) = -ref%density(r)*coriolis_term*sintheta(t)*buffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If(compute_quantity(Coriolis_pForce_r)) Then
            DO_PSI
                qty(PSI) = -ref%density(r)*coriolis_term*sintheta(t)*fbuffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If(compute_quantity(Coriolis_mForce_r)) Then
            DO_PSI
                qty(PSI) = -ref%density(r)*coriolis_term*sintheta(t)*m0_values(PSI2,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_Force_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*buffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_pForce_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*fbuffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_mForce_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*m0_values(PSI2,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(Coriolis_Force_Phi)) Then
            DO_PSI
                qty(PSI) = - coriolis_term*costheta(t)*buffer(PSI,vtheta) &
					       - coriolis_term*sintheta(t)*buffer(PSI,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_pForce_Phi)) Then
            DO_PSI
                qty(PSI) = - coriolis_term*costheta(t)*fbuffer(PSI,vtheta) &
					       - coriolis_term*sintheta(t)*fbuffer(PSI,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_mForce_Phi)) Then
            DO_PSI
                qty(PSI) = - coriolis_term*costheta(t)*m0_values(PSI2,vtheta) &
					       - coriolis_term*sintheta(t)*m0_values(PSI2,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Coriolis_Force



    Subroutine Compute_Viscous_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        ! Placeholder
        ! The radial force term can be grabbed easily enough following the solve.
        ! 
        ! The other two components are a bit complicated to calculate due to the higher order
        ! derivatives involved.  Will return to this later...
        ! It's  possible that this should be special routine, allowed to transform
        ! back to ell-space for the derivatives.
        If (compute_quantity(viscous_Force_r)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(viscous_pForce_r)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(viscous_mForce_r)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)
        Endif



        If (compute_quantity(viscous_Force_theta)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(viscous_pForce_theta)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(viscous_mForce_theta)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)
        Endif



        If (compute_quantity(viscous_Force_phi)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(viscous_pForce_phi)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(viscous_mForce_phi)) Then
            DO_PSI
                qty(PSI) = 0.0d0
            END_DO
            Call Add_Quantity(qty)

        Endif

    End Subroutine Compute_Viscous_Force

End Module Diagnostics_Linear_Forces
