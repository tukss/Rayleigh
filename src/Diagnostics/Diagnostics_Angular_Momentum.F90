#define DO_PSI Do t = my_theta%min, my_theta%max; Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max; Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

Module Diagnostics_Angular_Momentum
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_Angular_Momentum_Balance(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        If (compute_quantity(amom_fluct_r)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*radius(r)*sintheta(t) &
                    & *(fbuffer(PSI,vr)*fbuffer(PSI,vphi))
            END_DO

            Call Add_Quantity(qty)
        Endif 

        If (compute_quantity(amom_fluct_theta)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*radius(r)*sintheta(t) &
                    & *(fbuffer(PSI,vtheta)*fbuffer(PSI,vphi))
            END_DO

            Call Add_Quantity(qty)
        Endif 


        If (compute_quantity(amom_dr_r)) Then

            DO_PSI
                    qty(PSI) = ref%density(r)*radius(r)*sintheta(t) &
                        & *(m0_values(PSI2,vr)*m0_values(PSI2,vphi))
            END_DO

            Call Add_Quantity(qty)
        Endif 

        If (compute_quantity(amom_dr_theta)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*radius(r)*sintheta(t) &
                    & *(m0_values(PSI2,vtheta)*m0_values(PSI2,vphi))
            END_DO

            Call Add_Quantity(qty)
        Endif 

        !These need to be adjusted to handle non-dimensionalization
        If (compute_quantity(amom_mean_r)) Then
            
            DO_PSI
                qty(PSI) = ref%density(r)*((radius(r)*sintheta(t))**2) &
                    & *(m0_values(PSI2,vr)*Angular_Velocity)
            END_DO

            Call Add_Quantity(qty)
        Endif 

        If (compute_quantity(amom_mean_theta)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*((radius(r)*sintheta(t))**2) &
                    & *(m0_values(PSI2,vtheta)*Angular_Velocity)
            END_DO

            Call Add_Quantity(qty)
        Endif           

        If (magnetism) Then
            If (compute_quantity(maxwell_stress_r)) Then
                DO_PSI
                    qty(PSI) = -radius(r)*sintheta(t)*fbuffer(PSI,br)*fbuffer(PSI,bphi)*ref%Lorentz_Coeff
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(maxwell_stress_theta)) Then
                DO_PSI
                    qty(PSI) = -radius(r)*sintheta(t)*fbuffer(PSI,btheta)*fbuffer(PSI,bphi)*ref%Lorentz_Coeff
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(magnetic_torque_r)) Then
                DO_PSI
                    qty(PSI) = -radius(r)*sintheta(t)*m0_values(PSI2,br)*m0_values(PSI2,bphi)*ref%Lorentz_Coeff
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(magnetic_torque_theta)) Then
                DO_PSI
                    qty(PSI) = -radius(r)*sintheta(t)*m0_values(PSI2,btheta)*m0_values(PSI2,bphi)*ref%Lorentz_Coeff
                END_DO
                Call Add_Quantity(qty)
            Endif


        Endif
        
    End Subroutine Compute_Angular_Momentum_Balance

End Module Diagnostics_Angular_Momentum
