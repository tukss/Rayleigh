#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
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
                        & *(fbuffer(PSI,vtheta)*fbuffer(:,r,t,vphi))
                END_DO

                Call Add_Quantity(qty)
            Endif 


            If (compute_quantity(amom_dr_r)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*radius(r)*sintheta(t) &
                            & *(m0_values(r,t,vr)*m0_values(r,t,vphi))
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif 

            If (compute_quantity(amom_dr_theta)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*radius(r)*sintheta(t) &
                            & *(m0_values(r,t,vtheta)*m0_values(r,t,vphi))
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif 

            If (compute_quantity(amom_mean_r)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*((radius(r)*sintheta(t))**2) &
                            & *(m0_values(r,t,vr)*Angular_Velocity)
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif 

            If (compute_quantity(amom_mean_theta)) Then

                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        qty(:,r,t) = ref%density(r)*((radius(r)*sintheta(t))**2) &
                            & *(m0_values(r,t,vtheta)*Angular_Velocity)
                    EndDo
                EndDo

                Call Add_Quantity(qty)
            Endif           
        
    End Subroutine Compute_Angular_Momentum_Balance

End Module Diagnostics_Angular_Momentum
