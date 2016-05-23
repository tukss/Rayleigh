#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_CURRENT_DENSITY
!               This module computes the components of del x B. 
!               Zonal means and fluctuations about those means are 
!               also computed (if desired).
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Current_Density
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_J_Components(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        !/////////////////////////////////////////
        ! 1. terms involving J_r
        If (compute_quantity(j_r)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(j_r_sq) .or. compute_quantity(j_sq)) Then 
            DO_PSI
                qty(PSI) = buffer(PSI,jr)**2
            END_DO
            If (compute_quantity(j_r_sq)) Call Add_Quantity(qty)
            If (compute_quantity(j_sq)) tmp1 = qty
        Endif


        !/////////////////////////////////////////
        ! 2. terms involving J_theta
        If (compute_quantity(j_theta)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jtheta)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(j_theta_sq) .or. compute_quantity(j_sq)) Then 
            DO_PSI
                qty(PSI) = buffer(PSI,jtheta)**2
            END_DO
            If (compute_quantity(j_theta_sq)) Call Add_Quantity(qty)
            If (compute_quantity(j_sq)) tmp1 = tmp1+qty
        Endif




        !/////////////////////////////////////////
        ! 3. terms involving J_phi
        If (compute_quantity(j_phi)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jphi)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(j_phi_sq) .or. compute_quantity(j_sq)) Then 
            DO_PSI
                qty(PSI) = buffer(PSI,jphi)**2
            END_DO
            If (compute_quantity(j_phi_sq)) Call Add_Quantity(qty)
            If (compute_quantity(j_sq)) Then
                tmp1 = tmp1+qty
                Call Add_Quantity(tmp1)
            Endif 
        Endif


        !////////////////////////////////////////////////////
        !       Now the perturbation and mean terms
        !1.)  J_r
        If (compute_quantity(jp_r)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,jr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jp_r_sq) .or. compute_quantity(jp_sq)) Then 
            DO_PSI
                qty(PSI) = fbuffer(PSI,jr)**2
            END_DO
            If (compute_quantity(jp_r_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jp_sq)) tmp1= qty
        Endif

        If (compute_quantity(jm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,jr)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !2.) J_theta
        If (compute_quantity(jp_theta)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,jtheta)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jp_theta_sq) .or. compute_quantity(jp_sq)) Then 
            DO_PSI
                qty(PSI) = fbuffer(PSI,jtheta)**2
            END_DO
            If (compute_quantity(jp_theta_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jp_sq)) tmp1= tmp1+qty
        Endif

        !3.) J_phi
        If (compute_quantity(jm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,jtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        If (compute_quantity(jp_phi)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,jphi)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jp_phi_sq) .or. compute_quantity(jp_sq)) Then 
            DO_PSI
                qty(PSI) = fbuffer(PSI,jphi)**2
            END_DO
            If (compute_quantity(jp_phi_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jp_sq)) Then
                tmp1= tmp1+qty
                Call Add_Quantity(tmp1)
            Endif 
        Endif

        If (compute_quantity(jm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,jphi)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !///////////////////////////////////////////
        !  Also handle ohmic heating here.
        If (compute_quantity(ohmic_heat)) Then
            If (ohmic_heating) Then
                DO_PSI
                    qty(PSI) = ohmic_heating_coeff(r)*(buffer(PSI,jr)**2 + &
                               &   buffer(PSI,jtheta)**2 + &
                               &   buffer(PSI,jphi)**2)
                END_DO
            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(ohmic_heat_pp)) Then
            If (ohmic_heating) Then
                DO_PSI
                    qty(PSI) = ohmic_heating_coeff(r)*(fbuffer(PSI,jr)**2 + &
                               &   fbuffer(PSI,jtheta)**2 + &
                               &   fbuffer(PSI,jphi)**2)
                END_DO
            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif        

        If (compute_quantity(ohmic_heat_mm)) Then
            If (ohmic_heating) Then
                DO_PSI
                    qty(PSI) = ohmic_heating_coeff(r)*(m0_values(PSI2,jr)**2 + &
                               &   m0_values(PSI2,jtheta)**2 + &
                               &   m0_values(PSI2,jphi)**2)
                END_DO

            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif   

    End Subroutine Compute_J_Components

End Module Diagnostics_Current_Density
