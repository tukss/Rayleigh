#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_CURRENT_DENSITY
!               This module computes the components of the current density. 
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
        If (compute_quantity(jp_r)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,jr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,jr)
            END_DO
            Call Add_Quantity(qty)
        Endif	


        !/////////////////////////////////////////
        ! 2. terms involving J_theta
        If (compute_quantity(j_theta)) Then	
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jtheta)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jp_theta)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,jtheta)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,jtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif	


        !/////////////////////////////////////////
        ! 3. terms involving J_phi
        If (compute_quantity(j_phi)) Then	
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,jphi)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jp_phi)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,jphi)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,jphi)
            END_DO
            Call Add_Quantity(qty)
        Endif	


        
    End Subroutine Compute_J_Components

End Module Diagnostics_Current_Density
