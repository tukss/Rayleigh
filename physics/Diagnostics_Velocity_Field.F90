#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_VELOCITY_FIELD
!               This module computes the components of the velocity field 
!               and their derivatives.  Zonal means and fluctuations about
!               those means are also computed.
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Velocity_Field
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_Velocity_Components(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        !/////////////////////////////////////////
        ! 1. terms involving radial velocity
        If (compute_quantity(v_r)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(vp_r)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,vr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(vm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- radial derivatives of v_r
        If (compute_quantity(dv_r_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvrdr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_r_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvrdr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_r_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- theta derivatives of v_r
        If (compute_quantity(dv_r_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvrdt)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_r_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvrdt)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_r_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- phi derivatives of v_r
        If (compute_quantity(dv_r_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvrdp)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_r_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvrdp)
            Call Add_Quantity(qty)
        Endif		

        If (compute_quantity(dvm_r_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif		

        !-- -- {1/r d(v_r)/dtheta}
        If (compute_quantity(dv_r_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_r_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_r_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- {1/(r sin[theta]) d(v_r)/dphi
        If (compute_quantity(dv_r_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_r_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        If (compute_quantity(dvm_r_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif		

        !/////////////////////////////////////////
        ! 2. terms involving theta velocity
        If (compute_quantity(v_theta)) Then	
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vtheta)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(vp_theta)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,vtheta)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(vm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- radial derivatives of v_theta
        If (compute_quantity(dv_theta_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvtdr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_theta_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvtdr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_theta_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- theta derivatives of v_theta
        If (compute_quantity(dv_theta_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvtdt)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_theta_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvtdt)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_theta_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- phi derivatives of v_theta
        If (compute_quantity(dv_theta_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvtdp)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_theta_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvtdp)
            Call Add_Quantity(qty)
        Endif		

        If (compute_quantity(dvm_theta_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- {1/r d(v_theta)/dtheta}
        If (compute_quantity(dv_theta_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_theta_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_theta_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- {1/(r sin[theta]) d(v_theta)/dphi
        If (compute_quantity(dv_theta_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_theta_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif		

        If (compute_quantity(dvm_theta_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif	



        !/////////////////////////////////////////
        ! 3. terms involving phi velocity
        If (compute_quantity(v_phi)) Then	
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(vp_phi)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,vphi)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(vm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- radial derivatives of v_phi
        If (compute_quantity(dv_phi_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvpdr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_phi_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvpdr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_phi_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- theta derivatives of v_phi
        If (compute_quantity(dv_phi_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvpdt)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_phi_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvpdt)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_phi_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdt)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- phi derivatives of v_phi
        If (compute_quantity(dv_phi_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvpdp)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_phi_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvpdp)
            Call Add_Quantity(qty)
        Endif		

        If (compute_quantity(dvm_phi_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- {1/r d(v_phi)/dtheta}
        If (compute_quantity(dv_phi_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_phi_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_phi_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !-- -- {1/(r sin[theta]) d(v_phi)/dphi
        If (compute_quantity(dv_phi_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvp_phi_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(dvm_phi_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !/////////////////////////////////////////////////////////////
        ! Finally, if desired, we compute the mass flux.

        If (compute_quantity(rhov_r)) Then
            DO_PSI			
                qty(PSI) = buffer(PSI,vr)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        If (compute_quantity(rhovp_r)) Then
            DO_PSI			
                qty(PSI) = fbuffer(PSI,vr)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		

        If (compute_quantity(rhovm_r)) Then
            DO_PSI			
                qty(PSI) = m0_values(PSI2,vr)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif	


        If (compute_quantity(rhov_theta)) Then
            DO_PSI			
                qty(PSI) = buffer(PSI,vtheta)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		

        If (compute_quantity(rhovp_theta)) Then
            DO_PSI			
                qty(PSI) = fbuffer(PSI,vtheta)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		

        If (compute_quantity(rhovm_theta)) Then
            DO_PSI			
                qty(PSI) = m0_values(PSI2,vtheta)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif	


        If (compute_quantity(rhov_phi)) Then
            DO_PSI			
                qty(PSI) = buffer(PSI,vphi)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        If (compute_quantity(rhovp_phi)) Then
            DO_PSI			
                qty(PSI) = fbuffer(PSI,vphi)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif		

        If (compute_quantity(rhovm_phi)) Then
            DO_PSI			
                qty(PSI) = m0_values(PSI2,vphi)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif	

        
    End Subroutine Compute_Velocity_Components

End Module Diagnostics_Velocity_Field
