#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

Module Diagnostics_Advection

    Use Diagnostics_Base
Contains
    Subroutine Compute_Inertial_Terms(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8, Allocatable :: tbuffer(:,:,:,:)
        Integer :: vr_p     = 1, dvrdr_p = 2,dvrdt_p = 3
        Integer :: vtheta_p = 4, dvtdr_p = 5,dvtdt_p = 6
        Integer :: vphi_p   = 7, dvpdr_p = 8,dvpdt_p = 9
        Allocate(tbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:9))
        !The azimuthal means of v and its derivatives have already been computed, but we need 
        ! the perturbations as well - possibly have a flag for computing
        ! if needed. 
        ! Note that the azimuthal mean of d[]dphi is always zero, so we don't compute
        ! fluctuations of d[vr]dphi, d[vtheta]dphi, or d[vphi]dphi

        DO_PSI
            tbuffer(PSI,vr_p) = buffer(PSI,vr) - m0_values(PSI2,vr) 
        END_DO

        DO_PSI
            tbuffer(PSI,dvrdr_p) = buffer(PSI,dvrdr) - m0_values(PSI2,dvrdr) 
        END_DO

        DO_PSI
            tbuffer(PSI,dvrdt_p) = buffer(PSI,dvrdt) - m0_values(PSI2,dvrdt) 
        END_DO

        DO_PSI
            tbuffer(PSI,vtheta_p) = buffer(PSI,vtheta) - m0_values(PSI2,vtheta) 
        END_DO

        DO_PSI
            tbuffer(PSI,dvtdr_p) = buffer(PSI,dvtdr) - m0_values(PSI2,dvtdr) 
        END_DO

        DO_PSI
            tbuffer(PSI,dvtdt_p) = buffer(PSI,dvtdt) - m0_values(PSI2,dvtdt) 
        END_DO

        DO_PSI
            tbuffer(PSI,vphi_p) = buffer(PSI,vphi) - m0_values(PSI2,vphi) 
        END_DO        

        DO_PSI
            tbuffer(PSI,dvpdr_p) = buffer(PSI,dvpdr) - m0_values(PSI2,dvpdr) 
        END_DO

        DO_PSI
            tbuffer(PSI,dvpdt_p) = buffer(PSI,dvpdt) - m0_values(PSI2,dvpdt) 
        END_DO


        !//////////////////////////////////////////////////////////////////////////////////
        !/////////////// v dot grad v (full)//////////////////

        !--- Radial Component
        If (compute_quantity(v_grad_v_r)) Then
            DO_PSI
			    qty(PSI) = buffer(PSI,vr)*   buffer(PSI,dvrdr)                  &
				       + ( buffer(PSI,vtheta) * ( buffer(PSI,dvrdt)             &
                                                 -buffer(PSI,vtheta) )          &
				       +   buffer(PSI,vphi)   * ( buffer(PSI,dvrdp)*csctheta(t) &
                                                 -buffer(PSI,vphi) ) )          &
                       * one_over_r(r)  
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        
        !--- Theta Component
        If (compute_quantity(v_grad_v_theta)) Then

            DO_PSI
                qty(PSI) =  buffer(PSI,vr  ) *  buffer(PSI,dvtdr)             & ! vr d[vtheta]/dr
                         + (buffer(PSI,vtheta)*(buffer(PSI,dvtdt)             & ! vtheta/r d[vtheta]/dtheta
                         +         buffer(PSI,vr))                            & ! (vtheta vr )/r
                         +  buffer(PSI,vphi) *( buffer(PSI,dvtdp)*csctheta(t) & ! vphi/sintheta/r d[vtheta]/dphi
				         -         buffer(PSI,vphi )*cottheta(t) ) )          & ! vphi^2 cot(theta)/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        !--- Phi Component
        If (compute_quantity(v_grad_v_phi)) Then

            DO_PSI
                qty(PSI) =  buffer(PSI,vr  ) *  buffer(PSI,dvpdr)             & ! vr d[vphi]/dr
                         + (buffer(PSI,vtheta)*(buffer(PSI,dvpdt)             & ! vtheta/r d[vphi]/dtheta
                         +         buffer(PSI,vphi)*cottheta(t))              & ! {vtheta vphi cot(theta)}/r
                         +  buffer(PSI,vphi) *( buffer(PSI,dvpdp)*csctheta(t) & ! vphi/sintheta/r d[vphi]/dphi
				         +         buffer(PSI,vr )) )                         & ! {vphi vr}/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif


        !//////////////////////////////////////////////////////////////////////////////////////////
        !/////////////// v' dot grad <v> //////////////////
        !--- Radial Component
        If (compute_quantity(vp_grad_vm_r)) Then
            DO_PSI
			    qty(PSI) = tbuffer(PSI,vr_p)*m0_values(PSI2,dvrdr) &
				    + tbuffer(PSI,vtheta_p) * ( m0_values(PSI2,dvrdt)-m0_values(PSI2,vtheta) )*one_over_r(r)    &
				    + tbuffer(PSI,vphi_p)*(-m0_values(PSI2,vphi) )*one_over_r(r)  
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !--- Theta Component
        If (compute_quantity(vp_grad_vm_theta)) Then

            DO_PSI
                qty(PSI) =  tbuffer(PSI,vr_p ) *  m0_values(PSI2,dvtdr)  & ! vr' d<vtheta>/dr
                         + (tbuffer(PSI,vtheta_p)*(m0_values(PSI2,dvtdt) & ! vtheta'/r d<vtheta>/dtheta
                         +         m0_values(PSI2,vr))                   & ! (vtheta' <vr> )/r
                         +  tbuffer(PSI,vphi_p) *(                       & ! NO vphi/sintheta/r d<vtheta>/dphi term
				         -         m0_values(PSI2,vphi )*cottheta(t) ) ) & ! {vphi' <vphi> cot(theta)}/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        !--- Phi Component
        If (compute_quantity(vp_grad_vm_phi)) Then

            DO_PSI
                qty(PSI) =  tbuffer(PSI,vr_p  ) *  m0_values(PSI2,dvpdr) & ! vr' d[<vphi>]/dr
                         + (tbuffer(PSI,vtheta_p)*(m0_values(PSI2,dvpdt) & ! vtheta'/r d[<vphi>]/dtheta
                         +         m0_values(PSI2,vphi)*cottheta(t))     & ! {vtheta' <vphi> cot(theta)}/r
                         +  tbuffer(PSI,vphi_p) *(                       & ! NO vphi'/sintheta/r d[<vphi>]/dphi term
				         +         m0_values(PSI2,vr )) )                & ! {vphi' <vr>}/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        !////////////////////////////////////////////////////////////////////////////////////////
        !/////////////// <v> dot grad v' //////////////////
        If (compute_quantity(vm_grad_vp_r)) Then
            DO_PSI
			    qty(PSI) = m0_values(PSI2,vr)*tbuffer(PSI,dvrdr) &
				    + m0_values(PSI2,vtheta) * ( tbuffer(PSI,dvrdt)-tbuffer(PSI,vtheta) )*one_over_r(r)    &
				    + m0_values(PSI2,vphi)*(buffer(PSI,dvrdp)*csctheta(t)-tbuffer(PSI,vphi) )*one_over_r(r)  
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        
        !--- Theta Component
        If (compute_quantity(vm_grad_vp_theta)) Then

            DO_PSI
                qty(PSI) =  m0_values(PSI2,vr  ) *  tbuffer(PSI,dvtdr_p)          & ! <vr> d[vtheta']/dr
                         + (m0_values(PSI2,vtheta)*(tbuffer(PSI,dvtdt_p)          & ! <vtheta>/r d[vtheta']/dtheta
                         +         tbuffer(PSI,vr_p))                            & ! (<vtheta> vr' )/r
                         +  m0_values(PSI2,vphi) *( buffer(PSI,dvtdp)*csctheta(t) & ! <vphi>/sintheta/r d[vtheta']/dphi
				         -         tbuffer(PSI,vphi_p )*cottheta(t) ) )          & ! {<vphi> vphi' cot(theta)}/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        !--- Phi Component
        If (compute_quantity(vm_grad_vp_phi)) Then

            DO_PSI
                qty(PSI) =  m0_values(PSI2,vr  ) *  tbuffer(PSI,dvpdr_p)          & ! <vr> d[vphi']/dr
                         + (m0_values(PSI2,vtheta)*(tbuffer(PSI,dvpdt_p)          & ! <vtheta>/r d[vphi']/dtheta
                         +         tbuffer(PSI,vphi_p)*cottheta(t))               & ! {<vtheta> vphi' cot(theta)}/r
                         +  m0_values(PSI2,vphi) *( buffer(PSI,dvpdp)*csctheta(t) & ! <vphi>/sintheta/r d[vphi']/dphi
				         +         tbuffer(PSI,vr_p)) )                           & ! {<vphi> vr'}/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////////////////////////////
        !/////////////// v' dot grad v' //////////////////
        !--- Radial Component
        If (compute_quantity(vp_grad_vp_r)) Then
            DO_PSI
			    qty(PSI) = tbuffer(PSI,vr_p)*tbuffer(PSI,dvrdr_p) &
				    + tbuffer(PSI,vtheta_p) * ( buffer(PSI,dvrdt_p)-tbuffer(PSI,vtheta_p) )*one_over_r(r)    &
				    + buffer(PSI,vphi_p)*(buffer(PSI,dvrdp)*csctheta(t)-tbuffer(PSI,vphi_p) )*one_over_r(r)  
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !--- Theta Component
        If (compute_quantity(vp_grad_vp_theta)) Then

            DO_PSI
                qty(PSI) =  tbuffer(PSI,vr_p)    * tbuffer(PSI,dvtdr_p)        & ! vr' d[vtheta']/dr
                         + (tbuffer(PSI,vtheta_p)*(tbuffer(PSI,dvtdt_p)        & ! vtheta'/r d[vtheta']/dtheta
                         +         tbuffer(PSI,vr))                            & ! (vtheta' vr' )/r
                         +  tbuffer(PSI,vphi) *( buffer(PSI,dvtdp)*csctheta(t) & ! vphi'/sintheta/r d[vtheta']/dphi
				         -         tbuffer(PSI,vphi )*cottheta(t) ) )          & ! vphi'^2 cot(theta)/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        !--- Phi Component
        If (compute_quantity(vp_grad_vp_phi)) Then

            DO_PSI
                qty(PSI) =  tbuffer(PSI,vr_p  ) *  tbuffer(PSI,dvpdr_p)             & ! vr' d[vphi']/dr
                         + (tbuffer(PSI,vtheta_p)*(tbuffer(PSI,dvpdt_p)             & ! vtheta'/r d[vphi']/dtheta
                         +         tbuffer(PSI,vphi_p)*cottheta(t))               & ! {vtheta' vphi' cot(theta)}/r
                         +  tbuffer(PSI,vphi_p) *( buffer(PSI,dvpdp)*csctheta(t)  & ! vphi'/sintheta/r d[vphi']/dphi
				         +         tbuffer(PSI,vr_p)) )                          & ! {vphi' vr'}/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif
        !//////////////////////////////////////////////////////////////////////////////////////////////
        !/////////////// <v> dot grad <v> //////////////////
        !--- Radial Component
        If (compute_quantity(vm_grad_vm_r)) Then
            DO_PSI
			    qty(PSI) = m0_values(PSI2,vr)*m0_values(PSI2,dvrdr) &
				    + m0_values(PSI2,vtheta) * ( m0_values(PSI2,dvrdt)-m0_values(PSI2,vtheta) )*one_over_r(r)    &
				    + m0_values(PSI2,vphi)*(-m0_values(PSI2,vphi) )*one_over_r(r)  
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !--- Theta Component
        If (compute_quantity(vm_grad_vm_theta)) Then

            DO_PSI
                qty(PSI) =  m0_values(PSI2,vr  ) *  m0_values(PSI2,dvtdr) & ! vr d[vtheta]/dr
                         + (m0_values(PSI2,vtheta)*(m0_values(PSI2,dvtdt) & ! vtheta/r d[vtheta]/dtheta
                         +         m0_values(PSI2,vr))                   & ! (vtheta vr )/r
                         +  m0_values(PSI2,vphi) *(                      & ! NO vphi/sintheta/r d[vtheta]/dphi term
				         -         m0_values(PSI2,vphi )*cottheta(t) ) ) & ! vphi^2 cot(theta)/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        !--- Phi Component
        If (compute_quantity(vm_grad_vm_phi)) Then

            DO_PSI
                qty(PSI) =  m0_values(PSI2,vr_p  ) *  m0_values(PSI2,dvpdr) & ! <vr> d[<vphi>]/dr
                         + (m0_values(PSI2,vtheta_p)*(m0_values(PSI2,dvpdt) & ! <vtheta>/r d[<vphi>]/dtheta
                         +         m0_values(PSI2,vphi)*cottheta(t))        & ! {<vtheta> <vphi> cot(theta)}/r
                         +  m0_values(PSI2,vphi_p) *(                       & ! NO <vphi>/sintheta/r d[<vphi>]/dphi term
				         +         m0_values(PSI2,vr )) )                   & ! {<vphi> <vr>}/r
                         *  one_over_r(r)     
                qty(PSI) = qty(PSI)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        DeAllocate(tbuffer)
    End Subroutine Compute_Inertial_Terms
End Module Diagnostics_Advection
