#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
!////////////////////// Diagnostics Energy Flux ///////////////////////
!
!       This module handles calculation of terms comprising 
!       the RADIAL energy flux.  These are:
!       1.  Enthalpy Flux
!       2.  Kinetic Energy Flux
!       3.  Conductive Heat Flux
!       4.  Poynting Flux
!       5.  Viscous flux of KE
!
!///////////////////////////////////////////////////////////////////
Module Diagnostics_Energy_Flux
    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
    Implicit None

Contains

    Subroutine Compute_Energy_Flux(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        !First, the radial viscous flux of energy
        If (compute_quantity(visc_flux_r)) Then

            !Radial contribution (mod rho*nu)
            DO_PSI		
                qty(PSI) = buffer(PSI,vr)*ref%dlnrho(r)/3.0d0+buffer(PSI,dvrdr)
                qty(PSI) = qty(PSI)*buffer(PSI,vr)*2.0d0
            END_DO

            !Theta contribution (mod rho*nu)
            DO_PSI		
                tmp1(PSI) = (2.0d0/3.0d0)*buffer(PSI,vr)*ref%dlnrho(r)
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,dvtdr)-buffer(PSI,vtheta)/radius(r)
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,dvrdt)/radius(r)
                tmp1(PSI) = tmp1(PSI)*buffer(PSI,vtheta)
            END_DO            


            DO_PSI		
                qty(PSI) = qty(PSI)+tmp1(PSI)
            END_DO

            !phi contribution (mod rho*nu)
            DO_PSI			
                tmp1(PSI) = (2.0d0/3.0d0)*buffer(PSI,vr)*ref%dlnrho(r)
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,dvpdr)-buffer(PSI,vphi)/radius(r)
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,dvrdp)/radius(r)/sintheta(t)
                tmp1(PSI) = tmp1(PSI)*buffer(PSI,vphi)
            END_DO             

            DO_PSI		
                qty(PSI) = qty(PSI)+tmp1(PSI)
            END_DO

            !Multiply by rho and nu
            DO_PSI			
                qty(PSI) = qty(PSI)*nu(r)*ref%density(r)                            
            END_DO
            Call Add_Quantity(qty)
        Endif

        !Conductive Flux
        If (compute_quantity(cond_flux_r)) Then
            DO_PSI
                qty(PSI) = -ref%density(r) &
                  & *ref%temperature(r)*kappa(r)*buffer(PSI,dtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif	
        


        !/////////// Radial component of ExB
        !This is an unnormalized Poynting flux
        If (magnetism) Then
            If (compute_quantity(ecrossb_r)) Then
                ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
                DO_PSI
                    tmp1(PSI) = buffer(PSI,vphi)*buffer(PSI,br)- &
                                & buffer(PSI,vr)*buffer(PSI,bphi)
                    tmp1(PSI) = tmp1(PSI)-eta(r)*buffer(PSI,jtheta)
                END_DO
                DO_PSI
                    qty(PSI) = tmp1(PSI)*buffer(PSI,bphi) ! E_theta B_phi
                END_DO

                !Next, compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
                DO_PSI
                    tmp1(PSI) = buffer(PSI,vr)*buffer(PSI,btheta)- &
                                & buffer(PSI,vtheta)*buffer(PSI,br)
                    tmp1(PSI) = tmp1(PSI)-eta(r)*buffer(PSI,jphi)
                END_DO

                DO_PSI
                    qty(PSI) = qty(PSI)-tmp1(PSI)*buffer(PSI,btheta) ! E_phi B_theta
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        

    End Subroutine Compute_Energy_Flux

End Module Diagnostics_Energy_Flux
