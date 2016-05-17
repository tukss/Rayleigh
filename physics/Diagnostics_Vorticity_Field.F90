#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_CURRENT_DENSITY
!               This module computes the components of del x v and enstrophy. 
!               Zonal means and fluctuations about those means are 
!               also computed (if desired).
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Vorticity_Field
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_Vorticity_Field(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        !/////////////////////////////////////////
        ! 1. terms involving radial vorticity
        If (compute_quantity(vort_r)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( buffer(PSI,dvpdt)- &
                           csctheta(t)*buffer(PSI,dvtdp)  +&
                           cottheta(t)*buffer(PSI,vphi) )
            END_DO
            Call Add_Quantity(qty)
        Endif	

        If (compute_quantity(vortp_r)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( fbuffer(PSI,dvpdt)- &
                           csctheta(t)*fbuffer(PSI,dvtdp)  +&
                           cottheta(t)*fbuffer(PSI,vphi) )
            END_DO
            Call Add_Quantity(qty)
        Endif	

        If (compute_quantity(vortm_r)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( m0_values(PSI2,dvpdt)- &
                           csctheta(t)*m0_values(PSI2,dvtdp)  +&
                           cottheta(t)*m0_values(PSI2,vphi) )
            END_DO
            Call Add_Quantity(qty)
        Endif	

        !/////////////////////////////////////////////////
        ! 2. terms involving theta vorticity
        If (compute_quantity(vort_theta)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( csctheta(t)*buffer(PSI,dvrdp) - &
                           buffer(PSI,vphi) )-buffer(PSI,dvpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(vortp_theta)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( csctheta(t)*fbuffer(PSI,dvrdp) - &
                           fbuffer(PSI,vphi) )-fbuffer(PSI,dvpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(vortm_theta)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( csctheta(t)*m0_values(PSI2,dvrdp) - &
                           m0_values(PSI2,vphi) )-m0_values(PSI2,dvpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////
        ! 3. terms involving phi vorticity

    End Subroutine Compute_Vorticity_Field

End Module Diagnostics_Vorticity_Field
