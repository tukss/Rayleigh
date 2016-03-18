#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

Module Diagnostics_Lorentz_Forces

    Use Diagnostics_Base
Contains

    Subroutine Compute_Lorentz_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8, Allocatable :: tbuffer(:,:,:,:)
        Integer :: jrad_p = 4,jtheta_p = 5,jphi_p = 6, brad_p = 1,btheta_p = 2, bphi_p = 3 
        Allocate(tbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:6))
        !The azimuthal means of j and b have already been computed, but we need 
        ! the perturbations as well - possibly have a flag for computing
        ! if needed

        DO_PSI
            tbuffer(PSI,jrad_p) = buffer(PSI,jr) - m0_values(PSI2,jr) 
        END_DO

        DO_PSI
            tbuffer(PSI,jtheta_p) = buffer(PSI,jtheta) - m0_values(PSI2,jtheta) 
        END_DO

        DO_PSI
            tbuffer(PSI,jphi_p) = buffer(PSI,jphi) - m0_values(PSI2,jphi) 
        END_DO        

        DO_PSI
            tbuffer(PSI,brad_p) = buffer(PSI,br) - m0_values(PSI2,br) 
        END_DO

        DO_PSI
            tbuffer(PSI,btheta_p) = buffer(PSI,btheta) - m0_values(PSI2,btheta) 
        END_DO

        DO_PSI
            tbuffer(PSI,bphi_p) = buffer(PSI,bphi) - m0_values(PSI2,bphi) 
        END_DO     



        ! Full JxB terms
        If (compute_quantity(j_cross_b_r)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,jtheta)*buffer(PSI,bphi)- &
                         & buffer(PSI,btheta)*buffer(PSI,jphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(j_cross_b_theta)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,br)*buffer(PSI,jphi)- &
                         & buffer(PSI,jr)*buffer(PSI,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(j_cross_b_phi)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,jr)*buffer(PSI,btheta)- &
                           & buffer(PSI,br)*buffer(PSI,jtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////////
        !                   <J> x <B> terms

        If (compute_quantity(jm_cross_bm_r)) Then
            DO_PSI2
                qty(1:n_phi,PSI2) = m0_values(PSI2,jtheta)*m0_values(PSI2,bphi)- &
                                  & m0_values(PSI2,btheta)*m0_values(PSI2,jphi)
            END_DO2
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jm_cross_bm_theta)) Then
            DO_PSI2
                qty(1:n_phi,PSI2) = m0_values(PSI2,br)*m0_values(PSI2,jphi)- &
                                  & m0_values(PSI2,jr)*m0_values(PSI2,bphi)
            END_DO2
            Call Add_Quantity(qty)
        Endif       

        If (compute_quantity(jm_cross_bm_phi)) Then
            DO_PSI2
                qty(1:n_phi,PSI2) = m0_values(PSI2,jr)*m0_values(PSI2,btheta)- &
                                  & m0_values(PSI2,br)*m0_values(PSI2,jtheta)
            END_DO2
            Call Add_Quantity(qty)
        Endif


        !////////////////////////////////////////////////////////
        !  J' X B' terms
        If (compute_quantity(jp_cross_bp_r)) Then
            DO_PSI
                qty(PSI) = tbuffer(PSI,jtheta_p)*tbuffer(PSI,bphi_p)- &
                         & tbuffer(PSI,btheta_p)*tbuffer(PSI,jphi_p)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jp_cross_bp_theta)) Then
            DO_PSI
                qty(PSI) = tbuffer(PSI,brad_p)*tbuffer(PSI,jphi_p)- &
                         & tbuffer(PSI,jrad_p)*tbuffer(PSI,bphi_p)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jp_cross_bp_phi)) Then
            DO_PSI
                qty(PSI) = tbuffer(PSI,jrad_p)*tbuffer(PSI,btheta_p)- &
                           & tbuffer(PSI,brad_p)*tbuffer(PSI,jtheta_p)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !//////////////////////////////////////////////////////////
        !               J' x <B> terms
        If (compute_quantity(jp_cross_bm_r)) Then
            DO_PSI
                qty(PSI) = tbuffer(PSI,jtheta_p)*m0_values(PSI2,bphi)- &
                         & m0_values(PSI2,btheta)*tbuffer(PSI,jphi_p)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jp_cross_bm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,br)*tbuffer(PSI,jphi_p)- &
                         & tbuffer(PSI,jrad_p)*m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jp_cross_bm_phi)) Then
            DO_PSI
                qty(PSI) = tbuffer(PSI,jrad_p)*m0_values(PSI2,btheta)- &
                           & m0_values(PSI2,br)*tbuffer(PSI,jtheta_p)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////////////////////////
        !               <J> x B' terms
        !////////////////////////////////////////////////////////
        !  J' X B' terms
        If (compute_quantity(jm_cross_bp_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,jtheta)*tbuffer(PSI,bphi_p)- &
                         & tbuffer(PSI,btheta_p)*m0_values(PSI2,jphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jm_cross_bp_theta)) Then
            DO_PSI
                qty(PSI) = tbuffer(PSI,brad_p)*m0_values(PSI2,jphi)- &
                         & m0_values(PSI2,jr)*tbuffer(PSI,bphi_p)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jm_cross_bp_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,jr)*tbuffer(PSI,btheta_p)- &
                           & tbuffer(PSI,brad_p)*m0_values(PSI2,jtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif


        DeAllocate(tbuffer)
    End Subroutine Compute_Lorentz_Forces


End Module Diagnostics_Lorentz_Forces
