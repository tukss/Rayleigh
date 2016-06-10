#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

Module Diagnostics_Lorentz_Forces

    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_Lorentz_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        ! Full JxB terms
        If (compute_quantity(j_cross_b_r)) Then
            DO_PSI
                qty(PSI) = (buffer(PSI,jtheta)*buffer(PSI,bphi)- &
                         & buffer(PSI,btheta)*buffer(PSI,jphi) ) *ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(j_cross_b_theta)) Then
            DO_PSI
                qty(PSI) = ( buffer(PSI,br)*buffer(PSI,jphi)- &
                         & buffer(PSI,jr)*buffer(PSI,bphi) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(j_cross_b_phi)) Then
            DO_PSI
                qty(PSI) = ( buffer(PSI,jr)*buffer(PSI,btheta)- &
                           & buffer(PSI,br)*buffer(PSI,jtheta) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////////
        !                   <J> x <B> terms

        If (compute_quantity(jm_cross_bm_r)) Then
            DO_PSI2
                qty(1:n_phi,PSI2) = ( m0_values(PSI2,jtheta)*m0_values(PSI2,bphi)- &
                                  & m0_values(PSI2,btheta)*m0_values(PSI2,jphi) )*ref%Lorentz_Coeff
            END_DO2
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jm_cross_bm_theta)) Then
            DO_PSI2
                qty(1:n_phi,PSI2) = ( m0_values(PSI2,br)*m0_values(PSI2,jphi)- &
                                  & m0_values(PSI2,jr)*m0_values(PSI2,bphi) )*ref%Lorentz_Coeff
            END_DO2
            Call Add_Quantity(qty)
        Endif       

        If (compute_quantity(jm_cross_bm_phi)) Then
            DO_PSI2
                qty(1:n_phi,PSI2) = ( m0_values(PSI2,jr)*m0_values(PSI2,btheta)- &
                                  & m0_values(PSI2,br)*m0_values(PSI2,jtheta) )*ref%Lorentz_Coeff
            END_DO2
            Call Add_Quantity(qty)
        Endif


        !////////////////////////////////////////////////////////
        !  J' X B' terms
        If (compute_quantity(jp_cross_bp_r)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,jtheta)*fbuffer(PSI,bphi)- &
                         & fbuffer(PSI,btheta)*fbuffer(PSI,jphi) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jp_cross_bp_theta)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,br)*fbuffer(PSI,jphi)- &
                         & fbuffer(PSI,jr)*fbuffer(PSI,bphi) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jp_cross_bp_phi)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,jr)*fbuffer(PSI,btheta)- &
                           & fbuffer(PSI,br)*fbuffer(PSI,jtheta) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif


        !//////////////////////////////////////////////////////////
        !               J' x <B> terms
        If (compute_quantity(jp_cross_bm_r)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,jtheta)*m0_values(PSI2,bphi)- &
                         & m0_values(PSI2,btheta)*fbuffer(PSI,jphi) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jp_cross_bm_theta)) Then
            DO_PSI
                qty(PSI) = ( m0_values(PSI2,br)*fbuffer(PSI,jphi)- &
                         & fbuffer(PSI,jr)*m0_values(PSI2,bphi) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jp_cross_bm_phi)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,jr)*m0_values(PSI2,btheta)- &
                           & m0_values(PSI2,br)*fbuffer(PSI,jtheta) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////////////////////////
        !               <J> x B' terms
        !////////////////////////////////////////////////////////
        !  J' X B' terms
        If (compute_quantity(jm_cross_bp_r)) Then
            DO_PSI
                qty(PSI) = ( m0_values(PSI2,jtheta)*fbuffer(PSI,bphi)- &
                         & fbuffer(PSI,btheta)*m0_values(PSI2,jphi) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jm_cross_bp_theta)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,br)*m0_values(PSI2,jphi)- &
                         & m0_values(PSI2,jr)*fbuffer(PSI,bphi) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(jm_cross_bp_phi)) Then
            DO_PSI
                qty(PSI) = ( m0_values(PSI2,jr)*fbuffer(PSI,btheta)- &
                           & fbuffer(PSI,br)*m0_values(PSI2,jtheta) )*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Lorentz_Forces


End Module Diagnostics_Lorentz_Forces
