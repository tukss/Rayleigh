#define DO_PSI Do t = my_theta%min, my_theta%max; Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max; Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

Module Diagnostics_Miscellaneous
    Use Diagnostics_Base
    Implicit None

    ! Place diagnostics that don't clearly belong in another module here.
Contains

    Subroutine Compute_Misc_Diagnostics(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: mypi
        !////////////////////////////////////////////////////////
        ! Diagnostics for verifying output is working
        If (compute_quantity(diagnostic1)) Then
            mypi = acos(-1.0d0)
            DO_PSI
                qty(PSI) = sin(k*2.0d0*mypi/n_phi) &
                    & *(sintheta(t)**2)*radius(r)
            END_DO
            Write(6,*)'Diagnostic1!', my_rank
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(diagnostic2)) Then
            mypi = acos(-1.0d0)
            DO_PSI
                qty(PSI) = sin(k*4.0d0*mypi/n_phi) &
                    & *(sintheta(t)*costheta(t))*radius(r)**2
            END_DO
            Call Add_Quantity(qty)
        Endif
    End Subroutine Compute_Misc_Diagnostics

End Module Diagnostics_Miscellaneous
