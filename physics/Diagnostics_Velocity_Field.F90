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
        
    End Subroutine Compute_Velocity_Components

End Module Diagnostics_Velocity_Field
