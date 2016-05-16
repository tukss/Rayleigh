Module Diagnostics_Linear_Forces
    Use Diagnostics_Base
    Implicit None

Contains

    Subroutine Compute_Linear_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Call Compute_Buoyancy_Force(buffer)
        Call Compute_Coriolis_Force(buffer)
        Call Compute_Pressure_Force(buffer)
        Call Compute_Viscous_Force(buffer)
    End Subroutine Compute_Linear_Forces

    Subroutine Compute_Buoyancy_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)

    End Subroutine Compute_Buoyancy_Force

    Subroutine Compute_Coriolis_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)

    End Subroutine Compute_Coriolis_Force

    Subroutine Compute_Pressure_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)

    End Subroutine Compute_Pressure_Force

    Subroutine Compute_Viscous_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)

    End Subroutine Compute_Viscous_Force

End Module Diagnostics_Linear_Forces
