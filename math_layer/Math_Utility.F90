Module Math_Utility
    ! This module contains potentially useful stand alone routines
Contains

    Subroutine tanh_profile(x,y,flip)
        Real*8, Intent(In) :: x(1:)
        Real*8, Intent(InOut) :: y(1:)
        Logical, Intent(In), Optional :: flip
        Real*8 :: flip_factor = 1.0d0
        Integer :: xsize(1), nx
        xsize = size(x)
        nx = xsize(1)
        if (present(flip)) Then
            if(flip) flip_factor = -1.0d0
        endif

        Do i = 1, nx
            y(i) = 0.5d0*(1.0d0+flip_factor*tanh(x(i)))
        Enddo

    End Subroutine tanh_profile

End Module Math_Utility
