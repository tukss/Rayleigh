!Pseudo-code for performing tricubic interpolation in parallel.
!       Serial is easy.  We use the post-processing routine provided by 
!       Kyle Augustson.  There are advantages to either approach.
! Each process in column 0 writes some range in z of the cube to disk
!/////////////////////////////////////////////////
!  Module-wide variables for cubic interpolation
!
Integer :: my_zmin, my_zmax ! The z-bounds I hold
Integer :: delta_z, my_zdisp
Integer :: my_rmin_cube, my_rmax_cube
Integer :: ncube            ! dimension of the cube (ncube**3)

Subroutine Initialize_TriCubic_Interpolation()
    Implicit None
    Integer :: p, leftover
    !First we establish our local z-bounds within the cube.
    
    my_zmin = 1
    my_zmax = 0
    delta_z = nproc1/ncube
    leftover = Mod(ncube,nproc1)
    If (my_column_rank .lt. leftover) delta_z = delta_z+1
    my_zdisp = 0 
    Do p = 0, my_column_rank-1    
        my_zdisp = my_zdisp+delta_z
        if (p .lt. leftover) my_zdisp = my_zdisp+1
    Enddo
    my_zmin = my_zdisp+1
    my_zmax = my_zdisp+delta_z-1   

    ! Next, we need to establish the radial bounds and the theta bounds of the annulus
    ! needed in processor to perform our interpolation.
    Call find_
End Subroutine Initialize_TriCubic_Interpolation
