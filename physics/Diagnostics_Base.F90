
Module Diagnostics_Base
    !//////////////////////////////////////////////////////////
    ! This module holds common variables that may be accessed
    ! by any diagnostics routine.  These variables are primarily
    ! temporary, allocatable arrays and the output menu codes. 
    Use ProblemSize
    Use Spherical_IO
    Use Fields
    Use Math_Constants

    Implicit None

    !/////////////////////////////////////////////////////////
    !  Quantity Codes (some care must be taken to
    Integer, Parameter :: V_r = 1,   V_theta = 2, V_phi = 3
    Integer, Parameter :: Temperature = 4,    Pressure = 5

    Integer, Parameter :: v_sq = 6, kinetic_energy = 7
    Integer, Parameter :: gradt_r = 8, cond_flux_r = 9
    Integer, Parameter :: zonal_ke = 10, merid_ke = 11
    Integer, Parameter :: vol_heating = 12

    Integer, Parameter :: rhoV_r = 13,   rhoV_theta = 14, rhoV_phi = 15
    Integer, Parameter :: thermalE_flux_radial = 16, radial_ke = 17
    Integer, Parameter :: ke_flux_radial = 18, enth_flux_radial = 19
    Integer, Parameter :: buoyancy_work = 20

    Integer, Parameter :: vort_r = 21, vort_theta = 22, vort_phi = 23
    Integer, Parameter :: enstrophy = 24

    !Angular Momentum Transport Diagnostics
    Integer, Parameter :: amom_fluct_r = 25, amom_fluct_theta = 26, &
         amom_dr_r = 27, amom_dr_theta = 28, amom_mean_r = 29, amom_mean_theta = 30

    !Viscous Fluxes
    Integer, Parameter :: visc_flux_r = 31
         
    ! We have some "known" outputs as well that allow us to verify that
    ! the spherical_io interface is functional
    Integer, Parameter :: diagnostic1 = 99, diagnostic2 = 100
    ! We also have some comparison outputs for checking the moments
    Integer, Parameter :: vr2 = 101, vt2 = 102, vp2 = 103
    Integer, Parameter :: vr3 = 104, vt3 = 105, vp3 = 106

    !/////////// Magnetic Outputs.  Start at 200 to organization room for hydro
    Integer, Parameter :: B_r = 201, B_theta = 202, B_phi = 203
    Integer, Parameter :: J_r = 204, J_theta = 205, J_phi = 206
    Integer, Parameter :: B_sq = 207, magnetic_energy=208, zonal_me = 209
    Integer, Parameter :: merid_me = 210, b_r2 = 211, b_theta2 = 212, b_phi2 = 213

    !/////////////////////////// Lorentz Forces ///////////////////////////////
    ! "m" and "< >" denote the azimuthal mean.
    ! "p" and " ' " denote perturbations about the azimuthal mean 

    Integer, Parameter :: j_cross_b_r       = 220 ! radial component of j x B
    Integer, Parameter :: j_cross_b_theta   = 221 !  theta component of j x B
    Integer, Parameter :: j_cross_b_phi     = 222 !    phi component of j x B

    Integer, Parameter :: jp_cross_bm_r     = 223 ! radial component of j' x <B>  
    Integer, Parameter :: jp_cross_bm_theta = 224 !  theta component of j' x <B>
    Integer, Parameter :: jp_cross_bm_phi   = 225 !    phi component of j' x <B>

    Integer, Parameter :: jm_cross_bp_r     = 226 ! radial component of <j> x B'
    Integer, Parameter :: jm_cross_bp_theta = 227 !  theta component of <j> x B'
    Integer, Parameter :: jm_cross_bp_phi   = 228 !    phi component of <j> x B'

    Integer, Parameter :: jm_cross_bm_r     = 229 ! radial component of <j> x <B>
    Integer, Parameter :: jm_cross_bm_theta = 230 !  theta component of <j> x <B>  
    Integer, Parameter :: jm_cross_bm_phi   = 231 !    phi component of <j> x <B> 

    Integer, Parameter :: jp_cross_bp_r     = 232 ! radial component of j' x B'  
    Integer, Parameter :: jp_cross_bp_theta = 233 !  theta component of j' x B'
    Integer, Parameter :: jp_cross_bp_phi   = 234 !    phi component of j' x B'


    !///////////// Induction Terms ///////////////////////////
    ! "m," "< >", "  '  ", and "p" retain the same meaning as above
    !  These will be more involved
    ! Integer, Parameter :: dBr_p_vm_Bp, dBtheta_p_vm_Bp


    !///////////////////////////////////
    Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output
    Real*8, Allocatable :: tmp1(:,:,:)  ! A work array
    Real*8, Allocatable :: rweights(:), tweights(:)

    !//////////////////////////////////
    Real*8, Allocatable :: ell0_values(:,:), m0_values(:,:,:)		

    Type(SphericalBuffer), public :: add_fields

Contains

    Subroutine Generate_Diagnostic_Labels()
        ! Define labels for our quantity codes
        Write(6,*)'A line of code.'
        !Call Load_Label(v_r,'V_r')
        !Call Load_Label(v_theta,'V_theta')
        !Call Load_Label(v_phi, 'V_phi')
    End Subroutine Generate_Diagnostic_Labels

End Module Diagnostics_Base
