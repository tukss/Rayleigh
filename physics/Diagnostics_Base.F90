#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
Module Diagnostics_Base
    !//////////////////////////////////////////////////////////
    ! This module holds common variables that may be accessed
    ! by any diagnostics routine.  These variables are primarily
    ! temporary, allocatable arrays and the output menu codes. 
    Use ProblemSize
    Use Spherical_IO
    Use Fields
    Use Math_Constants
    Use ReferenceState
    Use TransportCoefficients

    Implicit None

    !/////////////////////////////////////////////////////////
    !  Quantity Codes (some care must be taken to

    ! Velocity field components.
    ! Fluctuations (denoted by "p") and azimuthal means
    ! (denoted by "m") may also be output.
    Integer, Parameter :: voffset = 0

    Integer, Parameter :: v_r  = voffset+1         ! Radial velocity
    Integer, Parameter :: vp_r = voffset+2    
    Integer, Parameter :: vm_r = voffset+3    

    Integer, Parameter :: dv_r_dr    = voffset+4   ! d{v_r}/dr
    Integer, Parameter :: dvp_r_dr   = voffset+5    
    Integer, Parameter :: dvm_r_dr   = voffset+6    

    Integer, Parameter :: dv_r_dt    = voffset+7   ! d{v_r}/dtheta
    Integer, Parameter :: dvp_r_dt   = voffset+8    
    Integer, Parameter :: dvm_r_dt   = voffset+9    

    Integer, Parameter :: dv_r_dp    = voffset+10  ! d{v_r}/dphi
    Integer, Parameter :: dvp_r_dp   = voffset+11    

    Integer, Parameter :: dv_r_dtr   = voffset+12  ! (1/r)d{v_r}/dtheta
    Integer, Parameter :: dvp_r_dtr  = voffset+13    
    Integer, Parameter :: dvm_r_dtr  = voffset+14    

    Integer, Parameter :: dv_r_dprs  = voffset+15  ! (1/r 1/sin(theta))d{v_r}/dphi
    Integer, Parameter :: dvp_r_dprs = voffset+16   

    Integer, Parameter :: v_theta  = voffset+17        ! Theta velocity
    Integer, Parameter :: vp_theta = voffset+18    
    Integer, Parameter :: vm_theta = voffset+19    

    Integer, Parameter :: dv_theta_dr    = voffset+20  ! d{v_theta}/dr
    Integer, Parameter :: dvp_theta_dr   = voffset+21    
    Integer, Parameter :: dvm_theta_dr   = voffset+22    

    Integer, Parameter :: dv_theta_dt    = voffset+23  ! d{v_theta}/dtheta
    Integer, Parameter :: dvp_theta_dt   = voffset+24    
    Integer, Parameter :: dvm_theta_dt   = voffset+25    

    Integer, Parameter :: dv_theta_dp    = voffset+26  ! d{v_theta}/dphi
    Integer, Parameter :: dvp_theta_dp   = voffset+27    

    Integer, Parameter :: dv_theta_dtr   = voffset+28  ! (1/r)d{v_theta}/dtheta
    Integer, Parameter :: dvp_theta_dtr  = voffset+29    
    Integer, Parameter :: dvm_theta_dtr  = voffset+30    

    Integer, Parameter :: dv_theta_dprs  = voffset+31  ! (1/r 1/sin(theta))d{v_theta}/dphi
    Integer, Parameter :: dvp_theta_dprs = voffset+32   


    Integer, Parameter :: v_phi  = voffset+33        ! Phi velocity
    Integer, Parameter :: vp_phi = voffset+34    
    Integer, Parameter :: vm_phi = voffset+35    

    Integer, Parameter :: dv_phi_dr    = voffset+36  ! d{v_phi}/dr
    Integer, Parameter :: dvp_phi_dr   = voffset+37    
    Integer, Parameter :: dvm_phi_dr   = voffset+38    

    Integer, Parameter :: dv_phi_dt    = voffset+39  ! d{v_phi}/dtheta
    Integer, Parameter :: dvp_phi_dt   = voffset+40    
    Integer, Parameter :: dvm_phi_dt   = voffset+41    

    Integer, Parameter :: dv_phi_dp    = voffset+42  ! d{v_phi}/dphi
    Integer, Parameter :: dvp_phi_dp   = voffset+43    

    Integer, Parameter :: dv_phi_dtr   = voffset+44  ! (1/r)d{v_phi}/dtheta
    Integer, Parameter :: dvp_phi_dtr  = voffset+45    
    Integer, Parameter :: dvm_phi_dtr  = voffset+46    

    Integer, Parameter :: dv_phi_dprs  = voffset+47  ! (1/r 1/sin(theta))d{v_phi}/dphi
    Integer, Parameter :: dvp_phi_dprs = voffset+48   


    Integer, Parameter :: rhov_r      = voffset+49
    Integer, Parameter :: rhovp_r     = voffset+50
    Integer, Parameter :: rhovm_r     = voffset+51

    Integer, Parameter :: rhov_theta  = voffset+52
    Integer, Parameter :: rhovp_theta = voffset+53
    Integer, Parameter :: rhovm_theta = voffset+54

    Integer, Parameter :: rhov_phi    = voffset+55
    Integer, Parameter :: rhovp_phi   = voffset+56
    Integer, Parameter :: rhovm_phi   = voffset+57

    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Vorticity Outputs
    Integer, Parameter :: vort_off = 57

    Integer, Parameter :: vort_r  = vort_off+1      ! Radial
    Integer, Parameter :: vortp_r = vort_off+2    
    Integer, Parameter :: vortm_r = vort_off+3 

    Integer, Parameter :: vort_theta  = vort_off+4  ! Theta Vorticity
    Integer, Parameter :: vortp_theta = vort_off+5    
    Integer, Parameter :: vortm_theta = vort_off+6 

    Integer, Parameter :: vort_phi  = vort_off+7    ! Phi Vorticity
    Integer, Parameter :: vortp_phi = vort_off+8    
    Integer, Parameter :: vortm_phi = vort_off+9 

    Integer, Parameter :: enstrophy   = vort_off+10 ! Enstrophy
    Integer, Parameter :: enstrophypm = vort_off+11
    Integer, Parameter :: enstrophymm = vort_off+12
    Integer, Parameter :: enstropypp  = vort_off+13

    Integer, Parameter :: tpoffset = 70

    Integer, Parameter :: Temperature = 4,    Pressure = 5

    Integer, Parameter :: v_sq = 6, kinetic_energy = 7
    Integer, Parameter :: gradt_r = 8, cond_flux_r = 9
    Integer, Parameter :: zonal_ke = 10, merid_ke = 11
    Integer, Parameter :: vol_heating = 12


    Integer, Parameter :: thermalE_flux_radial = 16, radial_ke = 17
    Integer, Parameter :: ke_flux_radial = 18, enth_flux_radial = 19
    Integer, Parameter :: buoyancy_work = 20



    !Angular Momentum Transport Diagnostics
    Integer, Parameter :: amom_fluct_r = 25, amom_fluct_theta = 26, &
         amom_dr_r = 27, amom_dr_theta = 28, amom_mean_r = 29, amom_mean_theta = 30

    !Viscous Fluxes
    Integer, Parameter :: visc_flux_r = 31

    !////////////////////////  Advection Terms ////////////////////
    ! Reynolds decomposition about the azimuthal mean may also be output
    ! "m" and "< >" denote the azimuthal mean.
    ! "p" and " ' " denote perturbations about the azimuthal mean 
    ! NOTE:  ADVECTION TERMS ARE SCALED BY DENSITY (so that they represent a force density)

    Integer, Parameter :: vgv = 40  ! Output offset for advection terms  
    Integer, Parameter :: v_grad_v_r       = vgv+1 ! radial component of v dot grad v
    Integer, Parameter :: v_grad_v_theta   = vgv+2 !  theta component of v dot grad v
    Integer, Parameter :: v_grad_v_phi     = vgv+3 !    phi component of v dot grad v

    Integer, Parameter :: vp_grad_vm_r     = vgv+4 ! radial component of v' dot grad <v>
    Integer, Parameter :: vp_grad_vm_theta = vgv+5 !  theta component of v' dot grad <v>
    Integer, Parameter :: vp_grad_vm_phi   = vgv+6 !    phi component of v' dot grad <v>

    Integer, Parameter :: vm_grad_vp_r     = vgv+7 ! radial component of <v> dot grad v'
    Integer, Parameter :: vm_grad_vp_theta = vgv+8 !  theta component of <v> dot grad v'
    Integer, Parameter :: vm_grad_vp_phi   = vgv+9 !    phi component of <v> dot grad v'

    Integer, Parameter :: vp_grad_vp_r     = vgv+10 ! radial component of v' dot grad v'
    Integer, Parameter :: vp_grad_vp_theta = vgv+11 !  theta component of v' dot grad v'
    Integer, Parameter :: vp_grad_vp_phi   = vgv+12 !    phi component of v' dot grad v'

    Integer, Parameter :: vm_grad_vm_r     = vgv+13 ! radial component of <v> dot grad <v>
    Integer, Parameter :: vm_grad_vm_theta = vgv+14 !  theta component of <v> dot grad <v>
    Integer, Parameter :: vm_grad_vm_phi   = vgv+15 !    phi component of <v> dot grad <v>
         
    ! We have some "known" outputs as well that allow us to verify that
    ! the spherical_io interface is functional
    Integer, Parameter :: diagnostic1 = 99, diagnostic2 = 100
    ! We also have some comparison outputs for checking the moments
    Integer, Parameter :: vr2 = 101, vt2 = 102, vp2 = 103
    Integer, Parameter :: vr3 = 104, vt3 = 105, vp3 = 106

    !/////////// Magnetic Outputs.  Start at 200 to organization room for hydro


    !//////////////////////////////////////////////////
    !               Magnetic field components.
    ! Fluctuations (denoted by "p") and azimuthal means
    ! (denoted by "m") may also be output.
    Integer, Parameter :: boffset = 200

    Integer, Parameter :: b_r  = boffset+1             ! Radial Bfield
    Integer, Parameter :: bp_r = boffset+2    
    Integer, Parameter :: bm_r = boffset+3    

    Integer, Parameter :: db_r_dr    = boffset+4       ! d{b_r}/dr
    Integer, Parameter :: dbp_r_dr   = boffset+5    
    Integer, Parameter :: dbm_r_dr   = boffset+6    

    Integer, Parameter :: db_r_dt    = boffset+7       ! d{b_r}/dtheta
    Integer, Parameter :: dbp_r_dt   = boffset+8    
    Integer, Parameter :: dbm_r_dt   = boffset+9    

    Integer, Parameter :: db_r_dp    = boffset+10      ! d{b_r}/dphi
    Integer, Parameter :: dbp_r_dp   = boffset+11    

    Integer, Parameter :: db_r_dtr   = boffset+12      ! (1/r)d{b_r}/dtheta
    Integer, Parameter :: dbp_r_dtr  = boffset+13    
    Integer, Parameter :: dbm_r_dtr  = boffset+14    

    Integer, Parameter :: db_r_dprs  = boffset+15      ! (1/r 1/sin(theta))d{b_r}/dphi
    Integer, Parameter :: dbp_r_dprs = boffset+16   

    Integer, Parameter :: b_theta  = boffset+17        ! Theta Bfield
    Integer, Parameter :: bp_theta = boffset+18    
    Integer, Parameter :: bm_theta = boffset+19    

    Integer, Parameter :: db_theta_dr    = boffset+20  ! d{b_theta}/dr
    Integer, Parameter :: dbp_theta_dr   = boffset+21    
    Integer, Parameter :: dbm_theta_dr   = boffset+22    

    Integer, Parameter :: db_theta_dt    = boffset+23  ! d{b_theta}/dtheta
    Integer, Parameter :: dbp_theta_dt   = boffset+24    
    Integer, Parameter :: dbm_theta_dt   = boffset+25    

    Integer, Parameter :: db_theta_dp    = boffset+26  ! d{b_theta}/dphi
    Integer, Parameter :: dbp_theta_dp   = boffset+27    

    Integer, Parameter :: db_theta_dtr   = boffset+28  ! (1/r)d{b_theta}/dtheta
    Integer, Parameter :: dbp_theta_dtr  = boffset+29    
    Integer, Parameter :: dbm_theta_dtr  = boffset+30    

    Integer, Parameter :: db_theta_dprs  = boffset+31  ! (1/r 1/sin(theta))d{b_theta}/dphi
    Integer, Parameter :: dbp_theta_dprs = boffset+32   

    Integer, Parameter :: b_phi  = boffset+33          ! Phi Bfield
    Integer, Parameter :: bp_phi = boffset+34    
    Integer, Parameter :: bm_phi = boffset+35    

    Integer, Parameter :: db_phi_dr    = boffset+36    ! d{b_phi}/dr
    Integer, Parameter :: dbp_phi_dr   = boffset+37    
    Integer, Parameter :: dbm_phi_dr   = boffset+38    

    Integer, Parameter :: db_phi_dt    = boffset+39    ! d{b_phi}/dtheta
    Integer, Parameter :: dbp_phi_dt   = boffset+40    
    Integer, Parameter :: dbm_phi_dt   = boffset+41    

    Integer, Parameter :: db_phi_dp    = boffset+42    ! d{b_phi}/dphi
    Integer, Parameter :: dbp_phi_dp   = boffset+43    

    Integer, Parameter :: db_phi_dtr   = boffset+44    ! (1/r)d{b_phi}/dtheta
    Integer, Parameter :: dbp_phi_dtr  = boffset+45    
    Integer, Parameter :: dbm_phi_dtr  = boffset+46    

    Integer, Parameter :: db_phi_dprs  = boffset+47    ! (1/r 1/sin(theta))d{b_phi}/dphi
    Integer, Parameter :: dbp_phi_dprs = boffset+48   


    !///////////////////////////////////////////////////
    !       Current Density Outputs (Including Ohmic Heating)
    Integer, Parameter :: joffset = 250

    Integer, Parameter :: j_r  = joffset+1      ! Radial Current Density
    Integer, Parameter :: jp_r = joffset+2    
    Integer, Parameter :: jm_r = joffset+3 

    Integer, Parameter :: j_theta  = joffset+4  ! Theta Current Density
    Integer, Parameter :: jp_theta = joffset+5    
    Integer, Parameter :: jm_theta = joffset+6 

    Integer, Parameter :: j_phi  = joffset+7    ! Phi Current Density
    Integer, Parameter :: jp_phi = joffset+8    
    Integer, Parameter :: jm_phi = joffset+9 

    Integer, Parameter :: j_r_sq      = joffset+10 ! (j_r)^2
    Integer, Parameter :: jp_r_sq     = joffset+11 ! (jp_r)^2
    Integer, Parameter :: j_theta_sq  = joffset+12 ! (j_theta)^2
    Integer, Parameter :: jp_theta_sq = joffset+13 ! (jp_theta)^2
    Integer, Parameter :: j_phi_sq    = joffset+14 ! (j_theta)^2
    Integer, Parameter :: jp_phi_sq   = joffset+15 ! (jp_theta)^2
    Integer, Parameter :: j_sq        = joffset+16 ! j dot j
    Integer, Parameter :: jp_sq       = joffset+17 ! j' dot j'
    Integer, Parameter :: ohmic_heat    = joffset+18 ! eta{  j  dot  j}
    Integer, Parameter :: ohmic_heat_pp = joffset+19 ! eta{  j' dot  j'}
    Integer, Parameter :: ohmic_heat_mm = joffset+20 ! eta{ <j> dot <j>}


    Integer, Parameter :: B_sq = 207, magnetic_energy=208, zonal_me = 209
    Integer, Parameter :: merid_me = 210, b_r2 = 211, b_theta2 = 212, b_phi2 = 213

    !/////////////////////////// Lorentz Forces ///////////////////////////////
    ! "m" and "< >" denote the azimuthal mean.
    ! "p" and " ' " denote perturbations about the azimuthal mean 
    Integer, Parameter :: loff = 220
    Integer, Parameter :: j_cross_b_r       = loff+1  ! radial component of j x B
    Integer, Parameter :: j_cross_b_theta   = loff+2  !  theta component of j x B
    Integer, Parameter :: j_cross_b_phi     = loff+3  !    phi component of j x B

    Integer, Parameter :: jp_cross_bm_r     = loff+4  ! radial component of j' x <B>  
    Integer, Parameter :: jp_cross_bm_theta = loff+5  !  theta component of j' x <B>
    Integer, Parameter :: jp_cross_bm_phi   = loff+6  !    phi component of j' x <B>

    Integer, Parameter :: jm_cross_bp_r     = loff+7  ! radial component of <j> x B'
    Integer, Parameter :: jm_cross_bp_theta = loff+8  !  theta component of <j> x B'
    Integer, Parameter :: jm_cross_bp_phi   = loff+9  !    phi component of <j> x B'

    Integer, Parameter :: jm_cross_bm_r     = loff+10 ! radial component of <j> x <B>
    Integer, Parameter :: jm_cross_bm_theta = loff+11 !  theta component of <j> x <B>  
    Integer, Parameter :: jm_cross_bm_phi   = loff+12 !    phi component of <j> x <B> 

    Integer, Parameter :: jp_cross_bp_r     = loff+13 ! radial component of j' x B'  
    Integer, Parameter :: jp_cross_bp_theta = loff+14 !  theta component of j' x B'
    Integer, Parameter :: jp_cross_bp_phi   = loff+15 !    phi component of j' x B'



    !////////////////////////////// Induction Terms ///////////////////////////
    ! "m," "< >", "  '  ", and "p" retain the same meaning as above

    Integer, Parameter :: indoff = 250
    !--------------- Terms involving v x B (full)
    Integer, Parameter :: induction_shear_r          = indoff+1  ! radial component of {B dot grad v}
    Integer, Parameter :: induction_comp_r           = indoff+2  ! radial component of -{div dot v}B
    Integer, Parameter :: induction_advec_r          = indoff+3  ! radial component of -{v dot grad B}
    Integer, Parameter :: induction_r                = indoff+4  ! radial component of del cros {v x B}
    Integer, Parameter :: induction_diff_r           = indoff+5  ! radial component of -del x (eta {del x B})

    Integer, Parameter :: induction_shear_theta      = indoff+6  ! theta component of {B dot grad v}
    Integer, Parameter :: induction_comp_theta       = indoff+7  ! theta component of -{div dot v}B
    Integer, Parameter :: induction_advec_theta      = indoff+8  ! theta component of -{v dot grad B}
    Integer, Parameter :: induction_theta            = indoff+9  ! theta component of del cros {v x B}
    Integer, Parameter :: induction_diff_theta       = indoff+10 ! theta component of -del x (eta {del x B})

    Integer, Parameter :: induction_shear_phi        = indoff+11 ! phi component of {B dot grad v}
    Integer, Parameter :: induction_comp_phi         = indoff+12 ! phi component of -{div dot v}B
    Integer, Parameter :: induction_advec_phi        = indoff+13 ! phi component of -{v dot grad B}
    Integer, Parameter :: induction_phi              = indoff+14 ! phi component of del cros {v x B}
    Integer, Parameter :: induction_diff_phi         = indoff+15 ! phi component of -del x (eta {del x B})

    !--------------- Terms involving <v> x <B> 
    Integer, Parameter :: induction_shear_vmbm_r     = indoff+16 ! radial component of {<B> dot grad <v>}
    Integer, Parameter :: induction_comp_vmbm_r      = indoff+17 ! radial component of -{div dot <v>}<B>
    Integer, Parameter :: induction_advec_vmbm_r     = indoff+18 ! radial component of -{<v> dot grad <B>}
    Integer, Parameter :: induction_vmbm_r           = indoff+19 ! radial component of del cros {<v> x <B>}
    Integer, Parameter :: induction_diff_bm_r        = indoff+20 ! radial component of -del x (eta {del x <B>})

    Integer, Parameter :: induction_shear_vmbm_theta = indoff+21 ! theta component of {<B> dot grad <v>}
    Integer, Parameter :: induction_comp_vmbm_theta  = indoff+22 ! theta component of -{div dot <v>}<B>
    Integer, Parameter :: induction_advec_vmbm_theta = indoff+23 ! theta component of -{<v> dot grad <B>}
    Integer, Parameter :: induction_vmbm_theta       = indoff+24 ! theta component of del cros {<v> x <B>}
    Integer, Parameter :: induction_diff_bm_theta    = indoff+25 ! theta component of -del x (eta {del x <B>})

    Integer, Parameter :: induction_shear_vmbm_phi   = indoff+26 ! phi component of {<B> dot grad <v>}
    Integer, Parameter :: induction_comp_vmbm_phi    = indoff+27 ! phi component of -{div dot <v>}<B>
    Integer, Parameter :: induction_advec_vmbm_phi   = indoff+28 ! phi component of -{<v> dot grad <B>}
    Integer, Parameter :: induction_vmbm_phi         = indoff+29 ! phi component of del cros {<v> x <B>}
    Integer, Parameter :: induction_diff_bm_phi      = indoff+30 ! phi component of -del x (eta {del x <B>})

    !--------------- Terms involving <v> x B' 
    Integer, Parameter :: induction_shear_vmbp_r     = indoff+31 ! radial component of {B' dot grad <v>}
    Integer, Parameter :: induction_comp_vmbp_r      = indoff+32 ! radial component of -{div dot <v>}B'
    Integer, Parameter :: induction_advec_vmbp_r     = indoff+33 ! radial component of -{<v> dot grad B'}
    Integer, Parameter :: induction_vmbp_r           = indoff+34 ! radial component of del cros {<v> x B'}
    Integer, Parameter :: induction_diff_bp_r        = indoff+35 ! radial component of -del x (eta {del x B'})

    Integer, Parameter :: induction_shear_vmbp_theta = indoff+36 ! theta component of {B' dot grad <v>}
    Integer, Parameter :: induction_comp_vmbp_theta  = indoff+37 ! theta component of -{div dot <v>}B'
    Integer, Parameter :: induction_advec_vmbp_theta = indoff+38 ! theta component of -{<v> dot grad B'}
    Integer, Parameter :: induction_vmbp_theta       = indoff+39 ! theta component of del cros {<v> x B'}
    Integer, Parameter :: induction_diff_bp_theta    = indoff+40 ! theta component of -del x (eta {del x B'})

    Integer, Parameter :: induction_shear_vmbp_phi   = indoff+41 ! phi component of {B' dot grad <v>}
    Integer, Parameter :: induction_comp_vmbp_phi    = indoff+42 ! phi component of -{div dot <v>}B'
    Integer, Parameter :: induction_advec_vmbp_phi   = indoff+43 ! phi component of -{<v> dot grad B'}
    Integer, Parameter :: induction_vmbp_phi         = indoff+44 ! phi component of del cros {<v> x B'}
    Integer, Parameter :: induction_diff_bp_phi      = indoff+45 ! phi component of -del x (eta {del x B'})

    !--------------- Terms involving v' x <B> 
    Integer, Parameter :: induction_shear_vpbm_r     = indoff+46 ! radial component of {<B> dot grad v'}
    Integer, Parameter :: induction_comp_vpbm_r      = indoff+47 ! radial component of -{div dot v'}<B>
    Integer, Parameter :: induction_advec_vpbm_r     = indoff+48 ! radial component of -{v' dot grad <B>}
    Integer, Parameter :: induction_vpbm_r           = indoff+49 ! radial component of del cros {v' x <B>}

    Integer, Parameter :: induction_shear_vpbm_theta = indoff+50 ! theta component of {<B> dot grad v'}
    Integer, Parameter :: induction_comp_vpbm_theta  = indoff+51 ! theta component of -{div dot v'}<B>
    Integer, Parameter :: induction_advec_vpbm_theta = indoff+52 ! theta component of -{v' dot grad <B>}
    Integer, Parameter :: induction_vpbm_theta       = indoff+53 ! theta component of del cros {v' x <B>}

    Integer, Parameter :: induction_shear_vpbm_phi   = indoff+54 ! phi component of {<B> dot grad v'}
    Integer, Parameter :: induction_comp_vpbm_phi    = indoff+55 ! phi component of -{div dot v'}<B>
    Integer, Parameter :: induction_advec_vpbm_phi   = indoff+56 ! phi component of -{v' dot grad <B>}
    Integer, Parameter :: induction_vpbm_phi         = indoff+57 ! phi component of del cros {v' x <B>}

    !--------------- Terms involving v' x B' 
    Integer, Parameter :: induction_shear_vpbp_r     = indoff+58 ! radial component of {B' dot grad v'}
    Integer, Parameter :: induction_comp_vpbp_r      = indoff+59 ! radial component of -{div dot v'}B'
    Integer, Parameter :: induction_advec_vpbp_r     = indoff+60 ! radial component of -{v' dot grad B'}
    Integer, Parameter :: induction_vpbp_r           = indoff+61 ! radial component of del cros {v' x B'}

    Integer, Parameter :: induction_shear_vpbp_theta = indoff+62 ! theta component of {B' dot grad v'}
    Integer, Parameter :: induction_comp_vpbp_theta  = indoff+63 ! theta component of -{div dot v'}B'
    Integer, Parameter :: induction_advec_vpbp_theta = indoff+64 ! theta component of -{v' dot grad B'}
    Integer, Parameter :: induction_vpbp_theta       = indoff+65 ! theta component of del cros {v' x B'}

    Integer, Parameter :: induction_shear_vpbp_phi   = indoff+66 ! phi component of {B' dot grad v'}
    Integer, Parameter :: induction_comp_vpbp_phi    = indoff+67 ! phi component of -{div dot v'}B'
    Integer, Parameter :: induction_advec_vpbp_phi   = indoff+68 ! phi component of -{v' dot grad B'}
    Integer, Parameter :: induction_vpbp_phi         = indoff+69 ! phi component of del cros {v' x B'}

    !///////////////////////////////////
    Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output
    Real*8, Allocatable :: tmp1(:,:,:)  ! A work array
    Real*8, Allocatable :: rweights(:), tweights(:)

    !//////////////////////////////////
    Real*8, Allocatable :: ell0_values(:,:), m0_values(:,:,:)

    ! This array will hold fluctuating quantities from the buffer { q - <q>}      
    Real*8, Allocatable :: fbuffer(:,:,:,:)	
    Integer :: vindex(1:12), bindex(1:12)  ! the indices within fbuffer corresponding to v,B, and their derivatives

    ! These indices are used for references fluctuations in vfluct and bfluct
    ! "p" denotes "prime" as in v'
    Integer :: vr_p     = 1, dvrdr_p = 4 ,dvrdt_p = 5 , dvrdp_p = 6
    Integer :: vtheta_p = 2, dvtdr_p = 7 ,dvtdt_p = 8 , dvtdp_p = 9
    Integer :: vphi_p   = 3, dvpdr_p = 10,dvpdt_p = 11, dvpdp_p = 12

    Integer :: br_p     = 1, dbrdr_p = 4 ,dbrdt_p = 5 , dbrdp_p = 6
    Integer :: btheta_p = 2, dbtdr_p = 7 ,dbtdt_p = 8 , dbtdp_p = 9
    Integer :: bphi_p   = 3, dbpdr_p = 10,dbpdt_p = 11, dbpdp_p = 12

    Real*8, Allocatable :: bvars(:,:,:,:)  ! Holds the components of b and their derivatives
    ! These indices are used to reference the values held in bvars
    Integer :: br_i     = 1, dbrdr_i = 4 ,dbrdt_i = 5  , dbrdp_i = 6
    Integer :: btheta_i = 2, dbtdr_i = 7 ,dbtdt_i = 8  , dbtdp_i = 9
    Integer :: bphi_i   = 3, dbpdr_i = 10,dbpdt_i = 11 , dbpdp_i = 12

    ! These indices are used to reference the values held in add_fields%p3a
    ! That buffer contains the derivatives of each component of b
    Integer :: dbrdr = 1 , dbrdt = 2  , dbrdp = 7
    Integer :: dbtdr = 3 , dbtdt = 4  , dbtdp = 8
    Integer :: dbpdr = 5 , dbpdt = 6 ,  dbpdp = 9    

    Type(SphericalBuffer), public :: add_fields
    !Integer :: dbrdt, dbrdr, dbtdr, dbpdr
Contains

    Subroutine Generate_Diagnostic_Labels()
        ! Define labels for our quantity codes
        Write(6,*)'A line of code.'
        !Call Load_Label(v_r,'V_r')
        !Call Load_Label(v_theta,'V_theta')
        !Call Load_Label(v_phi, 'V_phi')
    End Subroutine Generate_Diagnostic_Labels



    Subroutine Adjust_Bfield(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t, j
        !This routine will handle the sintheta divisions for the magnetic components
        

		!Call sintheta_divd(add_fieldss%p3a,dbtdr_ai)
		!Call sintheta_divd(add_fields%p3a,dbpdr_ai)
		!Call sintheta_divd(add_fields%p3a,dbtdt_ai)
		!Call sintheta_divd(add_fields%p3a,dbrdt_ai)
		!Call sintheta_divd(add_fields%p3a,dbpdp_ai)
		!Call sintheta_divd(add_fields%p3a,dbtdp_ai)
        
    End Subroutine Adjust_Bfield



    Subroutine Compute_Fluctuations(buffer)
        Implicit None
        Integer :: r,k, t, j,jmax
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        jmax = size(buffer,4)
        Allocate(fbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:jmax))
        

        Do j = 1, jmax
            DO_PSI
                fbuffer(PSI,j) = buffer(PSI,j) - m0_values(PSI2,j) 
            END_DO
        Enddo

    End Subroutine Compute_Fluctuations

    Subroutine DeAllocate_Fluctuations()
        Implicit None
        DeAllocate(fbuffer)
    End Subroutine DeAllocate_Fluctuations

    Subroutine Initialize_VBIndices()
        Implicit None
        vindex(1)  = vr
        vindex(2)  = vtheta
        vindex(3)  = vphi
        vindex(4)  = dvrdr
        vindex(5)  = dvrdt
        vindex(6)  = dvrdp
        vindex(7)  = dvtdr
        vindex(8)  = dvtdt
        vindex(9)  = dvtdp
        vindex(10) = dvpdr
        vindex(11) = dvpdt
        vindex(12) = dvpdp

        bindex(1)  = br
        bindex(2)  = btheta
        bindex(3)  = bphi
        bindex(4)  = dbrdr
        bindex(5)  = dbrdt
        bindex(6)  = dbrdp
        bindex(7)  = dbtdr
        bindex(8)  = dbtdt
        bindex(9)  = dbtdp
        bindex(10) = dbpdr
        bindex(11) = dbpdt
        bindex(12) = dbpdp
    End Subroutine Initialize_VBindices

End Module Diagnostics_Base
