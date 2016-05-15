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
    !  Diagnostic Quantity Codes 
    !  Reynolds decompositions are often used in the outputs.
    !  As a result, some shorthand is used as follows:
    !   "m" and "< >" denote the azimuthal OR spherical mean.
    !   "p" and " ' " denote perturbations about the azimuthal mean 

    !////////////////////////////////////////////////////////
    ! Velocity field components and it's derivatives

    Integer, Parameter :: voffset = 0

    !------------ Field Components ----------!
    Integer, Parameter :: v_r      = voffset+1 ! Full 
    Integer, Parameter :: v_theta  = voffset+2     
    Integer, Parameter :: v_phi    = voffset+3          

    Integer, Parameter :: vp_r     = voffset+4 ! Fluctuating 
    Integer, Parameter :: vp_theta = voffset+5     
    Integer, Parameter :: vp_phi   = voffset+6    


    Integer, Parameter :: vm_r     = voffset+7 ! Mean
    Integer, Parameter :: vm_theta = voffset+8     
    Integer, Parameter :: vm_phi   = voffset+9    

    !------------ Radial Derivatives -------------!
    Integer, Parameter :: dv_r_dr      = voffset+10 ! Full 
    Integer, Parameter :: dv_theta_dr  = voffset+11  
    Integer, Parameter :: dv_phi_dr    = voffset+12   

    Integer, Parameter :: dvp_r_dr     = voffset+13 ! Fluctuating 
    Integer, Parameter :: dvp_theta_dr = voffset+14  
    Integer, Parameter :: dvp_phi_dr   = voffset+15    
  
    Integer, Parameter :: dvm_r_dr     = voffset+16 ! Mean    
    Integer, Parameter :: dvm_theta_dr = voffset+17  
    Integer, Parameter :: dvm_phi_dr   = voffset+18    

    !------------ Theta Derivatives --------------!
    Integer, Parameter :: dv_r_dt      = voffset+19 ! Full 
    Integer, Parameter :: dv_theta_dt  = voffset+20  
    Integer, Parameter :: dv_phi_dt    = voffset+21   

    Integer, Parameter :: dvp_r_dt     = voffset+22 ! Fluctuating 
    Integer, Parameter :: dvp_theta_dt = voffset+23    
    Integer, Parameter :: dvp_phi_dt   = voffset+24    

    Integer, Parameter :: dvm_r_dt     = voffset+25 ! Mean 
    Integer, Parameter :: dvm_theta_dt = voffset+26   
    Integer, Parameter :: dvm_phi_dt   = voffset+27     

    !------------ Phi Derivatives ----------------!
    Integer, Parameter :: dv_r_dp      = voffset+28  ! Full
    Integer, Parameter :: dv_theta_dp  = voffset+29
    Integer, Parameter :: dv_phi_dp    = voffset+30

    Integer, Parameter :: dvp_r_dp     = voffset+31 ! Fluctuating
    Integer, Parameter :: dvp_theta_dp = voffset+32 ! (same as full if    
    Integer, Parameter :: dvp_phi_dp   = voffset+33 !  mean is azimuthal)

    Integer, Parameter :: dvm_r_dp     = voffset+34 ! Mean    
    Integer, Parameter :: dvm_theta_dp = voffset+35 ! (nonzero only when mean is 
    Integer, Parameter :: dvm_phi_dp   = voffset+36 !  spherical, not azimuthal)

    !------------ (1/r) * Theta Derivatives -------!
    Integer, Parameter :: dv_r_dtr      = voffset+37 ! Full 
    Integer, Parameter :: dv_theta_dtr  = voffset+38 
    Integer, Parameter :: dv_phi_dtr    = voffset+39   

    Integer, Parameter :: dvp_r_dtr     = voffset+40 ! Fluctuating
    Integer, Parameter :: dvp_theta_dtr = voffset+41
    Integer, Parameter :: dvp_phi_dtr   = voffset+42

    Integer, Parameter :: dvm_r_dtr     = voffset+43 ! Mean  
    Integer, Parameter :: dvm_theta_dtr = voffset+44
    Integer, Parameter :: dvm_phi_dtr   = voffset+45


    !------(1/{r sintheta})* Phi Derivatives ---!

    Integer, Parameter :: dv_r_dprs      = voffset+46 ! Full
    Integer, Parameter :: dv_theta_dprs  = voffset+47  
    Integer, Parameter :: dv_phi_dprs    = voffset+48

    Integer, Parameter :: dvp_r_dprs     = voffset+49 ! Fluctuating   
    Integer, Parameter :: dvp_theta_dprs = voffset+50
    Integer, Parameter :: dvp_phi_dprs   = voffset+51
 
    Integer, Parameter :: dvm_r_dprs     = voffset+52 ! Mean   
    Integer, Parameter :: dvm_theta_dprs = voffset+53
    Integer, Parameter :: dvm_phi_dprs   = voffset+54  

    !------------ Mass Flux ---------------------!
    Integer, Parameter :: rhov_r      = voffset+55 ! Full 
    Integer, Parameter :: rhov_theta  = voffset+56
    Integer, Parameter :: rhov_phi    = voffset+57

    Integer, Parameter :: rhovp_r     = voffset+58 ! Fluctuating 
    Integer, Parameter :: rhovp_theta = voffset+59
    Integer, Parameter :: rhovp_phi   = voffset+60

    Integer, Parameter :: rhovm_r     = voffset+61 ! Mean
    Integer, Parameter :: rhovm_theta = voffset+62
    Integer, Parameter :: rhovm_phi   = voffset+63


    !//////////////////////////////////////////////////////////////////
    !  Pressure, Entropy or Temperature, and Their Derivatives
    !  Note:  In Boussinesq Mode, Temperature is Output Instead of Entropy
    Integer, Parameter :: pt_off = voffset+63 ! = 63

    !------------ Fields ---------------------!    
    Integer, Parameter :: entropy    = pt_off+1 ! Full
    Integer, Parameter :: pressure   = pt_off+2

    Integer, Parameter :: entropy_p  = pt_off+3 ! Fluctuating
    Integer, Parameter :: pressure_p = pt_off+4

    Integer, Parameter :: entropy_m  = pt_off+5 ! Mean
    Integer, Parameter :: pressure_m = pt_off+6

    !------------ Radial Derivatives --------------!
    Integer, Parameter :: entropy_dr     = pt_off+7  ! Full
    Integer, Parameter :: pressure_dr    = pt_off+8

    Integer, Parameter :: entropy_p_dr   = pt_off+9  ! Fluctuating
    Integer, Parameter :: pressure_p_dr  = pt_off+10

    Integer, Parameter :: entropy_m_dr   = pt_off+11 ! Mean
    Integer, Parameter :: pressure__m_dr = pt_off+12

    !------------ Theta Derivatives ---------------!
    Integer, Parameter :: entropy_dt     = pt_off+13 ! Full
    Integer, Parameter :: pressure_dt    = pt_off+14

    Integer, Parameter :: entropy_p_dt   = pt_off+15 ! Fluctuating
    Integer, Parameter :: pressure_p_dt  = pt_off+16

    Integer, Parameter :: entropy_m_dt   = pt_off+17 ! Mean
    Integer, Parameter :: pressure__m_dt = pt_off+18

    !------------ Phi Derivatives -----------------!
    Integer, Parameter :: entropy_dp     = pt_off+13 ! Full
    Integer, Parameter :: pressure_dp    = pt_off+14

    Integer, Parameter :: entropy_p_dp   = pt_off+15 ! Fluctuating
    Integer, Parameter :: pressure_p_dp  = pt_off+16

    Integer, Parameter :: entropy_m_dp   = pt_off+17 ! Mean
    Integer, Parameter :: pressure__m_dp = pt_off+18

    !------------ (1/r) * Theta Derivatives --------!
    Integer, Parameter :: entropy_dtr     = pt_off+19 ! Full
    Integer, Parameter :: pressure_dtr    = pt_off+20

    Integer, Parameter :: entropy_p_dtr   = pt_off+21 ! Fluctuating
    Integer, Parameter :: pressure_p_dtr  = pt_off+22

    Integer, Parameter :: entropy_m_dtr   = pt_off+23 ! Mean
    Integer, Parameter :: pressure__m_dtr = pt_off+24

    !--- (1/{r sintheta}) * Phi Derivatives ---------!
    Integer, Parameter :: entropy_dprs     = pt_off+25 ! Full
    Integer, Parameter :: pressure_dprs    = pt_off+26

    Integer, Parameter :: entropy_p_dprs   = pt_off+27 ! Fluctuating
    Integer, Parameter :: pressure_p_dprs  = pt_off+28

    Integer, Parameter :: entropy_m_dprs   = pt_off+29 ! Mean
    Integer, Parameter :: pressure__m_dprs = pt_off+30

    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Vorticity Outputs
    Integer, Parameter :: vort_off = pt_off+30 ! = 93

    Integer, Parameter :: vort_r      = vort_off+1  ! Full
    Integer, Parameter :: vort_theta  = vort_off+2
    Integer, Parameter :: vort_phi    = vort_off+3

    Integer, Parameter :: vortp_r     = vort_off+4  ! Fluctuating   
    Integer, Parameter :: vortp_theta = vort_off+5
    Integer, Parameter :: vortp_phi   = vort_off+6

    Integer, Parameter :: vortm_r     = vort_off+7  ! Mean 
    Integer, Parameter :: vortm_theta = vort_off+8
    Integer, Parameter :: vortm_phi   = vort_off+9

    Integer, Parameter :: enstrophy   = vort_off+10 ! Enstrophy
    Integer, Parameter :: enstrophypm = vort_off+11 ! (fluctuating-mean)
    Integer, Parameter :: enstrophymm = vort_off+12 ! (mean-mean)
    Integer, Parameter :: enstropypp  = vort_off+13 ! (fluct-fluct)


    !//////////////////////////////////////////////////////////
    !               Radial Energy Fluxes
    Integer, Parameter :: eoffset = vort_off+ 13 ! =106
    Integer, Parameter :: ecrossb_r            = eoffset+1 ! [ExB]_r (un-normalized Poynting flux)
    Integer, Parameter :: ke_flux_radial       = eoffset+2 ! vr*KE
    Integer, Parameter :: thermalE_flux_radial = eoffset+3 ! vr*rho_bar*T_bar*S OR vr*T
    Integer, Parameter :: enth_flux_radial     = eoffset+4 ! vr*cp*rho_bar*T
    Integer, Parameter :: visc_flux_r          = eoffset+5 ! -[Div dot D]_r
    Integer, Parameter :: vol_heat_flux        = eoffset+6 ! "Flux" associated with volumetric heating
    Integer, Parameter :: cond_flux_r          = eoffset+7 ! Thermal conductive flux
    Integer, Parameter :: vol_heating          = eoffset+8 ! Volumetric Heating Function



    !///////////////////////////////////////////////////////////
    !           Kinetic Energies
    Integer, Parameter :: keoffset = eoffset+8 ! = 114
    Integer, Parameter :: kinetic_energy = keoffset+1   ! 1/2 rho_bar v^2
    Integer, Parameter :: radial_ke      = keoffset+2   ! 1/2 rho_bar {v_r}^2
    Integer, Parameter :: theta_ke       = keoffset+3   ! 1/2 rho_bar {v_theta}^2
    Integer, Parameter :: phi_ke         = keoffset+4   ! 1/2 rho_bar {v_phi}^2

    Integer, Parameter :: mkinetic_energy = keoffset+5 ! 1/2 rho_bar <v>^2
    Integer, Parameter :: radial_mke      = keoffset+6 ! 1/2 rho_bar <v_r>^2
    Integer, Parameter :: theta_mke       = keoffset+7 ! 1/2 rho_bar <v_theta>^2
    Integer, Parameter :: phi_mke         = keoffset+8 ! 1/2 rho_bar <v_phi>^2

    Integer, Parameter :: pkinetic_energy = keoffset+9  ! 1/2 rho_bar {v'}^2
    Integer, Parameter :: radial_pke      = keoffset+10 ! 1/2 rho_bar {v_r'}^2
    Integer, Parameter :: theta_pke       = keoffset+11 ! 1/2 rho_bar {v_theta'}^2
    Integer, Parameter :: phi_pke         = keoffset+12 ! 1/2 rho_bar {v_phi'}^2

    !--- Since density varies with radius, it may be useful to output the 
    !--- squared fields, sans density, as well.
    Integer, Parameter :: vsq         = keoffset+13   ! v^2
    Integer, Parameter :: radial_vsq  = keoffset+14   ! {v_r}^2
    Integer, Parameter :: theta_vsq   = keoffset+15   ! {v_theta}^2
    Integer, Parameter :: phi_vsq     = keoffset+16   ! {v_phi}^2

    Integer, Parameter :: mvsq        = keoffset+17   ! <v>^2
    Integer, Parameter :: radial_mvsq = keoffset+18   ! <v_r>^2
    Integer, Parameter :: theta_mvsq  = keoffset+19   ! <v_theta>^2
    Integer, Parameter :: phi_mvsq    = keoffset+20   ! <v_phi>^2

    Integer, Parameter :: pvsq        = keoffset+21   ! {v'}^2
    Integer, Parameter :: radial_pvsq = keoffset+22   ! {v_r'}^2
    Integer, Parameter :: theta_pvsq  = keoffset+23   ! {v_theta'}^2
    Integer, Parameter :: phi_pvsq    = keoffset+24   ! {v_phi'}^2




    !///////////////////////////////////////////////////////////////////////////////////
    !Angular Momentum Transport Diagnostics
    !  Reynolds decomposition of the azimuthally-averaged angular momentum fluxes.
    Integer, Parameter :: amoff = keoffset+ 24        ! = 138
    Integer, Parameter :: amom_fluct_r = amoff+1      ! rho_bar * r * sintheta * {v_r'} * {v_phi'}
    Integer, Parameter :: amom_fluct_theta = amoff+2  ! rho_bar * r * sintheta * {v_theta'} * {v_phi'}
    Integer, Parameter :: amom_dr_r = amoff+3         ! rho_bar * r * sintheta * <v_r> * <v_phi>
    Integer, Parameter :: amom_dr_theta = amoff+4     ! rho_bar * r * sintheta * <v_theta> * <v_phi>
    Integer, Parameter :: amom_mean_r = amoff+5       ! rho_bar * r^2 * sintheta^2 * <v_r> * Omega
    Integer, Parameter :: amom_mean_theta = amoff+6   ! rho_bar * r^2 * sintheta^2 * <v_theta> * Omega


    !/////////////////////////////////////////////////////////////
    !  Now we get into the hydro forces.  

    !////////////////////////////////////////
    ! Placeholder for buoyancy force terms  (pressure is already taken care of above)

    !////////////////////////////////////////
    ! Placeholder for Viscous Forces -- ugh


    !////////////////////////  Advection Terms ////////////////////
    ! Reynolds decomposition about the azimuthal mean may also be output
    ! NOTE:  ADVECTION TERMS ARE SCALED BY DENSITY (so that they represent a force density)

    Integer, Parameter :: vgv = amoff+6 ! = 134    ! Output offset for advection terms  
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
    Integer, Parameter :: dcheck_off = vgv+ 15 ! = 149
    Integer, Parameter :: diagnostic1 = dcheck_off+1, diagnostic2 = dcheck_off+2


    !//////////////////////////////////////////////////////////
    !/////////// Magnetic Outputs.  Start at 200 to additional room for hydro


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
    Integer, Parameter :: joffset = boffset+50 ! = 250

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

    !///////////////////////////////////////////////////////////
    !           Magnetic Energies
    Integer, Parameter :: meoffset = joffset+20 ! = 270

    Integer, Parameter :: magnetic_energy = meoffset+1 ! B^2
    Integer, Parameter :: radial_me       = meoffset+2 ! {B_r}^2
    Integer, Parameter :: theta_me        = meoffset+3 ! {B_theta}^2
    Integer, Parameter :: phi_me        = meoffset+4 ! {B_phi}^2

    Integer, Parameter :: mmagnetic_energy = meoffset+5 ! <B>^2
    Integer, Parameter :: radial_mme       = meoffset+6 ! <B_r>^2
    Integer, Parameter :: theta_mme        = meoffset+7 ! <B_theta>^2
    Integer, Parameter :: phi_mme        = meoffset+8 ! <B_phi>^2

    Integer, Parameter :: pmagnetic_energy = meoffset+9  ! {B'}^2
    Integer, Parameter :: radial_pme       = meoffset+10 ! {B_r'}^2
    Integer, Parameter :: theta_pme        = meoffset+11 ! {B_theta'}^2
    Integer, Parameter :: phi_pme        = meoffset+12 ! {B_phi'}^2



    !/////////////////////////// Lorentz Forces ///////////////////////////////
    Integer, Parameter :: loff = joffset+12 ! =282
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
    Integer, Parameter :: indoff = loff + 15 ! = 297
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


    ! END QUANTITY CODES
    ! END QUANTITY CODES
    ! END QUANTITY CODES
    !/////////////////////////////////////////////////////////////////////////////////////////////


    !///////////////////////////////////
    Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output
    Real*8, Allocatable :: tmp1(:,:,:)  ! A work array
    Real*8, Allocatable :: rweights(:), tweights(:), tmp1d(:)

    !//////////////////////////////////
    Real*8, Allocatable :: ell0_values(:,:), m0_values(:,:,:)

    ! This array will hold fluctuating quantities from the buffer { q - <q>}      
    Real*8, Allocatable :: fbuffer(:,:,:,:)	


Contains

    Subroutine Generate_Diagnostic_Labels()
        ! Define labels for our quantity codes
        Write(6,*)'A line of code.'
        !Call Load_Label(v_r,'V_r')
        !Call Load_Label(v_theta,'V_theta')
        !Call Load_Label(v_phi, 'V_phi')
    End Subroutine Generate_Diagnostic_Labels

    Subroutine Initialize_Diagnostics_Buffer()
        Logical :: dbtrans, dbconfig
        Logical :: test_reduce

        dbtrans = .false.
        dbconfig = .false.
        test_reduce = .false.

        Call cobuffer%init(field_count = cbfcount, config = 'p1a', &
            dynamic_transpose =dbtrans, dynamic_config = dbconfig, &
            hold_cargo = test_reduce, padding = pad_alltoall)        
    End Subroutine Initialize_Diagnostics_Buffer





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


End Module Diagnostics_Base
