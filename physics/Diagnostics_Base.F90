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
    Real*8, Allocatable :: ell0_values(:,:), m0_values(:,:,:), bm0_values(:,:,:)

    ! These arrays are used to hold fluctuating v- and b-related quantities      
    Real*8, Allocatable :: vfluct(:,:,:,:), bfluct(:,:,:,:)	
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
    Integer :: dbrdr_ai = 1 , dbrdt_ai = 2  , dbrdp_ai = 7
    Integer :: dbtdr_ai = 3 , dbtdt_ai = 4  , dbtdp_ai = 8
    Integer :: dbpdr_ai = 5 , dbpdt_ai = 6 ,  dbpdp_ai = 9    

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

    Subroutine Compute_Vfluctuations(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Allocate(vfluct(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:12))



        ! V_r terms
        DO_PSI
            vfluct(PSI,vr_p) = buffer(PSI,vr) - m0_values(PSI2,vr) 
        END_DO

        DO_PSI
            vfluct(PSI,dvrdr_p) = buffer(PSI,dvrdr) - m0_values(PSI2,dvrdr) 
        END_DO

        DO_PSI
            vfluct(PSI,dvrdt_p) = buffer(PSI,dvrdt) - m0_values(PSI2,dvrdt) 
        END_DO

        DO_PSI
            vfluct(PSI,dvrdp_p) = buffer(PSI,dvrdp) 
        END_DO

        ! V_theta terms
        DO_PSI
            vfluct(PSI,vtheta_p) = buffer(PSI,vtheta) - m0_values(PSI2,vtheta) 
        END_DO

        DO_PSI
            vfluct(PSI,dvtdr_p) = buffer(PSI,dvtdr) - m0_values(PSI2,dvtdr) 
        END_DO

        DO_PSI
            vfluct(PSI,dvtdt_p) = buffer(PSI,dvtdt) - m0_values(PSI2,dvtdt) 
        END_DO

        DO_PSI
            vfluct(PSI,dvtdp_p) = buffer(PSI,dvtdp)  
        END_DO

        ! V_phi terms
        DO_PSI
            vfluct(PSI,vphi_p) = buffer(PSI,vphi) - m0_values(PSI2,vphi) 
        END_DO        

        DO_PSI
            vfluct(PSI,dvpdr_p) = buffer(PSI,dvpdr) - m0_values(PSI2,dvpdr) 
        END_DO

        DO_PSI
            vfluct(PSI,dvpdt_p) = buffer(PSI,dvpdt) - m0_values(PSI2,dvpdt) 
        END_DO

        DO_PSI
            vfluct(PSI,dvpdp_p) = buffer(PSI,dvpdp) 
        END_DO
    End Subroutine Compute_Vfluctuations

    Subroutine Load_Bfield(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t, j
        Allocate(bvars(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:12))
        Allocate(bm0_values(my_r%min:my_r%max,my_theta%min:my_theta%max,1:12))
        !At this point, we'll need to reform the b buffer
        Call add_fields%reform() !moves from p2a to p3a
        !Call Phi_DerivativesD()
		!Call fft_to_physical(wsp%p3a,rsc = .true.)

        

		!Call sintheta_divd(add_fieldss%p3a,dbtdr_ai)
		!Call sintheta_divd(add_fields%p3a,dbpdr_ai)
		!Call sintheta_divd(add_fields%p3a,dbtdt_ai)
		!Call sintheta_divd(add_fields%p3a,dbrdt_ai)
		!Call sintheta_divd(add_fields%p3a,dbpdp_ai)
		!Call sintheta_divd(add_fields%p3a,dbtdp_ai)
        DO_PSI
            bvars(PSI,br_i)     = buffer(PSI,br)
        END_DO
        DO_PSI
            bvars(PSI,btheta_i) = buffer(PSI,btheta)
        END_DO
        DO_PSI
            bvars(PSI,bphi_i)   = buffer(PSI,bphi)
        END_DO

        DO_PSI
            bvars(PSI,dbrdr_i)  = add_fields%p3a(PSI,dbrdr_ai)
        END_DO
        DO_PSI
            bvars(PSI,dbrdt_i)  = add_fields%p3a(PSI,dbrdt_ai)
        END_DO
        DO_PSI
            bvars(PSI,dbrdp_i)  = add_fields%p3a(PSI,dbrdp_ai)
        END_DO

        DO_PSI
            bvars(PSI,dbtdr_i)  = add_fields%p3a(PSI,dbtdr_ai)
        END_DO
        DO_PSI
            bvars(PSI,dbtdt_i)  = add_fields%p3a(PSI,dbtdt_ai)
        END_DO
        DO_PSI
            bvars(PSI,dbtdp_i)  = add_fields%p3a(PSI,dbtdp_ai)
        END_DO

        DO_PSI
            bvars(PSI,dbpdr_i)  = add_fields%p3a(PSI,dbpdr_ai)
        END_DO
        DO_PSI
            bvars(PSI,dbpdt_i)  = add_fields%p3a(PSI,dbpdt_ai)
        END_DO
        DO_PSI
            bvars(PSI,dbpdp_i)  = add_fields%p3a(PSI,dbpdp_ai)
        END_DO



        
        Call ComputeM0(bvars,bm0_values)
    End Subroutine Load_Bfield

    Subroutine Compute_Bfluctuations()
        Implicit None
        Integer :: r,k, t, j
        Allocate(bfluct(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:12))
        ! B_r terms
        ! The bvars array is assumed to be ordered as expected by
        ! the a Dot Grad B routines.  Thus, we just loop over all 12 variables

        Do j = 1, 12
            DO_PSI
                bfluct(PSI,j) = bvars(PSI,j) - bm0_values(PSI2,j) 
            END_DO
        Enddo

    End Subroutine Compute_Bfluctuations

End Module Diagnostics_Base
