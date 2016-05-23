#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_INDUCTION
!               This module computes del x (vxB), its 
!               constituent terms (i.e., B dot grad v), 
!               and their Reynolds decomposition
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Induction
    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
    Implicit None
    Logical :: allocate_indr = .false.
    Logical :: allocate_indt = .false.
    Logical :: allocate_indp = .false.
    Logical :: compute_shear = .false.
    Logical :: compute_advec = .false.
    Logical :: compute_vmbm_shear = .false.
    Logical :: compute_vmbm_advec = .false.
    Logical :: compute_vmbp_shear = .false.
    Logical :: compute_vmbp_advec = .false.
    Logical :: compute_vpbm_shear = .false.
    Logical :: compute_vpbm_advec = .false.
    Logical :: compute_vpbp_shear = .false.
    Logical :: compute_vpbp_advec = .false.

Contains

    Subroutine Compute_Induction_Terms(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Allocatable :: ind_r(:,:,:), ind_theta(:,:,:), ind_phi(:,:,:)
        Real*8, Allocatable :: cbuffer(:,:,:,:)  
        Integer :: r,k, t
        Allocate(cbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Call Reset_Induction_Flags()

        If (allocate_indr) Then
            Allocate(    ind_r(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif
        If (allocate_indt) Then
            Allocate(ind_theta(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif
        If (allocate_indp) Then
            Allocate(  ind_phi(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif

        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 1.    Terms resulting full v cross full B.
        !
        !////////////////////////////////////////////////////////////////////////
        !1a.  B dot grad v 
        If (compute_shear) Then

            Call ADotGradB(buffer,buffer,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induction_shear_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)                  
            Endif
            If (compute_quantity(induction_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)                   
            Endif
            If (compute_quantity(induction_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)                
            Endif
        Endif

        !1b.  -v dot grad B 
        If (compute_advec) Then

            Call ADotGradB(buffer,buffer,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induction_advec_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_advec_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)                    
            Endif
            If (compute_quantity(induction_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)                    
            Endif
            If (compute_quantity(induction_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)                    
            Endif
        Endif

        !1c.  -B (div dot v)
        ! Take care with the logic here...

        If (compute_quantity(induction_comp_r) .or. compute_quantity(induction_r)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,br)*buffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_r)) Call Add_Quantity(qty)
            If (compute_quantity(induction_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induction_comp_theta) .or. compute_quantity(induction_theta)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,btheta)*buffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induction_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induction_comp_phi) .or. compute_quantity(induction_phi)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,bphi)*buffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induction_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_phi)
            Endif
        Endif

        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 2.    Terms resulting from <v> x B'.
        !
        !////////////////////////////////////////////////////////////////////////
        !2a.  B' dot grad <v> 
        If (compute_vmbp_shear) Then

            Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induction_shear_vmbp_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vmbp_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vmbp_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_vmbp_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)                  
            Endif
            If (compute_quantity(induction_vmbp_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)                   
            Endif
            If (compute_quantity(induction_vmbp_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)                
            Endif
        Endif

        !2b.  -<v> dot grad B' 
        If (compute_vmbp_advec) Then

            Call ADotGradB(m0_values,fbuffer,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induction_advec_vmbp_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_advec_vmbp_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vmbp_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_vmbp_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)                    
            Endif
            If (compute_quantity(induction_vmbp_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)                    
            Endif
            If (compute_quantity(induction_vmbp_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)                    
            Endif
        Endif

        !2c.  -B' (div dot <v>)
        ! Take care with the logic here...

        If (compute_quantity(induction_comp_vmbp_r) .or. compute_quantity(induction_vmbp_r)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,br)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vmbp_r)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vmbp_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induction_comp_vmbp_theta) .or. compute_quantity(induction_vmbp_theta)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,btheta)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vmbp_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vmbp_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induction_comp_vmbp_phi) .or. compute_quantity(induction_vmbp_phi)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,bphi)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vmbp_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vmbp_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_phi)
            Endif
        Endif

        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 3.    Terms resulting from v' x <B>.
        !
        !////////////////////////////////////////////////////////////////////////
        !3a.  <B> dot grad v' 
        If (compute_vpbm_shear) Then

            Call ADotGradB(m0_values,fbuffer,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induction_shear_vpbm_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vpbm_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vpbm_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_vpbm_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)                  
            Endif
            If (compute_quantity(induction_vpbm_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)                   
            Endif
            If (compute_quantity(induction_vpbm_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)                
            Endif
        Endif

        !3b.  -v' dot grad <B> 
        If (compute_vpbm_advec) Then

            Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induction_advec_vpbm_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_advec_vpbm_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vpbm_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_vpbm_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)                    
            Endif
            If (compute_quantity(induction_vpbm_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)                    
            Endif
            If (compute_quantity(induction_vpbm_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)                    
            Endif
        Endif

        !3c.  -<B> (div dot v')
        ! Take care with the logic here...

        If (compute_quantity(induction_comp_vpbm_r) .or. compute_quantity(induction_vpbm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,br)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vpbm_r)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vpbm_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induction_comp_vpbm_theta) .or. compute_quantity(induction_vpbm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,btheta)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vpbm_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vpbm_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induction_comp_vpbm_phi) .or. compute_quantity(induction_vpbm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,bphi)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vpbm_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vpbm_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_phi)
            Endif
        Endif


        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 4.    Terms resulting from <v> x <B>.
        !
        !////////////////////////////////////////////////////////////////////////
        !4a.  <B> dot grad <v> 
        If (compute_vmbm_shear) Then

            Call ADotGradB(m0_values,m0_values,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induction_shear_vmbm_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vmbm_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vmbm_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_vmbm_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)                  
            Endif
            If (compute_quantity(induction_vmbm_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)                   
            Endif
            If (compute_quantity(induction_vmbm_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)                
            Endif
        Endif

        !4b.  -<v> dot grad <B> 
        If (compute_vmbm_advec) Then

            Call ADotGradB(m0_values,m0_values,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induction_advec_vmbm_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_advec_vmbm_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vmbm_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_vmbm_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)                    
            Endif
            If (compute_quantity(induction_vmbm_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)                    
            Endif
            If (compute_quantity(induction_vmbm_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)                    
            Endif
        Endif

        !4c.  -<B> (div dot <v>)
        ! Take care with the logic here...

        If (compute_quantity(induction_comp_vmbm_r) .or. compute_quantity(induction_vmbm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,br)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vmbm_r)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vmbm_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induction_comp_vmbm_theta) .or. compute_quantity(induction_vmbm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,btheta)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vmbm_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vmbm_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induction_comp_vmbm_phi) .or. compute_quantity(induction_vmbm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,bphi)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vmbm_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vmbm_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_phi)
            Endif
        Endif


        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 5.    Terms resulting from v' x B'.
        !
        !////////////////////////////////////////////////////////////////////////
        !5a.  B' dot grad v' 
        If (compute_vpbp_shear) Then

            Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induction_shear_vpbp_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vpbp_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vpbp_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_vpbp_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)                  
            Endif
            If (compute_quantity(induction_vpbp_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)                   
            Endif
            If (compute_quantity(induction_vpbp_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)                
            Endif
        Endif

        !5b.  -v' dot grad B' 
        If (compute_vpbp_advec) Then

            Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induction_advec_vpbp_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_advec_vpbp_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induction_shear_vpbp_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induction_vpbp_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)                    
            Endif
            If (compute_quantity(induction_vpbp_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)                    
            Endif
            If (compute_quantity(induction_vpbp_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)                    
            Endif
        Endif

        !5c.  -B' (div dot v')
        ! Take care with the logic here...

        If (compute_quantity(induction_comp_vpbp_r) .or. compute_quantity(induction_vpbp_r)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,br)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vpbp_r)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vpbp_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induction_comp_vpbp_theta) .or. compute_quantity(induction_vpbp_theta)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,btheta)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vpbp_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vpbp_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induction_comp_vpbp_phi) .or. compute_quantity(induction_vpbp_phi)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,bphi)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induction_comp_vpbp_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induction_vpbp_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)                    
                 Call Add_Quantity(ind_phi)
            Endif
        Endif


        If (    allocated(ind_r)) DeAllocate(ind_r)
        If (allocated(ind_theta)) DeAllocate(ind_theta)
        If (  allocated(ind_phi)) DeAllocate(ind_phi)
        DeAllocate(cbuffer)
    End Subroutine Compute_Induction_Terms

    Subroutine Compute_Magnetic_Diffusion(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        !Diffusion terms will be added in the near-ish future
        qty(:,:,:) = 0.0d0

        If (compute_quantity(induction_diff_r))  Call Add_Quantity(qty)
        If (compute_quantity(induction_diff_theta))  Call Add_Quantity(qty)
        If (compute_quantity(induction_diff_phi))  Call Add_Quantity(qty)

        If (compute_quantity(induction_diff_bm_r))  Call Add_Quantity(qty)
        If (compute_quantity(induction_diff_bm_theta))  Call Add_Quantity(qty)
        If (compute_quantity(induction_diff_bm_phi))  Call Add_Quantity(qty)

        If (compute_quantity(induction_diff_bp_r))  Call Add_Quantity(qty)
        If (compute_quantity(induction_diff_bp_theta))  Call Add_Quantity(qty)
        If (compute_quantity(induction_diff_bp_phi))  Call Add_Quantity(qty)

    End Subroutine Compute_Magnetic_Diffusion

    Subroutine Reset_Induction_Flags()
        Implicit None
        !/////////////////////////////////////////////////
        !   Compute_Induction_Terms makes use of several logical flags.
        !   Before carrying out the computations in that subroutine,
        !   we set these flags to the appropriate value.

        !  First, reset ALL flags to false
        !  (This routine should be called at each output iteration)
        allocate_indr = .false.
        allocate_indt = .false.
        allocate_indp = .false.
        compute_shear = .false.
        compute_advec = .false.
        compute_vmbm_shear = .false.
        compute_vmbm_advec = .false.
        compute_vmbp_shear = .false.
        compute_vmbp_advec = .false.
        compute_vpbp_shear = .false.
        compute_vpbp_advec = .false.
        compute_vpbm_shear = .false.
        compute_vpbm_advec = .false.

        ! 1.  Decide if we need to allocate the induction arrays
        If (compute_quantity(induction_r     )) allocate_indr = .true.
        If (compute_quantity(induction_vmbm_r)) allocate_indr = .true.
        If (compute_quantity(induction_vmbp_r)) allocate_indr = .true.
        If (compute_quantity(induction_vpbm_r)) allocate_indr = .true.
        If (compute_quantity(induction_vpbp_r)) allocate_indr = .true.

        If (compute_quantity(induction_theta     )) allocate_indt = .true.
        If (compute_quantity(induction_vmbm_theta)) allocate_indt = .true.
        If (compute_quantity(induction_vmbp_theta)) allocate_indt = .true.
        If (compute_quantity(induction_vpbm_theta)) allocate_indt = .true.
        If (compute_quantity(induction_vpbp_theta)) allocate_indt = .true.

        If (compute_quantity(induction_phi     )) allocate_indp = .true.
        If (compute_quantity(induction_vmbm_phi)) allocate_indp = .true.
        If (compute_quantity(induction_vmbp_phi)) allocate_indp = .true.
        If (compute_quantity(induction_vpbm_phi)) allocate_indp = .true.
        If (compute_quantity(induction_vpbp_phi)) allocate_indp = .true.

        !2.  Set flags related to full v x B        
        If (compute_quantity(induction_r)) Then
            compute_shear = .true.
            compute_advec = .true.
        Endif

        If (compute_quantity(induction_shear_r)) compute_shear = .true.
        If (compute_quantity(induction_advec_r)) compute_advec = .true.

        If (compute_quantity(induction_theta)) Then
            compute_shear = .true.
            compute_advec = .true.
        Endif

        If (compute_quantity(induction_shear_theta)) compute_shear = .true.
        If (compute_quantity(induction_advec_theta)) compute_advec = .true.

        If (compute_quantity(induction_phi)) Then
            compute_shear = .true.
            compute_advec = .true.
        Endif

        If (compute_quantity(induction_shear_phi)) compute_shear = .true.
        If (compute_quantity(induction_advec_phi)) compute_advec = .true.

        !3.  Set flags related to <v> x <B>        
        If (compute_quantity(induction_vmbm_r)) Then
            compute_vmbm_shear = .true.
            compute_vmbm_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vmbm_r)) compute_vmbm_shear = .true.
        If (compute_quantity(induction_advec_vmbm_r)) compute_vmbm_advec = .true.

        If (compute_quantity(induction_vmbm_theta)) Then
            compute_vmbm_shear = .true.
            compute_vmbm_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vmbm_theta)) compute_vmbm_shear = .true.
        If (compute_quantity(induction_advec_vmbm_theta)) compute_vmbm_advec = .true.

        If (compute_quantity(induction_vmbm_phi)) Then
            compute_vmbm_shear = .true.
            compute_vmbm_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vmbm_phi)) compute_vmbm_shear = .true.
        If (compute_quantity(induction_advec_vmbm_phi)) compute_vmbm_advec = .true.


        !4. Set flags related to v' x <B>
        If (compute_quantity(induction_vpbm_r)) Then
            compute_vpbm_shear = .true.
            compute_vpbm_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vpbm_r)) compute_vpbm_shear = .true.
        If (compute_quantity(induction_advec_vpbm_r)) compute_vpbm_advec = .true.

        If (compute_quantity(induction_vpbm_theta)) Then
            compute_vpbm_shear = .true.
            compute_vpbm_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vpbm_theta)) compute_vpbm_shear = .true.
        If (compute_quantity(induction_advec_vpbm_theta)) compute_vpbm_advec = .true.

        If (compute_quantity(induction_vpbm_phi)) Then
            compute_vpbm_shear = .true.
            compute_vpbm_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vpbm_phi)) compute_vpbm_shear = .true.
        If (compute_quantity(induction_advec_vpbm_phi)) compute_vpbm_advec = .true.

        !5. Set flags related to v' x 'B'
        If (compute_quantity(induction_vpbp_r)) Then
            compute_vpbp_shear = .true.
            compute_vpbp_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vpbp_r)) compute_vpbp_shear = .true.
        If (compute_quantity(induction_advec_vpbp_r)) compute_vpbp_advec = .true.

        If (compute_quantity(induction_vpbp_theta)) Then
            compute_vpbp_shear = .true.
            compute_vpbp_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vpbp_theta)) compute_vpbp_shear = .true.
        If (compute_quantity(induction_advec_vpbp_theta)) compute_vpbp_advec = .true.

        If (compute_quantity(induction_vpbp_phi)) Then
            compute_vpbp_shear = .true.
            compute_vpbp_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vpbp_phi)) compute_vpbp_shear = .true.
        If (compute_quantity(induction_advec_vpbp_phi)) compute_vpbp_advec = .true.

        !6. Set flags related to <v> x 'B'
        If (compute_quantity(induction_vmbp_r)) Then
            compute_vmbp_shear = .true.
            compute_vmbp_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vmbp_r)) compute_vmbp_shear = .true.
        If (compute_quantity(induction_advec_vmbp_r)) compute_vmbp_advec = .true.

        If (compute_quantity(induction_vmbp_theta)) Then
            compute_vmbp_shear = .true.
            compute_vmbp_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vmbp_theta)) compute_vmbp_shear = .true.
        If (compute_quantity(induction_advec_vmbp_theta)) compute_vmbp_advec = .true.

        If (compute_quantity(induction_vmbp_phi)) Then
            compute_vmbp_shear = .true.
            compute_vmbp_advec = .true.
        Endif

        If (compute_quantity(induction_shear_vmbp_phi)) compute_vmbp_shear = .true.
        If (compute_quantity(induction_advec_vmbp_phi)) compute_vmbp_advec = .true.

    End Subroutine Reset_Induction_Flags

End Module Diagnostics_Induction
