#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

Module Diagnostics_Inertial_Forces

    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
Contains
    Subroutine Compute_Inertial_Terms(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8, Allocatable :: tbuffer(:,:,:,:), cbuffer(:,:,:,:)

        Integer :: binds(1:12)

        Integer :: vr_p     = 1, dvrdr_p = 4 ,dvrdt_p = 5 , dvrdp_p = 6
        Integer :: vtheta_p = 2, dvtdr_p = 7 ,dvtdt_p = 8 , dvtdp_p = 9
        Integer :: vphi_p   = 3, dvpdr_p = 10,dvpdt_p = 11, dvpdp_p = 12

        Logical :: compute_fluctuations = .false.
        Logical :: compute_full_full    = .false.
        Logical :: compute_fluct_mean   = .false.
        Logical :: compute_mean_fluct   = .false.
        Logical :: compute_fluct_fluct  = .false.
        Logical :: compute_mean_mean    = .false.

        !First, figure out what we need to bother computing

        !///-- Flags for full dot grad full terms
        If (compute_quantity(v_grad_v_r)) Then
            compute_full_full = .true.
        Endif
        If (compute_quantity(v_grad_v_theta)) Then
            compute_full_full = .true.
        Endif
        If (compute_quantity(v_grad_v_phi)) Then
            compute_full_full = .true.
        Endif

        !///-- Flags for mean dot grad mean terms
        If (compute_quantity(vm_grad_vm_r)) Then
            compute_mean_mean   = .true.
        Endif
        If (compute_quantity(vm_grad_vm_theta)) Then
            compute_mean_mean   = .true.
        Endif
        If (compute_quantity(vm_grad_vm_phi)) Then
            compute_mean_mean   = .true.
        Endif

        !///-- Flags for fluctuating dot grad mean terms
        If (compute_quantity(vp_grad_vm_r)) Then
            compute_fluctuations = .true.
            compute_fluct_mean   = .true.
        Endif
        If (compute_quantity(vp_grad_vm_theta)) Then
            compute_fluctuations = .true.
            compute_fluct_mean   = .true.
        Endif
        If (compute_quantity(vp_grad_vm_phi)) Then
            compute_fluctuations = .true.
            compute_fluct_mean   = .true.
        Endif

        !///-- Flags for mean dot grad fluctuating terms
        If (compute_quantity(vm_grad_vp_r)) Then
            compute_fluctuations = .true.
            compute_mean_fluct   = .true.
        Endif
        If (compute_quantity(vm_grad_vp_theta)) Then
            compute_fluctuations = .true.
            compute_mean_fluct   = .true.
        Endif
        If (compute_quantity(vm_grad_vp_phi)) Then
            compute_fluctuations = .true.
            compute_mean_fluct   = .true.
        Endif

        !///-- Flags for fluctuating dot grad fluctuating terms
        If (compute_quantity(vp_grad_vp_r)) Then
            compute_fluctuations = .true.
            compute_fluct_fluct  = .true.
        Endif
        If (compute_quantity(vp_grad_vp_theta)) Then
            compute_fluctuations = .true.
            compute_fluct_fluct  = .true.
        Endif
        If (compute_quantity(vp_grad_vp_phi)) Then
            compute_fluctuations = .true.
            compute_fluct_fluct  = .true.
        Endif

        Allocate(cbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        If (compute_fluctuations) Then
            Allocate(tbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:12))
            !tbuffer holds perturbations about the azimuthal mean

            ! V_r terms
            DO_PSI
                tbuffer(PSI,vr_p) = buffer(PSI,vr) - m0_values(PSI2,vr) 
            END_DO

            DO_PSI
                tbuffer(PSI,dvrdr_p) = buffer(PSI,dvrdr) - m0_values(PSI2,dvrdr) 
            END_DO

            DO_PSI
                tbuffer(PSI,dvrdt_p) = buffer(PSI,dvrdt) - m0_values(PSI2,dvrdt) 
            END_DO

            DO_PSI
                tbuffer(PSI,dvrdp_p) = buffer(PSI,dvrdp) 
            END_DO

            ! V_theta terms
            DO_PSI
                tbuffer(PSI,vtheta_p) = buffer(PSI,vtheta) - m0_values(PSI2,vtheta) 
            END_DO

            DO_PSI
                tbuffer(PSI,dvtdr_p) = buffer(PSI,dvtdr) - m0_values(PSI2,dvtdr) 
            END_DO

            DO_PSI
                tbuffer(PSI,dvtdt_p) = buffer(PSI,dvtdt) - m0_values(PSI2,dvtdt) 
            END_DO

            DO_PSI
                tbuffer(PSI,dvtdp_p) = buffer(PSI,dvtdp)  
            END_DO

            ! V_phi terms
            DO_PSI
                tbuffer(PSI,vphi_p) = buffer(PSI,vphi) - m0_values(PSI2,vphi) 
            END_DO        

            DO_PSI
                tbuffer(PSI,dvpdr_p) = buffer(PSI,dvpdr) - m0_values(PSI2,dvpdr) 
            END_DO

            DO_PSI
                tbuffer(PSI,dvpdt_p) = buffer(PSI,dvpdt) - m0_values(PSI2,dvpdt) 
            END_DO

            DO_PSI
                tbuffer(PSI,dvpdp_p) = buffer(PSI,dvpdp) 
            END_DO

        Endif

        binds(1)  = vr
        binds(2)  = vtheta
        binds(3)  = vphi
        binds(4)  = dvrdr
        binds(5)  = dvrdt
        binds(6)  = dvrdp
        binds(7)  = dvtdr
        binds(8)  = dvtdt
        binds(9)  = dvtdp
        binds(10) = dvpdr
        binds(11) = dvpdt
        binds(12) = dvpdp

        !//////////////////////////////////////////////////////////////////////////////////
        !/////////////// v dot grad v (full)//////////////////
        If (compute_full_full) Then
            Call ADotGradB(buffer,buffer,cbuffer,aindices=binds,bindices=binds)
            If (compute_quantity(v_grad_v_r)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(v_grad_v_theta)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(v_grad_v_phi)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_fluct_fluct) Then
            Call ADotGradB(tbuffer,tbuffer,cbuffer)
            If (compute_quantity(vp_grad_vp_r)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(vp_grad_vp_theta)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(vp_grad_vp_phi)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_mean_mean) Then
            Call ADotGradB(m0_values,m0_values,cbuffer,aindices=binds,bindices=binds)
            If (compute_quantity(vm_grad_vm_r)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(vm_grad_vm_theta)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(vm_grad_vm_phi)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_fluct_mean) Then
            Call ADotGradB(tbuffer,m0_values,cbuffer,bindices=binds)
            If (compute_quantity(vp_grad_vm_r)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(vp_grad_vm_theta)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(vp_grad_vm_phi)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_mean_fluct) Then
            Call ADotGradB(m0_values,tbuffer,cbuffer,aindices=binds)
            If (compute_quantity(vm_grad_vp_r)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(vm_grad_vp_theta)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(vm_grad_vp_phi)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        DeAllocate(cbuffer)
        If (compute_fluctuations) Then
            DeAllocate(tbuffer)
        Endif
    End Subroutine Compute_Inertial_Terms
End Module Diagnostics_Inertial_Forces
