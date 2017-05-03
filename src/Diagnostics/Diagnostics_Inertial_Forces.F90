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
        Real*8, Allocatable :: cbuffer(:,:,:,:)



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


        !//////////////////////////////////////////////////////////////////////////////////
        !/////////////// v dot grad v (full)//////////////////
        If (compute_full_full) Then
            Call ADotGradB(buffer,buffer,cbuffer,aindices=vindex,bindices=vindex)
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

        !/////////////// v' dot grad v' //////////////////
        If (compute_fluct_fluct) Then
            Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices=vindex,bindices=vindex)
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

        !/////////////// <v> dot grad <v> //////////////////
        If (compute_mean_mean) Then
            Call ADotGradB(m0_values,m0_values,cbuffer,aindices=vindex,bindices=vindex)
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

        !/////////////// v' dot grad <v> //////////////////
        If (compute_fluct_mean) Then
            Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = vindex, bindices=vindex)
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

        !/////////////// <v> dot grad v' //////////////////
        If (compute_mean_fluct) Then
            Call ADotGradB(m0_values,fbuffer,cbuffer,aindices=vindex,bindices=vindex)
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

    End Subroutine Compute_Inertial_Terms
End Module Diagnostics_Inertial_Forces
