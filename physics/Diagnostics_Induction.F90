Module Diagnostics_Induction

Contains

    Subroutine Compute_Induction_Terms(buffer)

        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 1.    Terms involving full v and full B.
        !
        !////////////////////////////////////////////////////////////////////////
        !1a.  B dot grad v 
        If (compute_full_full_shear) Then

            Call ADotGradB(bbuffer,buffer,cbuffer,bindices=binds)

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
        If (compute_full_full_advec) Then

            Call ADotGradB(bbuffer,buffer,cbuffer,bindices=binds)

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
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,2)                    
            Endif
        Endif

    End Subroutine Compute_Induction_Terms

End Module Diagnostics_Induction
