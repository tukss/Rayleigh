#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
Module Diagnostics_Thermodynamic_Gradients
    Use Diagnostics_Base
    Implicit None

Contains

Subroutine Compute_Thermodynamic_Gradients(buffer)
    Implicit None
    Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
    Integer :: r,k, t
    ! Here we compute the gradient of pressure and entropy/temperature


    !////////////////////////////////////////
    !       Entropy

    !  Entropy: field
    If (compute_quantity(entropy)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,tvar)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_p)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,tvar)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_m)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,tvar)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! Entropy:  radial derivatives
    If (compute_quantity(entropy_dr)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dtdr)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_p_dr)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dtdr)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_m_dr)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dtdr)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! Entropy:  theta derivatives
    If (compute_quantity(entropy_dtheta)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dtdt)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_p_dtheta)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dtdt)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_m_dtheta)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dtdt)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! Entropy:  phi derivatives
    If (compute_quantity(entropy_dphi)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dtdp)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_p_dphi)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dtdp)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_m_dphi)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dtdp)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! Entropy:  (1/r)*theta derivatives
    If (compute_quantity(entropy_dtr)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dtdt)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_p_dtr)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dtdt)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_m_dtr)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dtdt)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! Entropy:  (1/{r sintheta}) * phi derivatives
    If (compute_quantity(entropy_dprs)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dtdp)*csctheta(t)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_p_dprs)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dtdp)*csctheta(t)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(entropy_m_dprs)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dtdp)*csctheta(t)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif


    !////////////////////////////////////////
    !       Pressure

    !  pressure: field
    If (compute_quantity(pressure)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,pvar)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_p)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,pvar)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_m)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,pvar)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! pressure:  radial derivatives
    If (compute_quantity(pressure_dr)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dpdr)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_p_dr)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dpdr)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_m_dr)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dpdr)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! pressure:  theta derivatives
    If (compute_quantity(pressure_dtheta)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dpdt)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_p_dtheta)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dpdt)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_m_dtheta)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dpdt)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! pressure:  phi derivatives
    If (compute_quantity(pressure_dphi)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dpdp)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_p_dphi)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dpdp)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_m_dphi)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dpdp)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! pressure:  (1/r)*theta derivatives
    If (compute_quantity(pressure_dtr)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dpdt)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_p_dtr)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dpdt)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_m_dtr)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dpdt)*One_Over_R(r)
        END_DO
        Call Add_Quantity(qty)
    Endif

    ! pressure:  (1/{r sintheta}) * phi derivatives
    If (compute_quantity(pressure_dprs)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dpdp)*One_Over_R(r)*csctheta(t)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_p_dprs)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dpdp)*One_Over_R(r)*csctheta(t)
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(pressure_m_dprs)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dpdp)*One_Over_R(r)*csctheta(t)
        END_DO
        Call Add_Quantity(qty)
    Endif


    !Pressure:  d_by_dr(P/rho_bar)
    If (compute_quantity(rhopressure_dr)) Then
        DO_PSI
            qty(PSI) = buffer(PSI,dpdr)-buffer(PSI,pvar)*ref%dlnrho(r)       
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(rhopressurep_dr)) Then
        DO_PSI
            qty(PSI) = fbuffer(PSI,dpdr)-fbuffer(PSI,pvar)*ref%dlnrho(r)       
        END_DO
        Call Add_Quantity(qty)
    Endif

    If (compute_quantity(rhopressurem_dr)) Then
        DO_PSI
            qty(PSI) = m0_values(PSI2,dpdr)-m0_values(PSI2,pvar)*ref%dlnrho(r)       
        END_DO
        Call Add_Quantity(qty)
    Endif

End Subroutine Compute_Thermodynamic_Gradients

End Module Diagnostics_Thermodynamic_Gradients
