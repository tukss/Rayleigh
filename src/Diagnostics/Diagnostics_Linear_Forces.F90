#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
#define DDBUFF d2buffer%p3a
Module Diagnostics_Linear_Forces
    Use Diagnostics_Base
    Implicit None

Contains

    Subroutine Compute_Linear_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Call Compute_Buoyancy_Force(buffer)
        Call Compute_Coriolis_Force(buffer)
        Call Compute_Viscous_Force(buffer)
    End Subroutine Compute_Linear_Forces

    Subroutine Compute_Buoyancy_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        ! Buoyancy forces as they contribute to the ell .ne. 0 components of the momentum equation

        ! -- full buoyancy 
        If (compute_quantity(buoyancy_force)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*(buffer(PSI,tvar)-&
                           & ell0_values(r,tvar))
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! -- fluctuating buoyancy (ell = 0, m =0 already subtracted)
        If (compute_quantity(buoyancy_pforce)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! -- mean buoyancy
        If (compute_quantity(buoyancy_mforce)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*(m0_values(PSI2,tvar)-&
                           & ell0_values(r,tvar))
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Buoyancy_Force

    Subroutine Compute_Coriolis_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: coriolis_term
        coriolis_term = ref%Coriolis_Coeff

        If(compute_quantity(Coriolis_Force_r)) Then
            DO_PSI
                qty(PSI) = -ref%density(r)*coriolis_term*sintheta(t)*buffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If(compute_quantity(Coriolis_pForce_r)) Then
            DO_PSI
                qty(PSI) = -ref%density(r)*coriolis_term*sintheta(t)*fbuffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If(compute_quantity(Coriolis_mForce_r)) Then
            DO_PSI
                qty(PSI) = -ref%density(r)*coriolis_term*sintheta(t)*m0_values(PSI2,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_Force_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*buffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_pForce_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*fbuffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_mForce_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*m0_values(PSI2,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(Coriolis_Force_Phi)) Then
            DO_PSI
                qty(PSI) = - coriolis_term*costheta(t)*buffer(PSI,vtheta) &
					       - coriolis_term*sintheta(t)*buffer(PSI,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_pForce_Phi)) Then
            DO_PSI
                qty(PSI) = - coriolis_term*costheta(t)*fbuffer(PSI,vtheta) &
					       - coriolis_term*sintheta(t)*fbuffer(PSI,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_mForce_Phi)) Then
            DO_PSI
                qty(PSI) = - coriolis_term*costheta(t)*m0_values(PSI2,vtheta) &
					       - coriolis_term*sintheta(t)*m0_values(PSI2,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Coriolis_Force


    Subroutine Compute_Viscous_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: mypi, amp,del2u, estress
        Real*8, Allocatable :: mu_visc(:), dmudr(:), ovstheta(:), ovs2theta(:)

        Allocate(ovstheta(1:N_theta), ovs2theta(1:N_theta)) ! 1/sin; 1/sin^2
        ovstheta = 1.0d0/sintheta
        ovs2theta = 1.0d0/sin2theta

        !Compute the dynamic viscosity mu=rho*nu (nu is kinematic viscosity)
        Allocate(mu_visc(1:N_R), dmudr(1:N_R))

        mu_visc = ref%density*nu
        dmudr = mu_visc*ref%density*(ref%dlnrho+dlnu)


        !////////////////////////////////////////////////////////

        ! r-direction; Full
        If (compute_quantity(viscous_force_r)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
                del2u = DDBUFF(PSI,dvrdrdr)+Two_Over_R(r)*buffer(PSI,dvrdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvrdtdt)+cottheta(t)*buffer(PSI,dvrdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u-2.0d0*OneOverRsquared(r)*( &
                        buffer(PSI,vr) + & 
                        buffer(PSI,dvtdt)+buffer(PSI,vtheta)*cottheta(t) + &
                        ovstheta(t)*buffer(PSI,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*(buffer(PSI,dvrdr)*ref%dlnrho(r)+ &
                        buffer(PSI,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
                estress = buffer(PSI,dvrdr)-One_Third*buffer(PSI,vr)*ref%dlnrho(r)

                qty(PSI) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u


            END_DO

            Call Add_Quantity(qty)
        Endif

        ! r-direction; fluctuating
        If (compute_quantity(viscous_pforce_r)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
                del2u = d2_fbuffer(PSI,dvrdrdr)+Two_Over_R(r)*fbuffer(PSI,dvrdr)
                del2u = del2u+OneOverRSquared(r)*(d2_fbuffer(PSI,dvrdtdt)+cottheta(t)*fbuffer(PSI,dvrdt))
                del2u = del2u+OneOverRSquared(r)*d2_fbuffer(PSI,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u-2.0d0*OneOverRsquared(r)*( &
                        fbuffer(PSI,vr) + & 
                        fbuffer(PSI,dvtdt)+fbuffer(PSI,vtheta)*cottheta(t) + &
                        ovstheta(t)*fbuffer(PSI,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*(fbuffer(PSI,dvrdr)*ref%dlnrho(r)+ &
                        fbuffer(PSI,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
                estress = fbuffer(PSI,dvrdr)-One_Third*fbuffer(PSI,vr)*ref%dlnrho(r)

                qty(PSI) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u


            END_DO

            Call Add_Quantity(qty)
        Endif

        ! r-direction; mean
        If (compute_quantity(viscous_mforce_r)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
                del2u = d2_m0(PSI2,dvrdrdr)+Two_Over_R(r)*m0_values(PSI2,dvrdr)
                del2u = del2u+OneOverRSquared(r)*(d2_m0(PSI2,dvrdtdt)+cottheta(t)*m0_values(PSI2,dvrdt))
                del2u = del2u+OneOverRSquared(r)*d2_m0(PSI2,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u-2.0d0*OneOverRsquared(r)*( &
                        m0_values(PSI2,vr) + & 
                        m0_values(PSI2,dvtdt)+m0_values(PSI2,vtheta)*cottheta(t) + &
                        ovstheta(t)*m0_values(PSI2,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*(m0_values(PSI2,dvrdr)*ref%dlnrho(r)+ &
                        m0_values(PSI2,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
                estress = m0_values(PSI2,dvrdr)-One_Third*m0_values(PSI2,vr)*ref%dlnrho(r)

                qty(PSI) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u


            END_DO

            Call Add_Quantity(qty)
        Endif


        !Theta-direction; Full
        If (compute_quantity(viscous_force_theta)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_theta}
                del2u = DDBUFF(PSI,dvtdrdr)+Two_Over_R(r)*buffer(PSI,dvtdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvtdtdt)+cottheta(t)*buffer(PSI,dvtdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_theta
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*buffer(PSI,dvrdt) - &
                        ovs2theta(t)*(   buffer(PSI,vtheta) + &
                        2.0d0*costheta(t)*buffer(PSI,dvpdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*buffer(PSI,dvrdt)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(buffer(PSI,dvrdt)-buffer(PSI,vtheta) ) +buffer(PSI,dvtdr)

                qty(PSI) = dmudr(r)*estress +mu_visc(r)*del2u
            END_DO

            Call Add_Quantity(qty)
        Endif

        !Theta-direction; Fluctuating
        If (compute_quantity(viscous_force_theta)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_theta}
                del2u = d2_fbuffer(PSI,dvtdrdr)+Two_Over_R(r)*fbuffer(PSI,dvtdr)
                del2u = del2u+OneOverRSquared(r)*(d2_fbuffer(PSI,dvtdtdt)+cottheta(t)*fbuffer(PSI,dvtdt))
                del2u = del2u+OneOverRSquared(r)*d2_fbuffer(PSI,dvtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_theta
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*fbuffer(PSI,dvrdt) - &
                        ovs2theta(t)*(   fbuffer(PSI,vtheta) + &
                        2.0d0*costheta(t)*fbuffer(PSI,dvpdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*fbuffer(PSI,dvrdt)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(fbuffer(PSI,dvrdt)-fbuffer(PSI,vtheta) ) +fbuffer(PSI,dvtdr)

                qty(PSI) = dmudr(r)*estress +mu_visc(r)*del2u
            END_DO

            Call Add_Quantity(qty)
        Endif



        !Phi-direction
        If (compute_quantity(viscous_force_phi)) Then

            DO_PSI
                del2u = DDBUFF(PSI,dvpdrdr)+Two_Over_R(r)*buffer(PSI,dvpdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvpdtdt)+cottheta(t)*buffer(PSI,dvpdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvpdpdp)*ovs2theta(t)

                ! Add geometric terms here
                ! del2u - 
                !Add geometric terms to make this { Del^2{u} }_phi
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*buffer(PSI,dvrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   buffer(PSI,vphi) - &
                        2.0d0*costheta(t)*buffer(PSI,dvtdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*ovstheta(t)*buffer(PSI,dvrdp)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(ovstheta(t)*buffer(PSI,dvrdp)-buffer(PSI,vphi) )+ &
                          buffer(PSI,dvpdr)

                qty(PSI) =dmudr(r)*estress + mu_visc(r)*del2u
            END_DO

            Call Add_Quantity(qty)
        Endif


        DeAllocate(mu_visc, dmudr)
        DeAllocate(ovstheta,ovs2theta)
    End Subroutine Compute_Viscous_Force

End Module Diagnostics_Linear_Forces
