! REFERENCE STATE MODULE
! Contains routines for initializing the reference state structure
! Reference state structure contains all information related to background
! stratification.  It DOES NOT contain transport variable (e.g. nu, kappa) information



Module ReferenceState
	Use ProblemSize
    Use Controls
    Use Math_Constants
    Use Math_Utility
    Use General_MPI, Only : BCAST2D
	Implicit None
	Type ReferenceInfo
		Real*8, Allocatable :: Density(:)
		Real*8, Allocatable :: dlnrho(:)
		Real*8, Allocatable :: d2lnrho(:)

		Real*8, Allocatable :: Pressure(:)

		Real*8, Allocatable :: Temperature(:)
		Real*8, Allocatable :: dlnT(:)

		Real*8, Allocatable :: Entropy(:)
		Real*8, Allocatable :: dsdr(:)

		Real*8, Allocatable :: Gravity(:)
        Real*8, Allocatable :: Gravity_term_s(:)    ! -(gravity/rho)*drho_by_ds ..typically = gravity/cp
		Real*8 :: gamma
        Real*8, Allocatable :: heating(:)
		Real*8 :: rho_twiddle, g_twiddle, p_twiddle, s_twiddle, t_twiddle
	End Type ReferenceInfo
    Real*8, Allocatable :: s_conductive(:)

    Integer :: reference_type
    Integer :: heating_type = 0 ! 0 means no reference heating.  > 0 selects optional reference heating
    Real*8  :: Luminosity
    Real*8  :: heating_factor, heating_r0
    Type(ReferenceInfo) :: ref

    Real*8 :: pressure_specific_heat  ! CP (not CV)
    Real*8 :: poly_n
    Real*8 :: poly_Nrho
    Real*8 :: poly_mass
    Real*8 :: poly_rho_i
    Real*8 :: Gravitational_Constant = 6.67d-8
    Real*8 :: rho_twiddle, g_twiddle, p_twiddle, s_twiddle, t_twiddle, length_twiddle ! also in ref structure..

    Real*8 :: Angular_Velocity = 1.0d0


    
	!/////////////////////////////////////////////////////////////////////////////////////
	! Nondimensional Parameters
	Real*8 :: Rayleigh_Number         = 0.0d0
	Real*8 :: Ekman_Number            = 0.0d0
	Real*8 :: Prandtl_Number          = 1.0d0
	Real*8 :: Magnetic_Prandtl_Number = 1.0d0
    Real*8 :: gravity_power           = 0.0d0
    Logical :: Dimensional = .true.  ! By Default code is dimensional
    Character*120 :: custom_reference_file ='nothing'
	Namelist /Reference_Namelist/ reference_type,poly_n, poly_Nrho, poly_mass,poly_rho_i, &
            & pressure_specific_heat, heating_type, luminosity, Angular_Velocity, &
            & Rayleigh_Number, Ekman_Number, Prandtl_Number, Magnetic_Prandtl_Number, &
            & gravity_power, dimensional,heating_factor, heating_r0, custom_reference_file
Contains

	Subroutine Initialize_Reference()
		Implicit None
		Call Allocate_Reference_State()
		If (reference_type .eq. 1) Then
			Call Constant_Reference()
		Endif

		If (reference_type .eq. 2) Then
			Call Polytropic_Reference()
		Endif

        If (reference_type .eq. 3) Then
            Call Get_Custom_Reference()
        Endif

		Call Write_Reference()
	End Subroutine Initialize_Reference

	Subroutine Allocate_Reference_State
		Implicit None
		Allocate(ref%density(1:N_R))
		Allocate(ref%pressure(1:N_R))
		Allocate(ref%temperature(1:N_R))
		Allocate(ref%entropy(1:N_R))
		Allocate(ref%gravity(1:N_R))
		Allocate(ref%dlnrho(1:N_R))
		Allocate(ref%d2lnrho(1:N_R))
		Allocate(ref%dlnt(1:N_R))
		Allocate(ref%dsdr(1:N_R))
        Allocate(ref%gravity_term_s(1:N_R))
	End Subroutine Allocate_Reference_State

    Subroutine Polytropic_Reference()
        Real*8 :: zeta_0,  c0, c1, d
        Real*8 :: rho_c, P_c, T_c,denom
        Real*8 :: beta, Gas_Constant
        Real*8, Allocatable :: zeta(:)
        Real*8 :: One, ee
        Real*8 :: InnerRadius, OuterRadius
        Integer :: r
        If (my_rank .eq. 0) Call stdout%print('Initializing polytropic reference state.')
        ! Adiabatic, Polytropic Reference State (see, e.g., Jones et al. 2011)
        ! The following parameters are read from the input file.
        ! poly_n
        ! poly_Nrho
        ! poly_mass
        ! poly_rho_i

        ! Note that cp must also be specified.
        InnerRadius = Radius(N_r)
        OuterRadius = Radius(1)
		

        One = 1.0d0
        !-----------------------------------------------------------
        beta = InnerRadius/OuterRadius

        denom = beta * exp(poly_Nrho / poly_n) + 1.d0
        zeta_0 = (beta+1.d0)/denom

        c0 = (2.d0 * zeta_0 - beta - 1.d0) / (1.d0 - beta)

        denom = (1.d0 - beta)**2
        c1 = (1.d0+beta)*(1.d0-zeta_0)/denom

        !-----------------------------------------------------------
        ! allocate and define zeta
        ! also rho_c, T_c, P_c

        Allocate(zeta(N_R))

        d = OuterRadius - InnerRadius    

        zeta = c0 + c1 * d / Radius

        rho_c = poly_rho_i / zeta(N_R)**poly_n

        denom = (poly_n+1.d0) * d * c1
        P_c = Gravitational_Constant * poly_mass * rho_c / denom

        T_c = (poly_n+1.d0) * P_c / (Pressure_Specific_Heat * rho_c)

        !-----------------------------------------------------------
        ! Initialize reference structure 
        ref%gamma = (poly_n+1.0D0)/(poly_n)


        Gas_Constant = (ref%Gamma-one)*Pressure_Specific_Heat/ref%Gamma

        Ref%Gravity = Gravitational_Constant * poly_mass / Radius**2

        Ref%Density = rho_c * zeta**poly_n

        Ref%dlnrho = - poly_n * c1 * d / (zeta * Radius**2)
        Ref%d2lnrho = - Ref%dlnrho*(2.0d0/Radius-c1*d/zeta/Radius**2)

        Ref%Temperature = T_c * zeta
        Ref%dlnT = -(c1*d/Radius**2)/zeta

        Ref%Pressure = P_c * zeta**(poly_n+1)

        denom = P_c**(1.d0/ref%gamma)
        Ref%Entropy = Pressure_Specific_Heat * log(denom/rho_c)

        Ref%dsdr = 0.d0

        Ref%gravity_term_s = ref%gravity/Pressure_Specific_Heat*ref%density

        !We initialize s_conductive (modulo delta_s, specified by the boundary conditions)
        Allocate(s_conductive(1:N_R))
        s_conductive(:) = 0.0d0
        ee = -1.d0*poly_n
        denom = zeta(1)**ee - zeta(N_R)**ee
        Do r = 1, N_R
          s_conductive(r) = (zeta(1)**ee - zeta(r)**ee) / denom
        Enddo

        Deallocate(zeta)

        Call Initialize_Reference_Heating()

	End Subroutine Polytropic_Reference


    Subroutine Initialize_Reference_Heating()
        Implicit None
        ! This is where a volumetric heating function Phi(r) is computed
        ! This function appears in the entropy equation as
        ! dSdt = Phi(r)
        ! Phi(r) may represent internal heating of any type.  For stars, this heating would be
        ! associated with temperature diffusion of the reference state and/or nuclear burning.

        If (heating_type .gt. 0) Then
            If (.not. Allocated(ref%heating)) Allocate(ref%heating(1:N_R))
            ref%heating(:) = 0.0d0
        Endif

        If (heating_type .eq. 1) Then
            Call Constant_Reference_Heating()
        Endif

        If (heating_type .eq. 2) Then
            Call Tanh_Reference_Heating()
        Endif

	If (heating_type .eq. 3) Then
	    Call Bouss_Reference_Heating()
	Endif

    End Subroutine Initialize_Reference_Heating

    Subroutine Constant_Reference_Heating()
        Implicit None
        Real*8 :: integral, alpha
        Real*8, Allocatable :: temp(:)

        ! Luminosity is specified as an input
        ! Phi(r) is set to alpha such that 
        ! Integral_r=rinner_r=router (4*pi*alpha*rho(r)*T(r)*r^2 dr) = Luminosity
        Allocate(temp(1:N_R))

        temp = ref%density*ref%temperature
        Call Integrate_in_radius(temp,integral)
        integral = integral*4.0d0*pi
        alpha = Luminosity/integral
        ref%heating(:) = alpha
        DeAllocate(temp)
    End Subroutine Constant_Reference_Heating

    Subroutine Bouss_Reference_Heating()
	Implicit None
	ref%heating(:) = 1.0d0
    End Subroutine Bouss_Reference_Heating

    Subroutine Tanh_Reference_Heating()
        Implicit None
        Real*8 :: integral, alpha
        Real*8, Allocatable :: temp(:), x(:), temp2(:)
		  Integer :: i
        Character*120 :: heating_file

        ! Luminosity is specified as an input
        ! Heating is set so that temp * 4 pi r^2 integrates to one Lsun 
        ! Integral_r=rinner_r=router (4*pi*alpha*rho(r)*T(r)*r^2 dr) = Luminosity
        Allocate(temp(1:N_R))
		  Allocate(temp2(1:N_R))
        Allocate(x(1:N_R))
        x = heating_factor*(radius-heating_r0)/ (maxval(radius)-minval(radius)) ! x runs from zero to 1 if heating_r0 is min(radius)
        !Call tanh_profile(x,temp)
        Do i = 1, n_r
            temp2(i) = 0.5d0*(1.0d0-tanh(x(i))*tanh(x(i)))*heating_factor/(maxval(radius)-minval(radius))
        Enddo

		  !temp2 = heating_factor*(1-(temp*temp))/ (maxval(radius)-minval(radius))

		  temp2 = -temp2/(4*pi*radius*radius)

        Call Integrate_in_radius(temp2,integral)

        integral = integral*4.0d0*pi
        alpha = Luminosity/integral
        
		  ref%heating(:) = alpha*temp2/(ref%density*ref%temperature)
        if (my_rank .eq. 0) Then
            heating_file = Trim(my_path)//'reference_heating'
            Open(unit=15,file=heating_file,form='unformatted', status='replace',access='stream')
			Write(15)n_r
			Write(15)(radius(i),i=1,n_r)
			Write(15)(ref%heating(i), i = 1, n_r)
    		Write(15)(temp(i), i = 1, n_r)
			Write(15)(temp2(i), i = 1, n_r)
			Write(15)(ref%density(i), i = 1, n_r)
			Write(15)(ref%temperature(i), i = 1, n_r)
        Endif
		  
        DeAllocate(x,temp, temp2)
    End Subroutine Tanh_Reference_Heating


    Subroutine Integrate_in_radius(func,int_func)
        Implicit None
        Real*8, Intent(In) :: func(1:)
        Real*8, Intent(Out) :: int_func
        Integer :: i
        Real*8 :: delr, riweight
        !compute integrate_r=rmin_r=rmax func*r^2 dr
        int_func = 0.0d0
        Do i = 2, n_r-1
            delr = (radius(i-1)-radius(i+1))/2.0d0
            riweight = delr*radius(i)**2
            int_func = int_func+func(i)*riweight
        Enddo
        delr = (radius(1)-radius(2))/2.0d0
        riweight = delr*radius(1)**2
        int_func = int_func+riweight*func(1)

        delr = (radius(n_r-1)-radius(n_r))/2.0d0
        riweight = delr*radius(n_r)**2
        int_func = int_func+riweight*func(n_r)

    End Subroutine Integrate_in_radius

    Subroutine Indefinite_Integral(func,int_func)
        Implicit None
        Real*8, Intent(In) :: func(1:)
        Real*8, Intent(Out) :: int_func(1:)
        Integer :: i
        Real*8 :: delr
        !computes indefinite integral func dr
        int_func(1) = 0.0d0
        Do i = 2, n_r-1
            delr = (radius(i-1)-radius(i+1))/2.0d0
            int_func(i) = int_func(i-1)+func(i-1)*delr
        Enddo
        delr = (radius(n_r-1)-radius(n_r))/2.0d0

        int_func(n_r) = int_func(n_r-1)+delr*func(n_r-1)


    End Subroutine Indefinite_Integral

	Subroutine Write_Reference(filename)
		Implicit None
		Character*120, Optional, Intent(In) :: filename
		Character*120 :: ref_file
		Integer :: i,sig = 314
		if (present(filename)) then
			ref_file = Trim(my_path)//filename
		else
			ref_file = Trim(my_path)//'reference'
		endif

		If (my_rank .eq. 0) Then
			Open(unit=15,file=ref_file,form='unformatted', status='replace',access='stream')
            Write(15)sig
			Write(15)n_r
			Write(15)(radius(i),i=1,n_r)
			Write(15)(ref%density(i),i=1,n_r)
			Write(15)(ref%dlnrho(i),i=1,n_r)
			Write(15)(ref%d2lnrho(i),i=1,n_r)
			Write(15)(ref%pressure(i),i=1,n_r)
			Write(15)(ref%temperature(i),i=1,n_r)
			Write(15)(ref%dlnT(i),i=1,n_r)
			Write(15)(ref%dsdr(i),i=1,n_r)
			Write(15)(ref%entropy(i),i=1,n_r)			
			Write(15)(ref%gravity(i),i=1,n_r)
			Close(15)
		Endif
	End Subroutine Write_Reference

	Subroutine Write_Profile(arr,filename)
		Implicit None
		Character*120, Optional, Intent(In) :: filename
		Character*120 :: ref_file
		Integer :: i,j,nq,sig = 314
        Real*8, Intent(In) :: arr(1:,1:)
        nq = size(arr,2)
		If (my_rank .eq. 0) Then
			Open(unit=15,file=filename,form='unformatted', status='replace',access='stream')
            Write(15)sig
			Write(15)n_r
            Write(15)nq
			Write(15)((arr(i,j),i=1,n_r),j = 1, nq)

			Close(15)
		Endif
	End Subroutine Write_Profile

    Subroutine Get_Custom_Reference()
        Implicit None
        Real*8, Allocatable :: ref_arr(:,:), integrand(:)
        Allocate(ref_arr(1:n_r,1:9))
        If (my_rank .eq. 0) Then
            Call Read_Reference(custom_reference_file, ref_arr)
        Endif
        If (my_row_rank .eq. 0) Then
            ! Broadcast along the column
            Call BCAST2D(ref_arr,grp = pfi%ccomm)
        Endif
        Call BCAST2D(ref_arr,grp = pfi%rcomm)

        ref%density(:) = ref_arr(:,1)
        ref%dlnrho(:) = ref_arr(:,2) 
        ref%d2lnrho(:) = ref_arr(:,3)
        ref%pressure(:) = ref_arr(:,4)
        ref%temperature(:) = ref_arr(:,5)
        ref%dlnT(:) = ref_arr(:,6)
        ref%dsdr(:) = ref_arr(:,7)
        ref%entropy(:) = ref_arr(:,8)			
        ref%gravity(:) = ref_arr(:,9)
        ref%gravity_term_s(:) = -ref%temperature*ref%dlnT
        DeAllocate(ref_arr)

        ! This conductive profile is based on the assumption that kappa is constant
        ! And it is modulo kappa (s_conductive/kappa)
        Allocate(s_conductive(1:N_R))
        s_conductive(:) = 0.0d0
        Allocate(integrand(1:N_R))
        integrand = 1.0d0/(ref%density*ref%temperature)
        integrand = integrand*OneOverRSquared
        !The routine below integrates r^2 * integrand, hence the extra division by r^2
        Call Indefinite_Integral(integrand,s_conductive)
        s_conductive = s_conductive-s_conductive(1)
        s_conductive = s_conductive/s_conductive(N_R)
        DeAllocate(integrand)
    End Subroutine Get_Custom_Reference

    Subroutine Read_Reference(filename,ref_arr)
		Character*120, Intent(In), Optional :: filename
        Character*120 :: ref_file
        Integer :: pi_integer,nr_ref
        Integer :: nqvals =9 ! Reference state file contains nqvals quantities + radius
        Integer :: i, k
        Real*8, Allocatable :: ref_arr_old(:,:), rtmp(:), rtmp2(:)
        Real*8, Intent(InOut) :: ref_arr(:,:)
        Real*8, Allocatable :: old_radius(:)
        If (present(filename)) Then
            ref_file = Trim(my_path)//filename
        Else
            ref_file = 'reference'
        Endif
		Open(unit=15,file=ref_file,form='unformatted', status='old',access='stream')
        Read(15)pi_integer

        If (pi_integer .ne. 314) Then
            close(15)
            Open(unit=15,file=ref_file,form='unformatted', status='old', &
                 CONVERT = 'BIG_ENDIAN' , access='stream')
            Read(15)pi_integer
            If (pi_integer .ne. 314) Then
                Close(15)
                Open(unit=15,file=ref_file,form='unformatted', status='old', &
                 CONVERT = 'LITTLE_ENDIAN' , access='stream')
                Read(15)pi_integer
            Endif
        Endif
                
        If (pi_integer .eq. 314) Then 

		    Read(15)nr_ref
            Allocate(ref_arr_old(1:nr_ref,1:nqvals))  !10 quantities are stored in reference state
            Allocate(old_radius(1:nr_ref))

            Write(6,*)'nr_ref is: ', nr_ref

            Read(15)(old_radius(i),i=1,nr_ref)
            Do k = 1, nqvals
		        Read(15)(ref_arr_old(i,k) , i=1 , nr_ref)
            Enddo
   		    

            !Check to see if radius isn't reversed
            !If it is not, then reverse it
            If (old_radius(1) .lt. old_radius(nr_ref)) Then
                Write(6,*)'Reversing Radial Indices in Custom Ref File!'
                Allocate(rtmp(1:nr_ref))

                Do i = 1, nr_ref
                    old_radius(i) = rtmp(nr_ref-i+1)
                Enddo

                Do k = 1, nqvals
                    rtmp(:) = ref_arr_old(:,k)
                    Do i = 1, nr_ref
                        ref_arr_old(i,k) = rtmp(nr_ref-i+1)
                    Enddo
                Enddo

                DeAllocate(rtmp)

            Endif

            Close(15)


            If (nr_ref .ne. n_r) Then 
                !Interpolate onto the current radial grid if necessary
                !Note that the underlying assumption here is that same # of grid points
                ! means same grid - come back to this later for generality
                Allocate(rtmp2(1:n_r))
                Allocate(rtmp(1:nr_ref))
                
                Do k = 1, nqvals
 
                    rtmp(:) = ref_arr_old(:,k)
                    rtmp2(:) = 0.0d0
                    Call Spline_Interpolate(rtmp, old_radius, rtmp2, radius)
                    ref_arr(1:n_r,k) = rtmp2 
                Enddo

                DeAllocate(rtmp,rtmp2)
            Else
                ! Bit redundant here, but may want to do filtering on ref_arr array
                ref_arr(1:n_r,1:nqvals) = ref_arr_old(1:n_r,1:nqvals)

            Endif


            DeAllocate(ref_arr_old,old_radius)
        Endif
	End Subroutine Read_Reference

    

	Subroutine Constant_Reference()
            Implicit None
            Integer :: i
            Real*8 :: r_outer, r_inner, prefactor
			ref%density = 1.0d0
			ref%dlnrho = 0.0d0
			ref%d2lnrho = 0.0d0
			ref%pressure = 1.0d0
			ref%temperature = 1.0d0
			ref%dlnT = 0.0d0
			ref%dsdr = 0.0d0
			ref%pressure = 1.0d0
			ref%gravity = 0.0d0 ! Not used with constant reference right now
            ref%gravity_term_s = 0.0d0 ! Set to Ra later in equation_coefficients
			pressure_specific_heat = 1.0d0
			Call initialize_reference_heating()
            Allocate(s_conductive(1:N_R))
            r_outer = radius(1)
            r_inner = radius(N_R)
            prefactor = r_outer*r_inner/(r_inner-r_outer)
            Do i = 1, N_R
                s_conductive(i) = prefactor*(1.0d0/r_outer-1.0d0/radius(i))
            Enddo
	End Subroutine Constant_Reference

    Subroutine Read_Profile_File(filename,arr)
		Character*120, Intent(In) :: filename
        Character*120 :: ref_file, full_path
        Integer :: pi_integer,nr_ref, ncolumns
        Integer :: nqvals =9 ! Reference state file contains nqvals quantities + radius
        Integer :: i, k
        Real*8, Allocatable :: arr_old(:,:), rtmp(:), rtmp2(:)
        Real*8, Intent(InOut) :: arr(:,:)
        Real*8, Allocatable :: old_radius(:)
        full_path = Trim(my_path)//filename
        If (my_rank .eq. 0) Then
            !Only one processes actually opens the file
            !After that, the contents of the array are broadcast across columns and rows
		    Open(unit=15,file=full_path,form='unformatted', status='old',access='stream')
            Read(15)pi_integer

            If (pi_integer .ne. 314) Then
                close(15)
                Open(unit=15,file=full_path,form='unformatted', status='old', &
                     CONVERT = 'BIG_ENDIAN' , access='stream')
                Read(15)pi_integer
                If (pi_integer .ne. 314) Then
                    Close(15)
                    Open(unit=15,file=full_path,form='unformatted', status='old', &
                     CONVERT = 'LITTLE_ENDIAN' , access='stream')
                    Read(15)pi_integer
                Endif
            Endif
                    
            If (pi_integer .eq. 314) Then 

		        Read(15)nr_ref
                Read(15)ncolumns
                nqvals = ncolumns-1
                Allocate(arr_old(1:nr_ref,1:nqvals))  
                !10 quantities are stored in reference state

                Allocate(old_radius(1:nr_ref))

                Write(6,*)'nr_ref is: ', nr_ref

                Read(15)(old_radius(i),i=1,nr_ref)
                
                Do k = 1, nqvals
		            Read(15)(arr_old(i,k) , i=1 , nr_ref)
                Enddo
       		    

                !Check to see if radius isn't reversed
                !If it is not, then reverse it
                If (old_radius(1) .lt. old_radius(nr_ref)) Then
                    Write(6,*)'Reversing Radial Indices in Custom Ref File!'
                    Allocate(rtmp(1:nr_ref))
                    rtmp(:) = old_radius(:)
                    Do i = 1, nr_ref
                        old_radius(i) = rtmp(nr_ref-i+1)
                    Enddo

                    Do k = 1, nqvals
                        rtmp(:) = arr_old(:,k)
                        Do i = 1, nr_ref
                            arr_old(i,k) = rtmp(nr_ref-i+1)
                        Enddo
                    Enddo

                    DeAllocate(rtmp)

                Endif

                Close(15)


                If (nr_ref .ne. n_r) Then 
                    !Interpolate onto the current radial grid if necessary
                    !Note that the underlying assumption here is that same # of grid points
                    ! means same grid - come back to this later for generality
                    Allocate(rtmp2(1:n_r))
                    Allocate(rtmp(1:nr_ref))
                    
                    Do k = 1, nqvals
     
                        rtmp(:) = arr_old(:,k)
                        rtmp2(:) = 0.0d0
                        Call Spline_Interpolate(rtmp, old_radius, rtmp2, radius)
                        arr(1:n_r,k) = rtmp2 
                    Enddo

                    DeAllocate(rtmp,rtmp2)
                Else
                    ! Bit redundant here, but may want to do filtering on arr array
                    arr(1:n_r,1:nqvals) = arr_old(1:n_r,1:nqvals)

                Endif


                DeAllocate(arr_old,old_radius)
            Endif
        Endif


        If (my_row_rank .eq. 0) Then
            ! Broadcast along the column
            Call BCAST2D(arr,grp = pfi%ccomm)
        Endif
        Call BCAST2D(arr,grp = pfi%rcomm)

	End Subroutine Read_Profile_File

End Module ReferenceState
