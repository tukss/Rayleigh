Module Fields
	Use Linear_Solve
	Use Parallel_Framework
	Use ProblemSize
	Use Controls
	!///////////////////////////////////////////////////////////
	! The scalable framework of this program is centered around
	! a large buffer that holds multiple fields.
	! The buffer is reshaped during each transpose, and fields are
	! added and removed as the buffer moves amongst the different
	! configuration spaces.

	! The field type is meta data that allows us to know where
	! a field is in the buffer at a particular time, or if it's
	! even there at all
	Implicit None
	Type Field
		Integer :: i1a, i2a, i3a
		Integer :: i1b, i2b, i3b
		Integer :: ip3a = -1	 ! index of location in p3a
								    ! (physical space location where NL terms are computed)
		Integer :: is2a = -1  ! index of location in s2b (where theta derivatives are computed)
	End Type

	Type(Field), Allocatable :: vars(:)
	Integer :: c3a_counter = 0
	Integer :: c2a_counter = 0
	Integer :: c1a_counter = 0
	Integer :: c3b_counter = 0
	Integer :: c2b_counter = 0
	Integer :: c1b_counter = 0
	Integer :: global_counter = 0
	!////////////////////////////////////////////////////////////////////////////
	! Variable locations in the global field meta data buffer
	Integer :: wvar, pvar, tvar, zvar
	Integer :: dpdr
	Integer :: dwdr,  d3wdr3, dtdr, dzdr, d2zdr2, d2tdr2, d2wdr2
	Integer :: vr, vtheta, vphi,dtdt,dvrdt,dvtdr,dvpdr,dvrdr
	Integer :: dvrdp, dvtdp, dvpdp, dtdp

	Integer :: avar, dadr, d2adr2 ! Toroidal magnetic streamfunction and its radial derivatives
	Integer :: cvar, dcdr, d2cdr2 !Poloidal magnetic streamfunction and its radial derivatives

	Integer :: br,btheta,bphi, jr,jtheta,jphi	
	Integer :: emfr,emftheta,emfphi
	!///////////////////////////////////////////////////////////////////////////

	!///////////////////////////////////////////////////////////////////////////
	!  Indices for easy reference in physical space
	Integer :: vr_ps, vt_ps, vp_ps, t_ps ! etc.
	Integer :: nl_advec_r, nl_advec_t, nl_advec_p, temp_advec

	!==============================================================================
	!		Integer parameters used to reference specific equations and variables
	!				FOR LINEAR Solve
	!==============================================================================
	Integer, parameter :: weq = 1,  peq = 2,  teq = 3
	Integer, parameter :: zeq = 4,  ceq = 5,  aeq = 6


	Integer :: wsfcount(3,2)
	Real*8, Allocatable :: ABsave(:,:,:,:) !  Have to think about where to do this
	! Type(rmcontainer) :: ABSave(:)			! Might be most natural to keep this in rlm space
	Type(SphericalBuffer) :: wsp ! Primary workspace for the entire run

Contains

	Subroutine Add_Field(ind, arr,overwrite)
		Implicit None
		Integer, intent(InOut) :: ind
		Integer, Intent(in) :: arr(1:6)
		Integer, Intent(in), Optional :: overwrite


		If (present(overwrite)) Then
				!vars(ind)%i1a = vars(overwrite)%i1a
				!vars(ind)%i2a = vars(overwrite)%i2a
				!vars(ind)%i3a = vars(overwrite)%i3a
				!vars(ind)%i1b = vars(overwrite)%i1b
				!vars(ind)%i2b = vars(overwrite)%i2b
				!vars(ind)%i3b = vars(overwrite)%i3b
				ind = overwrite
		Else
			global_counter = global_counter+1
			ind = global_counter
			if (arr(1) .eq. 1) then 
				c1a_counter = c1a_counter+1
				vars(ind)%i1a = c1a_counter
			endif
			if (arr(2) .eq. 1) then 
				c2a_counter = c2a_counter+1
				vars(ind)%i2a = c2a_counter
			endif
			if (arr(3) .eq. 1) then
				c3a_counter = c3a_counter+1
				vars(ind)%i3a = c3a_counter
			endif
			if (arr(4) .eq. 1) then 
				c1b_counter = c1b_counter+1
				vars(ind)%i1b = c1b_counter
			endif
			if (arr(5) .eq. 1) then 
				c2b_counter = c2b_counter+1
				vars(ind)%i2b = c2b_counter
			endif
			if (arr(6) .eq. 1) then
				c3b_counter = c3b_counter+1
				vars(ind)%i3b = c3b_counter
			endif
		Endif
	End Subroutine Add_Field



	Subroutine Initialize_Field_Structure()
		Implicit None
		Integer :: arr(1:6)
		Integer :: max_fields = 30

		! Much of this is really waaay to convoluted.  We can do most of what we need
		! to do with simple integers and get rid of this Add_Field business.   TODO


		!Some care needs to be taken here in the order of the loading

		!////////////////////////////////////////////////////////////
		!	Fields and radial derivatives
		!  To start off with, I just want to get these out to physical space

		Allocate(vars(1:max_fields))


		!//////////////////////////////
		! It is rather important that the primary fields
		! we solve for are added first
		! Their numbers should be 1 through nfields and it is OK
		! if they are not used in all configurations.  We want the 
		! equation numbering and the field numbering to agree


		!//////////////////////////////////
		!  First, we need an accounting of all fields that will be stored
		!  (even temporarily) in the p1a buffer.  We will assume that
		!  These fields persist out to p3a for now
		arr(1) = 1		! clumsy right now, but let's just get this thing started
		arr(2) = 1
		arr(3) = 1		
		arr(4) = 0
		arr(5) = 0
		arr(6) = 0

		! Add the primary fields
		Call Add_Field(Wvar,  arr)
		Call Add_Field(Pvar,  arr)	! Pressure doesn't really need to be transposed, but doing so now for debugging
		Call Add_Field(Tvar,  arr)
		Call Add_Field(Zvar,  arr)
		if (magnetism) Then
		  Call Add_Field(cvar,arr)
		  Call Add_field(avar,arr)
		Endif
		

		Call Add_Field(d3Wdr3  ,arr)
		Call Add_Field(dPdr    ,arr)
		Call Add_Field(d2Tdr2  ,arr, overwrite = dpdr)	! replaces dpdr
		Call Add_field(d2Zdr2  ,arr)

		Call Add_Field(dWdr    ,arr, overwrite = d3wdr3)
		Call Add_Field(d2Wdr2  ,arr)
		Call Add_Field(dTdr    ,arr, overwrite = dpdr)	!replaces d2tdr2, which replaced dpdr
		Call Add_Field(dZdr    ,arr, overwrite = d2Zdr2)

		if (magnetism) Then
		  Call Add_field(dcdr,arr)
		  Call Add_field(dadr,arr)
		  Call Add_field(d2adr2,arr)
		  Call Add_field(d2cdr2,arr,overwrite = d2adr2)
		Endif


		!//////////////////////////////////////////////////////////
		! Next, we want to account for fields that we build in s2a (many are d by dtheta fields)
		arr(1) = 0
		Call Add_Field(vr,arr,overwrite = wvar) ! We only needed W to create Vr
		Call Add_Field(vtheta,arr)
		Call Add_Field(vphi,arr) 
		Call Add_Field(dvtdr,arr)
		Call Add_Field(dvpdr,arr)
		Call Add_Field(dvrdr,arr,overwrite = dwdr)
		Call Add_Field(dtdt,arr, overwrite = d2wdr2)
		Call Add_Field(dvrdt,arr, overwrite = dzdr)
		If (magnetism) Then
			Call Add_Field(Br     ,arr, overwrite = cvar)
			Call Add_Field(Btheta ,arr, overwrite = avar)
			Call Add_Field(Bphi   ,arr, overwrite = dcdr)

			Call Add_Field(Jr     ,arr)			! We need to make an extra slot.  Jr overwrites nothing.
			Call Add_Field(Jtheta ,arr, overwrite = d2cdr2)
			Call Add_Field(Jphi   ,arr, overwrite = dadr)
		Endif
		!///////////////////////////////////////////////////////
		!  Finally we have fields of the d_by_dphi variety
		!  that we add to the p3a buffer
		arr(2) = 0
		Call Add_Field(dvrdp,arr)
		Call Add_field(dvtdp,arr)
		Call Add_field(dvpdp,arr)
		Call Add_field(dtdp,arr,overwrite = tvar)		

		wsfcount(1,1) = c1a_counter
		wsfcount(2,1) = c2a_counter
		wsfcount(3,1) = c3a_counter

		if (.not. magnetism) then
			wsfcount(1,2) = 4		! four RHS's go back for the solve
			wsfcount(2,2) = 4
			wsfcount(3,2) = 4
		else
			emfr = avar
			emftheta = cvar
			emfphi = avar+1
			wsfcount(1,2) = 7		! seven RHS's go back for the solve (1 field is differentiated and combined at the end)
			wsfcount(2,2) = 7
			wsfcount(3,2) = 7
		endif
		!Write(6,*)'Fields initialized'
		!Write(6,*)'c1a_counter is: ', c1a_counter
		!Write(6,*)'c2a_counter is: ', c2a_counter
		!Write(6,*)'c3a_counter is: ', c3a_counter

	End Subroutine Initialize_Field_Structure


End Module Fields
