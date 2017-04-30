#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t

!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_SECOND_DERIVATIVES
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Second_Derivatives
    Use Diagnostics_Base
    Implicit None

    ! Support variables (These should probably be wrapped into
    Integer :: dd_vr, dd_vt, dd_vp, dd_t, dd_p
    Integer :: dd_dvrdr, dd_dvtdr, dd_dvpdr, dd_dtdr, dd_dpdr
    Integer :: dd_dvrdt
    Type(SphericalBuffer) :: supbuff
    Integer, Allocatable :: ddindmap(:,:)
    Integer :: nddfields
Contains

Subroutine Initialize_Second_Derivatives
    Type(Field_Indexer) :: dd_indices
    CHARACTER*3 :: config
    INTEGER :: ndind
    INTEGER :: ddfcount(3,2)
    Logical :: compute_velocity_dd = .true.  ! change this later
    !First, indentify the primitive variables that we want to send over to spectral space
    config = 'p3b'

    nddfields = 0
    ndind = 1
    
    IF (compute_velocity_dd) THEN
        nddfields = nddfields+3  
    ENDIF

    !Add other "if compute" logic here

    Allocate(ddindmap(nddfields,2))


    !////////////////////////////////////////
    ! Next, there are two branches of "if computes"
    ! First, we 'load' all the radial derivatives
    If (compute_velocity_dd) THEN
        dd_dvrdr = nind
        ddindmap(1,ndind) = dd_dvrdr
        ddindmap(2,ndind) = dvrdr
        ndind += 1

        dd_dvtdr = nind
        ddindmap(1,ndind) = dd_dvtdr
        ddindmap(2,ndind) = dvtdr       
        ndind += 1

        dd_dvpdt = nind
        ddindmap(1,ndind) = dd_dvpdr
        ddindmap(2,ndind) = dvpdr       
        ndind += 1
    ENDIF

    ! IF compute thermal --> radial

    ! IF Compute mag -->  radial


    !/////////////////////////////////////////
    ! N

    Call dd_indices%Add_Field(dd_dvrdr, config)


    Call dd_indices%Add_Field(dd_dvrdt, config)

    !Call dd_indices%Add_Field(dd_vp, config)

    !Call dd_indices%Add_Field(dd_t, config)
    !Call dd_indices%Add_Field(dd_p, config)

    ! We may want to send the theta derivatives back to s2a and leave them there.
    ! This will avoid unecessary complications related to the sin_theta factors


    !Next, identify those variables we want to bring back from spectral space
    config = 'p1a'
    Call dd_indices%Add_Field(dvrdrdr, config)

    !Call dd_indices%Add_Field(dvtdrdr, config)
    !Call dd_indices%Add_Field(dvpdrdr, config)

    !Call dd_indices%Add_Field(dd_dvrdr, config)  ! These will be needed for computing
    !Call dd_indices%Add_Field(dd_dvtdr, config)  ! cross derivatives (e.g., d2_by_drdt)
    !Call dd_indices%Add_Field(dd_dvpdr, config)

    !Call dd_indices%Add_Field(dtdrdr, config)
    !Call dd_indices%Add_Field(dpdrdr, config)

    !Call dd_indices%Add_Field(dd_dtdr, config)
    !Call dd_indices%Add_Field(dd_dpdr, config)


    config='p2a'
    Call dd_indices%Add_Field(dvrdrdt, config)
    Call dd_indices%Add_Field(dvrdtdt, config)

    !Call dd_indices%Add_Field(dvtdrdt, config)
    !Call dd_indices%Add_Field(dvpdrdt, config)

    !Call dd_indices%Add_Field(dvrdtdt, config)
    !Call dd_indices%Add_Field(dvtdtdt, config)
    !Call dd_indices%Add_Field(dvpdtdt, config)


    config = 'p3a'
    Call dd_indices%Add_Field(dvrdrdp,config)
    Call dd_indices%Add_Field(dvrdtdp,config)
    Call dd_indices%Add_Field(dvrdpdp,config)

	ddfcount(1,1) = dd_indices%c1a_counter
	ddfcount(2,1) = dd_indices%c2a_counter
	ddfcount(3,1) = dd_indices%c3a_counter
    ddfcount(3,2) = dd_indices%c3b_counter
    ddfcount(2,2) = dd_indices%c3b_counter
    ddfcount(1,2) = dd_indices%c3b_counter


    Call d2buffer%init(field_count = ddfcount, config = 'p1a')
    Call supbuff%init(field_count = ddfcount, config = 'p1b')
End Subroutine Initialize_Second_Derivatives

Subroutine Compute_Second_Derivatives(inbuffer)
    Implicit None
    Real*8, Intent(InOut) :: inbuffer(1:,my_r%min:,my_theta%min:,1:)
    Real*8, Allocatable :: ddwork(:,:,:,:)

    Write(6,*)'Constructing the buffer...', d2buffer%nf3b
    Call d2buffer%construct('p3b')
    Write(6,*)'Complete...'

    ! First, load radial and theta derivatives
    Do i = 1, nddfields*2
        supbuff%p3b(:,:,:,ddindmap(1,i)) = inbuffer(:,:,:,ddindamp(2,i)) 
    Enddo


    !FFT and move to s2b configuration
    Call fft_to_spectral(supbuff%p3b, rsc = .true.)
    Call supbuff%reform()                          
    Call supbuff%construct('s2b')
    Call Legendre_Transform(supbuff%p2b,d2buffer%s2b)
    Call supbuff%deconstruct('p2b')

    ! ------>  might copy out dxdtheta here
    

    !Move to p1b configuration
    Call d2buffer%reform()
    ! Left off here -- these allocation dimensions are wrong. Probably better to just use a work
    ! buffer and avoid the headache
    Allocate(ddwork(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,d2buffer%nf3b))
    Call gridcp%To_Spectral(d2buffer%p1b,ddwork)

End Subroutine Compute_Second_Derivatives

End Module Diagnostics_Second_Derivatives
