#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
#define DO_IDX2 Do mp = my_mp%min, my_mp%max; m = m_values(mp); Do imi = 1, 2; Do r = my_r%min, my_r%max
#define IDX2 m:l_max,r,imi
!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_SECOND_DERIVATIVES
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Second_Derivatives
    Use Diagnostics_Base
    Use Structures
    Use Spectral_Derivatives
    Implicit None


    Integer, Allocatable :: ddindmap(:,:)
    Integer :: nddfields
    Logical :: compute_vr_dd   = .false. 
    Logical :: compute_vt_dd   = .false. 
    Logical :: compute_vp_dd   = .false.

    Logical :: compute_pvar_dd = .false.
    Logical :: compute_tvar_dd = .false.

    Logical :: compute_br_dd   = .false. 
    Logical :: compute_bt_dd   = .false. 
    Logical :: compute_bp_dd   = .false.
Contains

    Subroutine Init_Derivative_Logic()
        IMPLICIT NONE

        ! Execute a lot of compute_q logic here to see if the
        ! different compute_xx_dd variables should be set to true.

        If (compute_vr_dd) need_second_derivatives = .true.
        If (compute_vt_dd) need_second_derivatives = .true.       
        If (compute_vp_dd) need_second_derivatives = .true. 

        If (compute_tvar_dd) need_second_derivatives = .true.       
        If (compute_pvar_dd) need_second_derivatives = .true. 

        If (compute_br_dd) need_second_derivatives = .true.
        If (compute_bt_dd) need_second_derivatives = .true.       
        If (compute_bp_dd) need_second_derivatives = .true.


    End Subroutine Init_Derivative_Logic

    Subroutine Initialize_Second_Derivatives()
        ! Initializes all the indexing related to computing and
        ! accessing second derivatives at output time.
        ! Most of the actual indexing is handled by Set_DD_Indices
        IMPLICIT NONE
        INTEGER :: ndind
        INTEGER :: ddfcount(3,2)


        Call Init_Derivative_Logic()


        nddfields = 0   ! Number of fields whose second derivatives we want
        ndind = 0       ! internal indexing variable

        IF (compute_vr_dd) nddfields = nddfields +1        
        IF (compute_vt_dd) nddfields = nddfields +1  
        IF (compute_vp_dd) nddfields = nddfields +1  

        IF (compute_tvar_dd) nddfields = nddfields +1
        IF (compute_pvar_dd) nddfields = nddfields +1

        IF (magnetism) THEN
            IF (compute_br_dd) nddfields = nddfields +1        
            IF (compute_bt_dd) nddfields = nddfields +1  
            IF (compute_bp_dd) nddfields = nddfields +1  
        ENDIF




        Allocate(ddindmap(4,nddfields*2))
        ddindmap(:,:) = -1


        If (compute_vr_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dvrdrdr, dvrdrdt, dvrdrdp, &
                                      dvrdtdt, dvrdtdp, dvrdpdp, &
                                      dvrdr  , dvrdt  , dvrdp)
        Endif

        IF (compute_vt_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dvtdrdr, dvtdrdt, dvtdrdp, &
                                      dvtdtdt, dvtdtdp, dvtdpdp, &
                                      dvtdr  , dvtdt  , dvtdp)
        ENDIF

        IF (compute_vp_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dvpdrdr, dvpdrdt, dvpdrdp, &
                                      dvpdtdt, dvpdtdp, dvpdpdp, &
                                      dvpdr  , dvpdt  , dvpdp)
        ENDIF

        IF (compute_tvar_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dtdrdr, dtdrdt, dtdrdp, &
                                      dtdtdt, dtdtdp, dtdpdp, &
                                      dtdr  , dtdt  , dtdp)
        ENDIF

        IF (compute_pvar_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dpdrdr, dpdrdt, dpdrdp, &
                                      dpdtdt, dpdtdp, dpdpdp, &
                                      dpdr  , dpdt  , dpdp)
        ENDIF

        If (magnetism) THEN
            If (compute_br_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, ddindmap, &
                                          dbrdrdr, dbrdrdt, dbrdrdp, &
                                          dbrdtdt, dbrdtdp, dbrdpdp, &
                                          dbrdr  , dbrdt  , dbrdp)
            Endif

            IF (compute_bt_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, ddindmap, &
                                          dbtdrdr, dbtdrdt, dbtdrdp, &
                                          dbtdtdt, dbtdtdp, dbtdpdp, &
                                          dbtdr  , dbtdt  , dbtdp)
            ENDIF

            IF (compute_bp_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, ddindmap,      &
                                          dbpdrdr, dbpdrdt, dbpdrdp, &
                                          dbpdtdt, dbpdtdp, dbpdpdp, &
                                          dbpdr  , dbpdt  , dbpdp)
            ENDIF
        ENDIF



	    ddfcount(1,1) = nddfields*4 ! config 1a
	    ddfcount(2,1) = nddfields*4 ! 2a
	    ddfcount(3,1) = nddfields*7 ! 3a
        ddfcount(3,2) = nddfields*2 ! 3b
        ddfcount(2,2) = nddfields*2 ! 2b
        ddfcount(1,2) = nddfields*2 ! 1b


        Call d2buffer%init(field_count = ddfcount, config = 'p3b')

        Call d2buffer%construct('p3a')
        WRITE(6,*)'BCHECK: ', shape(d2buffer%p3a)
        Call d2buffer%deconstruct('p3a')
    End Subroutine Initialize_Second_Derivatives

    Subroutine Compute_Second_Derivatives(inbuffer)
        Implicit None
        INTEGER :: i,j, imi, mp, m 
        INTEGER :: r,k, t
        Real*8, Intent(InOut) :: inbuffer(1:,my_r%min:,my_theta%min:,1:)
        Type(rmcontainer3D), Allocatable :: ddtemp(:)

        ! Here were compute all second derivatives for N variables
        ! Outline of the process:
        ! 1)  Intialize p3b space of d2buffer
        ! 2)  Load slots [1 : N] of d2buffer with d_by_dr for each variable
        ! 3)  Load slots [N+1 : 2N] with d_by_dtheta for each variable
        ! 4)  Move to p1b/p1a configuration
        ! 5)  Load slots [2N+1 : 3N] with d2_by_dr2 for each variable
        ! 6)  Load slots [3N+1 : 4N] with d2_by_drdt for each variable; move to s2a
        ! 7)  Move to s2a and load sintheta*dxdtdt into slots [N+1 : 2N] (ovewriting dxdt)
        ! 8)  Move to p3a and load d2_by_drdphi into slots [4N+1: 5N] for each variable
        ! 9)  Load d2_by_dphi2 into slots [5N+1: 6N] for each variable
        ! 10) Load d2_by_dtdphi into slots [6N+1 : 7N] for each variable
        ! 11) FFT, correct for sintheta factor, compute means and fluctuations

        ! When this routine is complete, the contents of d2buffer%p3a will be
        ! [   1 : N  ] -- workspace; 
        ! [ N+1 : 2N ] -- dxdtdt
        ! [2N+1 : 3N ] -- dxdrdr
        ! [3N+1 : 4N ] -- dxdrdt
        ! [4N+1 : 5N ] -- dxdrdp
        ! [5N+1 : 6N ] -- dxdpdp
        ! [6N+1 : 7N ] -- dxdtdp

        !///////////////////////////////////////////////////////////
        ! Step 1: Initialize p3b portion of d2buffer
        !Write(6,*)'Constructing the buffer...', d2buffer%nf3b

        Call d2buffer%construct('p3b')
        d2buffer%config = 'p3b'

        !Write(6,*)'Complete...'

        !///////////////////////////////////////////////////////////
        ! Steps 2-3:  Load radial and theta derivatives
        Do i = 1, nddfields*2
            d2buffer%p3b(:,:,:,ddindmap(1,i)) = inbuffer(:,:,:,ddindmap(2,i)) 
        Enddo

        !Write(6,*)'Steps 1-3 complete.'

        !////////////////////////////////////////////////////////////////
        ! Step 4:  Move to p1b/p1a configuration
        Call fft_to_spectral(d2buffer%p3b, rsc = .true.)
        Call d2buffer%reform()                          
        Call d2buffer%construct('s2b')
        Call Legendre_Transform(d2buffer%p2b,d2buffer%s2b)
        Call d2buffer%deconstruct('p2b')
        d2buffer%config ='s2b'



        !Move to p1b configuration
        Call d2buffer%reform()
        Call d2buffer%construct('p1a')
        Call gridcp%To_Spectral(d2buffer%p1b,d2buffer%p1a)
        Call d2buffer%deconstruct('p1b')  ! Note, we could use p1b as a buffer space to store derivatives
        d2buffer%config='p1a'

        !Write(6,*)'Steps 4 complete.'

        !////////////////////////////////////////////////////////////
        ! Steps 5-6:  Load d2_by_dr2 and d2_by_drdt into the buffer
        Do i = 1, nddfields*2
            j = i+nddfields*2
            Call gridcp%d_by_dr_cp(i,j,d2buffer%p1a,1)
        Enddo

        ! Ordering of fields in buffer is now dxdr, dxdt, dxdrdr, dxdrdt
        !Write(6,*)'Steps 5-6 complete.'

        !///////////////////////////////////////////////////////////////
        ! Step 7:  Move to s2a & overwrite dxdt with sintheta*{dxdtdt}
        Call d2buffer%reform()
        
        !Write(6,*)'Step 7 reformation complete'

        Call Allocate_rlm_Field(ddtemp)

        !Write(6,*)'Step 7 allocation complete', nddfields

        Do i = nddfields+1,nddfields*2
            ! We overwrite dxdt with sintheta* {dxdtdt}
            Call d_by_dtheta(d2buffer%s2a,i,ddtemp)
            DO_IDX2
                d2buffer%s2a(mp)%data(IDX2,i) = ddtemp(mp)%data(IDX2)
            END_DO
        Enddo

        !Write(6,*)'Step 7 derivatves complete'
        
        Call DeAllocate_rlm_Field(ddtemp)

        !Write(6,*)'Step 7 deallocation complete'

		Call d2buffer%construct('p2a')
		Call Legendre_Transform(d2buffer%s2a,d2buffer%p2a)
		Call d2buffer%deconstruct('s2a')
		d2buffer%config = 'p2a'	

        !Write(6,*)'Step 7 complete.'

        !/////////////////////////////////////////////////////////////////////
        !  Steps 8-10 : phi derivatives        
        Call d2buffer%reform() ! move to p3a

        ! Ordering of fields in buffer is now dxdr, sintheta*{dxdtdt}, dxdrdr, dxdrdt

        !Compute dxdrdp
        Do i = 1, nddfields
            j = i+nddfields*4
            Call d_by_dphi(d2buffer%p3a,i,j)
        Enddo

        !Copy dxdp into dxdr space, then calculate dxdpdp
        Do i = 1, nddfields
            j = i+nddfields*5
            d2buffer%p3a(:,:,:,i) = inbuffer(:,:,:,ddindmap(3,i)) 
            Call d_by_dphi(d2buffer%p3a,i,j)
        Enddo

        !Copy dxdt into dxdr space, then calculate dxdtdp
        Do i = 1, nddfields
            j = i+nddfields*6
            d2buffer%p3a(:,:,:,i) = inbuffer(:,:,:,ddindmap(4,i)) 
            Call d_by_dphi(d2buffer%p3a,i,j)
        Enddo

       ! Write(6,*)'Steps 8-10 complete.'

        !//////////////////////////////////////////
        ! Step 11:   Finalize
        ! FFT
        Call fft_to_physical(d2buffer%p3a,rsc = .true.)

        ! Convert sintheta*{dxdtdt} to dxdtdt
        Do i = nddfields*2+1,nddfields*3
		    DO_PSI
			    d2buffer%p3a(PSI,i) = d2buffer%p3a(PSI,i)*csctheta(t)	
		    END_DO
        Enddo
        !D2buffer is now initialized.  Ordering of fields is:
        ! [dxdt], dxdtdt,dxdrdr,dxdrdt, dxdrdp, dxdpdp, dxdtdp
        ! Note that we have one redundant derivative, dxdt, stored for each field

        ! Now compute the means and fluctuations
        Allocate(d2_ell0(my_r%min:my_r%max,1:nddfields*7))
        Allocate(d2_m0(my_r%min:my_r%max,my_theta%min:my_theta%max,1:nddfields*7))
        Allocate(d2_fbuffer(1:n_phi,my_r%min:my_r%max, &
            my_theta%min:my_theta%max,1:nddfields*7))

        Call ComputeEll0(d2buffer%p3a,d2_ell0)
        Call   ComputeM0(d2buffer%p3a,d2_m0)

        DO j = 1,nddfields*7
            DO_PSI
                d2_fbuffer(PSI,j) = d2buffer%p3a(PSI,j) - d2_m0(PSI2,j) 
            END_DO
        ENDDO

        !Write(6,*)'Step 11 complete.'
    End Subroutine Compute_Second_Derivatives

	Subroutine Allocate_rlm_Field(arr)
		Implicit None
		Type(rmcontainer3D), Intent(InOut), Allocatable :: arr(:)
		Integer :: mp,m


		Allocate(arr(my_mp%min:my_mp%max))
		Do mp = my_mp%min, my_mp%max
			m = m_values(mp)
			Allocate(arr(mp)%data(m:l_max,my_r%min:my_r%max,1:2))
			arr(mp)%data(:,:,:) = 0.0d0
		Enddo
	End Subroutine Allocate_rlm_Field

	Subroutine DeAllocate_rlm_Field(arr)
		Implicit None
		Type(rmcontainer3D), Intent(InOut), Allocatable :: arr(:)
		Integer :: mp
		Do mp = my_mp%min, my_mp%max
			DeAllocate(arr(mp)%data)
		Enddo
		DeAllocate(arr)
	End Subroutine DeAllocate_rlm_Field

    Subroutine Set_DD_Indices(iind, nskip, iindmap, & 
                                     dxdrdr, dxdrdt, dxdrdp, &
                                     dxdtdp, dxdtdt, dxdpdp, &
                                       dxdr,   dxdt, dxdp)
        ! Sets indexing within indmap and assigned values to
        ! dxdidj consistent with the logic used in 
        ! Compute_Second_Derivatives()
        ! [ N+1 : 2N ] -- dxdtdt
        ! [2N+1 : 3N ] -- dxdrdr
        ! [3N+1 : 4N ] -- dxdrdt
        ! [4N+1 : 5N ] -- dxdrdp
        ! [5N+1 : 6N ] -- dxdpdp
        ! [6N+1 : 7N ] -- dxdtdp                                           
        INTEGER, Intent(In)    :: iind, nskip
        INTEGER, INTENT(OUT)   :: dxdtdt, dxdrdr, dxdrdt
        INTEGER, INTENT(OUT)   :: dxdrdp, dxdpdp, dxdtdp
        INTEGER, INTENT(IN)    :: dxdr, dxdt, dxdp
        INTEGER, INTENT(INOUT) :: iindmap(:,:)
        dxdrdt = iind+nskip
        dxdrdr = iind+nskip*2
        dxdrdt = iind+nskip*3
        dxdrdp = iind+nskip*4
        dxdpdp = iind+nskip*5
        dxdtdp = iind+nskip*6

        iindmap(1,iind          ) = iind
        iindmap(2,iind          ) = dxdr
        iindmap(3,iind          ) = dxdp
        iindmap(4,iind          ) = dxdt
        iindmap(1,iind          ) = iind+nskip
        iindmap(1,iind+nddfields) = iind
        iindmap(2,iind+nddfields) = dxdt


    End Subroutine Set_DD_Indices
End Module Diagnostics_Second_Derivatives
