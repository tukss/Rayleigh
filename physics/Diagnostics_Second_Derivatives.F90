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
    Use Diagnostics_ADotGradB
    Implicit None

Contains

Subroutine Initialize_Second_Derivatives
    Type(Field_Indexer) :: dd_indices

    !First, indentify the primitive variables that we want to send over to spectral space
    config = 'p3b'

    Call dd_indices%Add_Field(dd_vr, config)
    Call dd_indices%Add_Field(dd_vt, config)
    Call dd_indices%Add_Field(dd_vp, config)

    Call dd_indices%Add_Field(dd_t, config)
    Call dd_indices%Add_Field(dd_p, config)

    ! We may want to send the theta derivatives back to s2a and leave them there.
    ! This will avoid unecessary complications related to the sin_theta factors


    !Next, identify those variables we want to bring back from spectral space
    config = 'p1a'
    Call dd_indices%Add_Field(dvrdrdr, config)
    Call dd_indices%Add_Field(dvtdrdr, config)
    Call dd_indices%Add_Field(dvpdrdr, config)

    Call dd_indices%Add_Field(dd_dvrdr, config)  ! These will be needed for computing
    Call dd_indices%Add_Field(dd_dvtdr, config)  ! cross derivatives (e.g., d2_by_drdt)
    Call dd_indices%Add_Field(dd_dvpdr, config)

    Call dd_indices%Add_Field(dtdrdr, config)
    Call dd_indices%Add_Field(dpdrdr, config)

    Call dd_indices%Add_Field(dd_dtdr, config)
    Call dd_indices%Add_Field(dd_dpdr, config)


    config='p2a'
    Call dd_indices%Add_Field(dvrdrdt, config)
    Call dd_indices%Add_Field(dvtdrdt, config)
    Call dd_indices%Add_Field(dvpdrdt, config)

    Call dd_indices%Add_Field(dvrdtdt, config)
    Call dd_indices%Add_Field(dvtdtdt, config)
    Call dd_indices%Add_Field(dvpdtdt, config)



    Call d2buffer%init(field_count = wsfcount, config = 'p1b')
End Subroutine Initialize_Second_Derivatives


End Module Diagnostics_Second_Derivatives
