*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                               config.f                               *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

C---> Declaration of configurable parameters (Do not modify)
      INTEGER mnx      ,  ! Maximum number of points in I-direction
     .        mny      ,  ! Maximum number of points in J-direction
     .        mgri     ,  ! Maximum number of regions in I-direction
     .        mgrj     ,  ! Maximum number of regions in J-direction
     .        mfnmlgth ,  ! Maximum length of filename strings
     .        mlinelgt ,  ! Maximum length of parsing line strings
     .        mntr     ,  ! Maximum number of trajectories
     .        maxtspts    ! Maximum number of time-series check points


C---: Definition of configurable parameters (Modify as required)
      parameter ( mnx      =    302 )
      parameter ( mny      =    302 )
      parameter ( mgri     =     20 )
      parameter ( mgrj     =     10 )
      parameter ( mfnmlgth =     32 )
      parameter ( mlinelgt =     80 )
      parameter ( mntr     =      1 )
      parameter ( maxtspts =     20 )

C-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
