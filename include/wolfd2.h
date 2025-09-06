/*
  Types definitions
*/

#define   REAL     real*8
#define   FLOAT    real*8
#define   INTEGER  integer*4
#define   LOGICAL  logical*4


/*
  Region and boundary condition definitions
*/

#define  RM_BLOCKG   0
#define  RM_INTERN   1
#define  RM_POROUS   2

#define  BM_INTERN   0
#define  BM_WALL1    1
#define  BM_WALL2    2
#define  BM_INLET    3
#define  BM_OUTLT1   4
#define  BM_OUTLT2   5

#define  RT_NOSRCE   0
#define  RT_HEATGN   1
#define  RT_TEMPER   2

#define  BT_INTERN   0
#define  BT_TEMPER   1
#define  BT_HTFLUX   2


/*
  Region wall locations
*/

#define  WEST        1
#define  EAST        2
#define  SOUTH       3
#define  NORTH       4


/*
  Variable name pseudonyms
*/

#define  _U_         1
#define  _V_         2
#define  _P_         3
#define  _T_         4

/*
  Grid direction pseudonyms
*/

#define  _I_         1
#define  _J_         2


/*
  Logical pseudonyms
*/

#define   LFALSE   .false.
#define   LTRUE    .true.


/*
  Standard output, error and log file unit definition
*/

#define   STDOUT     6
#define   STDERR     7
#define   STDLOG     8

/*
  File format type definitions
*/

#define  FT_FORMATTED      0
#define  FT_UNFORMATTED    1

