
/*
  Section and End keyword definitions
*/

#define    TSECTION             1
#define    TEND                 2


/* 
  Input parameters file:  Keyword definitions
*/

#define    TGRIDFILE            3

#define    TTHERMENERGY         4
#define    TBUOYANCY            5
#define    TSMALLSCALE          6
#define    TTRAJECTORIES        7

#define    TNTIMESTEPS          8
#define    TTIMESTEPSIZE        9

#define    TMAXQLITER          10
#define    TQLTOLERANCE        11

#define    TMAXMEITER          12
#define    TMETOLERANCE        13

#define    TFILTERU            14
#define    TFILTERV            15
#define    TFILTERT            16

#define    TPPESOLVER          17
#define    TMAXSORITER         18
#define    TSORTOLERANCE       19
#define    TSORRELAX           20

#define    TATDCU              21
#define    TATDTSC             22
#define    TATDHSC             23
#define    TATDBNCR            24
#define    TATDRMAX            25
#define    TATDRMEXP           26
#define    TATDTEM             27

#define    TATDPPESLVR         28
#define    TATDMAXSORIT        29
#define    TATDSORTOL          30
#define    TATDSORRELAX        31
#define    TATDFILTPARU        32
#define    TATDFILTPARV        33
#define    TATDFILTPART        34

#define    TREFLENGTH          35
#define    TREFVELOCITY        36
#define    TREFTEMPER          37
#define    TMAXREFTEMPER       38

#define    TFLUIDPROPFILE      39
#define    TFLUIDPROPTAB       40

#define    TRESTART            41
#define    TWRITERESTART       42

#define    TTIMEAVGMODE        43

#define    TOUTPUTFORMAT       44
#define    TOUTPUTPREFIX       45

#define    TSNAPSHOTSPREF      46
#define    TSNAPSHOTSFREQ      47
#define    TSNAPSHOTSFRST      48

#define    TTIMESRSPREFIX      49
#define    TTIMESRSFREQ        50
#define    TTIMESRSPOINT       51

#define    TPRINTDIFFFREQ      52
#define    TPRINTBANNFREQ      53

/* 
  Boundary conditions file:  Keyword definitions
*/

#define    TNUMOFREGIONS       3
#define    TIBORDERS           4
#define    TJBORDERS           5

#define    TBLOCKAGE           6
#define    TPOROUSREGION       7
#define    TWALL               8
#define    TINLET              9
#define    TOUTLET            10

#define    THEATGENREGION     11
#define    TTEMPERREGION      12
#define    TFIXEDTEMPWALL     13
#define    THEATFLUXWALL      14


/* 
  Note that the following compiler directive is required for jnum, 
  dnum, etc..
  This implies that the function as it stands now might only work
  with HP fortran.
*/

#ifdef _HP_

$HP9000_800 INTRINSICS

#endif _HP_

