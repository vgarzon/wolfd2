*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                                parse.F                               *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* lReadIPFile   -
* nFindTokBC    -
* lReadBCFile   -
* nFindTokIP    -
* lReadFldProps -
* InitBCFlags   -
* mCellLocStr   -
* lCheckIJreg   -
*
*----------------------------------------------------------------------*

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"
#include "parse.h"


*----------------------------------------------------------------------*
*                              lReadIPFile                             *
*----------------------------------------------------------------------*
*
      LOGICAL function
     .   lReadIPFile(inputFile, gridFile, flPropFn, flName,
     .               resReadFl, resWriteFl, 
     .               outPrefix, ssPrefix, TSPrefix,
     .               tAvgOutPref,
     .               lCartesGrid, lTmAvgFlg,
     .               nthermen, neqstate, nsmallscl, ntraject,
     .               nrestart, nreswrite, nsnapshots, nTSOutFlag, 
     .               lPrntDiff, lPrDfBann,
     .               nts, mqiter, nmeiter, nfiltu, nfiltv, nfiltt, 
     .               nPpeSolver, msorit, 
     .               nssPpeSlvr, mssSorIt,
     .               nPostProc, nOutFlForm, nWrtRstFq,
     .               isfreq, isfirst,
     .               nTSFreq, nTSPoints, nPrDfFreq, nPrDfBnFq,
     .               iTS, jTS,
     .               dk, qtol, dmeittol, fpu, fpv, fpt, 
     .               sortol, sorrel,
     .               ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .               ssBnCrit, ssRMpMax, ssRMpExp,
     .               ssSorTol, ssSorRel,
     .               ssFiltPar,
     .               dlref, uref, tref, tmax)

      implicit none

      include "config.f"

      character*(*) inputFile

C---: Returning variables
      character*(*) gridFile, flPropFn, flName, resReadFl, resWriteFl
      character*(*) outPrefix
      character*(*) ssPrefix, TSPrefix
      character*(*) tAvgOutPref

      LOGICAL lCartesGrid, lTmAvgFlg

      INTEGER nthermen, neqstate, nsmallscl, ntraject
      INTEGER nrestart, nreswrite, nsnapshots, nTSOutFlag,
     .        lPrntDiff, lPrDfBann
      INTEGER nts, mqiter, nmeiter, nfiltu, nfiltv, nfiltt, 
     .        nPpeSolver, msorit, nPostProc, nOutFlForm, nWrtRstFq,
     .        isfreq, isfirst,
     .        nTSFreq, nTSPoints, nPrDfFreq, nPrDfBnFq,
     .        nssPpeSlvr, mssSorIt

      INTEGER iTS(*), jTS(*)

      REAL    dk, qtol, dmeittol, fpu, fpv, fpt, sortol, sorrel,
     .        ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .        ssBnCrit, ssRMpMax, ssRMpExp,
     .        ssSorTol, ssSorRel,
     .        dlref, uref, tref, tmax

      REAL    ssFiltPar(*)

C---: Functions used locally
      LOGICAL lCheckFile, lStIsInt, lStIsReal
      INTEGER mTrimLeft, mExtractW
      INTEGER nFindTokIP

#ifndef _HP_

      INTEGER jnum
      REAL    dnum

#endif _HP_

C---: Local variables
      LOGICAL lInSection, lFoundSection

      character*(mlinelgt) cBuf, cw
      INTEGER i, nError, nl, nstat, nw, nc, nTok

C---: Default initial state
      lReadIPFile = .true.

C---: Check that input file is readable
      if(.not.lCheckFile(inputFile,nError)) then
        print*, 'Error code :', nError
        lReadIPFile = .false.
        return
      end if

C---: Initialize flags
      nthermen   = 0
      neqstate   = 0
      nsmallscl  = 0
      ntraject   = 0
      nrestart   = 0
      nreswrite  = 0
      nfiltu     = 0
      nfiltv     = 0
      nfiltt     = 0
      nTSOutFlag = 0
      lPrntDiff  = 0
      lPrDfBann  = 0

C---: Initialize integer variables
      nts        = 0
      mqiter     = 20
      nmeiter    = 1
      nPpeSolver = 1
      msorit     = 2000
      nPostProc  = 2
      nOutFlForm = FT_FORMATTED
      nTSFreq    = 100
      nTSPoints  = 0 
      nPrDfFreq  = 1
      nPrDfBnFq  = 20
      nWrtRstFq  = 100 
      nssPpeSlvr = 1
      mssSorIt   = 2000


C---: Initialize TSPoints arrays
      do i=1, maxtspts
        iTS(i) = 1
        iTS(i) = 1
      end do

C---: Initialize double precision variables
      qtol       = 1.d-4
      dmeittol   = 1.d-4
      fpu        = 5.d2
      fpv        = 5.d2
      fpt        = 5.d0
      sortol     = 1.d-8
      sorrel     = 1.d0
      ssCu0      = 0.d0
      ssTsCoef   = 1.d0
      ssHsCoef   = 1.d0
      ssTemCoef  = 1.d0
      ssBnCrit   = 2.d2 
      ssRMpMax   = 0.95d0
      ssRMpExp   = 5.d0
      dlref      = 1.d0
      uref       = 1.d0
      tref       = 2.8357d2
      tmax       = 1.d3
      ssSorTol   = 1.d-8
      ssSorRel   = 1.d0

      ssFiltPar(_U_) = 5.d2
      ssFiltPar(_V_) = 5.d2
      ssFiltPar(_T_) = 5.d2


C---: Initialize character variables
      flPropFn    = 'fluidprop.dat'
      flName      = 'air01'
      resReadFl   = 'restart.rst'
      resWriteFl  = 'restart.rst'
      outPrefix   = 'output'
      ssPrefix    = 'snapshot'
      TSPrefix    = 'timesrs'
      tAvgOutPref = 'timeavg'

C---: Initialize logical variables
      lCartesGrid = .false.
      lTmAvgFlg  = .false.


C---: Open boundary-conditions file
      open(10,file=inputFile,status='old')


C---: Initialize local flags
      lInSection    = .false.
      lFoundSection = .false.


C---: Initialize line counter
      nl = 0

C---: Start reading lines
      do
C-----: Increment line counter
        nl = nl + 1

C-----: Read line
        read(10,'(a)',iostat=nstat) cBuf
        if(nstat.lt.0) exit   !---> Reached end of file

C-----> Skip blank lines
        if(mTrimLeft(cBuf).lt.1) cycle

C-----> Comment lines:  Don't parse line if it starts with a '#'
        if(cBuf(1:1).eq.'#') cycle

C-----> Initialize word counter
        nw = 0

C-----> Start parsing words
        do
          if(mTrimLeft(cBuf).lt.1) exit
          nc = mExtractW(cBuf,cw)

C-------> Convert word to lower case
          call mToLowerC(cw)

C-------> End of line (empty word)
          if(nc.lt.1) exit

C-------> Increment word counter
          nw = nw + 1

C-------> Find token in hash table
          nTok = nFindTokIP(nc, cw)


C-------> If section flags have not been defined, don't parse 
C-------> the rest of the current line
          if(.not.lInSection.and.nTok.gt.TEND) then
            print '(a)', '* Section name not defined'
            call SyntaxError(nl,nw)
            lReadIPFile = .false.
            return
          end if

C-------> If not in the right section, don't parse the rest of the line
          if(lInSection.and.(.not.lFoundSection).and.
     .       nTok.lt.TSECTION) exit


C-------> Perform operation or assignment corresponding to token
          select case(nTok)

          case (TSECTION)
C---------: Check that keyword is not out of place
            if(lInSection) then
              print '(a)', '* SECTION keyword out of pace'
              call SyntaxError(nl,nw)
            end if
C---------: Activate flag
            lInSection = .true.
C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
C---------> Convert section name to lower case
            call mToLowerC(cw)
C---------: Compare current section name with expected
            if(cw(1:nc).eq.'input_parameters') then
              lFoundSection = .true.
            end if
            exit   !---> Don't fetch more words from this line


          case (TEND)
C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
C---------> Convert word to lower case
            call mToLowerC(cw)
            if(cw(1:nc).eq.'section') then
              if(lInSection) then
                lInSection    = .false.
                lFoundSection = .false.
              else
                print '(a)', '* END SECTION statement out of place'
                call SyntaxError(nl,nw)
              end if
C------------> If end of required section, return
               if(lFoundSection) return
           else
              print '(a)', '* SECTION keyword missing'
              call SyntaxError(nl,nw)
            end if
            exit   !---> Don't fetch more words from this line


          case (TGRIDFILE)
C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter

C---------: Assign value to corresponding variable
            gridFile = cw(1:nc)

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) exit
            nw = nw + 1   !---> Increment word counter

C---------> Convert word to lower case
            call mToLowerC(cw)

C-----------: `Cartesian grid' flag
            if(cw(1:nc).eq.'cartesian_grid') then
              lCartesGrid = .true.
            else
              print*, 'Unrecognized keyword: [',cw(1:nc),']'
              call SyntaxError(nl,nw)
            end if

            exit   !---> Don't fetch more words from this line


          case (TTHERMENERGY)
            nthermen = 1
            exit   !---> Don't fetch more words from this line

          case (TBUOYANCY)
            neqstate = 1
            exit   !---> Don't fetch more words from this line

          case (TSMALLSCALE)
            nsmallscl = 1
            exit   !---> Don't fetch more words from this line

          case (TTRAJECTORIES)
            ntraject = 1
            exit   !---> Don't fetch more words from this line

          case (TNTIMESTEPS)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            nts = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TTIMESTEPSIZE)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw) 
            dk = dnum(cw(1:nc))

C---> Note: The quantity expected here is the 'Dimensional' 
C---------: time-step size.  It is non-dimensionalized (i.e.,
C---------: scaled wrt. uref and lref) after returing to the
C---------: calling procedure (subroutine ReadParam)

            exit   !---> Don't fetch more words from this line

          case (TMAXQLITER)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            mqiter = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TQLTOLERANCE)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            qtol = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TMAXMEITER)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            nmeiter = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TMETOLERANCE)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            dmeittol = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TFILTERU)
C---------: Activate corresponding flag
            nfiltu = 1

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            fpu = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TFILTERV)
C---------: Activate corresponding flag
            nfiltv = 1

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            fpv = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TFILTERT)
C---------: Activate corresponding flag
            nfiltt = 1

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            fpt = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TPPESOLVER)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

C---------> Convert word to lower case
            call mToLowerC(cw)

C---------: Test parsed cw against possible options
            if(cw(1:nc).eq.'sor')              then
              nPpeSolver = 1
            else if(cw(1:nc).eq.'lsor')        then
              nPpeSolver = 2
            else if(cw(1:nc).eq.'rb_lsor')     then
              nPpeSolver = 3
            else if(cw(1:nc).eq.'par_rb_lsor') then
              nPpeSolver = 4
            else if(cw(1:nc).eq.'rb_sor')      then
              nPpeSolver = 5
            else if(cw(1:nc).eq.'par_rb_sor')  then
              nPpeSolver = 6
            else
              print*, 'Unrecognized keyword: [',cw(1:nc),']'
              call SyntaxError(nl,nw)
            end if
            exit   !---> Don't fetch more words from this line

          case (TMAXSORITER)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            msorit = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TSORTOLERANCE)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            sortol = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TSORRELAX)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            sorrel = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDCU)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssCu0 = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDTSC)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssTsCoef = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDHSC)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssHsCoef = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDBNCR)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssBnCrit = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDRMAX)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssRMpMax = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDRMEXP)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssRMpExp = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDTEM)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssTemCoef = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDPPESLVR)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

C---------> Convert word to lower case
            call mToLowerC(cw)

C---------: Test parsed cw against possible options
            if(cw(1:nc).eq.'sor')              then
              nssPpeSlvr = 1
            else if(cw(1:nc).eq.'lsor')        then
              nssPpeSlvr = 2
            else if(cw(1:nc).eq.'rb_lsor')     then
              nssPpeSlvr = 3
            else if(cw(1:nc).eq.'par_rb_lsor') then
              nssPpeSlvr = 4
            else if(cw(1:nc).eq.'rb_sor')      then
              nssPpeSlvr = 5
            else if(cw(1:nc).eq.'par_rb_sor')  then
              nssPpeSlvr = 6
            else
              print*, 'Unrecognized keyword: [',cw(1:nc),']'
              call SyntaxError(nl,nw)
            end if
            exit   !---> Don't fetch more words from this line

          case (TATDMAXSORIT)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            mssSorIt = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDSORTOL)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssSorTol = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDSORRELAX)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssSorRel = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDFILTPARU)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssFiltPar(_U_) = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDFILTPARV)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssFiltPar(_V_) = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TATDFILTPART)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            ssFiltPar(_T_) = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TREFLENGTH)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            dlref = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TREFVELOCITY)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            uref = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TREFTEMPER)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            tref = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TMAXREFTEMPER)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            tmax = dnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TFLUIDPROPFILE)
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter
            flPropFn = cw(1:nc)
            exit   !---> Don't fetch more words from this line

          case (TFLUIDPROPTAB)
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter
            flName = cw(1:nc)
            exit   !---> Don't fetch more words from this line

          case (TRESTART)
C---------: Activate corresponding flag
            nrestart = 1

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter
            resReadFl = cw(1:nc)
            exit   !---> Don't fetch more words from this line

          case (TWRITERESTART)
C---------: Activate corresponding flag
            nreswrite = 1

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter
            resWriteFl = cw(1:nc)

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) then
              nreswrite = 1
            else
              if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
              nreswrite = 2
              nWrtRstFq = jnum(cw(1:nc))
            end if

            exit   !---> Don't fetch more words from this line

          case (TTIMEAVGMODE)
C---------: Activate corresponding flag
            lTmAvgFlg  = .true.

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter
            tAvgOutPref = cw(1:nc)
            exit   !---> Don't fetch more words from this line

          case (TOUTPUTFORMAT)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

C---------> Convert word to lower case
            call mToLowerC(cw)

C---------: Test parsed cw against possible options
            if(cw(1:nc).eq.'plot3d')  then
              nPostProc = 1
            else if(cw(1:nc).eq.'fv') then
              nPostProc = 2
            else
              print*, 'Unrecognized keyword: [',cw(1:nc),']'
              call SyntaxError(nl,nw)
            end if

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) exit
            nw = nw + 1   !---> Increment word counter

C---------> Convert word to lower case
            call mToLowerC(cw)

C---------: Test parsed cw against possible options
            if(cw(1:nc).eq.'formatted')  then
              nOutFlForm = FT_FORMATTED
            else if(cw(1:nc).eq.'unformatted') then
              nOutFlForm = FT_UNFORMATTED
            else
              print*, 'Unrecognized keyword: [',cw(1:nc),']'
              call SyntaxError(nl,nw)
            end if

            exit   !---> Don't fetch more words from this line

          case (TOUTPUTPREFIX)
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter
            outPrefix = cw(1:nc)
            exit   !---> Don't fetch more words from this line

          case (TSNAPSHOTSPREF)
C---------: Activate corresponding flag
            nsnapshots = 1

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter
            ssPrefix = cw(1:nc)
            exit   !---> Don't fetch more words from this line

          case (TSNAPSHOTSFREQ)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
              isfreq = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TSNAPSHOTSFRST)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            isfirst = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TTIMESRSPREFIX)
C---------: Activate corresponding flag
            nTSOutFlag = 1

C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            if(nc.lt.1) call SyntaxError(nl,nw)
            nw = nw + 1   !---> Increment word counter
            TSPrefix = cw(1:nc)
            exit   !---> Don't fetch more words from this line

          case (TTIMESRSFREQ)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            nTSFreq = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TTIMESRSPOINT)
C---------: Check that we're not exceeding maxtspts
            if((nTSPoints+1).gt.maxtspts) then
            print*, 'Exceeded allowed number of time-series ',
     .           'monitor points: ', maxtspts
            print*, 'Ignoring last TSPoint definition.'

C-----------: Don't read the TSPoint indeces and read next line
              exit
            end if

C---------: Increment nTSPoints counter
            nTSPoints = nTSPoints + 1

C---------: Fetch first index
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            iTS(nTSPoints) = jnum(cw(1:nc))

C---------: Fetch second index
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jTS(nTSPoints) = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TPRINTDIFFFREQ)
C---------: Activate corresponding flag
            lPrntDiff = 1

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            nPrDfFreq = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

          case (TPRINTBANNFREQ)
C---------: Activate corresponding flag
            lPrDfBann = 1

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            nPrDfBnFq = jnum(cw(1:nc))
            exit   !---> Don't fetch more words from this line

C-------: If word was not recognized, print an error message and stop
          case default

            write(STDOUT, '(a/,3a)')
     .      '* Error in Input_Parameters section:',
     .      '* Unrecognized keyword: [',cw(1:nc),']'

            write(STDLOG,'(a/,3a)')
     .      '* Error in Input_Parameters section:',
     .      '* Unrecognized keyword: [',cw(1:nc),']'

            call SyntaxError(nl,nw)

          end select

        end do  !---> Stop parsing words

      end do   !---> Stop reading lines


C---: Close input file
      close(10)

      end



*----------------------------------------------------------------------*
*                               nFindTokIP                             *
*----------------------------------------------------------------------*
*
*
      INTEGER function nFindTokIP(nTestLen,str)

      implicit none
      character*(*) str

      INTEGER nKeyWords
      parameter ( nKeyWords =  53 )  !-> Note: The number of keywords
                                     !-> must be updated everytime a new
                                     !-> keyword definition is added.
      INTEGER nTestLen

C---: Note: Must change size of array
      character*20 cToken(nKeyWords)

      INTEGER nStrLen

      INTEGER nLen, i

      data  cToken(TGRIDFILE)      / 'grid_file'           /
     .      cToken(TTHERMENERGY)   / 'thermal_energy'      /
     .      cToken(TBUOYANCY)      / 'buoyancy'            /
     .      cToken(TSMALLSCALE)    / 'small_scale'         /
     .      cToken(TTRAJECTORIES)  / 'trajectories'        /
C----:
     .      cToken(TNTIMESTEPS)    / 'n_time_steps'        /
     .      cToken(TTIMESTEPSIZE)  / 'time_step_size'      /
C----:
     .      cToken(TMAXQLITER)     / 'max_ql_iter'         /
     .      cToken(TQLTOLERANCE)   / 'ql_tolerance'        /
C----:
     .      cToken(TMAXMEITER)     / 'max_me_iter'         /
     .      cToken(TMETOLERANCE)   / 'me_tolerance'        /
C----:
     .      cToken(TFILTERU)       / 'filter_u'            /
     .      cToken(TFILTERV)       / 'filter_v'            /
     .      cToken(TFILTERT)       / 'filter_t'            /
C----:
     .      cToken(TPPESOLVER)     / 'ppe_solver'          /
     .      cToken(TMAXSORITER)    / 'max_sor_iter'        /
     .      cToken(TSORTOLERANCE)  / 'sor_tolerance'       /
     .      cToken(TSORRELAX)      / 'sor_relaxation'      /
C----:
     .      cToken(TATDCU)         / 'atd_cu0'       /
     .      cToken(TATDTSC)        / 'atd_tscf'      /
     .      cToken(TATDHSC)        / 'atd_hscf'      /
     .      cToken(TATDBNCR)       / 'atd_bnumcrit'  /
     .      cToken(TATDRMAX)       / 'atd_rmapmax'   /
     .      cToken(TATDRMEXP)      / 'atd_rmapexp'   /
     .      cToken(TATDTEM)        / 'atd_temcf'     /
C----:
     .      cToken(TATDPPESLVR)   / 'atd_ppe_solver'     /
     .      cToken(TATDMAXSORIT)  / 'atd_max_sor_iter'   /
     .      cToken(TATDSORTOL)    / 'atd_sor_tolerance'  /
     .      cToken(TATDSORRELAX)  / 'atd_sor_relaxation' /
     .      cToken(TATDFILTPARU)  / 'atd_filter_u'       /
     .      cToken(TATDFILTPARV)  / 'atd_filter_v'       /
     .      cToken(TATDFILTPART)  / 'atd_filter_t'       /
C----:
     .      cToken(TREFLENGTH)     / 'ref_length'          /
     .      cToken(TREFVELOCITY)   / 'ref_velocity'        /
     .      cToken(TREFTEMPER)     / 'ref_temperature'     /
     .      cToken(TMAXREFTEMPER)  / 'max_ref_temperature' /
C----:
     .      cToken(TFLUIDPROPFILE) / 'fluid_prop_file'     /
     .      cToken(TFLUIDPROPTAB)  / 'fluid_prop_table'    /
C----:
     .      cToken(TRESTART)       / 'restart'             /
     .      cToken(TWRITERESTART)  / 'write_restart'       /
C----:
     .      cToken(TTIMEAVGMODE)   / 'time_average_mode'   /
C----:
     .      cToken(TOUTPUTFORMAT)  / 'output_format'       /
     .      cToken(TOUTPUTPREFIX)  / 'output_prefix'       /
C----:
     .      cToken(TSNAPSHOTSPREF) / 'snapshots_prefix'    /
     .      cToken(TSNAPSHOTSFREQ) / 'snapshots_freq'      /
     .      cToken(TSNAPSHOTSFRST) / 'snapshots_first'     /
C----:
     .      cToken(TTIMESRSPREFIX) / 'time_series_prefix'  /
     .      cToken(TTIMESRSFREQ)   / 'time_series_freq'    /
     .      cToken(TTIMESRSPOINT)  / 'time_series_point'   /
C----:
     .      cToken(TPRINTDIFFFREQ) / 'print_diff_freq'     /
     .      cToken(TPRINTBANNFREQ) / 'print_bann_freq'     /
C----:
     .      cToken(TSECTION)       / 'section'             /
     .      cToken(TEND)           / 'end'                 /


C---: Initial state
      nFindTokIP = 0

c      nTestLen =  nStrLen(str)

C---: Check all tokens (Note the top index has to be equal to the
C---: total number of tokens)
      do i=1, nKeyWords

        nLen = nStrLen(cToken(i))

        if(cToken(i)(1:nLen).eq.str(1:nTestLen)) then

          nFindTokIP = i
          exit

         end if

      end do

      return
      end


*----------------------------------------------------------------------*
*                              lReadBCFile                             *
*----------------------------------------------------------------------*

      LOGICAL function
     .   lReadBCFile(nx, ny, cBCFile,
     .               nthermen, neqstate, nsmallscl, ntraject,
     .               dlref, uref, tref, tmax, densref,
     .               nReg, nRegBrd,
     .               nRegType, nMomBdTp, nTRgType, nTemBdTp,
     .               nPRPermEq, dPRDimCoef,
     .               dTRgVal, dBCVal)

      implicit none

      include "config.f"

      INTEGER nx, ny

      character*(*) cBCFile

      INTEGER nthermen, neqstate, nsmallscl, ntraject
      REAL    dlref, uref, tref, tmax, densref

      character*(mlinelgt) cBuf, cw
      LOGICAL lCheckFile, lCheckIJreg, lStIsInt, lStIsReal
      INTEGER mTrimLeft, mExtractW, mCellLocStr, nFindTokBC

#ifndef _HP_

      INTEGER jnum
      REAL    dnum

#endif _HP_


C---: Number  of regions and region borders
      INTEGER nReg(2)
      INTEGER nRegBrd(mgri,mgrj,4)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Momentum boundary types
      INTEGER nMomBdTp(mgri,mgrj,4)

C---: Declare Temperature region values array
      REAL     dTRgVal(mgri,mgrj)

C---: Temperature boundary types
      INTEGER nTemBdTp(mgri,mgrj,4)

C---: Declare porous region coefficient arrays (temporary)
      INTEGER nPRPermEq(mgri,mgrj)
      REAL    dPRDimCoef(mgri,mgrj,3)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,4)

C---: Local variables
      LOGICAL lFndNRegs, lInSection, lFoundSection
      LOGICAL lDefIBorders, lDefJBorders
      INTEGER i, j, ireg, jreg, nl, nw, nstat, nc, nerror, nbcloc
      INTEGER ntok
      INTEGER ndefw
      REAL    dumT, dumV, dumP

C---: Default status
      lReadBCFile = .true.

C---: Check that boundary file is readable
      if(.not.lCheckFile(cBCFile,nError)) then
        print*, 'Error code :', nError
        return
      end if

C---: Open boundary-conditions file
      open(10,file=cBCFile,status='old')

C---: Initialize internal flags
      lFndNRegs     = .false.
      lDefIBorders  = .false.
      lDefJBorders  = .false.
      lInSection    = .false.
      lFoundSection = .false.


C---: Initialize line counter
      nl = 0

C---: Start reading lines
      do
C-----: Increment line counter
        nl = nl + 1

C-----: Read line
        read(10,'(a)',iostat=nstat) cBuf
        if(nstat.lt.0) exit   !---> Reached end of file

C-----> Skip blank lines
        if(mTrimLeft(cBuf).lt.1) cycle

C-----> Comment lines:  Don't parse line if it starts with a '#'
        if(cBuf(1:1).eq.'#') cycle

C-----> Convert entire line to lower case
        call mToLowerC(cBuf)

C-----: Initialize word counter
        nw = 0

C-----: Start parsing words
        do
          if(mTrimLeft(cBuf).lt.1) exit
          nc = mExtractW(cBuf,cw)

C-------> End of line (empty word)
          if(nc.lt.1) exit

C-------> Increment word counter
          nw = nw + 1   !---> Increment word counter

C-------: Find token in hash table
          nTok = nFindTokBC(nc, cw)


C-------> If section flags have not been defined, don't parse 
C-------> the rest of the current line
          if(.not.lInSection.and.nTok.gt.TEND) then
            print '(a)', '* Section name not defined'
            call SyntaxError(nl,nw)
            lReadBCFile = .false.
            return
          end if

C-------> If not in the right section, don't parse the rest of the line
          if(lInSection.and.(.not.lFoundSection).and.
     .       nTok.lt.TSECTION) exit


C-------: Check that no. of regs and interior borders have been defined
          if(.not.lDefIBorders.or.(.not.lDefJBorders)) then
            select case(ntok)

            case(TSECTION, TEND, TNUMOFREGIONS, TIBORDERS, TJBORDERS)

            case default
              write(STDOUT, '(2(a/),a)')
     .        '* Syntax error in boundary definition:', 
     .        '* Number of regions or interior borders not defined',
     .        '* prior to region and boundary declarations'

              write(STDLOG, '(2(a/),a)')
     .        '* Syntax error in boundary definition:', 
     .        '* Number of regions or interior borders not defined',
     .        '* prior to region and boundary declarations'

              lReadBCFile = .false.
              return
            end select
          end if

C-------> Perform operation or assignment corresponding to token
          select case(nTok)

          case (TSECTION)
C---------: Check that keyword is not out of place
            if(lInSection) then
              print '(a)', '* SECTION keyword out of pace'
              call SyntaxError(nl,nw)
            end if
C---------: Activate flag
            lInSection = .true.
C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
C---------> Convert section name to lower case
            call mToLowerC(cw)
C---------: Compare current section name with expected
            if(cw(1:nc).eq.'boundary_conditions') then
              lFoundSection = .true.
            end if
            exit   !---> Don't fetch more words from this line


          case (TEND)
C---------: Fetch next word
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
C---------> Convert word to lower case
            call mToLowerC(cw)
            if(cw(1:nc).eq.'section') then
              if(lInSection) then
                lInSection    = .false.
                lFoundSection = .false.
              else
                print '(a)', '* END SECTION statement out of place'
                call SyntaxError(nl,nw)
              end if
C------------> If end of required section, return
               if(lFoundSection) return
           else
              print '(a)', '* SECTION keyword missing'
              call SyntaxError(nl,nw)
            end if
            exit   !---> Don't fetch more words from this line


          case (TNUMOFREGIONS)
C---------: Read number of I-regions 
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            nReg(_I_) = jnum(cw(1:nc))

C---------: Read number of J-regions 
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            nReg(_J_) = jnum(cw(1:nc))

C---------: Check that nRegI and nRegJ are consistent
            if(.not.lCheckIJreg(nReg(_I_), nReg(_J_))) then
              lReadBCFile = .false.
              return
            end if

C---------: Turn 'found-number-of-regions' flag on
            lFndNRegs = .true.

C---------: Don't read more words from this line
            exit

          case (TIBORDERS)
C---------: If nReg's havent' been read yet
            if(.not.lFndNRegs) then
              print*, 'Error: Number_of_Regions not defined'
              call SyntaxError(nl,nw)
            end if

C---------: Read interior borders
            do i=2, nReg(_I_)
              nc = mExtractW(cBuf,cw)
              nw = nw + 1   !---> Increment word counter
              if(nc.lt.1) call SyntaxError(nl,nw)
              if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
              nRegBrd(i,1,WEST) = jnum(cw(1:nc))
            end do

C---------: Turn lDefIBorders flag on
            lDefIBorders = .true.

C---------: If I and J interior borders have been defined,
C---------: then initialize region and boundary flags
            if(lDefJBorders) 
     .        call InitBCFlags(nReg, nRegType, nMomBdTp,
     .                         nTRgType, nTemBdTp, dBCVal)

          case (TJBORDERS)
C---------: If nReg's havent' been read yet
            if(.not.lFndNRegs) then
              print*, 'Error: Number_of_Regions not defined'
              call SyntaxError(nl,nw)
            end if

C---------: Read interior borders
            do j = 2, nReg(_J_)
              nc = mExtractW(cBuf,cw)
              nw = nw + 1   !---> Increment word counter
              if(nc.lt.1) call SyntaxError(nl,nw)
              if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
              nRegBrd(1,j,SOUTH) = jnum(cw(1:nc))
            end do

C---------: Turn lDefJBorders flag on
            lDefJBorders = .true.

C---------: If I and J interior borders have been defined,
C---------: then initialize region and boundary flags
            if(lDefIBorders)
     .        call InitBCFlags(nReg, nRegType, nMomBdTp,
     .                         nTRgType, nTemBdTp, dBCVal)


C-----------------------------------------------------------------------
C------------: Keywords for special region types - Momentum :-----------
C-----------------------------------------------------------------------

C---------------------------[ Blockage ]

          case (TBLOCKAGE)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Read Region type
            nRegType(ireg,jreg) = RM_BLOCKG

C---------: For a blockage reg , set walls to be no-slip for momentum
C---------: and adiabatic for temperature
            do nbcloc= WEST, NORTH
              nMomBdTp(ireg,jreg,nbcloc) = BM_WALL1
              nTemBdTp(ireg,jreg,nbcloc) = BT_HTFLUX
            end do


C---------: Coincide corresponding walls on neighboring regions
            if(ireg.gt.1) then
              nMomBdTp(ireg-1,jreg,EAST) = BM_WALL1
              nTemBdTp(ireg-1,jreg,EAST) = BT_HTFLUX
            end if 
            if(ireg.lt.nReg(_I_)) then
              nMomBdTp(ireg+1,jreg,WEST) = BM_WALL1
              nTemBdTp(ireg+1,jreg,WEST) = BT_HTFLUX
            end if 
            if(jreg.gt.1) then
              nMomBdTp(ireg,jreg-1,NORTH) = BM_WALL1
              nTemBdTp(ireg,jreg-1,NORTH) = BT_HTFLUX
            end if
            if(jreg.lt.nReg(_J_)) then
              nMomBdTp(ireg,jreg+1,SOUTH) = BM_WALL1
              nTemBdTp(ireg,jreg+1,SOUTH) = BT_HTFLUX
            end if


C---------: Parse remaining of words in line
            do
              nc = mExtractW(cBuf,cw)
              nw = nw + 1   !---> Increment word counter

C-----------: Stop reading line when cBuf is empty
              if(nc.lt.1) exit

C-----------: Temperature
              if(cw(1:nc).eq.'temper'.or.cw(1:nc).eq.'t') then

C-------------: Assign region type
                nTRgType(ireg,jreg) = RT_TEMPER

C-------------: Read temperature value
                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. to reference values
                dumT = (dnum(cw(1:nc))-tref)/(tmax-tref)
                dTRgVal(ireg,jreg) = dumT

C-------------: Set fixed temperature on walls of fixed temp. region
                do nbcloc= WEST, NORTH
                  nTemBdTp(ireg,jreg,nbcloc)   = BT_TEMPER
                  dBCVal(ireg,jreg,nbcloc,_T_) = dumT
                end do

C-------------: Coincide corresponding walls on neighboring regions
                if(ireg.gt.1) then
                  nTemBdTp(ireg-1,jreg,EAST)   = BT_TEMPER
                  dBCVal(ireg-1,jreg,EAST,_T_) = dumT
                end if 
                if(ireg.lt.nReg(_I_)) then
                  nTemBdTp(ireg+1,jreg,WEST)   = BT_TEMPER
                  dBCVal(ireg+1,jreg,WEST,_T_) = dumT
                end if 
                if(jreg.gt.1) then
                  nTemBdTp(ireg,jreg-1,NORTH)   = BT_TEMPER
                  dBCVal(ireg,jreg-1,NORTH,_T_) = dumT
                end if
                if(jreg.lt.nReg(_J_)) then
                  nTemBdTp(ireg,jreg+1,SOUTH)   = BT_TEMPER
                  dBCVal(ireg,jreg+1,SOUTH,_T_) = dumT
                end if

              else
                 call SyntaxError(nl,nw)
              end if

            end do

C---------------------------[ Porous_Region]

          case (TPOROUSREGION)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Read Region type
            nRegType(ireg,jreg) = RM_POROUS


C---------: Also read the porous medium coefficients:
C---------: Which equation to use in calculation of permeability
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

C---------: Test parsed cw against possible options
            if(cw(1:nc).eq.'carman-kozeny')         then
              nPRPermEq(ireg,jreg) = 1
            else if(cw(1:nc).eq.'kaviany')          then
              nPRPermEq(ireg,jreg) = 2
            else if(cw(1:nc).eq.'happen-brenner')   then
              nPRPermEq(ireg,jreg) = 3
            else
              print*, 'Unrecognized keyword: [',cw(1:nc),']'
              call SyntaxError(nl,nw)
            end if

C---------: Porosity
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            dPRDimCoef(ireg,jreg,1) = dnum(cw(1:nc))

C---------: Diameter in permeability equation
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            dPRDimCoef(ireg,jreg,2) = dnum(cw(1:nc))

C---------: Nonlinear drag coefficient (Ergun coefficient)
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            dPRDimCoef(ireg,jreg,3) = dnum(cw(1:nc))

            exit   !---> Don't fetch more words from this line



C-----------------------------------------------------------------------
C-----------: Keywords for special boundary types - Momentum -----------
C-----------------------------------------------------------------------

C---------------------------[ Wall ]

          case (TWALL)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Read cell face location
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

            nbcloc = mCellLocStr(nc, cw)
            if (nbcloc.lt.1) call SyntaxError(nl,nw)

C---------: Assume default wall type (no_slip)
            nMomBdTp(ireg,jreg,nbcloc) = BM_WALL1

C---------: Parse remaining of words in line
            do
              nc = mExtractW(cBuf,cw)
              nw = nw + 1   !---> Increment word counter

C-----------: Stop reading line when cBuf is empty
              if(nc.lt.1) exit

C-----------: No-slip wall type
              if(cw(1:nc).eq.'no_slip') then
                nMomBdTp(ireg,jreg,nbcloc) = BM_WALL1

C-----------: No-stress wall type
              else if(cw(1:nc).eq.'no_stress') then
                nMomBdTp(ireg,jreg,nbcloc) = BM_WALL2

C-----------: Tangential flow speed
              else if(cw(1:nc).eq.'tangent_vel'.or.
     .          cw(1:nc).eq.'tv') then
                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. reference velocity
                dumV = dnum(cw(1:nc)) / uref

C-------------: Assign to velocity component depending on orientation
                select case(nbcloc)
                case (WEST,EAST)
                  dBCVal(ireg,jreg,nbcloc,_V_) = dumV
                case (SOUTH,NORTH)
                  dBCVal(ireg,jreg,nbcloc,_U_) = dumV
                end select

C-----------: Pressure
              else if(cw(1:nc).eq.'press'.or.cw(1:nc).eq.'p') then
                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. reference values
                dumP = dnum(cw(1:nc))/(densref*uref*uref)
                dBCVal(ireg,jreg,nbcloc,_P_) = dumP

C-----------: Temperature
              else if(cw(1:nc).eq.'temper'.or.cw(1:nc).eq.'t') then

C-------------: Assign boundary type flag
                nTemBdTp(ireg,jreg,nbcloc) = BT_TEMPER

                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. to reference values
                dumT = (dnum(cw(1:nc))-tref)/(tmax-tref)
                dBCVal(ireg,jreg,nbcloc,_T_) = dumT

              else
                 call SyntaxError(nl,nw)
              end if

            end do


C---------------------------[ Inlet ]

          case (TINLET)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Read cell face location
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

            nbcloc = mCellLocStr(nc, cw)
            if (nbcloc.lt.1) call SyntaxError(nl,nw)

C---------: Assign boundary type flag
            nMomBdTp(ireg,jreg,nbcloc) = BM_INLET

C---------: Parse remaining of words in line
            ndefw = 0
            do
              nc = mExtractW(cBuf,cw)
              nw = nw + 1   !---> Increment word counter

C-----------: Stop reading line when cBuf is empty
              if(nc.lt.1) then
                if(ndefw.lt.1) call SyntaxError(nl,nw)
                exit
              end if

              ndefw = ndefw + 1

C-----------: Normal flow speed
              if(cw(1:nc).eq.'normal_vel'.or.cw(1:nc).eq.'nv') then
                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. reference velocity
                dumV = dnum(cw(1:nc)) / uref

C-------------: Assign to velocity component depending on orientation
                select case(nbcloc)
                case (WEST,EAST)
                  dBCVal(ireg,jreg,nbcloc,_U_) = dumV
                case (SOUTH,NORTH)
                  dBCVal(ireg,jreg,nbcloc,_V_) = dumV
                end select

C-----------: Tangential flow speed
              else if(cw(1:nc).eq.'tangent_vel'.or.cw(1:nc).eq.'tv') 
     .          then
                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. reference velocity
                dumV = dnum(cw(1:nc)) / uref

C-------------: Assign to velocity component depending on orientation
                select case(nbcloc)
                case (WEST,EAST)
                  dBCVal(ireg,jreg,nbcloc,_V_) = dumV
                case (SOUTH,NORTH)
                  dBCVal(ireg,jreg,nbcloc,_U_) = dumV
                end select

C-----------: Temperature
              else if(cw(1:nc).eq.'temper'.or.cw(1:nc).eq.'t') then

C-------------: Assign boundary type flag
                nTemBdTp(ireg,jreg,nbcloc) = BT_TEMPER

                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. to reference values
                dumT = (dnum(cw(1:nc))-tref)/(tmax-tref)
                dBCVal(ireg,jreg,nbcloc,_T_) = dumT

              else
                 call SyntaxError(nl,nw)
              end if

            end do


C---------------------------[ Outlet ]

          case (TOUTLET)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Read cell face location
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

            nbcloc = mCellLocStr(nc, cw)
            if (nbcloc.lt.1) call SyntaxError(nl,nw)

C---------: Assume outlet type (Mass-conservation)
            nMomBdTp(ireg,jreg,nbcloc) = BM_OUTLT2

C---------: Parse remaining words in line
            do
              nc = mExtractW(cBuf,cw)
              nw = nw + 1   !---> Increment word counter

C-----------: Stop reading line when cBuf is empty
              if(nc.lt.1) exit

C-----------: Fully developed outlet type
              if(cw(1:nc).eq.'fully_dev'.or.cw(1:nc).eq.'fd') then
                nMomBdTp(ireg,jreg,nbcloc) = BM_OUTLT1

C-----------: Mass conservation outlet type
              else if(cw(1:nc).eq.'mass_cons'.or.
     .          cw(1:nc).eq.'mc') then
                nMomBdTp(ireg,jreg,nbcloc) = BM_OUTLT2

C-----------: Pressure
              else if(cw(1:nc).eq.'press'.or.cw(1:nc).eq.'p') then
                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. reference values
                dumP = dnum(cw(1:nc))/(densref*uref*uref)
                dBCVal(ireg,jreg,nbcloc,_P_) = dumP

C-----------: Temperature
              else if(cw(1:nc).eq.'temper'.or.cw(1:nc).eq.'t') then

C-------------: Assign boundary type flag
                nTemBdTp(ireg,jreg,nbcloc) = BT_TEMPER
                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(nc.lt.1) call SyntaxError(nl,nw)
                if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C-------------: Scale wrt. to reference values
                dumT = (dnum(cw(1:nc))-tref)/(tmax-tref)
                dBCVal(ireg,jreg,nbcloc,_T_) = dumT

              else
                call SyntaxError(nl,nw)
              end if

            end do

C-----------------------------------------------------------------------
C------------: Keywords for special region types - Thermal :------------
C-----------------------------------------------------------------------

C---------------------------[ Heat_Generation_Region ]

          case (THEATGENREGION)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Assign region type
            nTRgType(ireg,jreg) = RT_HEATGN

C---------: Read volumetric heat generation value
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

            dTRgVal(ireg,jreg) = dnum(cw(1:nc))

C---------: Note that the heat generation coefficient is 
C---------: non-dimensionalized in a separate procedure.


C---------------------------[ Fixed_Temperature_Region ]

          case (TTEMPERREGION)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Assign region type
            nTRgType(ireg,jreg) = RT_TEMPER

C---------: Read temperature value
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C---------: Scale wrt. to reference values
            dumT = (dnum(cw(1:nc))-tref)/(tmax-tref)
            dTRgVal(ireg,jreg) = dumT

C---------: Set fixed temperature on walls of fixed temp. region
            do nbcloc= WEST, NORTH
              nTemBdTp(ireg,jreg,nbcloc)   = BT_TEMPER
              dBCVal(ireg,jreg,nbcloc,_T_) = dumT
            end do


C-----------------------------------------------------------------------
C-----------: Keywords for special boundary types - Thermal ------------
C-----------------------------------------------------------------------

C---------------------------[ Fixed_Temperature_Wall ]

          case (TFIXEDTEMPWALL)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Read cell face location
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

            nbcloc = mCellLocStr(nc, cw)
            if (nbcloc.lt.1) call SyntaxError(nl,nw)

C---------: Assign boundary type flag
            nTemBdTp(ireg,jreg,nbcloc) = BT_TEMPER

C---------: Read temperature value
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)

C---------: Scale wrt. to reference values
            dumT = (dnum(cw(1:nc))-tref)/(tmax-tref)
            dBCVal(ireg,jreg,nbcloc,_T_) = dumT


C---------------------------[ Heat_Flux_Wall ]

          case (THEATFLUXWALL)
C---------: Read ireg and jreg indeces
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            ireg = jnum(cw(1:nc))

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)
            if(.not.lStIsInt(cw)) call SyntaxError(nl,nw)
            jreg = jnum(cw(1:nc))

C---------: Check that ireg and jreg are consistent
            if(.not.lCheckIJreg(ireg, jreg)) then
              lReadBCFile = .false.
              return
            end if

C---------: Read cell face location
            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

            nbcloc = mCellLocStr(nc, cw)
            if (nbcloc.lt.1) call SyntaxError(nl,nw)

C---------: Assign boundary type flag
            nTemBdTp(ireg,jreg,nbcloc) = BT_HTFLUX


C***** Note: For now only consider adiabatic walls

C-------: If word was not recognized, print an error message and stop
          case default
            write(STDOUT,'(a/,2a)')
     .      '* Error in Boundary_Conditions section:',
     .      '* Unrecognized keyword: [',cw(1:nc),']'

            write(STDLOG,'(a/,2a)')
     .      '* Error in Boundary_Conditions section:',
     .      '* Unrecognized keyword: [',cw(1:nc),']'

            call SyntaxError(nl,nw)

          end select

        end do  !---> Stop parsing words

      end do   !---> Stop reading lines

C---: Close input file
      close(10)

      return
      end


*----------------------------------------------------------------------*
*                               nFindTokBC                             *
*----------------------------------------------------------------------*
*
*
      INTEGER function nFindTokBC(nTestLen,str)

      implicit none

      character*(*) str

      INTEGER nKeyWords
      parameter ( nKeyWords = 14 ) !---> Note: The number of key words
                                   !---> must be updated everytime a new
                                   !---> keyword definition is added.
      INTEGER nTestLen

C---: Note: Must change size of array
      character*26 cToken(nKeyWords)

      INTEGER nStrLen

      INTEGER nLen, i

      data  cToken(TNUMOFREGIONS)   / 'number_of_regions'      /
     .      cToken(TIBORDERS)       / 'i_borders'              /
     .      cToken(TJBORDERS)       / 'j_borders'              /
C----:
     .      cToken(TBLOCKAGE)       / 'blockage'               /
     .      cToken(TPOROUSREGION)   / 'porous_region'          /
C----:
     .      cToken(TWALL)           / 'wall'                   /
     .      cToken(TINLET)          / 'inlet'                  /
     .      cToken(TOUTLET)         / 'outlet'                 /
C----:
     .      cToken(THEATGENREGION)  / 'heat_generation_region'   /
     .      cToken(TTEMPERREGION)   / 'fixed_temperature_region' /
C----:
     .      cToken(TFIXEDTEMPWALL)  / 'fixed_temperature_wall'   /
     .      cToken(THEATFLUXWALL)   / 'heat_flux_wall'           /
C----:
     .      cToken(TSECTION)        / 'section'                  /
     .      cToken(TEND)            / 'end'                      /


C---: Initial state
      nFindTokBC = 0

c      nTestLen =  nStrLen(str)

C---: Check all tokens (Note the top index has to be equal to the
C---: total number of tokens)
      do i=1, nKeyWords

        nLen = nStrLen(cToken(i))

        if(cToken(i)(1:nLen).eq.str(1:nTestLen)) then

          nFindTokBC = i
          exit

         end if

      end do

      return
      end


*----------------------------------------------------------------------*
*                             lReadFldProps                            *
*----------------------------------------------------------------------*

      LOGICAL function lReadFldProps(flPropFn,srchdFld,dSrchdTmp,
     .                              dRetFldProp)


      implicit none

      include "config.f"

      character*(*) flPropFn
      character*(*) srchdFld

      REAL dSrchdTmp
      REAL dRetFldProp(*)

      character*(mlinelgt) cBuf, cw
      LOGICAL lCheckFile, lStIsReal
      INTEGER mTrimLeft, mExtractW, nStrLen
      REAL dLinInterp

      LOGICAL lInsFlDef, lFoundSrFl, lFndInterv
      INTEGER nstat, nError, np, nFlDefLine, nl, nc, nw
      REAL dFldProp(6,2)

#ifndef _HP_

      REAL dnum

#endif _HP_


C---: Default status
      lReadFldProps = .false.

      if(.not.lCheckFile(flPropFn,nError)) then
        print*, 'Error code :', nError
        return
      end if

C---: Open the file
      open(10,file=flPropFn,status='old')

C---: Initialize logical flag for "found searched fluid" state
      lFoundSrFl = .false.

C---: Initialize logical flag for "found interval"
      lFndInterv = .false.

C---: Initialize fluid properties arrays
      do np=1,6
        dFldProp(np,1)  = 0.d0
        dFldProp(np,2)  = 0.d0
        dRetFldProp(np) = 0.d0
      end do

C---: Initialize line counters and cw counter
      nl = 0
      nFlDefLine = 0
      nw = 0

C---: Start reading lines
      lInsFlDef = .false.   !---> Inside fluid def. flag
      do
C-----: Increment line counter
        nl = nl + 1

C-----: Read line
        read(10,'(a)',iostat=nstat) cBuf
        if(nstat.lt.0) exit   !---> Reached end of file

C-----> Skip blank lines
        if(mTrimLeft(cBuf).lt.1) cycle

        if(cBuf(1:1).eq.'#') then
          if(lInsFlDef) then
            print*, 'Error: Comments not allowed inside Fluid Def.'
            return
          end if
          cycle
        end if

C-----> Convert entire line to lower case
        call mToLowerC(cBuf)   

C-----: Start parsing words
        nw = 0
        do
          if(mTrimLeft(cBuf).lt.1) exit
          nc = mExtractW(cBuf,cw)

C-------> End of line (empty cw)
          if(nc.lt.1) exit

C-------> Increment word counter
          nw = nw + 1

          if(cw(1:nc).eq.'fluid'.and.(.not.lInsFlDef)) then

            nc = mExtractW(cBuf,cw)
            nw = nw + 1   !---> Increment word counter
            if(nc.lt.1) call SyntaxError(nl,nw)

C---------: Initialize fluid definition line counter
            nFlDefLine = 0

            if(cw(1:nc).eq.srchdFld(1:nStrLen(srchdFld))) then
              lFoundSrFl = .true.

C-----------: Read next word and check if it is gas constant
              if(mTrimLeft(cBuf).gt.0) then
                nc = mExtractW(cBuf,cw)
                nw = nw + 1   !---> Increment word counter
                if(cw(1:nc).eq.'gasconstant') then
                  nc = mExtractW(cBuf,cw)
                  nw = nw + 1   !---> Increment word counter
                  if(nc.lt.1) call SyntaxError(nl,nw)
                  if(.not.lStIsReal(cw))
     .              call SyntaxError(nl,nw)
                  dRetFldProp(7) = dnum(cw(1:nc))
                end if
              end if

            end if

            lInsFlDef = .true.
            exit

          else if(cw(1:nc).eq.'end'.and.lInsFlDef) then
            lInsFlDef = .false.

C---------: If reached end of searched fluid, and found interval
            if(lFoundSrFl.and.lFndInterv) then

C-----------: Fill out returning fluid properties array
              do np=1,6
                dRetFldProp(np) = dFldProp(np,1)
              end do
C-----------: Change returning value of function 
              lReadFldProps = .true.
              return

C---------: If reached end of searched fluid, but did not reach
C---------: find interval, print error
            else if(lFoundSrFl.and.(.not.lFndInterv)) then
              print*, 'Did not find requested interval'
              return
            end if

          else if(lInsFlDef.and.lFoundSrFl) then

            nFlDefLine = nFlDefLine + 1   !---> Inc. fl def line ctr.

C---------: If already found interval, don't do anything else
            if(lFndInterv) exit

C---------: Copy current line of properties to next line
            do np=1,6
              dFldProp(np,2) = dFldProp(np,1)
            end do

C---------: Read first property (temperature)
            if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
            dFldProp(1,1) = dnum(cw(1:nc))

C---------: Check that wanted temp. value is not < first entry
            if(nFlDefLine.eq.1) then
              if(dFldProp(1,1).gt.dSrchdTmp) then
                print*, '* Error: Searched value (',dSrchdTmp,
     .            ') is lower than first entry in table'
                return
              end if 
            else
C-----------: Check if temp. in second line is > than in first line
              if(dFldProp(1,2).gt.dFldProp(1,1)) then
                print*, '* Error in entry ', nFlDefLine  
                print*, '* Entries should be arranged in ascending '
                print*, '* order with respect to temperature'
                return
              end if
            end if

C---------: Read the other properties
            do np = 2,6
              nc = mExtractW(cBuf,cw)
              nw = nw + 1   !---> Increment word counter
              if(nc.lt.1) call SyntaxError(nl,nw)
              if(.not.lStIsReal(cw)) call SyntaxError(nl,nw)
              dFldProp(np,1) = dnum(cw(1:nc))
            end do

C---------: Check if searched value lies in interval
            if (dSrchdTmp.le.dFldProp(1,1).and.
     .          dSrchdTmp.ge.dFldProp(1,2)) then

              lFndInterv = .true.

C-----------: Interpolate fluid properties
              do np=2,6
                dFldProp(np,1)=dLinInterp(dSrchdTmp,dFldProp(1,1),
     .            dFldProp(1,2),dFldProp(np,1),dFldProp(np,2))
              end do
              dFldProp(1,1) = dSrchdTmp

            end if

            exit   !---> Don't read any more words from this line

          else if(lInsFlDef.and.(.not.lFoundSrFl)) then
            nFlDefLine = nFlDefLine + 1   !---> Inc. fl def line ctr.
            exit

          else
            call SyntaxError(nl,nw)
          end if

        end do  !---> Stop parsing words

      end do   !---> Stop reading lines

C---: Close input file
      close(10)

      return
      end


*----------------------------------------------------------------------*
*                              InitBCFlags                             *
*----------------------------------------------------------------------*


      subroutine InitBCFlags(nReg, nRegType, nMomBdTp,
     .                       nTRgType, nTemBdTp, dBCVal)


      implicit none

      include "config.f"

C---: Arrays corresponding to subroutine argumengs
C---: Number  of regions in each direction
      INTEGER nReg(2)

C---: Grid regions: Types and indeces (temporary)
      INTEGER nRegType(mgri,mgrj)

C---: Boundary condition values (temporary)
      REAL dBCVal(mgri,mgrj,4,4)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Boundary type (temporary)
      INTEGER nMomBdTp(mgri,mgrj,4)

C---: Boundary type (temporary)
      INTEGER nTemBdTp(mgri,mgrj,4)

      INTEGER ireg, jreg, k, l

      REAL   d0

      d0 = 0.d0

C-----: Initialize boundary condition values to zero
        do l=1, 4   !--- For all variables 1: u, 2: v, 3: P, 4: T
          do k=1, 4   !--- For all cell faces: WEST, EAST, SOUTH, NORTH
            do jreg=1, nReg(_J_)   !--- For all regions in J direction
              do ireg=1, nReg(_I_)   !--- For all regions in I direction
                dBCVal(ireg,jreg,k,l) = d0
              end do
            end do
          end do
        end do


C-----: Initialize region types

C-----: The default is RM_INTERN: Internal region, unbstructed flow
C-----: (Unbstructed flow is the default region type.  Only define the
C-----: type if it is porous or a blockage.)
        do jreg=1, nReg(_J_)
          do ireg=1, nReg(_I_)
            nRegType(ireg,jreg) = RM_INTERN
          end do
        end do

C-----: Boundary definition: Type - Momentum
C-----: The default is interfacial boundary or wall, depending on
C-----: location

C-----: Set default: Internal-boundary conditions
        do jreg=1, nReg(_J_)
          do ireg=1, nReg(_I_)
            nMomBdTp(ireg,jreg,WEST)  = BM_INTERN
            nMomBdTp(ireg,jreg,EAST)  = BM_INTERN
            nMomBdTp(ireg,jreg,SOUTH) = BM_INTERN
            nMomBdTp(ireg,jreg,NORTH) = BM_INTERN
          end do
        end do

C-----: Initially, external boundaries have no-slip wall conditions
        do jreg=1, nReg(_J_)
          nMomBdTp(1, jreg,WEST)     = BM_WALL1
          nMomBdTp(nReg(_I_),jreg,EAST) = BM_WALL1
        end do

        do ireg=1, nReg(_I_)
          nMomBdTp(ireg,1,SOUTH)     = BM_WALL1
          nMomBdTp(ireg,nReg(_J_),NORTH) = BM_WALL1
        end do


C--------:
C--------: Region and Boundary definitions for thermal energy equation
C--------:

cC-----: Only if thermal energy equation flag is active
c        if (nthermen.eq.1) then    !---: Note: Start non-indented struct.

C-----: Region definition: Type - Temperature
C-----: The default is 0: No heat generation
        do jreg=1, nReg(_J_)
          do ireg=1, nReg(_I_)
            nTRgType(ireg,jreg) = RT_NOSRCE
          end do
        end do

C-----: Boundary declarations: Temperature
C-----: The default is interfacial boundary or adiabatic wall
C-----:

C-----: Set default: Internal-boundary conditions
        do jreg=1, nReg(_J_)
          do ireg=1, nReg(_I_)
            nTemBdTp(ireg,jreg,WEST) = BT_INTERN
            nTemBdTp(ireg,jreg,EAST) = BT_INTERN
            nTemBdTp(ireg,jreg,SOUTH) = BT_INTERN
            nTemBdTp(ireg,jreg,NORTH) = BT_INTERN
          end do
        end do

C-----: Set default of external boundaries to adiabatic wall conditions
        do jreg=1, nReg(_J_)
          nTemBdTp(1,jreg,WEST)     = BT_HTFLUX
          nTemBdTp(nReg(_I_),jreg,EAST) = BT_HTFLUX
        end do
        do ireg=1, nReg(_I_)
          nTemBdTp(ireg,1,SOUTH)      = BT_HTFLUX
          nTemBdTp(ireg,nReg(_J_),NORTH)  = BT_HTFLUX
        end do

      return
      end


*----------------------------------------------------------------------*
*                              mCellLocStr                             *
*----------------------------------------------------------------------*


      INTEGER function mCellLocStr(nc, word)

      implicit none
      INTEGER nc
      character*(*) word

      if(word(1:nc).eq.'w'.or.word(1:nc).eq.'west')       then
        mCellLocStr = WEST

      else if(word(1:nc).eq.'e'.or.word(1:nc).eq.'east')  then
        mCellLocStr = EAST

      else if(word(1:nc).eq.'s'.or.word(1:nc).eq.'south') then
        mCellLocStr = SOUTH

      else if(word(1:nc).eq.'n'.or.word(1:nc).eq.'north') then
        mCellLocStr = NORTH

      else
        print*, 'Unrecognized cell face designation: [',
     .    word(1:nc), ']'
        mCellLocStr = -1
      end if

      return
      end


*----------------------------------------------------------------------*
*                              lCheckIJreg                             *
*----------------------------------------------------------------------*


      LOGICAL function lCheckIJreg(ireg, jreg)

      implicit none

      include "config.f"

C---: Subroutine arguments
      INTEGER ireg, jreg

C---: Default status
      lCheckIJreg = .true.

C---: Test that ireg and jreg are consistent
      if (ireg.lt.1) then
        print*, 'Error: nReg(_I_) or ireg must be at least 1'
        lCheckIJreg = .false.
      else if (ireg.gt.mgri) then
        print*, 
     .    'Error: Requested nReg(_I_) or ireg is larger than mgri'
        lCheckIJreg = .false.
      else if (jreg.lt.1) then
        print*, 'Error: nReg(_J_) or jreg must be at least 1'
        lCheckIJreg = .false.
      else if (jreg.gt.mgrj) then
        print*, 
     .    'Error: Requested nReg(_J_) or jreg is larger than mgrj'
        lCheckIJreg = .false.
      end if

      return
      end

*-----------------------|---|---|---V---|---|---|----------------------*
