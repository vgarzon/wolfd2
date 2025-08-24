*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                                main.F                                *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
*                             w o l F D 2
*
*                           ( 2 D F l o w )
*
*  This program solves the incompressible  Navier-Stokes equations of
*  fluid flow in generalized  coordinate form.  The pressure-velocity
*  coupling  is  treated  via  Gresho's  'Projection  1' method.  The
*  spatial  derivatives  are  discretized  on finite  volumes using a
*  cell-based  indexing.  The  integration  of  time  derivatives  is
*  carried  out  by  applying  the  Douglas-Gunn  ADI  time-splitting
*  method.  In addition,  the thermal  energy  equation can be solved
*  and buoyancy terms can be included as source terms in the momentum
*  equations.  The modified Shuman filter is used to control aliasing
*  and cell-Reynolds (or cell-Peclet) number oscillations.
*
*  An  important  feature  of this  code  is the use of the  Additive
*  Turbulent  Decomposition (ATD) formalism of Hylin and McDonough to
*  model  turbulent  fluctuations.  The  version  of ATD  implemented
*  herein uses linear combinations of (properly-scaled)  chaotic maps
*  to simulate the small-scale part of the solution.
*
*  The  current  implementation  can also model flow  through  porous
*  media by solving the momentum and thermal energy equations modifed
*  accordingly.  The code also features particle trajectory  tracking
*  capabilities.
*
*
*                                                    Victor E. Garzon
*                                                <victor@ewl.uky.edu>
*
*                                                                 and
*
*                                                  James M. McDonough
*                                             <super180@ukcc.uky.edu>
*
*                             Computational Fluid Dynamics Laboratory
*                                     Dept. of Mechanical Engineering
*                                              University of Kentucky
*
*  Version:       0.3
*  Revision:        i
*  Date:     99.06.13
*
*----------------------------------------------------------------------*

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"

/*
  The following compiler directive is required for
  the secnds implicit timing function
*/
#ifdef _HP_

$NOSTANDARD SYSTEM

#endif


*----------------------------------------------------------------------*
*                                 main                                 *
*----------------------------------------------------------------------*

      program wolfd2

      implicit none

C---: Configurable parameters
      include "config.f"

C---: Declare basic file names
      character*(mfnmlgth)  inputFile, bndryFile,  gridFile, 
     .                      partFile,  resReadFl,  resWriteFl,
     .                      logFile
      character*(mfnmlgth-4) outPrefix

#ifdef _HP_
      character*(mfnmlgth) wrkStr
#endif
      character*(mfnmlgth)   ssFile  !---> Snapshots file
      character*(mfnmlgth-4) ssPrefix
      character*3 ssNumber  !---> Assume a max 999 snapshots
      character*(mfnmlgth-4) TSPrefix  !---> Time-series files prefix

C**** Temporary
      parameter (partFile='particles.dat')

C---: Grid information variables
      REAL    gx(0:mnx,0:mny),  gy(0:mnx,0:mny)

C---: Metric information
      REAL    ran(0:mnx,0:mny), rbn(0:mnx,0:mny), rgn(0:mnx,0:mny),
     .        rac(0:mnx,0:mny), rbc(0:mnx,0:mny), rgc(0:mnx,0:mny),
     .        rau(0:mnx,0:mny), rbu(0:mnx,0:mny),
     .        rbv(0:mnx,0:mny), rgv(0:mnx,0:mny),
     .        dju(0:mnx,0:mny), djv(0:mnx,0:mny),
     .        djc(0:mnx,0:mny), djn(0:mnx,0:mny)

      REAL    xzn(0:mnx,0:mny), xen(0:mnx,0:mny),
     .        yzn(0:mnx,0:mny), yen(0:mnx,0:mny),   
     .        xzc(0:mnx,0:mny), xec(0:mnx,0:mny),
     .        yzc(0:mnx,0:mny), yec(0:mnx,0:mny), 
     .        xzv(0:mnx,0:mny), xeu(0:mnx,0:mny),
     .        yzv(0:mnx,0:mny), yeu(0:mnx,0:mny), 
     .        xzu(0:mnx,0:mny), xev(0:mnx,0:mny),
     .        yzu(0:mnx,0:mny), yev(0:mnx,0:mny)

C---: Large-scale variables: Velocity components u and v,
C---: pressure p and temperature t at time-levels one and n
C---: respectively.
      REAL    u(0:mnx,0:mny),   v(0:mnx,0:mny),
     .        p(0:mnx,0:mny),   t(0:mnx,0:mny),
     .        un(0:mnx,0:mny),  vn(0:mnx,0:mny),
     .        pn(0:mnx,0:mny),  tn(0:mnx,0:mny)

C---: Small-scale variables
      REAL    uss(0:mnx,0:mny), vss(0:mnx,0:mny),
     .        pss(0:mnx,0:mny), tss(0:mnx,0:mny),
     .        usn(0:mnx,0:mny), vsn(0:mnx,0:mny),
     .        psn(0:mnx,0:mny), tsn(0:mnx,0:mny)

C---: Temporary variables used during momentum-thermal energy iter.
      REAL    us(0:mnx,0:mny),  vs(0:mnx,0:mny),
     .        ts(0:mnx,0:mny)

C---: Density at time-levels one and n
      REAL    d(0:mnx,0:mny), dn(0:mnx,0:mny)

      logical meitconv

      REAL    dif(4)

C---> Time-series pts. indeces
      REAL    iTS(maxtspts), jTS(maxtspts)

C---: Flags
      INTEGER nthermen, neqstate, nsmallscl, ntraject,
     .        nrestart, nreswrite, nsnapshots, nTSOutFlag,
     .        lPrntDiff, lPrDfBann
      INTEGER nts, mqiter, nmeiter, nfiltu, nfiltv, nfiltt 
      INTEGER nPpeSolver, msorit, nssPpeSlvr, mssSorIt
      INTEGER nPostProc, nOutFlForm, nWrtRstFq, isfreq, isfirst,
     .        nTSFreq, nTSPoints, nPrDfFreq, nPrDfBnFq

      REAL    dk, qtol, dmeittol, fpu, fpv, fpt
      REAL    sortol, sorrel
      REAL    ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .        ssBnCrit, ssRMpMax, ssRMpExp,
     .        ssSorTol, ssSorRel

C---: Reference quantities
      REAL    dlref, uref, tref, tmax,
     .        densref, rconst, dconduct

C---: Scaling numbers:  Reynolds, Peclet, Froude
      REAL    re, pe, fr

C---: Number of grid points in each direction
      INTEGER nx, ny

C---: Typed procedures
      INTEGER nStrLen
      INTEGER nAuxMomentum
      REAL    DiffMaxNorm

C---: Time-averaging mode flag
      LOGICAL lTmAvgFlg

C---: Cartesian grid flag
      LOGICAL lCartesGrid

/*
  Only declare time-averaged variables when _TIMEAVG_ is defined
*/

#ifdef _TIMEAVG_

C---: Time averaged quantities:
      REAL    ubar(0:mnx,0:mny),   vbar(0:mnx,0:mny),
     .        tbar(0:mnx,0:mny),   pbar(0:mnx,0:mny),
     .        upb(0:mnx,0:mny),    vpb(0:mnx,0:mny),
     .        tpb(0:mnx,0:mny),
     .        upupb(0:mnx,0:mny),  vpvpb(0:mnx,0:mny),
     .        uptpb(0:mnx,0:mny),  vptpb(0:mnx,0:mny),
     .        upvpb(0:mnx,0:mny),  trbke(0:mnx,0:mny),
     .        dssrt(0:mnx,0:mny),  dtdyb(0:mnx,0:mny),
     .        upxsb(0:mnx,0:mny),  upysb(0:mnx,0:mny),
     .        vpxsb(0:mnx,0:mny),  vpysb(0:mnx,0:mny)


#endif

C---: Time averaged output files prefix
      character*(mfnmlgth-4) tAvgOutPref

C---: Small-scale Shuman filter parameters
      REAL    ssFiltPar(4)

C---: Number of regions in each direction
      INTEGER nReg(2)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,4)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,4)

C---: Boundary type (temporary)
      INTEGER nTemBdTp(mgri,mgrj,4)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,4) 

C---: Porous region dimensionless coefficients and local porosity
      REAL    dPRporos(mgri,mgrj), dPRporc1(mgri,mgrj), 
     .        dPRporc2(mgri,mgrj) 

C---: Temperature region values
      REAL    dTRgVal(mgri,mgrj)

C---: Dimensionless heat generation source term values
      REAL    dHGSTval(mgri,mgrj)


#ifdef _TRAJECT_

C---: Variables used in particle trajectory calculations
      INTEGER ntr, itrajfreq, ntsubstp,
     .        nTrMethod, nTrCdEq, mTrHTmit

      INTEGER nTOutBnd(mntr)
      REAL    dTrHTtol, dTrHTdel
      REAL    cpartx(mntr), cparty(mntr), repc(mntr)
      REAL    xp(mntr), yp(mntr), up(mntr), vp(mntr)

#endif

C---: Counters used in solution snap-shots
      INTEGER isnum

C---: Counters used during time-averaging mode
      INTEGER nTimeAvg, nTmAvStart, nTmAvEnd

C---: Time-stepping variables
      INTEGER ks
      REAL    dtime, rdimtime, difmax

C---: PPE convergence return value
      INTEGER nSorConv

C---: Return value of auxiliary-momentum procedure
      INTEGER nQLiter

#ifdef _TIMEAVG_

C---: Temporary variables used in time-averaging mode
      REAL    dupij,  dvpij,  dtpij,  dnts,
     .        dupipj, dupimj, dupijp, dupijm, 
     .        dvpipj, dvpimj, dvpijp, dvpijm,
     .        upzi,   upet,   vpzi,   vpet,
     .        dupx,   dupy,   dvpx,   dvpy

C---: Temporary variables used in calculation of Stanton number
      REAL    tzi, tet

#endif

#ifdef _HP_
C---> Execution timing variables
      real*4 rtsec, rthrs, rtmns

C---> Variables used during parsing of command line args.
      INTEGER nargs, nSkipNext, nchar
#endif

C---: Local indicial variables
      INTEGER i, j, k, l

C---: Double precision constants
      REAL    dZero, dHalf, dQrtr
      parameter ( dZero = 0.00d+00 ,
     .            dHalf = 0.50d+00 ,
     .            dQrtr = 0.25d+00 )



#ifdef _HP_
C---: Initialize runtime variable
      rtsec = secnds(0.0)
#endif


C---: Default input and boundary filenames
      inputFile(1:nStrLen('input.dat')) = 'input.dat'
      bndryFile(1:nStrLen('bcond.dat')) = 'input.dat'

#ifdef _HP_

C---: Note: The following uses HP Fortran extensions to ansi f77:
C---: Number of command line arguments
      nargs = iargc()

      if (nargs.gt.0) then  !---> Get file names from argument list

C-----: Parse command line arguments
        nSkipNext = 0
        do i=1, nargs

          if(nSkipNext.eq.1) then
            nSkipNext = 0
            cycle
          end if
          nchar = igetarg(i,wrkStr,mfnmlgth)

          select case(wrkStr)
          case('-i')
            nchar = igetarg(i+1,wrkStr,mfnmlgth)
            if(nchar.le.0.or.wrkStr.eq.'-b') then
              print*, 'Usage: wolfd2 [-i inputFile] [-b bndryFile]'
              stop
            end if
C---------: Initialize working string
            call mInitStr(inputFile)
            inputFile(1:nchar) = wrkStr(1:nchar)

            call mInitStr(bndryFile)
            bndryFile(1:nchar) = wrkStr(1:nchar)
            nSkipNext = 1

          case('-b')
            nchar = igetarg(i+1,wrkStr,mfnmlgth)
            if(nchar.le.0.or.wrkStr.eq.'-i') then
              print*, 'Usage: wolfd2 [-i inputFile] [-b bndryFile]'
              stop
            end if
C---------: Initialize working string
            call mInitStr(bndryFile)
            bndryFile(1:nchar) = wrkStr(1:nchar)
            nSkipNext = 1

          case default
            print*, 'Usage: wolfd2 [-i inputFile] [-b bndryFile]'
            stop

          end select

        end do

        if(nstrlen(inputFile).lt.1.or.nstrlen(bndryFile).lt.1) then
          print*, 'Usage: wolfd2 [-i inputFile -b bndryFile]'
          stop
        end if

      end if

#endif _HP_

C---: Read input parameters
      call ReadParam(inputFile, gridFile,
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
     .               dlref, uref, tref, tmax,
     .               densref, rconst, dconduct,
     .               re, pe, fr)


cC---: Print title to standard output
c      call PrintTitle(STDOUT)


C---: Write input parameters and other info to standard output
      call LogInputPar(STDOUT, inputFile, bndryFile, gridFile,
     .       nthermen, neqstate, nsmallscl, ntraject,
     .       nts, mqiter, nmeiter, nfiltu, nfiltv, nfiltt,
     .       nPpeSolver, msorit,
     .       dk, qtol, dmeittol, fpu, fpv, fpt, sortol, sorrel,
     .       ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .       dlref, uref, tref, tmax, re, pe, fr)

C---: Open the log file
      logFile = outPrefix(1:nstrlen(outPrefix))//'.log'
      open(STDLOG, file=logFile, status='unknown')

C---: Print title to standard output
      call PrintTitle(STDLOG)

C---: Write input parameters and other info to log file
      call LogInputPar(STDLOG, inputFile, bndryFile, gridFile,
     .       nthermen, neqstate, nsmallscl, ntraject,
     .       nts, mqiter, nmeiter, nfiltu, nfiltv, nfiltt,
     .       nPpeSolver, msorit,
     .       dk, qtol, dmeittol, fpu, fpv, fpt, sortol, sorrel,
     .       ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .       dlref, uref, tref, tmax, re, pe, fr)



C---: Read grid from file and compute metric information
      call Grid(gridFile,
     .          nx, ny,
     .          dlref,
     .          gx, gy,
     .          rau, rbu, rbv, rgv,
     .          ran, rbn, rgn,
     .          rac, rbc, rgc,
     .          dju, djv, djc, djn,
     .          xen, yen, xzn, yzn,
     .          xec, yec, xzc, yzc,
     .          xeu, yeu, xzv, yzv,
     .          xzu, yzu, xev, yev)

C---: Set-up boundary conditions
      call SetUpBCs(nx, ny, bndryFile,
     .              nthermen, neqstate, nsmallscl, ntraject,
     .              nReg, nRegBrd,
     .              nRegType, nTRgType, nMomBdTp, nTemBdTp,
     .              dlref, uref, tref, tmax,
     .              densref, rconst, dconduct,
     .              re, pe,
     .              dPRporos, dPRporc1, dPRporc2,
     .              dTRgVal, dHGSTval,
     .              dBCVal)


#ifdef _TRAJECT_

C---: Initialize particle trajectory data (if requested)
      if(ntraject.eq.1)
     .   call InitTraject(partFile,
     .                    ntr, nts,
     .                    itrajfreq, nTSubStp,
     .                    nTOutBnd,
     .                    nTrMethod, nTrCdEq, mTrHTmit,
     .                    dlref, densref, re,
     .                    dTrHTtol, dTrHTdel,
     .                    xp, yp, up, vp,
     .                    cpartx, cparty, repc)

#else

      if(ntraject.eq.1) then

        ntraject = 0

        write(STDOUT,'(2a)') '* Warning: Executable not compiled ',
     .    'for trajectory calculations.'
        write(STDOUT,'(2a)') '* Recompile with directive -D_TRAJECT_'
        write(STDOUT,*)

        write(STDLOG,'(2a)') '* Warning: Executable not compiled ',
     .    'for trajectory calculations.'
        write(STDLOG,'(2a)') '* Recompile with directive -D_TRAJECT_'
        write(STDLOG,*)

      end if

#endif


C---: Initialize solution snapshots counter (if requested)
      if(nsnapshots.eq.1) isnum   = isfirst

C---: Initialize time series files (if requested)
      if(nTSOutFlag.eq.1) then
        rdimtime = dZero
        call SaveTimeSrs(nx, ny, 0,
     .                   TSPrefix,
     .                   nthermen, nsmallscl,
     .                   nTSPoints,
     .                   iTS, jTS,
     .                   rdimtime,
     .                   u, v, p, t, uss, vss, pss, tss)
      end if

#ifdef _TIMEAVG_

C---: If in time-averaging mode
      if(lTmAvgFlg) then

C-----: Time average loop indeces
        nTmAvStart = 1
        nTmAvEnd   = 2

C-----: Initialize time-average variables (if requested)
        do j=0, ny+1
          do i=0, nx+1
            ubar(i,j)  = dZero
            vbar(i,j)  = dZero
            tbar(i,j)  = dZero
            pbar(i,j)  = dZero
            upb(i,j)   = dZero
            vpb(i,j)   = dZero
            tpb(i,j)   = dZero
            upupb(i,j) = dZero
            vpvpb(i,j) = dZero
            uptpb(i,j) = dZero
            vptpb(i,j) = dZero
            upxsb(i,j) = dZero
            upysb(i,j) = dZero
            vpxsb(i,j) = dZero
            vpysb(i,j) = dZero
            upvpb(i,j) = dZero
            trbke(i,j) = dZero
            dssrt(i,j) = dZero
            dtdyb(i,j) = dZero
          end do
        end do

      else  !---> Normal mode

        nTmAvStart = 0
        nTmAvEnd   = 0

      end if

#else

      if(lTmAvgFlg) then

        lTmAvgFlg = LFALSE
        nTmAvStart = 0
        nTmAvEnd   = 0

        write(STDOUT,'(2a)') '* Warning: Executable not compiled ',
     .    'for time-averaging calculations.'
        write(STDOUT,'(2a)') '* Recompile with directive -D_TIMEAVG_'
        write(STDOUT,*)

        write(STDLOG,'(2a)') '* Warning: Executable not compiled ',
     .    'for time-averaging calculations.'
        write(STDLOG,'(2a)') '* Recompile with directive -D_TIMEAVG_'
        write(STDLOG,*)

      end if

#endif

C---: Start time-averaged loop (if requested) ---> Note: Not indented
      do nTimeAvg = nTmAvStart, nTmAvEnd

C---: On the second pass, disable unnecessary post-processing flags.
C---: (i.e., nreswrite, ntraject, nsnapshots, nTSOutFlag, lPrntDiff)
      if (nTimeAvg.gt.1) then
        nreswrite  = 0
        ntraject   = 0
        nsnapshots = 0
        nTSOutFlag = 0
        lPrntDiff  = 0

        write(STDOUT,*) 
        write(STDOUT,'(2a)')
     .    '* Starting second pass for calculating time-averaged ',
     .    'quantities...'
        write(STDOUT,*)

        write(STDLOG,*)
        write(STDLOG,'(2a)') 
     .    '* Starting second pass for calculating time-averaged ',
     .    'quantities...'
        write(STDLOG,*) 

      end if


C---: Initial conditions
      call InitialCond(nx, ny, nrestart, resReadFl, ks, dtime, 
     .   u, v, p, t, uss, vss, pss, tss)


C---: If cold start, enforce global mass conservation
      if(nrestart.eq.0) then

          call VelBoundCond(nx, ny, nReg, nRegBrd,
     .                      nMomBdTp, dBCVal,
     .                      u, v)

C-------: Solve Pressure Poisson Equation
          call Ppe(nx, ny,
     .             nReg, nRegBrd,
     .             nRegType,
     .             lCartesGrid,
     .             nPpeSolver, msorit, nSorConv, 
     .             dk, sortol, sorrel,
     .             rau, rbu, rbv, rgv,
     .             xeu, yeu, xzv, yzv,
     .             u, v, p)

          call PresBoundCond(nx, ny,
     .                       nReg, nRegBrd,
     .                       nRegType, nMomBdTp, dBCVal,
     .                       p)

C-------: Project auxiliar velocity onto a solenoidal space
          call Project(nx, ny,
     .                 nReg, nRegBrd,
     .                 nRegType, nMomBdTp,
     .                 dk,
     .                 dju, djv,
     .                 yeu, xzv, yzu, xev,
     .                 p, u, v)

          call VelBoundCond(nx, ny, nReg, nRegBrd,
     .                      nMomBdTp, dBCVal,
     .                      u, v)

      end if

#ifdef _SMALLSCL_

C---: Initialize small-scale components
      if(nsmallscl.eq.1)
     .  call SmallScale(nx, ny, 0, nthermen,
     .                  lCartesGrid,
     .                  nReg, nRegBrd,
     .                  nRegType, nTRgType, nMomBdTp, nTemBdTp,
     .                  nssPpeSlvr, mssSorIt,
     .                  dlref, uref,tref, tmax,
     .                  dk, re, pe,
     .                  ssSorTol, ssSorRel,
     .                  ssFiltPar,
     .                  ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .                  ssBnCrit, ssRMpMax, ssRMpExp,
     .                  dTRgVal, dBCVal,
     .                  rau, rbu, rbv, rgv,
     .                  dju, djv, djc,
     .                  xeu, yeu, xzv, yzv,
     .                  xzu, yzu, xev, yev,
     .                  xec, yec, xzc, yzc,
     .                  u, v, t, uss,
     .                  vss, pss, tss )

#else

C---: If not compiled for ATD calculations
      if(nsmallscl.eq.1) then
        nsmallscl = 0
        write(STDOUT,'(2a)') '* Warning: Executable not compiled ',
     .    'for small-scale calculations.'
        write(STDOUT,'(2a)') '* Recompile with directive -D_SMALLSCL_'
        write(STDOUT,*)

        write(STDLOG,'(2a)') '* Warning: Executable not compiled ',
     .    'for small-scale calculations.'
        write(STDLOG,'(2a)') '* Recompile with directive -D_SMALLSCL_'
        write(STDLOG,*)

        stop
      end if

#endif

C-/-/-/-/-/-/-/-/-/-/-: Begin time stepping loop :-\-\-\-\-\-\-\-\-\-\-

C---: Start time stepping
      do k= ks+1, ks+nts

C-----: Dimensionless time
        dtime = dtime + dk

C-----: Initialize large-scale time-level n variables
        do j=0, ny+1
          do i=0, nx+1
            pn(i,j)  = p(i,j)
            un(i,j)  = u(i,j)
            vn(i,j)  = v(i,j)
            tn(i,j)  = t(i,j)
            dn(i,j)  = d(i,j)
          end do
        end do

#ifdef _SMALLSCL_

C---: If also calculating small-scale fluctuations
      if(nsmallscl.eq.1) then

C-----: Initialize small-scale time-level n variables
        do j=0, ny+1
          do i=0, nx+1
            usn(i,j) = uss(i,j)
            vsn(i,j) = vss(i,j)
            tsn(i,j) = tss(i,j)
          end do
        end do

C-----> Add small-scale velocity components to last iterate
        do j=0, ny+1
          do i=0, nx+1
            un(i,j) = un(i,j) + uss(i,j)
            vn(i,j) = vn(i,j) + vss(i,j)
            tn(i,j) = tn(i,j) + tss(i,j)
          end do
        end do

      end if

#endif

C-----: Default state of momentum-energy iterations convergence flag
        meitconv = .false.

        if(nthermen.ne.1.and.nmeiter.gt.0) nmeiter=1
C-----: Start momentum/energy iterations
        do l=1, nmeiter

C-------: Initialize starred quantities (momentum/energy iterations)
          do j=0,ny+1
            do i=0,nx+1
              us(i,j) = u(i,j)
              vs(i,j) = v(i,j)
              ts(i,j) = t(i,j)
            end do
          end do

C-------: Compute auxiliary velocities (auxiliary momentum equations)
C-------: The returning value is the number of iterations required.
C-------: If the number of iterations was larger than mqiter, the
C-------: returning value is negative.
          nQLiter = nAuxMomentum(nx, ny,
     .                           mqiter,
     .                           nReg, nRegBrd,
     .                           nRegType, nMomBdTp,
     .                           dk, re, fr, qtol,
     .                           dPRporos, dPRporc1, dPRporc2,
     .                           dBCVal,
     .                           ran, rbn, rgn,
     .                           rac, rbc, rgc,
     .                           dju, djv,
     .                           xec, yec, xzn, yzn,
     .                           xen, yen, xzc, yzc,
     .                           xeu, yeu, xzu, yzu,
     .                           xev, yev, xzv, yzv,
     .                           d, dn,
     .                           un, vn, us, vs)

C-------: If q-l iterations did not converge, write a warning to the log
          if(nQLiter.lt.0) then
            write(STDLOG,'(a,i4,a)') '-> Warning:' 
            write(STDLOG,'(a,i3,a)')
     .        '-> Q-l precedure failed to converge after ',
     .        mqiter, ' iterations.'
            write(STDLOG,'(3x,2(a,i5))') ' Time Step:', k,
     .                                 ', ME Iter  :', l
            write(STDLOG,*)
          end if


C-------: Apply Shuman filter to auxiliary velocities
          if(nfiltu.eq.1) call Filter(nx, ny, _U_,
     .                                nReg, nRegBrd,
     .                                nRegType, nMomBdTp, nTRgType,
     .                                fpu, us)

          if(nfiltv.eq.1) call Filter(nx, ny, _V_,
     .                                nReg, nRegBrd,
     .                                nRegType, nMomBdTp, nTRgType,
     .                                fpv, vs)

          call VelBoundCond(nx, ny, nReg, nRegBrd,
     .                      nMomBdTp, dBCVal,
     .                      us, vs)

          call PresBoundCond(nx, ny,
     .                       nReg, nRegBrd,
     .                       nRegType, nMomBdTp, dBCVal,
     .                       p)

C-------: Solve Pressure Poisson Equation
          call Ppe(nx, ny, 
     .             nReg, nRegBrd,
     .             nRegType,
     .             lCartesGrid,
     .             nPpeSolver, msorit, nSorConv, 
     .             dk, sortol, sorrel,
     .             rau, rbu, rbv, rgv,
     .             xeu, yeu, xzv, yzv,
     .             us, vs, p)

          call PresBoundCond(nx, ny,
     .                       nReg, nRegBrd,
     .                       nRegType, nMomBdTp, dBCVal,
     .                       p)


C-------: Update velocities to n+1 time level
          call Project(nx, ny,
     .                 nReg, nRegBrd,
     .                 nRegType, nMomBdTp,
     .                 dk,
     .                 dju, djv,
     .                 yeu, xzv, yzu, xev,
     .                 p, us, vs)


          call VelBoundCond(nx, ny, nReg, nRegBrd,
     .                      nMomBdTp, dBCVal,
     .                      us, vs)

          call PresBoundCond(nx, ny,
     .                       nReg, nRegBrd,
     .                       nRegType, nMomBdTp, dBCVal,
     .                       p)


C-------: Solve thermal energy equation (if requested)
          if(nthermen.eq.1)
     .      call ThermEnergy(nx, ny,
     .                       nReg, nRegBrd, 
     .                       nTRgType, nTemBdTp,
     .                       dk, pe,
     .                       dTRgVal, dHGSTval, dBCVal,
     .                       rau, rbu, rbv, rgv, djc,
     .                       xeu, yeu, xzv, yzv,
     .                       xec, yec, xzc, yzc,
     .                       un, vn, us, vs, tn, ts)


C-------: Compute density (ideal gas equation of state)
          if(neqstate.eq.1) call EqState(nx, ny, 
     .                           uref, densref, tmax, tref, rconst,
     .                           p, ts, d)

           dif(1) = DiffMaxNorm(nx,ny,u,us)
           dif(2) = DiffMaxNorm(nx,ny,v,vs)
           dif(3) = DiffMaxNorm(nx,ny,t,ts)

c           print 203, (dif(i),i=1,3)

C-------: Update starred variables
          do j=0, ny+1
            do i=0, nx+1
              u(i,j) = us(i,j)
              v(i,j) = vs(i,j)
              t(i,j) = ts(i,j)
            end do
          end do

C-------: Test for convergence of momentum/energy iterations
          difmax = max(dif(1),dif(2),dif(3))
          if(l.gt.1.and.difmax.lt.dmeittol) then
            meitconv = .true.
            exit
          end if

        end do   !---> End momentum/energy iterations


cC-----: If requested convergence was not achieved in the 
cC-----: momentum/energy iterations, print a warning to the effect
c        if(.not.meitconv)
c     .    print*, 'Warning: M-E iterations did not converge after ', 
c     .    nmeiter, ' iterations'


C-----: Apply Shuman filter to temperature
        if(nthermen.eq.1.and.nfiltt.eq.1)
     .    call Filter(nx, ny, _T_,
     .                nReg, nRegBrd,
     .                nRegType, nMomBdTp, nTRgType,
     .                fpt, t)

#ifdef _SMALLSCL_

C-----: If all solving for small-scale variables
        if(nsmallscl.eq.1) then

C-------: Substract small-scale velocity at time level n 
C-------: from complete velocity field at time level n+1
          do j=0, ny+1
            do i=0, nx+1
              u(i,j) = u(i,j) - usn(i,j)
              v(i,j) = v(i,j) - vsn(i,j)
              t(i,j) = t(i,j) - tsn(i,j)
            end do
          end do

C-------: Compute small-scale fluctuations at time level n+1
          call SmallScale(nx, ny, 1, nthermen,
     .                    lCartesGrid,
     .                    nReg, nRegBrd,
     .                    nRegType, nTRgType, nMomBdTp, nTemBdTp,
     .                    nssPpeSlvr, mssSorIt,
     .                    dlref, uref,tref, tmax,
     .                    dk, re, pe,
     .                    ssSorTol, ssSorRel,
     .                    ssFiltPar,
     .                    ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .                    ssBnCrit, ssRMpMax, ssRMpExp,
     .                    dTRgVal, dBCVal,
     .                    rau, rbu, rbv, rgv,
     .                    dju, djv, djc,
     .                    xeu, yeu, xzv, yzv,
     .                    xzu, yzu, xev, yev,
     .                    xec, yec, xzc, yzc,
     .                    u, v, t, 
     .                    uss, vss, pss, tss)

C-------: Add small-scale components at time level n+1 to
C-------: large-scale velocity (complete velocity)
          do j=0,ny+1
            do i=0,nx+1
              u(i,j) = u(i,j) + uss(i,j)
              v(i,j) = v(i,j) + vss(i,j)
              t(i,j) = t(i,j) + tss(i,j)
            end do
          end do

        end if

#endif

        call VelBoundCond(nx, ny, nReg, nRegBrd,
     .                    nMomBdTp, dBCVal,
     .                    u, v)

        call PresBoundCond(nx, ny,
     .                     nReg, nRegBrd,
     .                     nRegType, nMomBdTp, dBCVal,
     .                     p)

        if(nthermen.eq.1) call TempBoundCond(nx, ny,
     .                                       nReg, nRegBrd,
     .                                       nTRgType, nTemBdTp,
     .                                       dTRgVal, dBCVal,
     .                                       t)

C-----: Compute max-norm of differences between time-levels 1 and n
        dif(1) = DiffMaxNorm(nx,ny,pn,p)
        dif(2) = DiffMaxNorm(nx,ny,un,u)
        dif(3) = DiffMaxNorm(nx,ny,vn,v)
        dif(4) = DiffMaxNorm(nx,ny,tn,t)
        difmax = max(dif(1),dif(2),dif(3),dif(4))

C-----: Check for divergence of solutions ***[needs revision]***
        if(difmax.gt.1.d12) then
          print*, '* Solution diverged. Please reduce CFL number.'
          stop
        end if

C-----: Print max-norm of diff. if requested
        if(lPrntDiff.eq.1) then
          rdimtime = dtime*dlref/uref !---> Dimensional time

          call PrintDiff(nPrDfFreq, lPrDfBann, nPrDfBnFq,
     .                   nthermen, nsmallscl, k, ks, nQLIter, 
     .                   (l-1), nSorConv, rdimtime, dif)
        end if

C-----: Save time series of small-scale fluctuations (if requested)
        if(nTSOutFlag.eq.1) then
          rdimtime = dtime*dlref/uref !---> Dimensional time

          if (mod((k-ks),nTSFreq).eq.0)
     .      call SaveTimeSrs(nx, ny, 1,
     .                     TSPrefix,
     .                     nThermEn, nSmallScl,
     .                     nTSPoints,
     .                     iTS, jTS,
     .                     rdimtime,
     .                     u, v, p, t, uss, vss, pss, tss)
        end if

#ifdef _TRAJECT_

C-----: Compute trajectories (if requested)
        if(ntraject.eq.1) then

C-------: Compute average values and density at natural grid points
C-------: and store them in starred arrays (to save storage)
          call VelAvg(nx, ny, nReg, nRegBrd, nRegType,
     .                u,  v,  us,  vs)
          call VelAvg(nx, ny, nReg, nRegBrd, nRegType,
     .                un, vn, usn, vsn)
C-------: Note that this is actually d avg!
          call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, d,  ts)
C-------: Note that this is actually dn avg!
          call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, dn, tsn)

C-------: Integrage particle trajectories
          call Traject(nx, ny,
     .                 ntr, ntsubstp,
     .                 nTrMethod, nTrCdEq, mTrHTmit,
     .                 nTOutBnd,
     .                 dk,
     .                 densref, fr,
     .                 dTrHTtol, dTrHTdel,
     .                 cpartx, cparty, repc,
     .                 gx, gy, us, vs, usn, vsn,
     .                 ts, tsn,
     .                 xp, yp, up, vp)

C-------: Save trajectories
          if(mod((k-ks),itrajfreq).eq.0)
     .      call SaveTraject(ntr, dtime, xp, yp)

        end if

#endif


C-----: Checkpointing:  Save restart file (w/freq nWrtRstFq)
        if(nreswrite.eq.2.and.mod((k-ks),nWrtRstFq).eq.0) then

C-------: Log message
          write(STDLOG,'(a)') '*'
          write(STDLOG,'(a,i8)')
     .     '* Checkpointing at time step ', k
          write(STDLOG,'(a)') '*'

          call SaveRestart(nx, ny, resWriteFl, k-1, dtime, 
     .                     u, v, p, t, uss, vss, pss, tss)
        end if


C-----: Save solutions snapshot (if requested)
        if(nsnapshots.eq.1.and.mod((k-ks),isfreq).eq.0) then

C-------: Construct filename
          write(ssNumber, '(i3.3)') isnum
          ssFile = ssPrefix(1:nStrLen(ssPrefix))//ssNumber

C-------: Log message
          write(STDLOG,'(a)') '*'
          write(STDLOG,'(a,i8)')
     .     '* Solution snapshot at time step ', k
          write(STDLOG,'(a)') '*'

C-------: Compute average values at natural grid points
C-------: and store them in starred arrays (or level n for p)
          call VelAvg(nx, ny, nReg, nRegBrd, nRegType,
     .                u, v, us, vs)

          call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, p, pn)

c          call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, t, ts)
          call TAveraged(nx, ny, 0,
     .                   nReg, nRegBrd, nTRgType,
     .                   dTRgVal,
     .                   t, ts)

C-------: Compute small-scale average values at natural grid
C-------: locations and store them in time level n arrays
          call VelAvg(nx, ny, nReg, nRegBrd, nRegType,
     .                uss, vss, usn, vsn)

          call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, pss, psn)

c          call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, tss, tsn)
          call TAveraged(nx, ny, 1,
     .                   nReg, nRegBrd, nTRgType,
     .                   dTRgVal,
     .                   tss, tsn)

C-------: Store variables of interest (large and small scales)

          call SaveStdVarsP3D(nx, ny, nOutFlForm,
     .                       ssFile,
     .                       nthermen, nsmallscl,
     .                       LFALSE, LFALSE,
     .                       gx, gy,
     .                       us, vs, pn, ts, usn, vsn, psn, tsn)

          isnum = isnum + 1

          if(isnum.gt.999) then
            write(STDLOG,*)
     .       'Exceeding max number of snapshots (999)'
            nsnapshots = 0
          end if

        end if

#ifdef _TIMEAVG_

C-----: Accumulate variables required for time-averaging calculations
        select case (nTimeAvg)

C-----: 1. On first pass accumulate complete time-averaged data:
        case (1)

C-------: Compute average values at natural grid points
C-------: and store them in starred arrays (or level n for p)
          call VelAvg(nx, ny, nReg, nRegBrd, nRegType,
     .                u, v, us, vs)
          call TAveraged(nx, ny, 0,
     .                   nReg, nRegBrd, nTRgType,
     .                   dTRgVal,
     .                   t, ts)
          call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, p, pn)

          do j=1, ny
            do i=1, nx
              ubar(i,j) = ubar(i,j) + us(i,j)
              vbar(i,j) = vbar(i,j) + vs(i,j)
              tbar(i,j) = tbar(i,j) + ts(i,j)
              pbar(i,j) = pbar(i,j) + pn(i,j)
            end do
          end do
   
C-----: 2. On second pass, accumulate small-scale fluctuations
        case (2)

          call VelAvg(nx, ny, nReg, nRegBrd, nRegType,
     .                u, v, us, vs)
          call TAveraged(nx, ny, 0,
     .                   nReg, nRegBrd, nTRgType,
     .                   dTRgVal,
     .                   t, ts)

          do j=1, ny
            do i=1, nx
C-----------: Fluctuations: u', v' and T'
              dupij = us(i,j) - ubar(i,j)
              dvpij = vs(i,j) - vbar(i,j)
              dtpij = ts(i,j) - tbar(i,j)

              upb(i,j) = upb(i,j) + dupij
              vpb(i,j) = vpb(i,j) + dvpij
              tpb(i,j) = tpb(i,j) + dtpij

              upupb(i,j) = upupb(i,j) + dupij * dupij
              vpvpb(i,j) = vpvpb(i,j) + dvpij * dvpij
              upvpb(i,j) = upvpb(i,j) + dupij * dvpij
              uptpb(i,j) = uptpb(i,j) + dupij * dtpij
              vptpb(i,j) = vptpb(i,j) + dvpij * dtpij

C-----------: Gradients used to compute dissipation rate function
              dupipj = (u(i+1,j+1) + u(i+1,j) + u(i,j+1) + u(i,j)
     .                - ubar(i+1,j+1) - ubar(i+1,j) - ubar(i,j+1)
     .                - ubar(i,j)) * dQrtr
              dupimj = (u(i,j+1) + u(i,j) + u(i-1,j+1) + u(i-1,j)
     .                - ubar(i,j+1) - ubar(i,j) - ubar(i-1,j+1)
     .                - ubar(i-1,j)) * dQrtr
              dupijp = u(i,j+1) - ubar(i,j+1)
              dupijm = u(i,j)   - ubar(i,j)

              dvpipj = v(i+1,j) - vbar(i+1,j)
              dvpimj = v(i,j)   - vbar(i,j)
              dvpijp = (v(i+1,j+1) + v(i+1,j) + v(i,j+1) + v(i,j)
     .                - vbar(i+1,j+1) - vbar(i+1,j) - vbar(i,j+1)
     .                - vbar(i,j)) * dQrtr
              dvpijm = (v(i+1,j) + v(i,j) + v(i+1,j-1) + v(i,j-1)
     .                - vbar(i+1,j) - vbar(i,j) - vbar(i+1,j-1) 
     .                - vbar(i,j-1)) * dQrtr

C-----------: Here as elsewhere, we assume Delta_xi = Delta_eta = 1
              upzi = dupipj - dupimj
              upet = dupijp - dupijm
              vpzi = dvpipj - dvpimj
              vpet = dvpijp - dvpijm

C-----------: Compute gradients and store them in temporary variable
              dupx = djn(i,j)*( yen(i,j)*upzi - yzn(i,j)*upet)
              dupy = djn(i,j)*(-xen(i,j)*upzi + xzn(i,j)*upet)
              dvpx = djn(i,j)*( yen(i,j)*vpzi - yzn(i,j)*vpet)
              dvpy = djn(i,j)*(-xen(i,j)*vpzi + xzn(i,j)*vpet)

C-----------: Accumulate squares of gradients
              upxsb(i,j) = upxsb(i,j) + dupx * dupx
              upysb(i,j) = upysb(i,j) + dupy * dupy
              vpxsb(i,j) = vpxsb(i,j) + dvpx * dvpx
              vpysb(i,j) = vpysb(i,j) + dvpy * dvpy

C-----------: Gradient (y-deriv.) of complete temperature.
C-----------: Used to estimate the Stanton number.
              tzi = dHalf*( t(i+1,j+1) + t(i+1,j) - t(i,j+1) - t(i,j) )
              tet = dHalf*( t(i+1,j+1) + t(i,j+1) - t(i+1,j) - t(i,j) )
              dtdyb(i,j) = dtdyb(i,j) + djn(i,j)*(-xen(i,j)*tzi 
     .                   + xzn(i,j)*tet)
 
            end do
          end do

        end select

#endif

      end do   !---> End time stepping

C-/-/-/-/-/-/-/-/-/-/-/: End time stepping loop :\-\-\-\-\-\-\-\-\-\-\-


C---: Close time series files (if requested)
      if(nTSOutFlag.eq.1) then
        rdimtime = dZero
        call SaveTimeSrs(nx, ny, 2,
     .                   TSPrefix,
     .                   nthermen, nsmallscl,
     .                   nTSPoints,
     .                   iTS, jTS,
     .                   rdimtime,
     .                   u, v, p, t, uss, vss, pss, tss)
      end if

#ifdef _TRAJECT_

cC---: Save trajectories
c      if(ntraject.eq.1) call SaveTraject(ntr, xp, yp)

C---: Close trajectories output file
      if(ntraject.eq.1) close(13)

#endif

#ifdef _TIMEAVG_

      select case (nTimeAvg)

C---: On first time-averaging pass
      case (1)

C-----: Number of time-steps taken (casted in double precision form)
        dnts = dble(nts)

C-----: Divide time-averaged complete variables by number of time-steps
        do j=1, ny
          do i=1, nx
            ubar(i,j) = ubar(i,j) / dnts
            vbar(i,j) = vbar(i,j) / dnts
            tbar(i,j) = tbar(i,j) / dnts
            pbar(i,j) = pbar(i,j) / dnts
          end do
        end do

C---: On second time-averaging pass:
      case (2)

C-----: Number of time-steps taken (casted in double precision form)
        dnts = dble(nts)

C-----: Divide time-averaged small-scale variables by number of time-steps
        do j=1, ny
          do i=1, nx
            upb(i,j)   = upb(i,j)   / dnts
            vpb(i,j)   = vpb(i,j)   / dnts
            tpb(i,j)   = tpb(i,j)   / dnts
            upupb(i,j) = upupb(i,j) / dnts
            vpvpb(i,j) = vpvpb(i,j) / dnts
            upvpb(i,j) = upvpb(i,j) / dnts  !---> Reynolds stress (in 2D)
            uptpb(i,j) = uptpb(i,j) / dnts
            vptpb(i,j) = vptpb(i,j) / dnts

C---------: Averaged squares of gradients of velocity fluctuations
            upxsb(i,j) = upxsb(i,j) / dnts
            upysb(i,j) = upysb(i,j) / dnts
            vpxsb(i,j) = vpxsb(i,j) / dnts
            vpysb(i,j) = vpysb(i,j) / dnts

C---------: Dissipation rate function
            dssrt(i,j) = (upxsb(i,j) + upysb(i,j) + vpxsb(i,j)
     .                  + vpysb(i,j)) * uref * dlref / re

C---------: Specific turbulence kinetic energy (tke per unit mass)
            trbke(i,j) = (upupb(i,j) + vpvpb(i,j)) * dHalf

C---------: Averaged gradient of temperature (y-deriv.)
            dtdyb(i,j) = dtdyb(i,j) / dnts

          end do
        end do

        call SaveTmAvgP3D(nx, ny, nts, nOutFlForm, tAvgOutPref,
     .                    ubar, vbar, tbar, pbar,
     .                    upb, vpb, tpb,
     .                    upupb, vpvpb, upvpb, uptpb, vptpb,
     .                    trbke, dssrt, dtdyb)

      end select

#endif

      end do  !---> End time-averaging loop  --- End unindented loop



C---: Save restart information (at final time step)
      if(nreswrite.eq.1) 
     .  call SaveRestart(nx, ny, resWriteFl, k-1, dtime, 
     .                   u, v, p, t, uss, vss, pss, tss)


C---: Compute average values at natural grid points
C---: and store them in time level n arrays
      call VelAvg(nx, ny,  nReg, nRegBrd, nRegType,
     .            u, v, un, vn)
      call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, p, pn)
      call TAveraged(nx, ny, 0,
     .               nReg, nRegBrd, nTRgType,
     .               dTRgVal,
     .               t, tn)

C---: Compute small-scale average values at natural grid points
C---: and store them in time level n arrays
      call VelAvg(nx, ny, nReg, nRegBrd, nRegType,
     .            uss, vss, usn, vsn)
      call PTDAvg(nx, ny, nReg, nRegBrd, nRegType, pss, psn)
      call TAveraged(nx, ny, 1,
     .               nReg, nRegBrd, nTRgType,
     .               dTRgVal,
     .               tss, tsn)

C---: Store variables of interest (large and small scales)
       call SaveStdVarsP3D(nx, ny, nOutFlForm,
     .                     outPrefix,
     .                     nthermen, nsmallscl,
     .                     LTRUE, LTRUE,
     .                     gx, gy,
     .                     un, vn, pn, tn, usn, vsn, psn, tsn)


cC---: Compute average value of density at natural grid points
cC---: and store it in time level n array
c      call EqState(nx, ny, uref, densref, tmax, tref, rconst,
c     .              p, ts, d)
c      call PTDAvg(nx, ny,  nReg, nRegBrd, nRegType, d, dn)
c      call SaveFunc(nx, ny, nPostProc, 'den.qqq', dn)


#ifdef _HP_

C---: Calculate runtime and write it to standard output and log file

      rtsec = secnds(rtsec)
      rthrs = rtsec/3600.0
      rtmns = rtsec/60.0

      write (STDOUT,*)
      write (STDOUT,'(a,i4,a1,i2.2,a1,f5.2,6x,a,f10.2,a)') 
     .  '--> Execution time: ', 
     .  int(rthrs), ':', int(rtmns)-int(rthrs)*60,
     .  ':', amod(rtsec,60.0), '( ', rtsec, ' sec)'
      write (STDOUT,*)

      write (STDLOG,*)
      write (STDLOG,'(a,i4,a1,i2.2,a1,f5.2,6x,a,f10.2,a)') 
     .  '--> Execution time: ', 
     .  int(rthrs), ':', int(rtmns)-int(rthrs)*60,
     .  ':', amod(rtsec,60.0), '( ', rtsec, ' sec)'
      write (STDLOG,*)

#endif


C---: Close the log file
      close(STDLOG)


      stop
      end


*-----------------------|---|---|---V---|---|---|----------------------*
