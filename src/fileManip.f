*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                             fileManip.F                              *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*                                                                      *
*
* ReadParam     - Read program parameters from file
* SaveSolution  - Save pressure, velocity and temperature grid functions
* SaveGrid2DP3D - Save 2D grid information to file using Plot3D format
* SaveRestart   - Save restart information to file
* ReadRestart   - Read variables to restart computation
* SaveFunc      - Save a single grid function to file (Plot3D format)
* SaveTraject   - Save computed particle trajectories
* SaveFullSol   -
* SaveTimeSrs   -
* lCheckFile    -
* SaveTmAvg     -
*
*----------------------------------------------------------------------*

/*
  Preprocessor header file(s) 
*/

#include "wolfd2.h"

*----------------------------------------------------------------------*
*                              ReadParam                               *
*----------------------------------------------------------------------*

      subroutine ReadParam(inputFile, gridFile,
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

      implicit none

C---: Read configurable parameters
      include "config.f"

C---: Arguments
      character*(*) inputFile, gridFile, outPrefix
      character*(*) resReadFl, resWriteFl
      character*(*) ssPrefix, TSPrefix
      character*(*) tAvgOutPref

      LOGICAL lCartesGrid, lTmAvgFlg

      INTEGER nthermen, neqstate, nsmallscl, ntraject,
     .        nrestart, nreswrite, nsnapshots, nTSOutFlag, 
     .        lPrntDiff, lPrDfBann,
     .        nts, mqiter, nmeiter, nfiltu, nfiltv, nfiltt, 
     .        nPpeSolver, msorit, nPostProc, nOutFlForm, nWrtRstFq,
     .        isfreq, isfirst,
     .        nTSFreq, nTSPoints, nPrDfFreq, nPrDfBnFq,
     .        nssPpeSlvr, mssSorIt

      INTEGER iTS(*), jTS(*)   !---> T.S. pts. indeces

      REAL    dk, qtol, dmeittol, fpu, fpv, fpt, sortol, sorrel,
     .        ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .        ssBnCrit, ssRMpMax, ssRMpExp,
     .        ssSorTol, ssSorRel,
     .        dlref, uref, tref, tmax

      REAL    ssFiltPar(*)

      REAL    densref, rconst, dconduct, re, pe, fr


C---: Function declarations
      LOGICAL lReadFldProps, lReadIPFile

C---: Local definitions
C---: Fluid properties filename and table name
      character*(mfnmlgth) flPropFn, flName

C---: Fluid properties array
      REAL  dFldProp(7)

      REAL  dkinvisc, dprandtl, dgrav

C---: Parse input file
      if (.not.lReadIPFile(inputFile, gridFile, flPropFn, flName, 
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
     .               dlref, uref, tref, tmax)) then

        print*, 'Error parsing input file.  Aborting.'
        stop
      end if

C---: Scale physical time-step size wrt. reference velocity and length
      dk = dk * uref / dlref

C---: Read fluid themophysical properties from file
      if (.not.lReadFldProps(flPropFn,flName,tref,dFldProp)) then
        print*, '* Error reading thermophysical properties file.'
        stop
      end if

C---: Assign fluid properties from array
      tref     = dFldProp(1)
      densref  = dFldProp(2)
      dkinvisc = dFldProp(4) / densref
      dconduct = dFldProp(5)
      dprandtl = dFldProp(6)
      rconst   = dFldProp(7)


C***Temporary:
C***For now assume gravitational acceleration to be 9.81 m/s
      dgrav = 9.810d+00


C---: Calculate Reynolds, Peclet and Froude numbers
      re = dlref*uref/dkinvisc
      pe = re*dprandtl
      fr = uref**2.00d0/(dlref*dgrav)


      return
      end

*----------------------------------------------------------------------*
*                             SaveSolution                             *
*----------------------------------------------------------------------*
* Save pressure, velocity and temperature grid functions to file
* npostp -> Postprocessor flag:
*   1: Plain Plot3D format (2D: nx,ny)
*   2: FieldView 5.5 Plot3D format (2D: nx,ny,nvar)
*   3: FieldView 5.5 Unstructured format (3D)
*   4: Snapshot in Plot3D format (only solution file, no grid)
*
      subroutine SaveSolution(nx, ny, outPrefix, npostp,
     .                        x, y, u, v, p, t)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, npostp

      character*(*) outPrefix

      REAL    x(0:mnx,0:mny), y(0:mnx,0:mny)
      REAL    u(0:mnx,0:mny), v(0:mnx,0:mny),
     .        p(0:mnx,0:mny), t(0:mnx,0:mny)

      REAL    to(0:mnx,0:mny)

c      real*4 func(0:mnx,0:mny,10)

C---: External procedures
      INTEGER nstrlen

C---: Local variables
      INTEGER i, j, k, ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8
      character*(mfnmlgth) grdFile, solFile

C---: Double precision constants
      REAL    dZero
      parameter ( dZero = 0.00d+00 )


C---: Construct filenames
      grdFile = outPrefix(1:nstrlen(outPrefix))//'.xyz'
      solFile = outPrefix(1:nstrlen(outPrefix))//'.qqq'

C---: Chech for excessively small numbers (needs revision)
C????

C*** Temporary:  (yuk!)
 205  format(5(e14.6))


      select case(npostp)
      case(1)   !---> Plain Plot3D format (nx,ny)
c        call SaveGrid2D(nx, ny, grdFile, x, y)   !---> Save grid
        open(11,file=solFile,status='unknown')
          write(11,*) nx, ny
          write(11,205) ((p(i,j),i=1,nx),j=1,ny),
     .                  ((u(i,j),i=1,nx),j=1,ny),
     .                  ((v(i,j),i=1,nx),j=1,ny),
     .                  ((t(i,j),i=1,nx),j=1,ny)
        close(11)

      case(2)   !---> FieldView 5.5 Plot3D format (nx,ny,nvar)
c        call SaveGrid2D(nx, ny, grdFile, x, y)    !---> Save grid
        open(11,file=solFile,status='unknown')
          write(11,*) nx,ny,5
          write(11,205) ((p(i,j),i=1,nx),j=1,ny),
     .                  ((u(i,j),i=1,nx),j=1,ny),
     .                  ((v(i,j),i=1,nx),j=1,ny),
     .                  ((dZero ,i=1,nx),j=1,ny),
     .                  ((t(i,j),i=1,nx),j=1,ny)
        close(11)
        open(12,file='n.nam',status='unknown')
           write(12,*) 'Pressure'
           write(12,*) 'U-vel ; Velocity'
           write(12,*) 'V-vel'
           write(12,*) 'W-vel'
           write(12,*) 'Temperature'
           write(12,*)
        close(12)

      case(3)   !---> FieldView 5.5 unstructured format
        solFile = outPrefix(1:nstrlen(outPrefix))//'.fv'
        open(11, file=solFile, status='unknown')
          write(11,'(a)') 'FIELDVIEW 1 1'
          write(11,'(a)') 'Boundary Table'
          write(11,'(a)') '0'
          write(11,'(a)') 'Variable Names'
          write(11,*) '5'
          write(11,*) 'Pressure'
          write(11,*) 'U-vel ; Velocity'
          write(11,*) 'V-vel'
          write(11,*) 'W-vel'
          write(11,*) 'Temperature'
          write(11,'(a)') 'Nodes'
          write(11,*) (2*nx*ny)
          do k=0,1
            do j=1,ny
              do i=1,nx
                write(11,205) x(i,j),y(i,j),(0.1d0)*dble(k)
              end do
            end do
          end do
          write(11,'(a)') 'Boundary Faces'
          write(11,*) '0'
          write(11,'(a)') 'Elements'
          do j=1,ny-1
            do i=1,nx-1
              ne1 = (j-1)*nx+i
              ne2 = ne1+1
              ne3 = nx*ny+ne1
              ne4 = ne3+1
              ne5 = ne1+nx
              ne6 = ne5+1
              ne7 = nx*ny+ne5
              ne8 = ne7+1
              write(11,*) 2,1,ne1,ne2,ne3,ne4,ne5,ne6,ne7,ne8
            end do
          end do
          write(11,'(a)') 'Variables'
          write(11,205) (((p(i,j),i=1,nx),j=1,ny),k=1,2),
     .                  (((u(i,j),i=1,nx),j=1,ny),k=1,2),
     .                  (((v(i,j),i=1,nx),j=1,ny),k=1,2),
     .                  (((dZero ,i=1,nx),j=1,ny),k=1,2),
     .                  (((t(i,j),i=1,nx),j=1,ny),k=1,2)
        close(11)

      case(4)
        print*, '# Solution snapshot [',solFile,']'
        open(11,file=solFile,status='unknown')
          write(11,*) nx,ny,5
          write(11,205) ((p(i,j),i=1,nx),j=1,ny),
     .                  ((u(i,j),i=1,nx),j=1,ny),
     .                  ((v(i,j),i=1,nx),j=1,ny),
     .                  ((dZero ,i=1,nx),j=1,ny),
     .                  ((to(i,j),i=1,nx),j=1,ny)
        close(11)
      case default
        print*, 'Error: Wrong postproc. flag passed to SaveSolution: ', 
     .    npostp
        stop
      end select

      return
      end
 
*----------------------------------------------------------------------*
*                             SaveGrid2DP3D                            *
*----------------------------------------------------------------------*
*
* Save grid information in Plot3D format (2D)
*
      subroutine SaveGrid2DP3D(nx, ny, nForm, gridFile, x, y)

      implicit none

      include "config.f"

      INTEGER nx, ny, nForm
      character*(*) gridFile
      REAL    x(0:mnx,0:mny), y(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j

      select case(nForm)

      case (FT_FORMATTED)
        open(22, file=gridFile, status='unknown')
          write(22,*) nx, ny
          write(22,'(5(e14.6))') ((x(i,j), i=1,nx), j=1,ny),
     .                           ((y(i,j), i=1,nx), j=1,ny)
        close(22)

      case (FT_UNFORMATTED)
        open(22, file=gridFile, form='unformatted', status='unknown')
          write(22) nx, ny
          write(22) ((sngl(x(i,j)), i=1,nx), j=1,ny),
     .              ((sngl(y(i,j)), i=1,nx), j=1,ny)
        close(22)

      case default
        write(STDERR,'(a)')
     .    '* Wrong nForm flag passed to SaveGrid2DP3D'
        stop

      end select

      return
      end

*----------------------------------------------------------------------*
*                             SaveRestart                              *
*----------------------------------------------------------------------*
*
* Save restart information
*
      subroutine SaveRestart(nx, ny, 
     .                       resFile,
     .                       k,
     .                       dtime, 
     .                       u, v, p, t, uss, vss, pss, tss)

      implicit none

      include "config.f"

C---: Subroutine arguments
      INTEGER nx, ny
      character*(*) resFile
      INTEGER k
      REAL    dtime
      REAL    u(0:mnx,0:mny),   v(0:mnx,0:mny), 
     .        p(0:mnx,0:mny),   t(0:mnx,0:mny),
     .        uss(0:mnx,0:mny), vss(0:mnx,0:mny), 
     .        pss(0:mnx,0:mny), tss(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j

      open(12,file=resFile,form='unformatted',status='unknown')
        write(12) k, dtime
        write(12) nx, ny
        write(12) ((p(i,j),i=0,nx+1),j=0,ny+1),
     .            ((u(i,j),i=0,nx+1),j=0,ny+1),
     .            ((v(i,j),i=0,nx+1),j=0,ny+1),
     .            ((t(i,j),i=0,nx+1),j=0,ny+1),
     .            ((pss(i,j),i=0,nx+1),j=0,ny+1),
     .            ((uss(i,j),i=0,nx+1),j=0,ny+1),
     .            ((vss(i,j),i=0,nx+1),j=0,ny+1),
     .            ((tss(i,j),i=0,nx+1),j=0,ny+1)
      close(12)

      return
      end

*----------------------------------------------------------------------*
*                             ReadRestart                              *
*----------------------------------------------------------------------*
*
* Read restart information
*
      subroutine ReadRestart(nx, ny, 
     .                       resFile,
     .                       k,
     .                       dtime, 
     .                       u, v, p, t, uss, vss, pss, tss)

      implicit none

      include "config.f"

C---: Subroutine arguments
      INTEGER nx, ny
      character*(*) resFile
      INTEGER k
      REAL    dtime
      REAL    u(0:mnx,0:mny),   v(0:mnx,0:mny), 
     .        p(0:mnx,0:mny),   t(0:mnx,0:mny),
     .        uss(0:mnx,0:mny), vss(0:mnx,0:mny), 
     .        pss(0:mnx,0:mny), tss(0:mnx,0:mny)

C---: External procedures
      LOGICAL lCheckFile

C---: Local variables
      INTEGER i, j, nxfile, nyfile, nError

C---: Check that restart file is readable
      if(.not.lCheckFile(resFile, nError)) then
        print*, 'Error code :', nError
        stop
      end if

C---: Open file
      open(23,file=resFile,form='unformatted',status='unknown')
        read(23) k,dtime
        read(23) nxfile,nyfile

        if((nxfile.ne.nx).or.(nyfile.ne.ny)) then
          print*, 'Error: Index mismatch in restart file: ', 
     .            nx, ny, nxfile, nyfile
          stop
        end if

        read(23) ((p(i,j),i=0,nx+1),j=0,ny+1),
     .           ((u(i,j),i=0,nx+1),j=0,ny+1),
     .           ((v(i,j),i=0,nx+1),j=0,ny+1),
     .           ((t(i,j),i=0,nx+1),j=0,ny+1),
     .           ((pss(i,j),i=0,nx+1),j=0,ny+1),
     .           ((uss(i,j),i=0,nx+1),j=0,ny+1),
     .           ((vss(i,j),i=0,nx+1),j=0,ny+1),
     .           ((tss(i,j),i=0,nx+1),j=0,ny+1)
      close(23)

      return
      end

*----------------------------------------------------------------------*
*                               SaveFunc                               *
*----------------------------------------------------------------------*
*
* Save single grid function to file
* npostp -> Postprocessor flag:
*   1: Plain Plot3D format (nx,ny)
*   2: FieldView 5.5 Plot3D format (nx,ny,nvar)
*
      subroutine SaveFunc(nx, ny, npostp, solFile, f)

      implicit none

      include "config.f"

      INTEGER nx, ny, npostp
      character*(*) solFile
      REAL    f(0:mnx,0:mny)

      INTEGER i, j

      open(12, file=solFile, status='unknown')

      select case(npostp)

      case(1,2)
        write(12,*) nx, ny, 1

      case(3)
        print*, 'Warning: Flag npostp=3 not implemented in SaveFunc.'
        print*, '         Using npostp=2 instead.'

        write(12,*) nx, ny, 1

      case default
        print*, 'Error: Wrong npostp flag passed to SaveFunc: ', 
     .    npostp
        stop

      end select

      write(12,'(5e13.6)') ((f(i,j),i=1,nx),j=1,ny)

      close(12)

      return
      end
 
*----------------------------------------------------------------------*
*                             SaveTraject                              *
*----------------------------------------------------------------------*

      subroutine SaveTraject(ntr, dtime, xp, yp)

      implicit none

      include "config.f"

      INTEGER ntr
      REAL    dtime
      REAL    xp(mntr), yp(mntr)

      INTEGER i

      write(13,'(e13.6)') dtime

      do i=1, ntr
        write(13,'(2(1x,e13.6))') xp(i), yp(i)
      end do

      return
      end


*----------------------------------------------------------------------*
*                             SaveFullSol                              *
*----------------------------------------------------------------------*
*  nSaveGrid flag:
*     0 -> Don't save grid
*     1 -> Save grid as well
*

      subroutine SaveFullSol(nx, ny, nForm,
     .                       outPrefix, 
     .                       nThermEn, nSmallScl,
     .                       lSaveGrid, lSaveNamFl,
     .                       gx, gy,
     .                       u, v, p, t, us, vs, ps, ts)

      implicit none

      include "config.f"

C---: Subroutine arguments
      INTEGER nx, ny, nForm
      character*(*) outPrefix
      INTEGER nThermEn, nSmallScl
      LOGICAL lSaveNamFl, lSaveGrid
      REAL    gx(0:mnx,0:mny), gy(0:mnx,0:mny)
      REAL    u(0:mnx,0:mny),  v(0:mnx,0:mny),
     .        p(0:mnx,0:mny),  t(0:mnx,0:mny),
     .        us(0:mnx,0:mny), vs(0:mnx,0:mny), 
     .        ps(0:mnx,0:mny), ts(0:mnx,0:mny)
 
C---: External procedures
      INTEGER nStrLen

C---: Local variables
      INTEGER i, j, l, nVars, nLamV
      character*(mfnmlgth) grdFile, solFile, namFile
      FLOAT   f(0:mnx,0:mny,10)

C---: Double precision constants
      REAL    dZero
      parameter ( dZero = 0.00d+00 )

C---: Construct filenames
      solFile = outPrefix(1:nStrLen(outPrefix))//'.qqq'
      grdFile = outPrefix(1:nStrLen(outPrefix))//'.xyz'
      namFile = outPrefix(1:nStrLen(outPrefix))//'.nam'

C---: Save grid to file if requested
      if(lSaveGrid) call SaveGrid2DP3D(nx, ny, nForm, grdFile, gx, gy)

C---: Supress excessively small numbers (for reading into FieldView)
C---: Select variables for output according to flags
C---: Default number of variables
      nLamV = 4    !---> Laminar variables
      nVars = 4    !---> Total number of variables

      if(nThermEn.eq.1)  then
        nVars = nVars + 1    !---> +T
        nLamV = nLamV + 1
      end if

      if(nSmallScl.eq.1) nVars = 2 * nVars    !---> + Small-scale vars

C---: Default variables
      do j=1,ny
        do i=1,nx
          f(i,j,1) = p(i,j)
          f(i,j,2) = u(i,j)
          f(i,j,3) = v(i,j)
          f(i,j,4) = dZero      !---> 2D solution files
        end do
      end do

C---: Thermal Energy (laminar)
      if(nThermEn.eq.1) then
        do j=1,ny
          do i=1,nx
            f(i,j,5) = t(i,j)
          end do
        end do
      end if

C---: Small-scale variables
      if(nSmallScl.eq.1) then
        do j=1,ny
          do i=1,nx
            f(i,j,nLamV + 1) = ps(i,j)
            f(i,j,nLamV + 2) = us(i,j)
            f(i,j,nLamV + 3) = vs(i,j)
            f(i,j,nLamV + 4) = dZero      !---> 2D solution files
          end do
        end do

C-----: Thermal Energy (small-scale)
        if(nThermEn.eq.1) then
          do j=1,ny
            do i=1,nx
              f(i,j,nLamV + 5) = ts(i,j)
            end do
          end do
        end if
      end if


C---: FieldView 5.5 two-dimensional Plot3D format (nx,ny,nVars)
      open(15,file=solFile,status='unknown')
      write(15,*) nx, ny, nVars
C-----: Write to file using default format to save space
        write(15,*) (((f(i,j,l),i=1,nx),j=1,ny),l=1,nVars)
      close(15)

C---: Save names file (for FieldView) if requested
      if(lSaveNamFl) then
        open(16,file=namFile,status='unknown')
C-------: Default variable names
          write(16,*) 'Complete Pressure'
          write(16,*) 'Complete U ; Complete Velocity'
          write(16,*) 'Complete V'
          write(16,*) 'Complete W'

C-------: If solving thermal energy 
          if(nThermEn.eq.1) write(16,*) 'Complete Temperature'

C-------: If solving small-scale variables
          if(nSmallScl.eq.1) then
            write(16,*) 'Small-Scale Pressure'
            write(16,*) 'Small-Scale U ; Small-Scale Velocity'
            write(16,*) 'Small-Scale V'
            write(16,*) 'Small-Scale W'

C---------: If solving thermal energy 
            if(nThermEn.eq.1) write(16,*) 'Small-Scale Temperature'

          end if

        close(16)

      end if

      return
      end


*----------------------------------------------------------------------*
*                             SaveTimeSrs                              *
*----------------------------------------------------------------------*
*
      subroutine SaveTimeSrs(nx, ny, nCase,
     .                       TSPrefix,
     .                       nThermEn, nSmallScl,
     .                       nTSPoints,
     .                       iTS, jTS,
     .                       dtime,
     .                       u, v, p, t, us, vs, ps, ts)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, nCase
      INTEGER nThermEn, nSmallScl
      INTEGER nTSPoints

      INTEGER iTS(*), jTS(*)   !---> T.S. pts. indeces

      character*(*) TSPrefix

      REAL    dtime

      REAL    u(0:mnx,0:mny),  v(0:mnx,0:mny),
     .        p(0:mnx,0:mny),  t(0:mnx,0:mny)
      REAL    us(0:mnx,0:mny), vs(0:mnx,0:mny),
     .        ps(0:mnx,0:mny), ts(0:mnx,0:mny)

C---: External procedures
      INTEGER nStrLen

C---: Local variables
      INTEGER j, nts, nTSU
      INTEGER nLamV, nVars

      INTEGER nUnitMask

C---: Set a number above which to declare file units
      parameter(nUnitMask = 50)

      character*3 TSNumber  !---> Assume abs max 999 TS Points
      character*(mfnmlgth) TSFile   !---> Time-series output filename

      FLOAT dTSv(maxtspts,8)   !---> Working t.s. array (8 vars.)

      save nLamV, nVars
 

C---: Select type of operation
      select case(nCase)

      case(0)   !---> Initialize files and counter

C-----: Select variables for output according to flags
C-----: Default number of variables
        nLamV = 3    !---> Laminar variables
        nVars = 3    !---> Total number of variables

        if(nThermEn.eq.1)  then
          nLamV = nLamV + 1
          nVars = nVars + 1    !---> +T
        end if

        if(nSmallScl.eq.1) nVars = 2 * nVars    !---> + Small-scale vars

C-----: Log message
        write(STDLOG,'(a)') '* Initializing time-series files:'
        write(STDLOG,'(7x,a4,15x,a4,2x,2(a4,x))')
     .     'File', 'Unit', '-I-', '-J-'

C-----: Construct filenames and open files
        do nts=1, nTSPoints

C-------: Check that monitor points are inside domain
          if((iTS(nts).lt.1.or.iTS(nts).gt.nx).or.
     .       (jTS(nts).lt.1.or.jTS(nts).gt.ny)) then
            print*, 'Warning: Requested t.s. indeces are not in ',
     .        'domain. Using (1,1).'
            iTS(nts) = 1
            jTS(nts) = 1
          end if

C-------: Write suffix number to internal file
          write(TSNumber, '(i3.3)') nts

          TSFile = TSPrefix(1:nStrLen(TSPrefix))//TSNumber//'.ts'

          write(STDLOG,'(4x,a20,1x,i4,2x,2(1x,i4))')
     .      TSFile, nUnitMask+nts, iTS(nts), jTS(nts)

          open(nUnitMask+nts,file=TSFile,status='unknown')

C-------: Write a comment in header of file
          nTSU = nUnitMask+nts
          write(nTSU,'(a)')     '#'
          write(nTSU,'(a,i4)')  '# Time-series No. ', nts
          write(nTSU,'(a,2i4)') '# Location: ', iTS(nts), jTS(nts)
          write(nTSU,'(a)')     '#'

          if(nSmallScl.eq.1) then
            if(nThermEn.eq.1) then
              write(nTSU,'(a,8a14)') '#    Time','u','v','P','T',
     .          'u*','v*','P*','T*'
            else
              write(nTSU,'(a,8a14)') '#    Time','u','v','P',
     .          'u*','v*','P*'
            end if
          else
            if(nThermEn.eq.1) then
              write(nTSU,'(a,8a14)') '#    Time','u','v','P','T'
            else
              write(nTSU,'(a,8a14)') '#    Time','u','v','P'
            end if
          end if

        end do

        write(STDLOG,*)


      case(1)   !---> Output time-series (Normal case)

C-----: For all time-series monitor points
        do nts=1,nTSPoints
C-------: Gather default grid function values
          dTSv(nts,1) = u(iTS(nts),jTS(nts))
          dTSv(nts,2) = v(iTS(nts),jTS(nts))
          dTSv(nts,3) = p(iTS(nts),jTS(nts))

C-------: Thermal Energy (laminar)
          if(nThermEn.eq.1) dTSv(nts,4) = t(iTS(nts),jTS(nts))

C-------: Small-scale default variables
          if(nSmallScl.eq.1) then
            dTSv(nts,nLamV+1) = us(iTS(nts),jTS(nts))
            dTSv(nts,nLamV+2) = vs(iTS(nts),jTS(nts))
            dTSv(nts,nLamV+3) = ps(iTS(nts),jTS(nts))

C---------: Thermal Energy (small-scale)
            if(nThermEn.eq.1) 
     .        dTSv(nts,nLamV+4) = ts(iTS(nts),jTS(nts))

          end if

C-------: and write them to file
          write(nUnitMask+nts,'(17(e14.6))') dtime, 
     .      (dTSv(nts,j),j=1,nVars)

        end do


      case(2)   !---> Close files

C---: Close files
        write(STDLOG,*)
        write(STDLOG,'(a)') '* Closing time-series files.'

        do nts=1, nTSPoints
          close(nUnitMask+nts)
        end do

      case default
        print*, 'Wrong nCase flag passed to SaveTimeSrs'
        stop

      end select

      return
      end


*----------------------------------------------------------------------*
*                             lCheckFile                               *
*----------------------------------------------------------------------*

      LOGICAL function lCheckFile(fname,nError)

      implicit none

      character*(*) fname
      INTEGER nStrLen
      INTEGER nError

      logical*4 lex   !---> Note, this is required for inquire
                      !---> 'exist' clause

C---: Default value
      lCheckFile = .true.

      inquire(file=fname,iostat=nError,exist=lex)

      if(.not.lex.or.nError.ne.0) then
        lCheckFile = .false.
        write(STDOUT,'(3a)')
     .    '-> Error: Could not open file [',fname(1:nStrLen(fname)),
     .    '].'

C***Note:  Log file needs to be declared before using it at this point.
c        write(STDLOG,'(3a)')
c     .    '-> Error: Could not open file [',fname(1:nStrLen(fname)),
c     .    '].'
c        return

      end if

      return
      end


*----------------------------------------------------------------------*
*                             SaveTmAvgP3D                             *
*----------------------------------------------------------------------*

      subroutine SaveTmAvgP3D(nx, ny, nts, nForm, tAvgFile,
     .                       ubar, vbar, tbar, pbar,
     .                       upb, vpb, tpb,
     .                       upupb, vpvpb, upvpb, uptpb, vptpb,
     .                       tke, diss, dtdyb)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, nts, nForm
      character*(*) tAvgFile
      character*(mfnmlgth) tAvgQqq, tAvgNam

      REAL ubar(0:mnx,0:mny) , vbar(0:mnx,0:mny) ,
     .     tbar(0:mnx,0:mny) , pbar(0:mnx,0:mny) ,
     .     upb(0:mnx,0:mny)  , vpb(0:mnx,0:mny)  , tpb(0:mnx,0:mny) ,
     .     upupb(0:mnx,0:mny), vpvpb(0:mnx,0:mny), upvpb(0:mnx,0:mny),
     .     uptpb(0:mnx,0:mny), vptpb(0:mnx,0:mny),
     .     tke(0:mnx,0:mny)  , diss(0:mnx,0:mny) , dtdyb(0:mnx,0:mny)

      INTEGER nVars
      parameter ( nVars =   16 )  !--> Number of variables

C---: External procedures
      INTEGER nStrLen

C---: Local variables
      INTEGER i, j

C---: Single precision constants
      FLOAT   sZero
      parameter ( sZero = 0.00e+00 )

C---: Construct filenames
      tAvgQqq = tAvgFile(1:nStrLen(tAvgFile))//'.qqq'
      tAvgNam = tAvgFile(1:nStrLen(tAvgFile))//'.nam'

      select case(nForm)

      case (FT_FORMATTED)
C-----: Open time-average output file (formatted)
        open(22, file=tAvgQqq, status='unknown')
        write(22,*) nx, ny, nVars
        write(22,*) ((sngl( ubar(i,j)),  i=1,nx), j=1,ny),
     .              ((sngl( vbar(i,j)),  i=1,nx), j=1,ny),
     .              ((      sZero,       i=1,nx), j=1,ny),
     .              ((sngl( tbar(i,j)),  i=1,nx), j=1,ny),
     .              ((sngl( pbar(i,j)),  i=1,nx), j=1,ny),
     .              ((sngl( upb(i,j)  ), i=1,nx), j=1,ny),
     .              ((sngl( vpb(i,j)  ), i=1,nx), j=1,ny),
     .              ((sngl( tpb(i,j)  ), i=1,nx), j=1,ny),
     .              ((sngl( upupb(i,j)), i=1,nx), j=1,ny),
     .              ((sngl( vpvpb(i,j)), i=1,nx), j=1,ny),
     .              ((sngl( upvpb(i,j)), i=1,nx), j=1,ny),
     .              ((sngl( uptpb(i,j)), i=1,nx), j=1,ny),
     .              ((sngl( vptpb(i,j)), i=1,nx), j=1,ny),
     .              ((sngl( tke(i,j)  ), i=1,nx), j=1,ny),
     .              ((sngl( diss(i,j) ), i=1,nx), j=1,ny),
     .              ((sngl( dtdyb(i,j)), i=1,nx), j=1,ny)
        close(22)

      case (FT_UNFORMATTED)
C-----: Open time-average output file (unformatted)
        open(22, file=tAvgQqq, form='unformatted', status='unknown')
        write(22) nx, ny, nVars
        write(22) ((sngl( ubar(i,j)),  i=1,nx), j=1,ny),
     .            ((sngl( vbar(i,j)),  i=1,nx), j=1,ny),
     .            ((      sZero,       i=1,nx), j=1,ny),
     .            ((sngl( tbar(i,j)),  i=1,nx), j=1,ny),
     .            ((sngl( pbar(i,j)),  i=1,nx), j=1,ny),
     .            ((sngl( upb(i,j)  ), i=1,nx), j=1,ny),
     .            ((sngl( vpb(i,j)  ), i=1,nx), j=1,ny),
     .            ((sngl( tpb(i,j)  ), i=1,nx), j=1,ny),
     .            ((sngl( upupb(i,j)), i=1,nx), j=1,ny),
     .            ((sngl( vpvpb(i,j)), i=1,nx), j=1,ny),
     .            ((sngl( upvpb(i,j)), i=1,nx), j=1,ny),
     .            ((sngl( uptpb(i,j)), i=1,nx), j=1,ny),
     .            ((sngl( vptpb(i,j)), i=1,nx), j=1,ny),
     .            ((sngl( tke(i,j)  ), i=1,nx), j=1,ny),
     .            ((sngl( diss(i,j) ), i=1,nx), j=1,ny),
     .            ((sngl( dtdyb(i,j)), i=1,nx), j=1,ny)
        close(22)

      case default
        write(STDERR,'(a)')
     .    '* Wrong nForm flag passed to SaveTmAvgP3D'
        stop
      end select

C-----: Open time-average function names file
        open(23, file=tAvgNam, status='unknown')

C-----: Read number of grid points and time-steps
        write(23,*) 'Avgd. U-vel ; Avgd. Velocity'
        write(23,*) 'Avgd. V-vel'
        write(23,*) 'Avgd. W-vel (null)'
        write(23,*) 'Avgd. T'
        write(23,*) 'Avgd. P'
        write(23,*) 'Avgd. up'
        write(23,*) 'Avgd. vp'
        write(23,*) 'Avgd. tp'
        write(23,*) 'Avgd. up * up'
        write(23,*) 'Avgd. vp * vp'
        write(23,*) 'Avgd. up * vp'
        write(23,*) 'Avgd. up * tp'
        write(23,*) 'Avgd. vp * tp'
        write(23,*) 'Avgd. Turb. Kinetic En.'
        write(23,*) 'Avgd. Dissipation Rate '
        write(23,*) 'Avgd. dT*/dy*'

        close(23)

      return
      end


*----------------------------------------------------------------------*
*                           SaveStdVarsP3D                             *
*----------------------------------------------------------------------*
*  nSaveGrid flag:
*     0 -> Don't save grid
*     1 -> Save grid as well
*
      subroutine SaveStdVarsP3D(nx, ny, nForm, outPrefix,
     .                         nThermEn, nSmallScl,
     .                         lSaveGrid, lSaveNamFl,
     .                         gx, gy,
     .                         u, v, p, t, us, vs, ps, ts)

      implicit none

      include "config.f"

C---: Subroutine arguments
      INTEGER nx, ny, nForm
      character*(*) outPrefix
      INTEGER nThermEn, nSmallScl
      LOGICAL lSaveNamFl, lSaveGrid
      REAL    gx(0:mnx,0:mny), gy(0:mnx,0:mny)
      REAL    u(0:mnx,0:mny),  v(0:mnx,0:mny), 
     .        p(0:mnx,0:mny),  t(0:mnx,0:mny),
     .        us(0:mnx,0:mny), vs(0:mnx,0:mny),
     .        ps(0:mnx,0:mny), ts(0:mnx,0:mny)

C---: External procedures
      INTEGER nStrLen
 
C---: Local variables
      INTEGER i, j, l, nLamV, nVars
      character*(mfnmlgth) grdFile, solFile, namFile
      FLOAT f(0:mnx,0:mny,10)

C---: Double precision constants
      REAL    dZero
      parameter ( dZero = 0.00d+00 )


C---: Construct filenames
      solFile = outPrefix(1:nStrLen(outPrefix))//'.qqq'
      grdFile = outPrefix(1:nStrLen(outPrefix))//'.xyz'
      namFile = outPrefix(1:nStrLen(outPrefix))//'.nam'

C---: Save grid to file if requested
      if(lSaveGrid) call SaveGrid2DP3D(nx, ny, nForm, grdFile, gx, gy)

C---: Supress excessively small numbers (for reading into FieldView)
C---: Select variables for output according to flags
C---: Default number of variables
      nLamV = 4    !---> Laminar variables
      nVars = 4    !---> Total number of variables

      if(nThermEn.eq.1)  then
        nVars = nVars + 1    !---> +T
        nLamV = nLamV + 1
      end if

      if(nSmallScl.eq.1) nVars = 2 * nVars    !---> + Small-scale vars

C---: Default variables
      do j=1,ny
        do i=1,nx
          f(i,j,1) = p(i,j)
          f(i,j,2) = u(i,j)
          f(i,j,3) = v(i,j)
          f(i,j,4) = dZero      !---> 2D solution files
        end do
      end do

C---: Thermal Energy (laminar)
      if(nThermEn.eq.1) then
        do j=1,ny
          do i=1,nx
            f(i,j,5) = t(i,j)
          end do
        end do
      end if

C---: Small-scale variables
      if(nSmallScl.eq.1) then
        do j=1,ny
          do i=1,nx
            f(i,j,nLamV + 1) = ps(i,j)
            f(i,j,nLamV + 2) = us(i,j)
            f(i,j,nLamV + 3) = vs(i,j)
            f(i,j,nLamV + 4) = dZero      !---> 2D solution files
          end do
        end do

C-----: Thermal Energy (small-scale)
        if(nThermEn.eq.1) then
          do j=1,ny
            do i=1,nx
              f(i,j,nLamV + 5) = ts(i,j)
            end do
          end do
        end if
      end if


C---: Two-dimensional Plot3D format (nx, ny, nVars)
      select case(nForm)

      case (FT_FORMATTED)
        open(15, file=solFile, status='unknown')
        write(15,*) nx, ny, nVars
C-----: Write to file using default single prec. format to save space
        write(15,*) (((f(i,j,l), i=1,nx), j=1,ny), l=1,nVars)
        close(15)

      case (FT_UNFORMATTED)
        open(15, file=solFile, form='unformatted', status='unknown')
        write(15) nx, ny, nVars
        write(15) (((f(i,j,l), i=1,nx), j=1,ny), l=1,nVars)
        close(15)

      case default
        write(STDERR,'(a)')
     .    '* Wrong nForm flag passed to SaveStdVarsP3D'
        stop
      end select

C---: Save names file (for FieldView) if requested
      if(lSaveNamFl) then
        open(16,file=namFile,status='unknown')
C-------: Default variable names
          write(16,*) 'Complete Pressure'
          write(16,*) 'Complete U ; Complete Velocity'
          write(16,*) 'Complete V'
          write(16,*) 'Complete W'

C-------: If solving thermal energy 
          if(nThermEn.eq.1) write(16,*) 'Complete Temperature'

C-------: If solving small-scale variables
          if(nSmallScl.eq.1) then
            write(16,*) 'Small-Scale Pressure'
            write(16,*) 'Small-Scale U ; Small-Scale Velocity'
            write(16,*) 'Small-Scale V'
            write(16,*) 'Small-Scale W'

C---------: If solving thermal energy 
            if(nThermEn.eq.1) write(16,*) 'Small-Scale Temperature'

          end if

        close(16)

      end if

      return
      end


*----------------------------------------------------------------------*
*                            SaveStdVarsFV                             *
*----------------------------------------------------------------------*
*
* Save standard variables (complete and small-scale pressure, velocity
* and temperature) in FieldView 5.1 format.
*
      subroutine SaveStdVarsFV(nx, ny, outPrefix, npostp,
     .                         x, y, u, v, p, t)

      implicit none

C---: Configurable parameters
      include "config.f"

C---: Arguments
      INTEGER nx, ny, npostp

      character*(*) outPrefix

      REAL    x(0:mnx,0:mny), y(0:mnx,0:mny)
      REAL    u(0:mnx,0:mny), v(0:mnx,0:mny),
     .        p(0:mnx,0:mny), t(0:mnx,0:mny)

c      real*4 func(0:mnx,0:mny,10)

C---: External procedures
      INTEGER nstrlen

C---: Local variables
      INTEGER i, j, k, ne1, ne2, ne3, ne4, ne5, ne6, ne7, ne8
      character*(mfnmlgth) solFile

C---: Double precision constants
      REAL    dZero
      parameter ( dZero = 0.00d+00 )


C---: Construct filenames
      solFile = outPrefix(1:nStrLen(outPrefix))//'.fv'

C---: Chech for excessively small numbers (needs revision)
C????

C*** Temporary:  (yuk!)
 205  format(5(e14.6))

C---: Open file
      open(11, file=solFile, status='unknown')

C---: Header
      write(11,'(a)') 'FIELDVIEW 1 1'
      write(11,'(a)') 'Boundary Table'
      write(11,'(a)') '0'

C---: Variable names
      write(11,'(a)') 'Variable Names'
      write(11,*) '5'
      write(11,*) 'Pressure'
      write(11,*) 'U-vel ; Velocity'
      write(11,*) 'V-vel'
      write(11,*) 'W-vel'
      write(11,*) 'Temperature'

C---: Geometry:  Nodes
      write(11,'(a)') 'Nodes'
      write(11,*) (2*nx*ny)
      do k=0, 1
        do j=1, ny
          do i=1, nx
            write(11,205) x(i,j),y(i,j),(0.1d0)*dble(k)
          end do
        end do
      end do

C---: Boundary conditions
      write(11,'(a)') 'Boundary Faces'
      write(11,*) '0'

C---: Geometry:  Cells
      write(11,'(a)') 'Elements'
      do j=1,ny-1
        do i=1,nx-1
          ne1 = (j-1)*nx+i
          ne2 = ne1+1
          ne3 = nx*ny+ne1
          ne4 = ne3+1
          ne5 = ne1+nx
          ne6 = ne5+1
          ne7 = nx*ny+ne5
          ne8 = ne7+1
          write(11,*) 2,1,ne1,ne2,ne3,ne4,ne5,ne6,ne7,ne8
        end do
      end do

C---: Grid functions
      write(11,'(a)') 'Variables'
      write(11,205) (((p(i,j), i=1,nx), j=1,ny), k=1,2),
     .              (((u(i,j), i=1,nx), j=1,ny), k=1,2),
     .              (((v(i,j), i=1,nx), j=1,ny), k=1,2),
     .              (((dZero , i=1,nx), j=1,ny), k=1,2),
     .              (((t(i,j), i=1,nx), j=1,ny), k=1,2)
      close(11)

      return
      end


*----------------------------------------------------------------------*
*                              SaveVarP3D                              *
*----------------------------------------------------------------------*
*
* Save single grid function to file using 2D Plot3D format
*
      subroutine SaveVarP3D(nx, ny, fnPrefix, varName, var)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny
      character*(*) fnPrefix
      character*(*) varName

      REAL    var(0:mnx,0:mny)

C---: External procedures
      INTEGER nStrLen

C---: Local variables
      INTEGER i, j
      character*(mfnmlgth) qqqFile, namFile

C---: Construct filenames
      qqqFile = fnPrefix(1:nStrLen(fnPrefix))//'.qqq'
      namFile = fnPrefix(1:nStrLen(fnPrefix))//'.nam'

C---: Open function file
      open(12, file=qqqFile, status='unknown')

C---: Write indeces and number of variables
      write(12,*) nx, ny, 1

C---: Output function values in single precision format
      write(12,*) ((sngl(var(i,j)), i=1,nx), j=1,ny)

      close(12)

C---: Open names file
      open(13, file=namFile, status='unknown')

C---: Write function name
      write(13,*) varName

      close(13)

      return
      end

*-----------------------|---|---|---V---|---|---|----------------------*
