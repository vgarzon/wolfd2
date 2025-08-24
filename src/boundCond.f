*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                              boundCond.F                             *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* SetUpBCs      - Set up regions, region types, boundary types and
*               - boundary condition values
* InitialCond   - Initial conditions
* VelBoundCond  - Boundary conditions and image points for velocity 
*                 components
* PresBoundCond - Boundary conditions (image points) for pressure
* TempBoundCond - Boundary conditions (image points) for temperature
* SmlSclBC      - BCs (using image pts.) for small-scale variables
* PorRegCnst    - Porous region dimensionless constants
* NDHGenCoef    -
* VelOutflowBCs -
*
*----------------------------------------------------------------------*

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"


*----------------------------------------------------------------------*
*                               SetUpBCs                               *
*----------------------------------------------------------------------*

      subroutine SetUpBCs(nx, ny, bndryFile,
     .                    nthermen, neqstate, nsmallscl, ntraject,
     .                    nReg, nRegBrd,
     .                    nRegType, nTRgType, nMomBdTp, nTemBdTp,
     .                    dlref, uref, tref, tmax,
     .                    densref, rconst, dconduct,
     .                    re, pe,
     .                    dPRporos, dPRporc1, dPRporc2,
     .                    dTRgVal, dHGSTval,
     .                    dBCVal)


      implicit none

      include "config.f"

C---: Arguments
C---: Number of grid points in each direction
      INTEGER nx, ny

C---: Boundary conditions file
      character*(*) bndryFile

C---: Flags
      INTEGER nthermen, neqstate, nsmallscl, ntraject

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

      REAL    dlref, uref, tref, tmax,
     .        densref, rconst, dconduct,
     .        re, pe

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,4)

C---: Boundary type (temporary)
      INTEGER nTemBdTp(mgri,mgrj,4)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL    dBCVal(mgri,mgrj,4,4) 

C---: Porous region dimensionless coefficients and local porosity
      REAL    dPRporos(mgri,mgrj), dPRporc1(mgri,mgrj), 
     .        dPRporc2(mgri,mgrj) 

C---: Temperature region value
      REAL    dTRgVal(mgri,mgrj)

C---: Dimensionless heat generation source term values
      REAL    dHGSTval(mgri,mgrj)


C---: External procedures
      LOGICAL lReadBCFile


C---: Local variables
      INTEGER i, l, ireg, jreg

C---: Short temporary strings
      character*6 tmpStr(4)

C---: Porous region coefficient arrays
      INTEGER  nPRPermEq(mgri,mgrj)
      REAL     dPRDimCoef(mgri,mgrj,3)

C---: Parse boundary conditions file
      if (.not.lReadBCFile(nx, ny, bndryFile,
     .                     nthermen, neqstate, nsmallscl, ntraject,
     .                     dlref, uref, tref, tmax, densref,
     .                     nReg, nRegBrd,
     .                     nRegType, nMomBdTp, nTRgType, nTemBdTp,
     .                     nPRPermEq, dPRDimCoef,
     .                     dTRgVal, dBCVal)) then
        print*, 'Error parsing boundary file.  Aborting.'
        stop
      end if


C---: Simple consistency check of interior border definitions
      do ireg=2, nReg(_I_)
C-----: Negative numbers and ends
        if(nRegBrd(ireg,1,WEST).lt.2.or.
     .     nRegBrd(ireg,1,WEST).gt.nx-2) then
          print*, 'Bad iRegW index: ', nRegBrd(ireg,1,WEST), 
     .     '. Must be >= 2 and <= nx-2.'
          stop
        end if
      end do

      do ireg=2, nReg(_I_)-1
C-----: Not sequentially increasing
        if((nRegBrd(ireg+1,1,WEST)-nRegBrd(ireg,1,WEST)).lt.2) then
          print*, 'Bad iRegW index: ', nRegBrd(ireg+1,1,WEST) 
          print*, 'Indeces must be sequentially increasing and at ',
     .      'least 2 indeces apart'
          stop
        end if
      end do

      do jreg=2, nReg(_J_)
C-----: Negative numbers and ends
        if(nRegBrd(1,jreg,SOUTH).lt.2.or.
     .     nRegBrd(1,jreg,SOUTH).gt.ny-2) then
          print*, 'Bad jRegS index: ', nRegBrd(1,jreg,SOUTH), 
     .      '. Must be >= 2 and <= ny-2'
          stop
        end if
      end do

      do jreg=2, nReg(_J_)-1
C-----: Not sequentially increasing
        if((nRegBrd(1,jreg+1,SOUTH)-nRegBrd(1,jreg,SOUTH)).lt.2) then
          print*, 'Bad jRegS index: ', nRegBrd(1,jreg,SOUTH)
          print*, 'Indeces must be sequentially increasing and at ',
     .      'least 2 indeces apart'
          stop
        end if
      end do


C-----: Assign border locations to the other regions
C-----: Ends
        do jreg=1,nReg(_J_)
          nRegBrd(1,jreg,WEST)     = 1
          nRegBrd(nReg(_I_),jreg,EAST) = nx
        end do
        do ireg=1,nReg(_I_)
          nRegBrd(ireg,1,SOUTH)     = 1
          nRegBrd(ireg,nReg(_J_),NORTH) = ny
        end do
C-----: Interior
        do jreg=2,nReg(_J_)
          do ireg=2,nReg(_I_)
            nRegBrd(ireg,jreg,WEST)  = nRegBrd(ireg,1,WEST)
            nRegBrd(ireg,jreg,SOUTH) = nRegBrd(1,jreg,SOUTH)
          end do
        end do
C-----: Match W and E
        do jreg=1,nReg(_J_)
          do ireg=1,nReg(_I_)-1
            nRegBrd(ireg,jreg,EAST) = nRegBrd(ireg+1,1,WEST)
          end do
        end do

        do jreg=1,nReg(_J_)-1
          do ireg=1,nReg(_I_)
            nRegBrd(ireg,jreg,NORTH) = nRegBrd(1,jreg+1,SOUTH)
          end do
        end do

C-----: Simple boundary type consistency check
        do jreg=1, nReg(_J_)
          do ireg=1, nReg(_I_)-1
C---------: Inlets (also used to fix velocity at internal walls)
            if (nMomBdTp(ireg,jreg,EAST).eq.  BM_INLET) 
     .          nMomBdTp(ireg+1,jreg,WEST)  = BM_INLET
            if (nMomBdTp(ireg+1,jreg,WEST).eq.BM_INLET) 
     .          nMomBdTp(ireg,jreg,EAST)    = BM_INLET
C---------: Walls
            if (nMomBdTp(ireg,jreg,EAST).eq.  BM_WALL1) 
     .          nMomBdTp(ireg+1,jreg,WEST)  = BM_WALL1
            if (nMomBdTp(ireg+1,jreg,WEST).eq.BM_WALL1) 
     .          nMomBdTp(ireg,jreg,EAST)    = BM_WALL1
         end do
        end do
        do jreg=1, nReg(_J_)-1
          do ireg=1, nReg(_I_)
C---------: Inlets (also used to fix velocity at internal walls)
            if (nMomBdTp(ireg,jreg+1,SOUTH).eq.BM_INLET) 
     .          nMomBdTp(ireg,jreg,NORTH)    = BM_INLET
            if (nMomBdTp(ireg,jreg,NORTH).eq.  BM_INLET) 
     .          nMomBdTp(ireg,jreg+1,SOUTH)  = BM_INLET
C---------: Walls
            if (nMomBdTp(ireg,jreg+1,SOUTH).eq.BM_WALL1) 
     .          nMomBdTp(ireg,jreg,NORTH)    = BM_WALL1
            if (nMomBdTp(ireg,jreg,NORTH).eq.  BM_WALL1) 
     .          nMomBdTp(ireg,jreg+1,SOUTH)  = BM_WALL1
          end do
        end do


C-----: Simple boundary type consistency check -Temperature-
        do jreg=1, nReg(_J_)
          do ireg=1, nReg(_I_)-1
C---------: Fixed temperature
            if (nTemBdTp(ireg,jreg,EAST).eq.  BT_TEMPER)
     .          nTemBdTp(ireg+1,jreg,WEST)  = BT_TEMPER
            if (nTemBdTp(ireg+1,jreg,WEST).eq.BT_TEMPER)
     .          nTemBdTp(ireg,jreg,EAST)    = BT_TEMPER
C---------: Temperature flux
            if (nTemBdTp(ireg,jreg,EAST).eq.  BT_HTFLUX)
     .          nTemBdTp(ireg+1,jreg,WEST)  = BT_HTFLUX
            if (nTemBdTp(ireg+1,jreg,WEST).eq.BT_HTFLUX)
     .          nTemBdTp(ireg,jreg,EAST)    = BT_HTFLUX
          end do
        end do

        do jreg=1, nReg(_J_)-1
          do ireg=1, nReg(_I_)
C---------: Fixed temperature
            if (nTemBdTp(ireg,jreg+1,SOUTH).eq.BT_TEMPER)
     .          nTemBdTp(ireg,jreg,NORTH)    = BT_TEMPER
            if (nTemBdTp(ireg,jreg,NORTH).eq.  BT_TEMPER)
     .          nTemBdTp(ireg,jreg+1,SOUTH)    = BT_TEMPER
C---------: Temperature flux
            if (nTemBdTp(ireg,jreg+1,SOUTH).eq.BT_HTFLUX)
     .          nTemBdTp(ireg,jreg,NORTH)    = BT_HTFLUX
            if (nTemBdTp(ireg,jreg,NORTH).eq.  BT_HTFLUX)
     .          nTemBdTp(ireg,jreg+1,SOUTH)  = BT_HTFLUX
          end do
        end do


C-----: End thermal energy equation part.


C---: Check if any region was declared as a porous medium and
C---: compute the dimensionless porosity constant for that region.
       call PorRegCnst(nReg, nRegType,
     .                 nPRPermEq,
     .                 dlref, re,     
     .                 dPRDimCoef,
     .                 dPRporos, dPRporc1, dPRporc2)


C---: Check to see if any region was declared as a heat-generation
C---: block and compute the dimensionless heat generation src term.
C---: (Only if thermal energy equation flag is active)
      if (nthermen.eq.1)
     .  call NDHGenCoef(nReg, nTRgType,
     .                  dlref, tref, tmax, dconduct, pe,
     .                  dTRgVal,
     .                  dHGSTval)


C---: Write boundary condition information to the log file

      write(STDLOG,'(a)') '* Boundary conditions:'
      write(STDLOG,'(2x,a,2i4)') '- Number of regions:', (nReg(i),i=1,2)
      write(STDLOG,'(2x,a)') '- Region indeces:'
      write(STDLOG,'(2x,6(2x,a5))') 'ireg', 'jreg',
     .  'iRegW', 'iRegE', 'iRegS', 'iRegN'     
      do jreg=1, nReg(_J_)
        do ireg=1, nReg(_I_)
          write(STDLOG,'(1x,6(3x,i4))') ireg, jreg, 
     .     (nRegBrd(ireg,jreg,i),i=1,4)
        end do
      end do

      write(STDLOG,*) 
      write(STDLOG,'(2x,a)') '- Momentum region types:'
      write(STDLOG,'(2x,3(3x,a4))') 'ireg', 'jreg', 'Type'    
      do jreg=1, nReg(_J_)
        do ireg=1, nReg(_I_)

          select case(nRegType(ireg,jreg))
          case (RM_BLOCKG)
            tmpStr(1) = 'BLOCKG'
          case (RM_INTERN)
            tmpStr(1) = 'INTERN'
          case (RM_POROUS)
            tmpStr(1) = 'POROUS'
          case default
            tmpStr(1) = 'UNKNWN'
          end select

c          write(STDLOG,'(1x,3(3x,i4))') ireg, jreg, nRegType(ireg,jreg)
          write(STDLOG,'(1x,2(3x,i4),3x,a6)') ireg, jreg, tmpStr(1)

        end do 
      end do

      write(STDLOG,*)
      write(STDLOG,'(2x,a)') '- Momentum boundary types:'
c      write(STDLOG,'(2x,6(3x,a4))') 'ireg', 'jreg',
c     .  'nMBW', 'nMBE', 'nMBS', 'nMBN'
      write(STDLOG,'(2x,2(3x,a4),3x,4(a5,4x))') 'ireg', 'jreg',
     .  'West ', 'East ', 'South', 'North'
      do jreg=1, nReg(_J_)
        do ireg=1, nReg(_I_)

          do l=WEST, NORTH
            select case(nMomBdTp(ireg,jreg,l))
            case(BM_INTERN)
              tmpStr(l) = 'INTERN'
            case(BM_WALL1)
              tmpStr(l) = 'WALL1'
            case(BM_WALL2)
              tmpStr(l) = 'WALL2'
            case(BM_INLET)
              tmpStr(l) = 'INLET'
            case(BM_OUTLT1)
              tmpStr(l) = 'OUTLT1'
            case(BM_OUTLT2)
              tmpStr(l) = 'OUTLT2'
            case default
              tmpStr(l) = 'UNKNWN'
            end select
          end do

c          write(STDLOG,'(1x,6(3x,i4))') ireg, jreg,
c     .      (nMomBdTp(ireg,jreg,l), l=WEST,NORTH)

          write(STDLOG,'(1x,2(3x,i4),4(3x,a6))') ireg, jreg,
     .      (tmpStr(l), l=WEST,NORTH)

        end do
      end do

      if (nthermen.eq.1) then
        write(STDLOG,*)
        write(STDLOG,'(2x,a)') '- Temperature region types:'
        write(STDLOG,'(2x,3(3x,a4))') 'ireg', 'jreg', 'Type'    
        do jreg=1, nReg(_J_)
          do ireg=1, nReg(_I_)

            select case(nTRgType(ireg,jreg))
            case (RT_NOSRCE)
              tmpStr(1) = 'NOSRCE'
            case (RT_HEATGN)
              tmpStr(1) = 'HEATGN'
            case (RT_TEMPER)
              tmpStr(1) = 'TEMPER'
            case default
              tmpStr(1) = 'UNKNWN'
            end select

c          write(STDLOG,'(1x,3(3x,i4))') ireg, jreg, nTRgType(ireg,jreg)
            write(STDLOG,'(1x,2(3x,i4),3x,a6)') ireg, jreg, tmpStr(1)

          end do
        end do

        write(STDLOG,*)
        write(STDLOG,'(2x,a)') '- Temperature boundary types:'
c        write(STDLOG,'(2x,6(3x,a4))') 'ireg', 'jreg',
c       .  'nTBW', 'nTBE', 'nTBS', 'nTBN'

        write(STDLOG,'(2x,2(3x,a4),3x,4(a5,4x))') 'ireg', 'jreg',
     .    'West ', 'East ', 'South', 'North'

        do jreg=1, nReg(_J_)
          do ireg=1, nReg(_I_)

            do l=WEST, NORTH
              select case(nTemBdTp(ireg,jreg,l))
              case(BT_INTERN)
                tmpStr(l) = 'INTERN'
              case(BT_TEMPER)
                tmpStr(l) = 'TEMPER'
              case(BT_HTFLUX)
                tmpStr(l) = 'HTFLUX'
              case default
                tmpStr(l) = 'UNKNWN'
              end select
            end do

c            write(STDLOG,'(1x,6(3x,i4))') ireg, jreg,
c     .        (nTemBdTp(ireg,jreg,l), l=WEST,NORTH)

            write(STDLOG,'(1x,2(3x,i4),4(3x,a6))') ireg, jreg,
     .        (tmpStr(l), l=WEST,NORTH)

          end do
        end do

      end if

      write(STDLOG,*)
      write(STDLOG,'(2x,a)') '- Boundary condition values:'
      write(STDLOG,'(2x,2(3x,a4),4x,a3,4(a10,2x))') 
     .  'ireg', 'jreg', 'loc',
     .  '-U-', '-V-', '-P-', '-T-'
      do jreg=1, nReg(_J_)
        do ireg=1, nReg(_I_)
          write(STDLOG,'(4(1x,2(3x,i4),6x,a,2x,4(e12.5)/))')
     .      ireg, jreg, 'W', (dBCVal(ireg,jreg,WEST,i),i=1,4), 
     .      ireg, jreg, 'E', (dBCVal(ireg,jreg,EAST,i),i=1,4),
     .      ireg, jreg, 'S', (dBCVal(ireg,jreg,SOUTH,i),i=1,4),
     .      ireg, jreg, 'N', (dBCVal(ireg,jreg,NORTH,i),i=1,4)
        end do
      end do

      write(STDLOG,*)
      write(STDLOG,'(2x,a)') 
     .  '- Porous region dimensionless coefficients:'
      write(STDLOG,'(2x,2(3x,a4),3(a11,1x))') 
     .  'ireg', 'jreg', 'Poros', 'PorC1', 'PorC2'
      do jreg=1, nReg(_J_)   !--- For all regions in J direction
        do ireg=1, nReg(_I_)   !--- For all regions in I direction
          write(STDLOG,'(1x,2(3x,i4),2x,3(e12.5))')
     .      ireg, jreg, dPRporos(ireg,jreg), 
     .      dPRporc1(ireg,jreg), dPRporc2(ireg,jreg)
        end do   !---: End J-region sweep
      end do   !---: End I-region sweep
      write(STDLOG,*)

      return
      end


*----------------------------------------------------------------------*
*                             InitialCond                              *
*----------------------------------------------------------------------*

      subroutine InitialCond(nx, ny, nrestart, resFile, ks, dtime, 
     .                       u, v, p, t, uss,vss,pss,tss)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, nrestart, ks

      character*(*) resFile

      REAL    dtime

      REAL    u(0:mnx,0:mny),   v(0:mnx,0:mny),
     .        p(0:mnx,0:mny),   t(0:mnx,0:mny),
     .        uss(0:mnx,0:mny), vss(0:mnx,0:mny),
     .        pss(0:mnx,0:mny), tss(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j

C---: Double precision constants
      REAL    dZero
      parameter ( dZero = 0.00d+00 )


C---: If restarting from file
      if(nrestart.ne.0) then

         call ReadRestart(nx, ny, resFile, ks, dtime,
     .                    u, v, p, t, uss, vss, pss, tss)
      else

        ks    = 0
        dtime = dZero

C-----: Cold start:  Initialize the variables with quiescent conditions
        do j=0, ny+1
          do i=0, nx+1
            u(i,j)   = dZero
            v(i,j)   = dZero
            p(i,j)   = dZero
            t(i,j)   = dZero
            uss(i,j) = dZero
            vss(i,j) = dZero
            pss(i,j) = dZero
            tss(i,j) = dZero
          end do
        end do

      end if

      return
      end


*----------------------------------------------------------------------*
*                             VelBoundCond                             *
*----------------------------------------------------------------------*

      subroutine VelBoundCond(nx, ny, nReg, nRegBrd,
     .                        nMomBdTp, dBCVal,
     .                        u, v)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,4) 

      REAL    u(0:mnx,0:mny), v(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN

C---: Double precision constants
      REAL    dZero, dTwo, dThree, dFour, dFive, dEight
      parameter ( dZero  = 0.00d+00 ,
     .            dTwo   = 2.00d+00 ,
     .            dThree = 3.00d+00 ,
     .            dFour  = 4.00d+00 ,
     .            dFive  = 5.00d+00 ,
     .            dEight = 8.00d+00 )


C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-------: W E S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,WEST))
          case(BM_INTERN)   !---> Internal (interface)
                            !---> Don't do anything

          case(BM_WALL1)   !---> No-slip wall
C---[ U ]
            do j=jS, jN
              u(iW,j) = dZero
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = dTwo*dBCVal(ireg,jreg,WEST,_V_) - v(iW+1,j)
            end do

          case(BM_WALL2)   !---> No-stress wall
C---[ U ]
            do j=jS, jN
              u(iW,j) = dZero
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = v(iW+1,j)
            end do

          case(BM_INLET)   !---> Inlet
C---[ U ]
            do j=jS, jN
              u(iW,j) = dBCVal(ireg,jreg,WEST,_U_)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = dTwo*dBCVal(ireg,jreg,WEST,_V_) - v(iW+1,j)
            end do

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do j=jS, jN
              u(iW-1,j) = dBCVal(ireg,jreg,WEST,_U_) + u(iW,j)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = -v(iW+1,j)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation
C---[ U ]
            do j=jS+1, jN
              u(iW,j) = u(iW+1,j) + v(iW+1,j) - v(iW+1,j-1)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = - v(iW,j-1) + dFive*(v(iW+1,j) - v(iW+1,j-1))
     .          + dEight*(u(iW+1,j) - u(iW,j))
            end do

          case default
            print*, 'Wrong nBdTypeW flag in region ',ireg,',',jreg
            stop
          end select


C-------: E A S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,EAST))
          case(BM_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BM_WALL1)   !---> No-slip wall
C---[ U ]
            do j=jS, jN
              u(iE,j) = dZero
            end do
C---[ V ]
            do j=jS+1, jN
              v(iE+1,j) = dTwo*dBCVal(ireg,jreg,EAST,_V_) - v(iE,j)
            end do

          case(BM_WALL2)   !---> No-stress wall
C---[ U ]
            do j=jS, jN
              u(iE,j) = dZero
            end do
C---[ V ]
            do j=jS+1, jN
              v(iE+1,j) = v(iE,j)
            end do

          case(BM_INLET)   !---> Inlet
C---[ U ]
            do j=jS, jN
              u(iE,j) = dBCVal(ireg,jreg,EAST,_U_)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iE+1,j) = dTwo*dBCVal(ireg,jreg,EAST,_V_) - v(iE,j)
            end do

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do j=jS, jN
              u(iE+1,j) = dBCVal(ireg,jreg,EAST,_U_) + u(iE,j)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iE+1,j) = -v(iE,j)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation

cC---[ U ]
c            do j=jS+1, jN
c              u(iE,j) = u(iE-1,j) - v(iE,j) + v(iE,j-1)
c            end do
cC---[ V ]
c            do j=jS+1, jN
c              v(iE+1,j) = v(iE+1,j-1) + dFive*(v(iE,j) - v(iE,j-1))
c     .          + dEight*(u(iE,j) - u(iE-1,j))
c            end do

C---[ U ]
            do j=jS+1, jN
              u(iE,j) = u(iE-1,j) - ( v(iE,j) - v(iE,j-1) )
            end do
C---[ V ]
            do j=jS+1, jN-1
              v(iE+1,j) = v(iE+1,j-1) + dThree*(v(iE,j-1) - v(iE,j))
     .          - dFour*(u(iE,j) - u(iE-1,j))
            end do


          case default
            print*, 'Wrong nBdTypeE flag in region ',ireg,',',jreg
            stop
          end select


C-------: S O U T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,SOUTH))
          case(BM_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BM_WALL1)   !---> No-slip wall
C---[ U ]
            do i=iW+1, iE
              u(i,jS) = dTwo*dBCVal(ireg,jreg,SOUTH,_U_) - u(i,jS+1)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS) = dZero
            end do

          case(BM_WALL2)   !---> No-stress wall
C---[ U ]
            do i=iW+1, iE
              u(i,jS) = u(i,jS+1)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS) = dZero
            end do

          case(BM_INLET)   !---> Inlet
C---[ U ]
            do i=iW+1, iE
              u(i,jS) = dTwo*dBCVal(ireg,jreg,SOUTH,_U_) - u(i,jS+1)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS) = dBCVal(ireg,jreg,SOUTH,_V_)
            end do

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do i=iW+1, iE
              u(i,jS) = -u(i,jS+1)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS) = dBCVal(ireg,jreg,SOUTH,_V_) + v(i,jS)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation

cC---[ V ]
c            do i=iW+1, iE
c              v(i,jS) =  v(i,jS+1) + u(i,jS+1) - u(i-1,jS+1)
c            end do
cC---[ U ]
c            do i=iW+1, iE
c              u(i,jS) = -u(i-1,jS) + dFive*(u(i,jS+1) - u(i-1,jS+1))
c     .          + dEight*(v(i,jS+1) - v(i,jS))
c            end do

C---[ V ]
            do i=iW+1, iE
              v(i,jS) =  v(i,jS+1) + ( u(i,jS+1) - u(i-1,jS+1) )
            end do
C---[ U ]
            do i=iW+1, iE-1
              u(i,jS) = u(i-1,jS) + dThree*(u(i-1,jS+1) - u(i,jS+1))
     .          - dFour*(v(i,jS+1) - v(i,jS))
            end do

          case default
            print*, 'Wrong nBdTypeS flag in region ',ireg,',',jreg
            stop
          end select


C-------: N O R T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,NORTH))
          case(BM_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BM_WALL1)   !---> No-slip wall
C---[ U ]
            do i=iW+1, iE
              u(i,jN+1) = dTwo*dBCVal(ireg,jreg,NORTH,_U_) - u(i,jN)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jN) = dZero
            end do

          case(BM_WALL2)   !---> No-stress wall
C---[ U ]
            do i=iW+1, iE
              u(i,jN+1) = u(i,jN)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jN) = dZero
            end do

          case(BM_INLET)   !---> Inlet
C---[ U ]
            do i=iW+1, iE
              u(i,jN+1) = dTwo*dBCVal(ireg,jreg,NORTH,_U_) - u(i,jN)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jN) = dBCVal(ireg,jreg,NORTH,_V_)
            end do

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do i=iW+1, iE
              u(i,jN+1) = -u(i,jN)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jN+1) = dBCVal(ireg,jreg,NORTH,_V_) + v(i,jN)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation

cC---[ V ]
c            do i=iW+1, iE
c              v(i,jN) = v(i,jN-1) - u(i,jN) + u(i-1,jN)
c            end do
cC---[ U ]
c            do i=iW+1, iE
c              u(i,jN+1) = u(i-1,jN+1) + dFive*(u(i,jN) - u(i-1,jN))
c     .          + dEight*(v(i,jN) - v(i,jN-1))
c            end do

C---[ V ]
            do i=iW, iE
              v(i,jN) = v(i,jN-1) - ( u(i,jN) - u(i-1,jN) )
            end do
C---[ U ]
            do i=iW+1, iE-1
              u(i,jN+1) = u(i-1,jN+1) + dThree*(u(i-1,jN) - u(i,jN))
     .          - dFour*(v(i,jN) - v(i,jN-1))
            end do


          case default
            print*, 'Wrong nBdTypeN flag in region ',ireg,',',jreg
            stop
          end select

        end do
      end do


      return
      end

*----------------------------------------------------------------------*
*                             PresBoundCond                            *
*----------------------------------------------------------------------*

      subroutine PresBoundCond(nx, ny,
     .                         nReg, nRegBrd,
     .                         nRegType, nMomBdTp, dBCVal,
     .                         p)

      implicit none

      include "config.f"

C---: Arguments
C---: Grid dimensions
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,4) 

      REAL    p(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN

C---: Double precision constants
      REAL    dZero, dTwo
      parameter ( dZero = 0.00d+00 ,
     .            dTwo  = 2.00d+00 )

C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-------: If current region is a blockage
          if(nRegType(ireg,jreg).eq.RM_BLOCKG) then

C---------: Internal points
            do j=jS+1,jN
              do i=iW+1,iE
                p(i,j) = dZero
              end do
            end do

C---------: Assume zero normal gradient on blockage walls
C---------: W E S T   Boundary :-------]
            do j=jS+1, jN
              p(iW+1,j) = dBCVal(ireg,jreg,WEST,_P_) + p(iW,j)
            end do

C---------: E A S T   Boundary :-------]
            do j=jS+1, jN
              p(iE,j) = dBCVal(ireg,jreg,EAST,_P_) + p(iE+1,j)
            end do

C---------: S O U T H   Boundary :-------]
            do i=iW+1, iE
              p(i,jS+1) = dBCVal(ireg,jreg,SOUTH,_P_) + p(i,jS)
            end do

C---------: N O R T H   Boundary :-------]
            do i=iW+1, iE
              p(i,jN) = dBCVal(ireg,jreg,NORTH,_P_) + p(i,jN+1)
            end do

C---------: Goto to next region
            cycle

          end if

C-------> If region is either flow-through or porous medium:

C-------: W E S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,WEST))
          case(BM_INTERN)

          case(BM_WALL1, BM_WALL2, BM_INLET)
            do j=jS+1, jN
              p(iW,j) = dBCVal(ireg,jreg,WEST,_P_) + p(iW+1,j)
            end do

          case(BM_OUTLT1, BM_OUTLT2)
            do j=jS+1, jN
              p(iW,j) = dTwo*dBCVal(ireg,jreg,WEST,_P_) - p(iW+1,j)
            end do

          case default
            print*, 'Wrong nBdTypeW flag in region ',ireg,',',jreg
            stop
          end select


C-------: E A S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,EAST))
          case(BM_INTERN)

          case(BM_WALL1, BM_WALL2, BM_INLET)
            do j=jS+1, jN
              p(iE+1,j) = dBCVal(ireg,jreg,EAST,_P_) + p(iE,j)
            end do

          case(BM_OUTLT1, BM_OUTLT2) 
            do j=jS+1, jN
              p(iE+1,j) = dTwo*dBCVal(ireg,jreg,EAST,_P_) - p(iE,j)
            end do

          case default
            print*, 'Wrong nBdTypeE flag in region ',ireg,',',jreg
            stop
          end select


C-------: S O U T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,SOUTH))
          case(BM_INTERN)

          case(BM_WALL1, BM_WALL2, BM_INLET)
            do i=iW+1, iE
              p(i,jS) = dBCVal(ireg,jreg,SOUTH,_P_) + p(i,jS+1)
            end do

          case(BM_OUTLT1, BM_OUTLT2)
            do i=iW+1, iE
              p(i,jS) = dTwo*dBCVal(ireg,jreg,SOUTH,_P_) - p(i,jS+1)
            end do

          case default
            print*, 'Wrong nBdTypeS flag in region ',ireg,',',jreg
            stop
          end select


C-------: N O R T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,NORTH))
          case(BM_INTERN)   !---> Internal (interface)
                      !---> Don't do anything

          case(BM_WALL1, BM_WALL2, BM_INLET)
            do i=iW+1, iE
              p(i,jN+1) = dBCVal(ireg,jreg,NORTH,_P_) + p(i,jN)
            end do

          case(BM_OUTLT1, BM_OUTLT2) 
            do i=iW+1, iE
              p(i,jN+1) = dTwo*dBCVal(ireg,jreg,NORTH,_P_) - p(i,jN)
            end do

          case default
            print*, 'Wrong nBdTypeN flag in region ',ireg,',',jreg
          stop
          end select

        end do
      end do

      return
      end

*----------------------------------------------------------------------*
*                             TempBoundCond                            *
*----------------------------------------------------------------------*

      subroutine TempBoundCond(nx, ny,
     .                         nReg, nRegBrd,
     .                         nTRgType, nTemBdTp,
     .                         dTRgVal, dBCVal,
     .                         t)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Temperature Boundary type
      INTEGER nTemBdTp(mgri,mgrj,*)

C---: Temperature region values
      REAL    dTRgVal(mgri,mgrj)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,*)

      REAL t(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN

C---: Double precision constants
      REAL    dZero, dTwo
      parameter ( dZero = 0.00d+00 ,
     .            dTwo  = 2.00d+00 )

C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)


C-------: If a fixed-temperature region
          if(nTRgType(ireg,jreg).eq.RT_TEMPER) then
            do j=jS+1,jN
              do i=iW+1,iE
                T(i,j) = dTRgVal(ireg,jreg)
              end do
            end do

C---------: W E S T   Boundary :-------]
            do j=jS+1, jN
              T(iW+1,j) = dTwo*dBCVal(ireg,jreg,WEST,_T_) - T(iW,j)
            end do

C---------: E A S T   Boundary :-------]
            do j=jS+1, jN
              T(iE,j) = dTwo*dBCVal(ireg,jreg,EAST,_T_) - T(iE+1,j)
            end do

C---------: S O U T H   Boundary :-------]
            do i=iW+1, iE
              T(i,jS+1) = dTwo*dBCVal(ireg,jreg,SOUTH,_T_) - T(i,jS)
            end do

C---------: N O R T H   Boundary :-------]
            do i=iW+1, iE
              T(i,jN) = dTwo*dBCVal(ireg,jreg,NORTH,_T_) - T(i,jN+1)
            end do

C---------: Next region
            cycle

          end if

C-------> If region is not of fixed-temperature type

C-------: W E S T   Boundary :-------]
          select case(nTemBdTp(ireg,jreg,WEST))
          case(BT_INTERN)   !---> Internal (interface)

          case(BT_TEMPER)   !---> Fixed temperature
            do j=jS+1, jN
              T(iW,j) = dTwo*dBCVal(ireg,jreg,WEST,_T_) - T(iW+1,j)
            end do

          case(BT_HTFLUX)   !---> Temperature flux
            do j=jS+1, jN
             T(iW,j) = dBCVal(ireg,jreg,WEST,_T_) + T(iW+1,j)
            end do

          case default
            print*, 'Wrong nTemBdTp W flag in region ',ireg,',',jreg
            stop
          end select

C-------: E A S T   Boundary :-------]
          select case(nTemBdTp(ireg,jreg,EAST))
          case(BT_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BT_TEMPER)   !---> Fixed temperature
            do j=jS+1, jN
              T(iE+1,j) = dTwo*dBCVal(ireg,jreg,EAST,_T_) - T(iE,j)
            end do

          case(BT_HTFLUX)   !---> Temperature flux
            do j=jS+1, jN
              T(iE+1,j) = dBCVal(ireg,jreg,EAST,_T_) + T(iE,j)
            end do

          case default
            print*, 'Wrong nTemBdTp E flag in region ',ireg,',',jreg
            stop
          end select


C-------: S O U T H   Boundary :-------]
          select case(nTemBdTp(ireg,jreg,SOUTH))
          case(BT_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BT_TEMPER)   !---> Fixed temperature
            do i=iW+1, iE
              T(i,jS) = dTwo*dBCVal(ireg,jreg,SOUTH,_T_) - T(i,jS+1)
            end do

          case(BT_HTFLUX)   !---> Temperature flux
            do i=iW+1, iE
              T(i,jS) = dBCVal(ireg,jreg,SOUTH,_T_) + T(i,jS+1)
            end do

          case default
            print*, 'Wrong nTemBdTp S flag in region ',ireg,',',jreg
            stop
          end select


C-------: N O R T H   Boundary :-------]
          select case(nTemBdTp(ireg,jreg,NORTH))
          case(BT_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BT_TEMPER)   !---> Fixed temperature
            do i=iW+1, iE
              T(i,jN+1) = dTwo*dBCVal(ireg,jreg,NORTH,_T_) - T(i,jN)
            end do

          case(BT_HTFLUX)   !---> Temperature flux
            do i=iW+1, iE
              T(i,jN+1) = dBCVal(ireg,jreg,NORTH,_T_) + T(i,jN)
            end do

          case default
            print*, 'Wrong nTemBdTp N flag in region ',ireg,',',jreg
            stop
          end select

        end do
      end do


      return
      end


*----------------------------------------------------------------------*
*                               SmlSclBC                               *
*----------------------------------------------------------------------*

      subroutine SmlSclBC(nx, ny,
     .                    nReg, nRegBrd,
     .                    nRegType, nMomBdTp, nTRgType, nTemBdTp,
     .                    dBCVal,
     .                    u, v, p, t)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Momentum Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

C---: Momentum Boundary type
      INTEGER nTemBdTp(mgri,mgrj,*)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,*)

      REAL    u(0:mnx,0:mny), v(0:mnx,0:mny), 
     .        p(0:mnx,0:mny), t(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN

C---: Double precision constants
      REAL    dZero, dTwo, dThree, dFour, dFive, dEight
      parameter ( dZero  = 0.00d+00 ,
     .            dTwo   = 2.00d+00 ,
     .            dThree = 3.00d+00 ,
     .            dFour  = 4.00d+00 ,
     .            dFive  = 5.00d+00 ,
     .            dEight = 8.00d+00 )


C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-------: W E S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,WEST))
          case(BM_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BM_WALL1, BM_INLET)   !---> No-slip wall, inlet
C---[ U ]
            do j=jS, jN
              u(iW,j) = dZero
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = - v(iW+1,j)
            end do

          case(BM_WALL2)   !---> No-stress wall
C---[ U ]
            do j=jS, jN
              u(iW,j) = dZero
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = v(iW+1,j)
            end do

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do j=jS, jN
              u(iW-1,j) = +u(iW,j)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = -v(iW+1,j)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation
C---[ U ]
            do j=jS, jN
              u(iW,j) = u(iW+1,j) - v(iW,j) + v(iW,j-1)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = - v(iW,j-1) + dFive*(v(iW+1,j) - v(iW+1,j-1))
     .          + dEight*(u(iW+1,j) - u(iW,j))
            end do

          case default
            print*, 'Wrong nBdTypeW flag in region ',ireg,',',jreg
            stop
          end select


C-------: E A S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,EAST))
          case(BM_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BM_WALL1, BM_INLET)   !---> No-slip wall, inlet
C---[ U ]
            do j=jS, jN
              u(iE,j) = dZero
            end do
C---[ V ]
            do j=jS+1, jN
              v(iE+1,j) = - v(iE,j)
            end do

          case(BM_WALL2)   !---> No-stress wall
C---[ U ]
            do j=jS, jN
              u(iE,j) = dZero
            end do
C---[ V ]
            do j=jS+1, jN
              v(iE+1,j) = v(iE,j)
            end do

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do j=jS, jN
              u(iE+1,j) = +u(iE+1,j)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iE+1,j) = -v(iE,j)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation

cC---[ U ]
c            do j=jS, jN
c              u(iE,j) = u(iE-1,j) - v(iE,j) + v(iE,j-1)
c            end do
cC---[ V ]
c            do j=jS+1, jN
c              v(iE+1,j) = v(iE+1,j-1) + dFive*(v(iE,j) - v(iE,j-1))
c     .          + dEight*(u(iE,j) - u(iE-1,j))
c            end do

C---[ U ]
            do j=jS+1, jN
              u(iE,j) = u(iE-1,j) - ( v(iE,j) - v(iE,j-1) )
            end do
C---[ V ]
            do j=jS+1, jN-1
              v(iE+1,j) = v(iE+1,j-1) + dThree*(v(iE,j-1) - v(iE,j))
     .          - dFour*(u(iE,j) - u(iE-1,j))
            end do

          case default
            print*, 'Wrong nBdTypeE flag in region ',ireg,',',jreg
            stop
          end select


C-------: S O U T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,SOUTH))
          case(BM_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BM_WALL1, BM_INLET)   !---> No-slip wall, inlet
C---[ U ]
            do i=iW+1, iE
              u(i,jS) = - u(i,jS+1)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS) = dZero
            end do

          case(BM_WALL2)   !---> No-stress
C---[ U ]
            do i=iW+1, iE
              u(i,jS) = u(i,jS+1)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS) = dZero
            end do

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do i=iW+1, iE
              u(i,jS) = -u(i,jS+1)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS-1) = +v(i,jS)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation
C---[ U ]
            do i=iW+1, iE
              u(i,jS-1) = u(i-1,jS-1) + dFive*(u(i,jS) - u(i-1,jS))
     .          + dEight*(v(i,jS) - v(i-1,jS-1))
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS-1) =  v(i,jS) - u(i,jS) - u(i-1,jS)
            end do

          case default
            print*, 'Wrong nBdTypeS flag in region ',ireg,',',jreg
            stop
          end select


C-------: N O R T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,NORTH))
          case(BM_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BM_WALL1, BM_INLET)   !---> No-slip wall, inlet
C---[ U ]
            do i=iW+1, iE
              u(i,jN+1) = -u(i,jN)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jN) = dZero
            end do

          case(BM_WALL2)   !---> No-stress wall
C---[ U ]
            do i=iW+1, iE
              u(i,jN+1) = u(i,jN)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jN) = dZero
            end do

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do i=iW+1, iE
              u(i,jN+1) = -u(i,jN)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jN+1) = +v(i,jN)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation

cC---[ U ]
c            do i=iW+1, iE
c              u(i,jN+1) = u(i-1,jN+1) + dFive*(u(i,jN) - u(i-1,jN))
c     .          + dEight*(v(i,jN) - v(i-1,jN-1))
c            end do
cC---[ V ]
c            do i=iW, iE
c              v(i,jN+1) = v(i,jN) - u(i,jN) - u(i-1,jN)
c            end do

C---[ V ]
            do i=iW, iE
              v(i,jN) = v(i,jN-1) - ( u(i,jN) - u(i-1,jN) )
            end do
C---[ U ]
            do i=iW+1, iE-1
              u(i,jN+1) = u(i-1,jN+1) + dThree*(u(i-1,jN) - u(i,jN))
     .          - dFour*(v(i,jN) - v(i,jN-1))
            end do

          case default
            print*, 'Wrong nBdTypeN flag in region ',ireg,',',jreg
            stop
          end select

        end do
      end do



C*** This needs revision:
C---: For now try same boundary conditions as large-scale variables
      call PresBoundCond(nx, ny,
     .                   nReg, nRegBrd,
     .                   nRegType, nMomBdTp, dBCVal,
     .                   p)


C--------------------- Small-scale Temperature BCs

C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)


C-------: If a fixed-temperature region
          if(nTRgType(ireg,jreg).eq.RT_TEMPER) then
            do j=jS+1,jN
              do i=iW+1,iE
                T(i,j) = dZero
              end do
            end do

C---------: W E S T   Boundary :-------]
            do j=jS+1, jN
              T(iW+1,j) =  - T(iW,j)
            end do

C---------: E A S T   Boundary :-------]
            do j=jS+1, jN
              T(iE,j) = - T(iE+1,j)
            end do

C---------: S O U T H   Boundary :-------]
            do i=iW+1, iE
              T(i,jS+1) = - T(i,jS)
            end do

C---------: N O R T H   Boundary :-------]
            do i=iW+1, iE
              T(i,jN) = - T(i,jN+1)
            end do

C---------: Next region
            cycle

          end if

C-------> If region is not of fixed-temperature type

C-------: W E S T   Boundary :-------]
          select case(nTemBdTp(ireg,jreg,WEST))
          case(BT_INTERN)   !---> Internal (interface)

          case(BT_TEMPER)   !---> Fixed temperature
            do j=jS+1, jN
              T(iW,j) = - T(iW+1,j)
            end do

          case(BT_HTFLUX)   !---> Temperature flux
            do j=jS+1, jN
             T(iW,j) = T(iW+1,j)
            end do

          case default
            print*, 'Wrong nTemBdTp W flag in region ',ireg,',',jreg
            stop
          end select

C-------: E A S T   Boundary :-------]
          select case(nTemBdTp(ireg,jreg,EAST))
          case(BT_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BT_TEMPER)   !---> Fixed temperature
            do j=jS+1, jN
              T(iE+1,j) = - T(iE,j)
            end do

          case(BT_HTFLUX)   !---> Temperature flux
            do j=jS+1, jN
              T(iE+1,j) = T(iE,j)
            end do

          case default
            print*, 'Wrong nTemBdTp E flag in region ',ireg,',',jreg
            stop
          end select


C-------: S O U T H   Boundary :-------]
          select case(nTemBdTp(ireg,jreg,SOUTH))
          case(BT_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BT_TEMPER)   !---> Fixed temperature
            do i=iW+1, iE
              T(i,jS) = - T(i,jS+1)
            end do

          case(BT_HTFLUX)   !---> Temperature flux
            do i=iW+1, iE
              T(i,jS) = T(i,jS+1)
            end do

          case default
            print*, 'Wrong nTemBdTp S flag in region ',ireg,',',jreg
            stop
          end select


C-------: N O R T H   Boundary :-------]
          select case(nTemBdTp(ireg,jreg,NORTH))
          case(BT_INTERN)   !---> Internal (interface)
                    !---> Don't do anything

          case(BT_TEMPER)   !---> Fixed temperature
            do i=iW+1, iE
              T(i,jN+1) = - T(i,jN)
            end do

          case(BT_HTFLUX)   !---> Temperature flux
            do i=iW+1, iE
              T(i,jN+1) =  T(i,jN)
            end do

          case default
            print*, 'Wrong nTemBdTp N flag in region ',ireg,',',jreg
            stop
          end select

        end do
      end do


      return
      end


*----------------------------------------------------------------------*
*                            VelOutflowBCs                             *
*----------------------------------------------------------------------*

      subroutine VelOutflowBCs(nx, ny,
     .                         nReg, nRegBrd, nMomBdTp,
     .                         dBCVal,
     .                         u, v)

      implicit none

      include "config.f"

C---: Arguments
	   INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,*)

      REAL    u(0:mnx,0:mny), v(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN

C---: Double precision constants
      REAL    dZero, dTwo, dThree, dFour, dFive, dEight
      parameter ( dZero  = 0.00d+00 ,
     .            dTwo   = 2.00d+00 ,
     .            dThree = 3.00d+00 ,
     .            dFour  = 4.00d+00 ,
     .            dFive  = 5.00d+00 ,
     .            dEight = 8.00d+00 )


C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-------: W E S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,WEST))

          case(BM_INTERN, BM_WALL1, BM_WALL2, BM_INLET)
          !---> Do nothing

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do j=jS, jN
              u(iW-1,j) = dBCVal(ireg,jreg,WEST,_U_) + u(iW,j)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = -v(iW+1,j)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation
C---[ U ]
            do j=jS+1, jN
              u(iW,j) = u(iW+1,j) - v(iW+1,j) + v(iW+1,j-1)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iW,j) = - v(iW,j-1) + dFive*(v(iW+1,j) - v(iW+1,j-1))
     .          + dEight*(u(iW+1,j) - u(iW,j))
            end do

          case default
            print*, 'Wrong nBdTypeW flag in region ',ireg,',',jreg
            stop
          end select


C-------: E A S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,EAST))

          case(BM_INTERN, BM_WALL1, BM_WALL2, BM_INLET)
          !---> Do nothing

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do j=jS, jN
              u(iE+1,j) = dBCVal(ireg,jreg,EAST,_U_) + u(iE,j)
            end do
C---[ V ]
            do j=jS+1, jN
              v(iE+1,j) = -v(iE,j)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation

cC---[ U ]
c            do j=jS+1, jN
c              u(iE,j) = u(iE-1,j) - v(iE,j) + v(iE,j-1)
c            end do
cC---[ V ]
c            do j=jS+1, jN
c              v(iE+1,j) = v(iE+1,j-1) + dFive*(v(iE,j) - v(iE,j-1))
c     .          + dEight*(u(iE,j) - u(iE-1,j))
c            end do

C---[ U ]
            do j=jS+1, jN
              u(iE,j) = u(iE-1,j) - ( v(iE,j) - v(iE,j-1) )
            end do
C---[ V ]
            do j=jS+1, jN-1
              v(iE+1,j) = v(iE+1,j-1) + dThree*(v(iE,j-1) - v(iE,j))
     .          - dFour*(u(iE,j) - u(iE-1,j))
            end do

          case default
            print*, 'Wrong nBdTypeE flag in region ',ireg,',',jreg
            stop
          end select


C-------: S O U T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,SOUTH))

          case(BM_INTERN, BM_WALL1, BM_WALL2, BM_INLET)
          !---> Do nothing

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do i=iW+1, iE
              u(i,jS) = -u(i,jS+1)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jS) = dBCVal(ireg,jreg,SOUTH,_V_) + v(i,jS)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation

cC---[ V ]
c            do i=iW+1, iE
c              v(i,jS) =  v(i,jS+1) - u(i,jS+1) + u(i-1,jS+1)
c            end do
cC---[ U ]
c            do i=iW+1, iE
c              u(i,jS) = - u(i-1,jS) + dFive*(u(i,jS+1) - u(i-1,jS+1))
c     .          + dEight*(v(i,jS+1) - v(i,jS))
c            end do

C---[ V ]
            do i=iW+1, iE
              v(i,jS) =  v(i,jS+1) + ( u(i,jS+1) - u(i-1,jS+1) )
            end do
C---[ U ]
            do i=iW+1, iE-1
              u(i,jS) = u(i-1,jS) + dThree*(u(i-1,jS+1) - u(i,jS+1))
     .          - dFour*(v(i,jS+1) - v(i,jS))
            end do

          case default
            print*, 'Wrong nBdTypeS flag in region ',ireg,',',jreg
            stop
          end select


C-------: N O R T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,NORTH))

          case(BM_INTERN, BM_WALL1, BM_WALL2, BM_INLET)
          !---> Do nothing

          case(BM_OUTLT1)   !---> Outlet: fully developed flow conditions
C---[ U ]
            do i=iW+1, iE
              u(i,jN+1) = -u(i,jN)
            end do
C---[ V ]
            do i=iW, iE
              v(i,jN+1) = dBCVal(ireg,jreg,NORTH,_V_) + v(i,jN)
            end do

          case(BM_OUTLT2)   !---> Outlet: mass conservation

cC---[ V ]
c            do i=iW+1, iE
c              v(i,jN) = v(i,jN-1) - u(i,jN) + u(i-1,jN)
c            end do
cC---[ U ]
c            do i=iW+1, iE
c              u(i,jN+1) = u(i-1,jN+1) + dFive*(u(i,jN) - u(i-1,jN))
c     .          + dEight*(v(i,jN) - v(i,jN-1))
c            end do

C---[ V ]
            do i=iW, iE
              v(i,jN) = v(i,jN-1) - ( u(i,jN) - u(i-1,jN) )
            end do
C---[ U ]
            do i=iW+1, iE-1
              u(i,jN+1) = u(i-1,jN+1) + dThree*(u(i-1,jN) - u(i,jN))
     .          - dFour*(v(i,jN) - v(i,jN-1))
            end do

          case default
            print*, 'Wrong nBdTypeN flag in region ',ireg,',',jreg
            stop
          end select

        end do
      end do

      return
      end


*----------------------------------------------------------------------*
*                               PorRegCnst                             *
*----------------------------------------------------------------------*
*
* Calculate dimensionless porosity coefficients in a given region
*
      subroutine PorRegCnst(nReg, nRegType,
     .                      nPRPermEq,
     .                      dlref, re,
     .                      dPRDimCoef,
     .                      dPRporos, dPRporc1, dPRporc2)

      implicit none

      include "config.f"

C---: Number of regions
      INTEGER  nReg(*)

      INTEGER  nRegType(mgri,mgrj)
      INTEGER  nPRPermEq(mgri,mgrj)

      REAL     dlref, re

      REAL     dPRDimCoef(mgri,mgrj,3)

C---: Porous region dimensionless coefficients and local porosity 
C---: value (temporary)
      REAL     dPRporos(mgri,mgrj), dPRporc1(mgri,mgrj), 
     .         dPRporc2(mgri,mgrj)


C---: Local variables
      INTEGER  ireg, jreg
      INTEGER  nPermEq
      REAL     dPorosity, dPermDia, dPorCf, dPermeab,
     .         por, dia, pors, pors2, dKozeny
      REAL     dpi4


C---: Pi/4
      dpi4 = dacos(-1.00d0)/4.00d0

      do jreg=1, nReg(_J_)   !--- For all regions in J direction
        do ireg=1, nReg(_I_)   !--- For all regions in I direction

          if(nRegType(ireg,jreg).eq.RM_POROUS) then  !---: If a porous region

C***: Temporary:
C***: The following is done only to allow for backward compatibility
C***: with previous subroutines (for debugging purposes)
            nPermEq   = nPRPermEq(ireg,jreg)
            dPorosity = dPRDimCoef(ireg,jreg,1)
            dPermDia  = dPRDimCoef(ireg,jreg,2)
            dPorcf    = dPRDimCoef(ireg,jreg,3)

C---: Check for zero in denominator of permeability formulas
            if((1.d0-dPorosity).lt.1.d-4) then
              print*,'Warning: (1-Porosity) < Tol, using Porosity=0.98'
              dPorosity = 0.98d0
            end if

C---: Calculate permeability from porosity
            por = dPorosity
            dia = dPermDia

            select case(nPermEq)
            case(1)   !---> Carman-Kozeny (packed bed of spheres)
              dPermeab = (dia**2*por**3)/(1.8d2*(1.00d0 - por)**2)

            case(2)   !---> Kaviany (over cylinders, numerical)
              dPermeab = 0.0606d0*dpi4*(por**5.10d0)*dia**2/
     .                   (1.00d0 - por)

            case(3)   !---> Happen & Brenner (over cylinders, perpendic)
C-----------: Kozeny constant
              pors = 1.00d0 - por
              pors2 = pors**2
              dKozeny = (2.00d0*por**3)/pors/(dlog(1.00d0/pors) 
     .          - (1.00d0 - pors2)/(1.00d0 + pors2))
              dPermeab = dia**2*por**3/(16.00d0*dkozeny*pors2)

            case default
              print*, 'Error: Wrong nPermEq flag passed to ReadParam.'
              stop
            end select

C---: Dimensionless porosity coefficients and local porosity
            dPRporos(ireg,jreg) = dPRDimCoef(ireg,jreg,1)
            dPRporc1(ireg,jreg) = (dporosity*dlref**2)/
     .                            (re*dpermeab)
            dPRporc2(ireg,jreg) 
     .        = (dlref*dPorCf*dPorosity)/dsqrt(dPermeab)

          else 

C---------: If not a porous region, set the porosity to 1.d0 and the 
C---------: dimensionless coefficients to 0.d0

            dPRporos(ireg,jreg) = 1.00d0
            dPRporc1(ireg,jreg) = 0.00d0
            dPRporc2(ireg,jreg) = 0.00d0

          end if

        end do   !---: End J-region sweep
      end do   !---: End I-region sweep

      return
      end


*----------------------------------------------------------------------*
*                               NDHGenCoef                             *
*----------------------------------------------------------------------*
*
* Calculate dimensionless heat generation coefficients for any 
* required given region
*
      subroutine NDHGenCoef(nReg, nTRgType,
     .                      dlref, tref, tmax, dconduct, pe,
     .                      dTRgVal,
     .                      dHGSTval)

      implicit none

      include "config.f"

C---: Number of regions
      INTEGER nReg(*)

C---: Region types: Thermal energy
      INTEGER   nTRgType(mgri,mgrj)

C---: Declare temperature-region value
      REAL   dlref, tref, tmax, dconduct, pe
      REAL   dTRgVal(mgri,mgrj)
      REAL   dHGSTval(mgri,mgrj)

C---: Local variables
      INTEGER ireg, jreg

      do jreg=1, nReg(_J_)   !--- For all regions in J direction
        do ireg=1, nReg(_I_)   !--- For all regions in I direction

          if(nTRgType(ireg,jreg).eq.RT_HEATGN) then  !---: If a heat gen. reg

            dHGSTval(ireg,jreg) = (dTRgVal(ireg,jreg)*dlref**2)/
     .        (pe*dconduct*(tmax-tref))

          end if

        end do   !---: End J-region sweep
      end do   !---: End I-region sweep

      return
      end



*-----------------------|---|---|---V---|---|---|----------------------*
