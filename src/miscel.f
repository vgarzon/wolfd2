*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                               miscel.F                               *
*
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* Filter       - Shuman filter
* Project      - Project aux-vel onto a divergence-free space
* DiffMaxNorm  - Maximum difference (in abs. value) between two arrays
* DMaxNorm     - Maximum component (in absolute value) of an array
* PTDAvg       - Average pressure (also t and d) at natural grid points
* VelAvg       - Average velocity components at natural grid points
* dLinInterp   -
* TAveraged    -
*
*----------------------------------------------------------------------*

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"


*----------------------------------------------------------------------*
*                                Filter                                *
*----------------------------------------------------------------------*
*
* ncomp -> Component (variable) flag:
*
      subroutine Filter(nx, ny, ncomp,
     .                  nReg, nRegBrd,
     .                  nRegType, nMomBdTp, nTRgType,
     .                  fp, qu)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, ncomp
      REAL    qu(0:mnx,0:mny)

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

      INTEGER nRegType(mgri,mgrj)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,4)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Filter parameter
      REAL fp 

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN
      REAL    qh(0:mnx,0:mny)

C---: Double precision constants
      REAL    dFour
      parameter ( dFour = 4.00d+00 )

C---: Load function to be filtered into a temporary array
      do j=0, ny+1
        do i=0, nx+1
          qh(i,j) = qu(i,j)
        end do
      end do

      select case(ncomp)
      case(_U_)   !------[ u: x-component of velocity ]------!

C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

          if(nRegType(ireg,jreg).eq.RM_BLOCKG) cycle  !---: If a blockage

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-------: Filter interior points
          do j=jS+1, jN
            do i=iW+1, iE-1
              qh(i,j) = (qu(i,j-1)+qu(i-1,j)+qu(i,j+1)+qu(i+1,j)
     .                +fp*qu(i,j))/( fp + dFour)
            end do
          end do

C-------: Filter boundaries

C-------: W E S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,WEST))

          case(BM_OUTLT1)   !---> Outlet: fully developed flow
            do j=jS+1, jN
              qh(iW,j) = (qu(iW,j-1)+qu(iW-1,j)+qu(iW,j+1)+qu(iW+1,j)
     .          +fp*qu(iW,j)) / (fp + dFour)
            end do

          case(BM_INTERN, BM_WALL1, BM_WALL2, BM_INLET, BM_OUTLT2)
          !---> Internal boundary, walls, inlet, outlets: Don't filter
          !---> Note that internal boundaries are filtered only at
          !---> East and North walls to avoid filtering twice.
 
          case default
            print*, 'Wrong nBdTypeW flag in region ',ireg,',',jreg
          end select

C-------: E A S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,EAST))

          case(BM_INTERN, BM_OUTLT1)   !---> Internal (interface)
                                       !---> Outlet: fully developed flow
            do j=jS+1, jN
              qh(iE,j) = (qu(iE,j-1)+qu(iE-1,j)+qu(iE,j+1)+qu(iE+1,j)
     .          +fp*qu(iE,j)) / (fp + dFour)
            end do

          case(BM_WALL1, BM_WALL2, BM_INLET, BM_OUTLT2)
          !---> Walls, inlet, outlets: Don't filter

          case default
            print*, 'Wrong nBdTypeE flag in region ',ireg,',',jreg
          end select

        end do
      end do

      case(_V_)   !-----[ v: y-component of velocity ]------!

C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_) 
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

          if(nRegType(ireg,jreg).eq.RM_BLOCKG) cycle  !---: If a blockage

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-------: Filter interior points
          do j=jS+1, jN-1
            do i=iW+1, iE
              qh(i,j) = (qu(i,j-1)+qu(i-1,j)+qu(i,j+1)+qu(i+1,j)
     .                +fp*qu(i,j)) / (fp + dFour)
            end do
          end do

C-------: Filter boundaries

C-------: S O U T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,SOUTH))

          case(BM_OUTLT1)   !---> Outlet: fully developed flow
            do i=iW+1, iE
              qh(i,jS) = (qu(i,jS-1)+qu(i-1,jS)+qu(i,jS+1)+qu(i+1,jS)
     .          +fp*qu(i,jS)) / (fp + dFour)
            end do

          case(BM_INTERN, BM_WALL1, BM_WALL2, BM_INLET, BM_OUTLT2)
          !---> Internal, walls, inlet, outlets: Don't filter
          !---> Note that internal boundaries are filtered only at
          !---> East and North walls to avoid filtering twice.

          case default
            print*, 'Wrong nBdTypeS flag in region ',ireg,',',jreg
          end select


C-------: N O R T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,NORTH))

          case(BM_INTERN, BM_OUTLT1)   !---> Internal (interface)
                                       !---> Outlet: fully developed flow
            do i=iW+1, iE
              qh(i,jN) = (qu(i,jN-1)+qu(i-1,jN)+qu(i,jN+1)+qu(i+1,jN)
     .          +fp*qu(i,jN)) / (fp + dFour)
            end do

          case(BM_WALL1, BM_WALL2, BM_INLET, BM_OUTLT2)
          !---> Walls, inlet, outlets: Don't filter

          case default
            print*, 'Wrong nBdTypeN flag in region ',ireg,',',jreg
          end select

        end do
      end do

      case(_T_)   !---> T: Temperature, density or other scalars

C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

          if(nTRgType(ireg,jreg).eq.BT_TEMPER) cycle  !---: If fixed T reg

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-------: Filter interior points
          do j=jS+1, jN
            do i=iW+1, iE
              qh(i,j) = (qu(i,j-1)+qu(i-1,j)+qu(i,j+1)+qu(i+1,j)
     .                +fp*qu(i,j)) / (fp + dFour)
            end do
          end do

        end do
      end do

      case default
        print*, 'Wrong ncomp flag passed to Filter'
        stop
      end select


C---: Re-load filtered values
      do j=0, ny+1
        do i=0, nx+1
          qu(i,j) = qh(i,j)
        end do
      end do


      return
      end

*----------------------------------------------------------------------*
*                               Project                                *
*----------------------------------------------------------------------*

      subroutine Project(nx, ny,
     .                   nReg, nRegBrd,
     .                   nRegType, nMomBdTp,
     .                   dk,
     .                   dju, djv,
     .                   yeu, xzv, yzu, xev,
     .                   p, u, v)

      implicit none

      include "config.f"

      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Time-step size
      REAL    dk

C---: Metric information
      REAL    dju(0:mnx,0:mny),  djv(0:mnx,0:mny),
     .        yzu(0:mnx,0:mny),  xev(0:mnx,0:mny),
     .        yeu(0:mnx,0:mny),  xzv(0:mnx,0:mny)

      REAL    p(0:mnx,0:mny)
      REAL    u(0:mnx,0:mny), v(0:mnx,0:mny)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN
      REAL    pzi, pet, djk

C---: Double precision constants
      REAL    dFour
      parameter ( dFour = 4.00d+00 )

C---------------------------[ Update U ]---------------------------

C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

          if(nRegType(ireg,jreg).eq.RM_BLOCKG) cycle  !---: If a blockage

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-----:  Interior points
          do j=jS+1, jN
            do i=iW+1, iE-1
              pzi = p(i+1,j) - p(i,j)
              pet = (p(i+1,j+1)+p(i,j+1)-p(i+1,j-1)-p(i,j-1))/dFour
              djk = dju(i,j)*dk
              u(i,j) = u(i,j) - djk*(yeu(i,j)*pzi-yzu(i,j)*pet)
            end do
          end do

C-------: W E S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,WEST))

          case(BM_OUTLT1)   !---> Internal (interface)

            do j=jS+1,jN
              pzi = p(iW+1,j) - p(iW,j)
              pet = (p(iW+1,j+1)+p(iW,j+1)-p(iW+1,j-1)-p(iW,j-1))/dFour
              djk = dju(iW,j)*dk
              u(iW,j) = u(iW,j) - djk*(yeu(iW,j)*pzi-yzu(iW,j)*pet)
            end do

          case(BM_INTERN, BM_WALL1, BM_WALL2, BM_INLET, BM_OUTLT2)
          !---> Don't update internal, walls, inlet, outlet boundaries.
          !---> Note that internal boundaries are updated only at
          !---> East and North walls to avoid filtering twice.

          case default
            print*, 'Wrong nBdTypeW flag in region ',ireg,',',jreg
          end select

C-------: E A S T   Boundary :-------]

          select case(nMomBdTp(ireg,jreg,EAST))
          case(BM_INTERN, BM_OUTLT1)   !---> Internal (interface)
                                       !---> Outlet: fully developed flow
            do j=jS+1, jN
              pzi = p(iE+1,j) - p(iE,j)
              pet = (p(iE+1,j+1)+p(iE,j+1)-p(iE+1,j-1)-p(iE,j-1))/dFour
              djk = dju(iE,j)*dk
              u(iE,j) = u(iE,j) - djk*(yeu(iE,j)*pzi-yzu(iE,j)*pet)
             end do

          case(BM_WALL1, BM_WALL2, BM_INLET, BM_OUTLT2)
          !---> Don't update internal, walls, inlet, outlet boundaries.
          !---> Note that internal boundaries are updated only at
          !---> East and North walls to avoid filtering twice.

          case default
            print*, 'Wrong nBdTypeE flag in region ',ireg,',',jreg
          end select

        end do
      end do

C---------------------------[ Update V ]---------------------------

C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

          if(nRegType(ireg,jreg).ne.RM_BLOCKG) then  !---: If not a blockage

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

C-------: Interior points
          do j=jS+1, jN-1
            do i=iW+1, iE
              pzi = (p(i+1,j+1)+p(i+1,j)-p(i-1,j+1)-p(i-1,j))/dFour
              pet = p(i,j+1) - p(i,j)
              djk = djv(i,j)*dk
              v(i,j) = v(i,j) - djk*(-xev(i,j)*pzi+xzv(i,j)*pet)
            end do
          end do

C-------: S O U T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,SOUTH))
          case(BM_OUTLT1)   !---> Internal (interface)
                            !---> Outlet: fully developed flow
            do i=iW+1, iE
              pzi = (p(i+1,jS+1)+p(i+1,jS)-p(i-1,jS+1)-p(i-1,jS))/dFour
              pet = p(i,jS+1) - p(i,jS)
              djk = djv(i,jS)*dk
              v(i,jS) = v(i,jS) - djk*(-xev(i,jS)*pzi+xzv(i,jS)*pet)
            end do

          case(BM_INTERN, BM_WALL1, BM_WALL2, BM_INLET, BM_OUTLT2)
          !---> Don't update internal, walls, inlet, outlet boundaries.
          !---> Note that internal boundaries are updated only at
          !---> East and North walls to avoid filtering twice.

          case default
            print*, 'Wrong nBdTypeS flag in region ',ireg,',',jreg
          end select


C-------: N O R T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,NORTH))
          case(BM_INTERN, BM_OUTLT1)   !---> Internal (interface)
                                       !---> Outlet: fully developed flow
            do i=iW+1,iE
              pzi = (p(i+1,jN+1)+p(i+1,jN)-p(i-1,jN+1)-p(i-1,jN))/dFour
              pet = p(i,jN+1) - p(i,jN)
              djk = djv(i,jN)*dk
              v(i,jN) = v(i,jN) - djk*(-xev(i,jN)*pzi+xzv(i,jN)*pet)
            end do

          case(BM_WALL1, BM_WALL2, BM_INLET, BM_OUTLT2)
          !---> Don't update internal, walls, inlet, outlet boundaries.
          !---> Note that internal boundaries are updated only at
          !---> East and North walls to avoid filtering twice.

          case default
            print*, 'Wrong nBdTypeN flag in region ',ireg,',',jreg
          end select

          end if

        end do
      end do

      return
      end

*----------------------------------------------------------------------*
*                             DiffMaxNorm                              *
*----------------------------------------------------------------------*

      REAL function DiffMaxNorm(nx, ny, un, u)

      implicit none

      include "config.f"

      INTEGER nx, ny
      REAL    u(0:mnx,0:mny), un(0:mnx,0:mny)

      INTEGER i, j
      REAL    diff

c      diff = dabs(un(1,1)-u(1,1))
      diff = dabs(un(2,2)-u(2,2))

c***note: needs revision
c      do j=0,ny+1
c        do i=0,nx+1
      do j=2,ny-1
        do i=2,nx-1
          diff = max(diff,dabs(un(i,j)-u(i,j)))
        end do
      end do

      DiffMaxNorm = diff

      return
      end

*----------------------------------------------------------------------*
*                              DMaxNorm                                *
*----------------------------------------------------------------------*

      REAL function DMaxNorm(nx, ny, u)

      implicit none

      include "config.f"

      INTEGER nx, ny
      REAL    u(0:mnx,0:mny)

      INTEGER i, j
      REAL    diff

c      diff = dabs(u(1,1))
c      diff = dabs(u(2,2))
      diff = dabs(u(5,5))

c***note: needs revision
c      do j=1,ny-1
c        do i=1,nx-1
      do j=2,ny-1
        do i=2,nx-1
          diff = max(diff, dabs(u(i,j)))
        end do
      end do

      DMaxNorm = diff
    
      return
      end

*----------------------------------------------------------------------*
*                                PTDAvg                                *
*----------------------------------------------------------------------*

      subroutine PTDAvg(nx, ny, nReg, nRegBrd, nRegType, p, pav)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region type
      INTEGER nRegType(mgri,mgrj)

      REAL    p(0:mnx,0:mny), pav(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN
      REAL    sum

 
C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

          if(nRegType(ireg,jreg).eq.RM_BLOCKG) then

            do j=jS+1, jN-1
              do i=iW+1, iE-1
                pav(i,j) = 0.00d0
              end do
            end do

C---------: Next blockage
            cycle

          end if

          do j=jS, jN
            do i=iW, iE
              sum = (p(i,j) + p(i,j+1) + p(i+1,j+1) + p(i+1,j))/4.00d0
              if(dabs(sum).lt.1.d-20) sum = 0.00d0
              pav(i,j) = sum
            end do
          end do

        end do   !---: End I sweep
      end do   !---: End J sweep

      return
      end

*----------------------------------------------------------------------*
*                                VelAvg                                *
*----------------------------------------------------------------------*

      subroutine VelAvg(nx, ny, nReg, nRegBrd, nRegType,
     .                  u, v, util, vbar)

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

      REAL    u(0:mnx,0:mny),    v(0:mnx,0:mny)
      REAL    util(0:mnx,0:mny), vbar(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN


C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

          if(nRegType(ireg,jreg).eq.RM_BLOCKG) then  !---: If a blockage
            do j=jS,jN
              do i=iW,iE
                util(i,j) = 0.00d0
                vbar(i,j) = 0.00d0
              end do
            end do

C---------: Next region
            cycle

          end if

C-------: If not a blockage, then average velocities to natural grid
C-------: points depending on the boundary types

C-------: Velocities at natural grid points - Interior
c          do j=jS+1,jN-1
c            do i=iW+1,iE-1
          do j=jS,jN
            do i=iW,iE
              util(i,j) = (u(i,j) + u(i,j+1))/2.00d0
              vbar(i,j) = (v(i,j) + v(i+1,j))/2.00d0
            end do
          end do

        end do   !---: End I sweep
      end do   !---: End J sweep
 
      return
      end


*----------------------------------------------------------------------*
*                             dLinInterp                               *
*----------------------------------------------------------------------*

      REAL function dLinInterp(xs, x1, x2, y1, y2)

      implicit none
      REAL xs, x1, x2, y1, y2

      dLinInterp = y1 + (y2 - y1) * (xs - x1) / (x2 - x1)

      return
      end

*----------------------------------------------------------------------*
*                               TAveraged                              *
*----------------------------------------------------------------------*

      subroutine TAveraged(nx, ny, nScale,
     .                    nReg, nRegBrd, nTRgType,
     .                    dTRgVal,
     .                    t, tav)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, nScale

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Temperature region values
      REAL    dTRgVal(mgri,mgrj)

      REAL t(0:mnx,0:mny), tav(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN

 
C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

          if(nTRgType(ireg,jreg).eq.RT_TEMPER) then

            do j=jS, jN
              do i=iW, iE

                if (nScale.eq.0) then   !---> Large-scale
                  tav(i,j) = dTRgVal(ireg,jreg)
                else                    !---> Small-scale
                  tav(i,j) = 0.d0
                end if

              end do
            end do

C---------: Go to the next region
            cycle

          end if

          do j=jS, jN
            do i=iW, iE
              tav(i,j) = (t(i,j) + t(i,j+1) + t(i+1,j+1)
     .                  + t(i+1,j))/4.d0
            end do
          end do

        end do   !---: End I sweep
      end do   !---: End J sweep


      return
      end

*-----------------------|---|---|---V---|---|---|----------------------*
