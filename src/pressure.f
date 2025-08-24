*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                             pressure.F                               *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* Ppe         - Pseudo-pressure Poisson equation
* Sor         - Solve Poisson problem using Successive Overrelaxation
* Slor        - Successive Line Overrelaxation (regular ordering)
* SlorRB      - Red-black ordered Successive Line Overrelaxation
* SlorRBP     - Parallelizable Red-black Line SOR
* Divergence  - Divergence of Velocity vector
* RhsPpe      - Right-hand side of pseudo-pressure Poisson equation
*
*----------------------------------------------------------------------*

/*
  Preprocessor header file(s) 
*/

#include "wolfd2.h"


*----------------------------------------------------------------------*
*                                  Ppe                                 *
*----------------------------------------------------------------------*

      subroutine Ppe(nx, ny, 
     .               nReg, nRegBrd,
     .               nRegType,
     .               lCartesGrid,
     .               nPpeSolver, msorit, nSorConv, 
     .               dk, sortol, sorrel,
     .               rau, rbu, rbv, rgv,
     .               xeu, yeu, xzv, yzv,
     .               u, v, p)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn    ! Maximum number of total points
      parameter ( mn  = mnx*mny )

C---: Procedure arguments
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Cartesian grid flag
      LOGICAL lCartesGrid

      INTEGER nPpeSolver, msorit, nSorConv

      REAL    dk, sortol, sorrel

      REAL    u(0:mnx,0:mny), v(0:mnx,0:mny), p(0:mnx,0:mny)

      REAL    rau(0:mnx,0:mny), rbu(0:mnx,0:mny),
     .        rbv(0:mnx,0:mny), rgv(0:mnx,0:mny),
     .        xeu(0:mnx,0:mny), yeu(0:mnx,0:mny),
     .        xzv(0:mnx,0:mny), yzv(0:mnx,0:mny)

      REAL    div(0:mnx,0:mny)
      REAL    a(mn,5), b(mn)

C---: Local variables
      INTEGER i, j, ind, ireg, jreg, iW, iE, jS, jN

C---: Double precision constants
      REAL    dZero, dOne
      parameter ( dZero = 0.00d+00 ,
     .            dOne  = 1.00d+00 )


C---: Note:  Boundary conditions for pressure must be enforced
C---: prior to calling this subroutine.


C---: Compute divergence of auxiliary velocity
      call Divergence(nx, ny, 1,
     .                xeu, yeu, xzv, yzv,
     .                u, v, div)

C---: Construct discrete left-hand side of ppe
      do j=2, ny
        do i=2, nx
          ind = (j-2)*(nx-1)+i-1
          a(ind,1) = rgv(i,j-1)
          a(ind,2) = rau(i-1,j)
          a(ind,3) = -rau(i,j)-rau(i-1,j)-rgv(i,j)-rgv(i,j-1)
          a(ind,4) = rau(i,j)
          a(ind,5) = rgv(i,j)
        end do
      end do

C---: On blockages assign a value to the pressure at internal points
      do jreg=1, nReg(_J_)  !---: Scan regions: Vertical sweep
        do ireg=1, nReg(_I_)  !-----: Scan regions: Horizontal sweep

          if(nRegType(ireg,jreg).eq.RM_BLOCKG) then  !---: If a blockage

C---------: Extract boundary indeces
            iW = nRegBrd(ireg,jreg,WEST)
            iE = nRegBrd(ireg,jreg,EAST)
            jS = nRegBrd(ireg,jreg,SOUTH)
            jN = nRegBrd(ireg,jreg,NORTH)

C********** [Needs revision:  Indeces]
c            do j=jS+1, jN
c              do i=iW+1, iE
            do j=jS+2, jN-1
              do i=iW+2, iE-1

                ind = (j-2)*(nx-1)+i-1
                a(ind,1) = dZero
                a(ind,2) = dZero
                a(ind,3) = dOne
                a(ind,4) = dZero
                a(ind,5) = dZero

C************** [Needs revision]
                div(i,j) = dZero

              end do
            end do

C---------: Also if neighboring regions are also blockages, zero out
C---------: common boundaries.


            if(ireg.gt.1.and.nRegType(ireg-1,jreg).eq.RM_BLOCKG) then
              do j=jS+1, jN
                ind = (j-2)*(nx-1)+iW
                a(ind,1) = dZero
                a(ind,2) = dZero
                a(ind,3) = dOne
                a(ind,4) = dZero
                a(ind,5) = dZero
                div(iW+1,j) = dZero
              end do
            end if

            if(ireg.lt.nReg(_I_).and.nRegType(ireg+1,jreg).eq.
     .        RM_BLOCKG) then
              do j=jS+1, jN
                ind = (j-2)*(nx-1)+iE-1
                a(ind,1) = dZero
                a(ind,2) = dZero
                a(ind,3) = dOne
                a(ind,4) = dZero
                a(ind,5) = dZero
                div(iE,j) = dZero
              end do
            end if

            if(jreg.gt.1.and.nRegType(ireg,jreg-1).eq.RM_BLOCKG) then
              do i=iW+1, iE
                ind = (jS-1)*(nx-1)+i-1
                a(ind,1) = dZero
                a(ind,2) = dZero
                a(ind,3) = dOne
                a(ind,4) = dZero
                a(ind,5) = dZero
                div(i,jS+1) = dZero
              end do
            end if

            if(jreg.lt.nReg(_J_).and.nRegType(ireg,jreg+1).eq.
     .        RM_BLOCKG) then
              do i=iW+1, iE
                ind = (jN-2)*(nx-1)+i-1
                a(ind,1) = dZero
                a(ind,2) = dZero
                a(ind,3) = dOne
                a(ind,4) = dZero
                a(ind,5) = dZero
                div(i,jN) = dZero
              end do
            end if
       
          end if

        end do  !---: End vertical sweep
      end do  !---: End horizonta sweep


      nSorConv = 0

      select case(nPpeSolver)
      case(1)   !---> Point successive over-relaxation
        call Sor(nx, ny, lCartesGrid, 
     .           msorit, nSorConv, 
     .           dk, sortol, sorrel,
     .           rbu, rbv, a, b, div, p)

      case(2)   !---> Line successive over-relaxation
        call Slor(nx, ny, lCartesGrid, 1, 
     .            msorit, nSorConv, 
     .            dk, sortol, sorrel,
     .            rbu, rbv, a, b, div, p)

      case(3)   !---> Line SOR with Red/Black ordering of lines 
        call SlorRB(nx, ny, lCartesGrid, 1, 
     .              msorit, nSorConv, 
     .              dk, sortol, sorrel,
     .              rbu, rbv, a, b, div, p)

      case(4)   !---> Parallel version of Line SOR with R/B ordering
        call SlorRBP(nx, ny, lCartesGrid, 
     .               msorit, nSorConv, 
     .               dk, sortol, sorrel,
     .               rbu, rbv, a, b, div, p)

      case(5)   !---> Red-Black ordered point SOR
        call SorRB(nx, ny, lCartesGrid, 
     .             msorit, nSorConv, 
     .             dk, sortol, sorrel,
     .             rbu, rbv, a, b, div, p)

      case(6)   !---> Parallel Red-Black ordered point SOR
        call SorRBP(nx, ny, lCartesGrid, 
     .              msorit, nSorConv, 
     .              dk, sortol, sorrel,
     .              rbu, rbv, a, b, div, p)

      case default
        print*, 'Wrong nPpeSolver flag passed to Ppe'
      end select

      if(nSorConv.lt.1) then
        print*, 'Warning: SOR iterations did not converge after ',
     .    msorit, ' iterations.'
        nSorConv = msorit
      end if

      return
      end

*----------------------------------------------------------------------*
*                              Divergence                              *
*----------------------------------------------------------------------*
*
* nloc -> Location flag:
*   1: Center of Pressure control volume
*   2: Natural grid points
* 
* The corresponding metric information to be passed for each case is:
*   nloc   xet   yet   xzi   yzi
*   ----   ---   ---   ---   ---
*     1    xeu   yeu   xzv   yzv
*     2    xev   yev   xzu   yzu 
*
      subroutine Divergence(nx, ny, nloc,
     .                      xet, yet, xzi, yzi,
     .                      u, v, div)

      implicit none

C---: Configurable parameters
      include "config.f"

C---: Arguments
      INTEGER nx, ny, nloc
      REAL    xet(0:mnx,0:mny), yet(0:mnx,0:mny),
     .        xzi(0:mnx,0:mny), yzi(0:mnx,0:mny)

      REAL    u(0:mnx,0:mny), v(0:mnx,0:mny), div(0:mnx,0:mny)

C---: Local declarations
      INTEGER i, j
      REAL    uci1j, ucij, vcij1, vcij

      select case(nloc)
      case(1)   !---> Divergence at center of Pressure C.V.
        do j=1,ny
          do i=1,nx
            ucij  = yet(i,j)*u(i,j) - xet(i,j)*(v(i+1,j) + v(i,j)
     .            + v(i+1,j-1) + v(i,j-1))/4.d0
            uci1j = yet(i-1,j)*u(i-1,j) - xet(i-1,j)*(v(i,j) 
     .            + v(i-1,j) + v(i,j-1) + v(i-1,j-1))/4.d0
            vcij  = xzi(i,j)*v(i,j) - yzi(i,j)*(u(i,j+1) + u(i-1,j+1)
     .            + u(i,j) + u(i-1,j))/4.d0
            vcij1 = xzi(i,j-1)*v(i,j-1) - yzi(i,j-1)*(u(i,j)
     .            + u(i-1,j) + u(i,j-1) + u(i-1,j-1))/4.d0
            div(i,j) = ucij - uci1j + vcij - vcij1
          end do
        end do

      case(2)   !---> Divergence at natural grid points
        do j=1,ny
          do i=1,nx
            uci1j = yet(i+1,j)*(u(i,j+1) + u(i+1,j+1) + u(i,j)
     .            + u(i+1,j))/4.d0 - xet(i+1,j)*v(i+1,j)
            ucij  = yet(i,j)*(u(i-1,j+1) + u(i,j+1) + u(i-1,j)
     .            + u(i,j))/4.d0 - xet(i,j)*v(i,j)
            vcij1 = xzi(i,j+1)*(v(i,j+1) + v(i+1,j+1) + v(i,j)
     .            + v(i+1,j))/4.d0 - yzi(i,j+1)*u(i,j+1)
            vcij  = xzi(i,j)*(v(i,j) + v(i+1,j) + v(i,j-1) 
     .            + v(i+1,j-1))/4.d0 -yzi(i,j)*u(i,j)
            div(i,j) = uci1j - ucij + vcij1 - vcij
          end do
        end do

      case default
        print*, 'Error: Wrong location flag passed to Divergence: ', 
     .    nloc
        stop
      end select

      return
      end

*----------------------------------------------------------------------*
*                                RhsPpe                                *
*----------------------------------------------------------------------*

      subroutine RhsPpe(nx, ny, lCartesGrid, dk, 
     .                  rbu, rbv, div, p, b)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn    ! Maximum number of total points
      parameter ( mn  = mnx*mny )

C---: Arguments
      INTEGER nx, ny
      REAL    dk
      REAL    rbu(0:mnx,0:mny), rbv(0:mnx,0:mny)
      REAL    p(0:mnx,0:mny), div(0:mnx,0:mny)
      REAL    b(mn)

C---: Cartesian grid flag
      LOGICAL lCartesGrid

C---: Local declarations
      INTEGER i, j, ind

      do j=2,ny
        do i=2,nx
          ind = (j-2)*(nx-1)+i-1
          b(ind) = div(i,j)/dk
        end do
      end do

C---: If grid is not cartesian, add mixed derivative terms from LHS
      if (.not.lCartesGrid) then

        do j=2,ny
          do i=2,nx
            ind = (j-2)*(nx-1)+i-1
            b(ind) = b(ind)
     .           -(rbu(i,j)*(p(i+1,j+1)+p(i,j+1)-p(i+1,j-1)-p(i,j-1))
     .           -rbu(i-1,j)*(p(i,j+1)+p(i-1,j+1)-p(i,j-1)-p(i-1,j-1))
     .           +rbv(i,j)*(p(i+1,j+1)+p(i+1,j)-p(i-1,j+1)-p(i-1,j))
     .           -rbv(i,j-1)*(p(i+1,j)+p(i+1,j-1)-p(i-1,j)-p(i-1,j-1)))
          end do
        end do

      end if

      return
      end

*----------------------------------------------------------------------*
*                                  Sor                                 *
*----------------------------------------------------------------------*

      subroutine Sor(nx, ny, lCartesGrid,
     .               msorit, nConv, 
     .               dk, sortol, sorrel, 
     .               rbu, rbv, a, b, div, p)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn    ! Maximum number of total points
      parameter ( mn  = mnx*mny )

C---: Procedure arguments
      INTEGER nx, ny, msorit, nConv
      REAL    dk, sortol, sorrel
      REAL    rbu(0:mnx,0:mny), rbv(0:mnx,0:mny)
      REAL    a(mn,5), b(mn)
      REAL    div(0:mnx,0:mny), p(0:mnx,0:mny)

C---: Cartesian grid flag
      LOGICAL lCartesGrid

C---: Local declarations
      INTEGER m, i, j, ind
      REAL    dif, sum

      nConv = 0

C---: If grid is cartesian, evaluate right-hand side of ppe
C---: only once:  before startin SOR iterations
      if (lCartesGrid)
     .  call RhsPpe(nx, ny, lCartesGrid, dk,
     .              rbu, rbv, div, p, b)

C---: Start SOR iterations
      do m=1, msorit

C-----: If grid is not cartesian, evaluate right-hand side of ppe
C-----: at every SOR iteration
        if (.not.lCartesGrid)
     .    call RhsPpe(nx, ny, lCartesGrid, dk,
     .                rbu, rbv, div, p, b)

        dif = 0.d0

        do j=2,ny
          do i=2,nx
            ind = (j-2)*(nx-1)+i-1
            sum = b(ind) - a(ind,1)*p(i,j-1) - a(ind,2)*p(i-1,j)
     .          - a(ind,4)*p(i+1,j) - a(ind,5)*p(i,j+1)
            sum = sum/a(ind,3) - p(i,j)
            p(i,j) = p(i,j) + sorrel*sum
            sum = dabs(sum)
            dif = dmax1(dif, sum)
          end do
        end do

        if(m.gt.1.and.dif.lt.sortol) then
          nConv = m
          return 
        end if        
      end do   !---> End SOR iterations

      return
      end


*----------------------------------------------------------------------*
*                                 SorRB                                *
*----------------------------------------------------------------------*

      subroutine SorRB(nx, ny, lCartesGrid,
     .                 msorit, nConv, 
     .                 dk, sortol, sorrel, 
     .                 rbu, rbv, a, b, div, p)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn    ! Maximum number of total points
      parameter ( mn  = mnx*mny )

C---: Procedure arguments
      INTEGER nx, ny, msorit, nConv
      REAL    dk, sortol, sorrel
      REAL    rbu(0:mnx,0:mny), rbv(0:mnx,0:mny)
      REAL    a(mn,5), b(mn)
      REAL    div(0:mnx,0:mny), p(0:mnx,0:mny)

C---: Cartesian grid flag
      LOGICAL lCartesGrid

C---: Local declarations
      INTEGER m, i, j, ind
      REAL    dif, sum

      nConv = 0

C---: If grid is cartesian, evaluate right-hand side of ppe
C---: only once:  before startin SOR iterations
      if (lCartesGrid)
     .  call RhsPpe(nx, ny, lCartesGrid, dk,
     .              rbu, rbv, div, p, b)

C---: Start SOR iterations
      do m=1, msorit

C-----: If grid is not cartesian, evaluate right-hand side of ppe
C-----: at every SOR iteration
        if (.not.lCartesGrid)
     .    call RhsPpe(nx, ny, lCartesGrid, dk,
     .                rbu, rbv, div, p, b)

        dif = 0.d0

C-----: Relax `Black' points
        do j=2, ny
          do i=(2+mod(j,2)), nx, 2
            ind = (j-2)*(nx-1)+i-1

            sum = b(ind) - a(ind,1)*p(i,j-1) - a(ind,2)*p(i-1,j)
     .          - a(ind,4)*p(i+1,j) - a(ind,5)*p(i,j+1)
            sum = sum/a(ind,3) - p(i,j)
            p(i,j) = p(i,j) + sorrel*sum
            sum = dabs(sum)
            dif = dmax1(dif, sum)

          end do
        end do

C-----: Relax `Red' points
        do j=2, ny
          do i=(2+mod(j+1,2)), nx, 2
            ind = (j-2)*(nx-1)+i-1

            sum = b(ind) - a(ind,1)*p(i,j-1) - a(ind,2)*p(i-1,j)
     .          - a(ind,4)*p(i+1,j) - a(ind,5)*p(i,j+1)
            sum = sum/a(ind,3) - p(i,j)
            p(i,j) = p(i,j) + sorrel*sum
            sum = dabs(sum)
            dif = dmax1(dif, sum)

          end do
        end do

        if(m.gt.1.and.dif.lt.sortol) then
          nConv = m
          return 
        end if        
      end do   !---> End SOR iterations

      return
      end


*----------------------------------------------------------------------*
*                                 SorRBP                               *
*----------------------------------------------------------------------*

      subroutine SorRBP(nx, ny, lCartesGrid,
     .                  msorit, nConv, 
     .                  dk, sortol, sorrel, 
     .                  rbu, rbv, a, b, div, p)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn    ! Maximum number of total points
      parameter ( mn  = mnx*mny )

C---: Procedure arguments
      INTEGER nx, ny, msorit, nConv
      REAL    dk, sortol, sorrel
      REAL    rbu(0:mnx,0:mny), rbv(0:mnx,0:mny)
      REAL    a(mn,5), b(mn)
      REAL    div(0:mnx,0:mny), p(0:mnx,0:mny)

C---: Cartesian grid flag
      LOGICAL lCartesGrid

C---: Local declarations
      INTEGER m, i, j, ind
      REAL    dif(0:nx+1,0:ny+1)
      REAL    difmax, sum

      nConv = 0

C---: If grid is cartesian, evaluate right-hand side of ppe
C---: only once:  before startin SOR iterations
      if (lCartesGrid)
     .  call RhsPpe(nx, ny, lCartesGrid, dk,
     .              rbu, rbv, div, p, b)

C---: Start SOR iterations
      do m=1, msorit

C-----: If grid is not cartesian, evaluate right-hand side of ppe
C-----: at every SOR iteration
        if (.not.lCartesGrid)
     .    call RhsPpe(nx, ny, lCartesGrid, dk,
     .                rbu, rbv, div, p, b)

C-----: Initialize difference vector
        do j=1, ny
          do i=1, nx
            dif(i,j) = 0.00d0
          end do
        end do

C-----: Relax `Black' points

C$DIR loop_parallel(ivar=j,threads)
C$DIR loop_private(i,ind)
C$DIR loop_private(sum)

        do j=2, ny
          do i=(2+mod(j,2)), nx, 2
            ind = (j-2)*(nx-1)+i-1

            sum = b(ind) - a(ind,1)*p(i,j-1) - a(ind,2)*p(i-1,j)
     .          - a(ind,4)*p(i+1,j) - a(ind,5)*p(i,j+1)
            sum = sum/a(ind,3) - p(i,j)
            p(i,j) = p(i,j) + sorrel*sum

            dif(i,j) = dmax1(dif(i,j), dabs(sum))

          end do
        end do

C-----: Relax `Red' points

C$DIR loop_parallel(ivar=j,threads)
C$DIR loop_private(i,ind)
C$DIR loop_private(sum)

        do j=2, ny
          do i=(2+mod(j+1,2)), nx, 2
            ind = (j-2)*(nx-1)+i-1

            sum = b(ind) - a(ind,1)*p(i,j-1) - a(ind,2)*p(i-1,j)
     .          - a(ind,4)*p(i+1,j) - a(ind,5)*p(i,j+1)
            sum = sum/a(ind,3) - p(i,j)
            p(i,j) = p(i,j) + sorrel*sum

            dif(i,j) = dmax1(dif(i,j), dabs(sum))

          end do
        end do

C-----: Max-norm of dif
        difmax = 0.00d0
        do j=2, ny
          do i=2, nx
            difmax = dmax1(difmax, dif(i,j))
          end do
        end do

        if(m.gt.1.and.difmax.lt.sortol) then
          nConv = m
          return 
        end if        
      end do   !---> End SOR iterations

      return
      end


*----------------------------------------------------------------------*
*                                 Slor                                 *
*----------------------------------------------------------------------*
*
* Succesive Line Over-relaxation solver (regular ordering)
* 98.04.04
* vegarzon
*
* ndir     -> ndirection flag, 0 for i=const, 1 for j=const
* nlines   -> number of lines to solve
* ncomps   -> number of components in each line
* indal    -> matrix a index corresponding to aline(?,1)
* indar    -> matrix a index corresponding to aline(?,3)

      subroutine Slor(nx, ny, lCartesGrid, ndir, 
     .                msorit, nConv,
     .                dk, sortol, sorrel,
     .                rbu, rbv, a, b, div, p)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn  ,  ! Maximum number of total points
     .        mnl    ! Maximum number of elemenents in line-solver arrays 
      parameter ( mn  = mnx*mny        )
      parameter ( mnl = max0(mnx, mny) )

C---: Procedure arguments
      INTEGER nx, ny, ndir, msorit, nConv
      REAL    dk, sortol, sorrel
      REAL    rbu(0:mnx,0:mny), rbv(0:mnx,0:mny)
      REAL    a(mn,5), b(mn)
      REAL    div(0:mnx,0:mny), p(0:mnx,0:mny)

C---: Cartesian grid flag
      LOGICAL lCartesGrid

C---: Local declarations
      REAL    bline(mnl), plold(mnl), aline(3,mnl)
      INTEGER nlines, ncomps, indal, indar, indbl, indbr
      INTEGER ind, m, nc, nl
      REAL    plprev, plnext, sum, dif

      nConv = 0

      if (ndir.eq.0) then
        nlines = nx-1
        ncomps = ny-1
        indal = 1
        indar = 5
        indbl = 2
        indbr = 4
      else if (ndir.eq.1) then
        nlines = ny-1
        ncomps = nx-1
        indal = 2
        indar = 4
        indbl = 1
        indbr = 5
      else
        print *, 'Error in direction flag passed to SLor'
        return
      end if


C---: If grid is cartesian, evaluate right-hand side of ppe
C---: only once:  before startin SOR iterations
      if (lCartesGrid)
     .  call RhsPpe(nx, ny, lCartesGrid, dk,
     .              rbu, rbv, div, p, b)


C---: Start Sor iterations
      do m=1, msorit   

        dif = 0.0d0

C-----: If grid is not cartesian, evaluate right-hand side of ppe
C-----: at every SOR iteration
        if (.not.lCartesGrid)
     .    call RhsPpe(nx, ny, lCartesGrid, dk,
     .                rbu, rbv, div, p, b)


C-----: Solve for each line
        do nl=2, nlines+1

C-------: Fill out matrix aline and vector bline
          do nc=2, ncomps+1

            if(ndir.eq.0) then
              ind = (nc-2)*nlines+nl-1
              plold(nc-1) = p(nl,nc)
              plprev = p(nl-1,nc)
              plnext = p(nl+1,nc)

            else if(ndir.eq.1) then 
              ind = (nl-2)*ncomps+nc-1
              plold(nc-1) = p(nc,nl)
              plprev = p(nc,nl-1)
              plnext = p(nc,nl+1)

            end if

            aline(1,nc-1) = a(ind,indal)
            aline(2,nc-1) = a(ind,3)
            aline(3,nc-1) = a(ind,indar)
            bline(nc-1)   = b(ind) 
     .        -(a(ind,indbl)*plprev + a(ind,indbr)*plnext)

          end do

C---------: Solve the tri-diagonal system via LU factorization
C---------: and load returning vector onto array bline
            call AltTridLU(ncomps, aline, bline)

C-------: Overrelaxation
          do nc=1, ncomps
            sum = bline(nc)-plold(nc)
            bline(nc) = plold(nc) + sorrel*sum
            sum = dabs(sum)
            dif = dmax1(dif, sum)
          end do

C-------: Load line results into full matrix
          do nc=2, ncomps+1
            if(ndir.eq.0) p(nl,nc) = bline(nc-1)
            if(ndir.eq.1) p(nc,nl) = bline(nc-1)
          end do

        end do

        if ((m.gt.1).and.(dif.lt.sortol)) then
          nConv = m
          return 
        end if

      end do   !---> End SOR iterations

      return
      end


*----------------------------------------------------------------------*
*                                SlorRB                                *
*----------------------------------------------------------------------*
*
* Succesive Line Over-relaxation solver with red-black ordering
* 98.02.04
* vegarzon
*
* ndir     -> ndirection flag, 0 for i=const, 1 for j=const
* nlines   -> number of lines to solve
* ncomps   -> number of components in each line
* indal    -> matrix a index corresponding to aline(?,1)
* indar    -> matrix a index corresponding to aline(?,3)

      subroutine SlorRB(nx, ny, lCartesGrid, ndir, 
     .                  msorit, nConv, 
     .                  dk, sortol, sorrel,
     .                  rbu, rbv, a, b, div, p)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn  ,  ! Maximum number of total points
     .        mnl    ! Maximum number of elemenents in line-solver arrays 
      parameter ( mn  = mnx*mny        )
      parameter ( mnl = max0(mnx, mny) )

C---: Arguments
      INTEGER nx, ny, ndir, msorit, nConv
      REAL    dk, sortol, sorrel
      REAL    rbu(0:mnx,0:mny), rbv(0:mnx,0:mny)
      REAL    a(mn,5), b(mn)
      REAL    div(0:mnx,0:mny), p(0:mnx,0:mny)

C---: Cartesian grid flag
      LOGICAL lCartesGrid

C---: Local declarations
      INTEGER nlines, ncomps, indal, indar, indbl, indbr
      INTEGER ind, m, k, nc, nl
      REAL    plprev, plnext, sum, dif

      REAL    bline(mnl), plold(mnl), aline(3,mnl)

C---: Default state of returning iteration counter
      nConv = 0

      if (ndir.eq.0) then
        nlines = nx-1
        ncomps = ny-1
        indal = 1
        indar = 5
        indbl = 2
        indbr = 4
      else if (ndir.eq.1) then
        nlines = ny-1
        ncomps = nx-1
        indal = 2
        indar = 4
        indbl = 1
        indbr = 5
      else
        print *, 'Error in direction flag passed to SLor'
        return
      end if


C---: If grid is cartesian, evaluate right-hand side of ppe
C---: only once:  before startin SOR iterations
      if (lCartesGrid)
     .  call RhsPpe(nx, ny, lCartesGrid, dk,
     .              rbu, rbv, div, p, b)

C---: Start Sor iterations
      do m=1, msorit

cC-------: If grid is not cartesian, evaluate right-hand side of ppe
cC-------: at every SOR iteration
c          if (.not.lCartesGrid)
c     .      call RhsPpe(nx, ny, lCartesGrid, dk,
c     .                  rbu, rbv, div, p, b)

        dif = 0.00d0

C-----: Red solve: k=2, black solve: k=3
        do k=2, 3

C-------: If grid is not cartesian, evaluate right-hand side of ppe
C-------: at every SOR iteration
          if (.not.lCartesGrid)
     .      call RhsPpe(nx, ny, lCartesGrid, dk,
     .                  rbu, rbv, div, p, b)

C-------: Solve for each line
          do nl=k, nlines+1, 2

C---------: Fill out matrix aline and vector bline
            do nc=2, ncomps+1

              if(ndir.eq.0) then
                ind = (nc-2)*nlines+nl-1
                plold(nc-1) = p(nl,nc)
                plprev      = p(nl-1,nc)
                plnext      = p(nl+1,nc)

              else if(ndir.eq.1) then 
                ind         = (nl-2)*ncomps+nc-1
                plold(nc-1) = p(nc,nl)
                plprev      = p(nc,nl-1)
                plnext      = p(nc,nl+1)
              end if

C-----------: Load LHS matrix and RHS vector for tri-diagonal system
              aline(1,nc-1) = a(ind,indal)
              aline(2,nc-1) = a(ind,3)
              aline(3,nc-1) = a(ind,indar)
              bline(nc-1)   = b(ind) - (a(ind,indbl)*plprev
     .                      + a(ind,indbr)*plnext)

            end do

C---------: Solve the tri-diagonal system via LU factorization
C---------: and load returning vector onto array bline
            call AltTridLU(ncomps, aline, bline)

C---------: Overrelaxation
            do nc=1, ncomps
              sum = bline(nc) - plold(nc)
              bline(nc) = plold(nc) + sorrel*sum
              sum = dabs(sum)
              dif = dmax1(dif, sum)
            end do

C---------: Load line results into full matrix
            do nc=2, ncomps+1
              if(ndir.eq.0) p(nl,nc) = bline(nc-1)
              if(ndir.eq.1) p(nc,nl) = bline(nc-1)
            end do

          end do

        end do

        if ((m.gt.1).and.(dif.lt.sortol)) then
          nConv = m
          return 
        end if

      end do   !---> End SOR iterations

      return
      end


*----------------------------------------------------------------------*
*                               SlorRBP                                *
*----------------------------------------------------------------------*
*
* Parallelizable Line SOR solver with red-black ordering
* 98.02.03
* vegarzon
*
* ndir     -> Assume j=const
* nlines   -> number of lines to solve
* ncomps   -> number of components in each line
* indal    -> matrix a index corresponding to aline(?,1)
* indar    -> matrix a index corresponding to aline(?,3)

      subroutine SlorRBP(nx, ny, lCartesGrid, 
     .                   msorit, nConv,
     .                   dk, sortol, sorrel,
     .                   rbu, rbv, a, b, div, p)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn  ,  ! Maximum number of total points
     .        mnl    ! Maximum number of elemenents in line-solver arrays 
      parameter ( mn  = mnx*mny        )
      parameter ( mnl = max0(mnx, mny) )

C---: Arguments
      INTEGER nx, ny, msorit, nConv
      REAL    dk, sortol, sorrel
      REAL    rbu(0:mnx,0:mny), rbv(0:mnx,0:mny)
      REAL    a(mn,5), b(mn)
      REAL    div(0:mnx,0:mny), p(0:mnx,0:mny)

C---: Cartesian grid flag
      LOGICAL lCartesGrid


C---: Local declarations
      INTEGER nlines, ncomps, indal, indar, indbl, indbr
      INTEGER ind, m, nc, nl
      REAL    sum, dif

      REAL    pn(0:mnx,0:mny)
      REAL    aline(3,mnl), bline(mnl)

      nlines = ny-1
      ncomps = nx-1
      indal = 2
      indar = 4
      indbl = 1
      indbr = 5

C---: Default state of returning iteration counter
      nConv = 0

C---: If grid is cartesian, evaluate right-hand side of ppe
C---: only once:  before startin SOR iterations
      if (lCartesGrid)
     .  call RhsPpe(nx, ny, lCartesGrid, dk,
     .              rbu, rbv, div, p, b)

C---: Start Sor iterations
      do m=1, msorit

C-----: Initialize array pn (for overrelaxation and conv. test)
        do nl=1, nlines+2
          do nc=1, ncomps+2
            pn(nc,nl) = p(nc,nl)
          end do
        end do

C-----: If grid is not cartesian, evaluate right-hand side of ppe
C-----: at every SOR iteration
        if (.not.lCartesGrid)
     .    call RhsPpe(nx, ny, lCartesGrid, dk,
     .                rbu, rbv, div, p, b)


C-----: Solve for each line (red solve)

C$DIR loop_parallel(ivar=nl,threads)
C$DIR loop_private(ind,nc)
C$DIR loop_private(aline,bline)

        do nl=2, nlines+1, 2

C-------: Fill out tri-diagonal matrix aline and vector bline
          do nc=2, ncomps+1
            ind = (nl-2)*ncomps+nc-1
            aline(1,nc-1) = a(ind,indal)
            aline(2,nc-1) = a(ind,3)
            aline(3,nc-1) = a(ind,indar)
            bline(nc-1)   = b(ind) - (a(ind,indbl)*p(nc,nl-1) 
     .                    + a(ind,indbr)*p(nc,nl+1))
          end do

C-------: Solve the tri-diagonal system via LU factorization
C-------: and store returning vector in array bline
          call AltTridLU(ncomps, aline, bline)

C-------: Load line results into full matrix
          do nc=2, ncomps+1
            p(nc,nl) = bline(nc-1)
          end do
        end do

C-----: Overrelaxation
        do nl=2, nlines+1, 2
          do nc=2, ncomps+1
            p(nc,nl) = pn(nc,nl) + sorrel*(p(nc,nl)-pn(nc,nl))
          end do
        end do

C-----: If grid is not cartesian, evaluate right-hand side of ppe
C-----: at every SOR iteration
        if (.not.lCartesGrid)
     .    call RhsPpe(nx, ny, lCartesGrid, dk,
     .                rbu, rbv, div, p, b)


C-----: Solve for each line (black solve)

C$DIR loop_parallel(ivar=nl,threads)
C$DIR loop_private(ind,nc)
C$DIR loop_private(aline,bline)

        do nl=3, nlines+1, 2

C-------: Fill out tri-diagonal matrix aline and vector bline
          do nc=2, ncomps+1
            ind = (nl-2)*ncomps+nc-1
            aline(1,nc-1) = a(ind,indal)
            aline(2,nc-1) = a(ind,3)
            aline(3,nc-1) = a(ind,indar)
            bline(nc-1) = b(ind) - (a(ind,indbl)*p(nc,nl-1) 
     .        + a(ind,indbr)*p(nc,nl+1))
          end do

C-------: Solve the tri-diagonal system via LU factorization
C-------: and store returning vector in array bline
          call AltTridLU(ncomps, aline, bline)

C-------: Load line results into full matrix
          do nc=2, ncomps+1
            p(nc,nl) = bline(nc-1)
          end do
        end do

C-----: Overrelaxation
        do nl=3, nlines+1, 2
          do nc=2, ncomps+1
            p(nc,nl) = pn(nc,nl) + sorrel*(p(nc,nl)-pn(nc,nl))
          end do
        end do

C-----: Convergence test
        dif = 0.00d0
        do nl=2, nlines+1
          do nc=2, ncomps+1
            sum = dabs(pn(nc,nl)-p(nc,nl))
            dif = dmax1(dif, sum)
          end do
        end do

        if ((m.gt.1).and.(dif.lt.sortol)) then
          nConv = m
          return 
        end if

      end do   !---> End SOR iterations

      return
      end


*-----------------------|---|---|---V---|---|---|----------------------*
