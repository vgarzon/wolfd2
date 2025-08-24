*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                             momentum.F                               *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* nAuxMomentum - Driver function for Newton-Kantarovich iterations
* XMomentum    - Delta-form quasilinear auxiliary x-momentum equation
* YMomentum    - Delta-form quasilinear auxiliary y-momentum equation
* ConvCoef     - Discretization of convective terms (u,v,T)
* DConvU       - Discretization of x-momentum convective terms
* DDiffU       - Discretization of x-momentum diffusive terms
* DConvV       - Discretization of y-momentum convective terms
* DDiffV       - Discretization of y-momentum diffusive terms
* PorosCoef    - Discretization of Dupuit-Forchheimer terms
* TridLU       - Tridiagonal LU decomposition (Thomas algorithm)
* AltTridLU    - Alternate LU factorization
*
*----------------------------------------------------------------------*

/*
  Preprocessor header file(s) 
*/

#include "wolfd2.h"

*----------------------------------------------------------------------*
*                           nAuxMomentum                               *
*----------------------------------------------------------------------*

      INTEGER function nAuxMomentum(nx, ny, 
     .                              mqiter,
     .                              nReg, nRegBrd,
     .                              nRegType, nMomBdTp,
     .                              dk, re, fr, qtol,
     .                              dPRporos, dPRporc1, dPRporc2,
     .                              dBCVal,
     .                              ran, rbn, rgn,
     .                              rac, rbc, rgc,
     .                              dju, djv,
     .                              xec, yec, xzn, yzn,
     .                              xen, yen, xzc, yzc,
     .                              xeu, yeu, xzu, yzu,
     .                              xev, yev, xzv, yzv,
     .                              d, dn,
     .                              un, vn, us, vs)

      implicit none

C---: Configurable parameters
      include "config.f"

C---: Subroutine arguments
      INTEGER nx, ny, mqiter

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

      REAL    dk, re, fr, qtol

      REAL    dPRporos(mgri,mgrj), dPRporc1(mgri,mgrj), 
     .        dPRporc2(mgri,mgrj)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,*)

      REAL    ran(0:mnx,0:mny), rbn(0:mnx,0:mny), rgn(0:mnx,0:mny),
     .        rac(0:mnx,0:mny), rbc(0:mnx,0:mny), rgc(0:mnx,0:mny),
     .        dju(0:mnx,0:mny), djv(0:mnx,0:mny)

      REAL    xec(0:mnx,0:mny), yec(0:mnx,0:mny),
     .        xzn(0:mnx,0:mny), yzn(0:mnx,0:mny),
     .        xen(0:mnx,0:mny), yen(0:mnx,0:mny),
     .        xzc(0:mnx,0:mny), yzc(0:mnx,0:mny),
     .        xeu(0:mnx,0:mny), yeu(0:mnx,0:mny),
     .        xzu(0:mnx,0:mny), yzu(0:mnx,0:mny),
     .        xev(0:mnx,0:mny), yev(0:mnx,0:mny),
     .        xzv(0:mnx,0:mny), yzv(0:mnx,0:mny)

      REAL    d(0:mnx,0:mny),  dn(0:mnx,0:mny)    !---> Density
      REAL    un(0:mnx,0:mny), vn(0:mnx,0:mny),
     .        us(0:mnx,0:mny), vs(0:mnx,0:mny)

C---: External procedures
      REAL    DMaxNorm
      REAL    difmax

C---: Local variables
      INTEGER m, i, j
      REAL    dus(0:mnx,0:mny), dvs(0:mnx,0:mny)
      REAL    dif(2)

C---: Double precision constants
      REAL   dZero
      parameter ( dZero = 0.00d+00 )

C---: Default return value (indicating that the quasilinearization
C---: iterations did not converge to the requested tolerance after
C---: mqiter iterations.)
      nAuxMomentum = -1

C---: Initialize 'star' time-level arrays
      do j=1, ny+1
        do i=1, nx+1
          us(i,j) = un(i,j)
          vs(i,j) = vn(i,j)
        end do
      end do

C---> Start Newton-Kantorovich iterations
      do m=1, mqiter

C---? Comment:  Enforcing velocity boundary conditions on us and vs
C---? at this point seems to make the q-l iterations more difficult
C---? to converge.
C---?         call VelBoundCond(nx, ny, nReg, nRegBrd,
C---?     .                     nMomBdTp, dBCVal,
C---?     .                     us, vs)

C------: Intead, only enforce velocity outflow conditions at every 
C------: q-l iteration
         call VelOutflowBCs(nx, ny,
     .                      nReg, nRegBrd, nMomBdTp,
     .                      dBCVal,
     .                      us, vs)

C------- Initialize delta-u and delta-v
         do j=0,ny+1
           do i=0,nx+1
             dus(i,j) = dZero
             dvs(i,j) = dZero
           end do
         end do

C------: Solve for delta-u and delta-v using Douglas-Gunn ADI 
         call XMomentum(nx, ny,
     .                  nReg, nRegBrd,
     .                  nRegType, nMomBdTp,
     .                  dk, re,
     .                  dPRporos, dPRporc1, dPRporc2,
     .                  rbn, rgn, rac, rbc, dju,
     .                  xec, yec, xzn, yzn,
     .                  xeu, yeu, xzu, yzu,
     .                  us, vs, un, vn,
     .                  dus)

         call YMomentum(nx, ny, 
     .                  nReg, nRegBrd,
     .                  nRegType, nMomBdTp,
     .                  dk, re, fr, 
     .                  dPRporos, dPRporc1, dPRporc2,
     .                  ran, rbn, rbc, rgc, djv,
     .                  xen, yen, xzc, yzc,
     .                  xev, yev, xzv, yzv,
     .                  d, dn,
     .                  us, vs, un, vn,
     .                  dvs)

C------: Update us and vs
         do j=1,ny
           do i=1,nx
             us(i,j) = us(i,j) + dus(i,j)
             vs(i,j) = vs(i,j) + dvs(i,j)
          end do
         end do

C------: Calculate max-norm of increments
         dif(1) = DMaxNorm(nx, ny, dus)
         dif(2) = DMaxNorm(nx, ny, dvs)

         difmax = max(dif(1),dif(2))

C------: Check convergence
         if(difmax.le.qtol) then
           nAuxMomentum = m
           return
         end if

      end do   !-----> End Q-L iterations

      return
      end

*----------------------------------------------------------------------*
*                              XMomentum                               *
*----------------------------------------------------------------------*

      subroutine XMomentum(nx, ny,
     .                     nReg, nRegBrd,
     .                     nRegType, nMomBdTp,
     .                     dk, re,
     .                     dPRporos, dPRporc1, dPRporc2,
     .                     rbn, rgn, rac, rbc, dju,
     .                     xec, yec, xzn, yzn,
     .                     xeu, yeu, xzu, yzu,
     .                     us, vs, un, vn,
     .                     dus)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn    ! Maximum number of total points
      parameter ( mn  = mnx*mny )

C---: Subroutine arguments
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

      REAL    dk, re

C---: Porous region dimensionless coefficients and local porosity
      REAL    dPRporos(mgri,mgrj), dPRporc1(mgri,mgrj), 
     .        dPRporc2(mgri,mgrj) 

      REAL    un(0:mnx,0:mny), vn(0:mnx,0:mny),
     .        us(0:mnx,0:mny), vs(0:mnx,0:mny),
     .        dus(0:mnx,0:mny)

      REAL    rbn(0:mnx,0:mny), rgn(0:mnx,0:mny),
     .        rac(0:mnx,0:mny), rbc(0:mnx,0:mny),
     .        dju(0:mnx,0:mny)

      REAL    xec(0:mnx,0:mny), yec(0:mnx,0:mny),
     .        xzn(0:mnx,0:mny), yzn(0:mnx,0:mny),
     .        xeu(0:mnx,0:mny), yeu(0:mnx,0:mny),
     .        xzu(0:mnx,0:mny), yzu(0:mnx,0:mny)


C---: Local variables
      INTEGER i, j, ind, ireg, jreg, iW, iE, jS, jN
      REAL    re1, dk2, rkj, dLocPoros

      REAL    cj1(0:mnx,0:mny),  cj2(0:mnx,0:mny), 
     .        c1s(0:mnx,0:mny),  c2s(0:mnx,0:mny),  
     .        c1n(0:mnx,0:mny),  c2n(0:mnx,0:mny)
      REAL    cps(0:mnx,0:mny),  cpn(0:mnx,0:mny),
     .        cpj(0:mnx,0:mny)
      REAL    cnvs(0:mnx,0:mny), difs(0:mnx,0:mny),
     .        cnvn(0:mnx,0:mny), difn(0:mnx,0:mny)

      REAL    a(3,mn), b(mn)

C---: Double precision constants
      REAL   dZero, dOne, dTwo, dHalf
      parameter ( dZero = 0.00d+00 ,
     .            dOne  = 1.00d+00 ,
     .            dTwo  = 2.00d+00 ,
     .            dHalf = 0.50d+00 ) 


      re1 = dOne/re
      dk2 = dk*dHalf

c      call ConvCoef(nx, ny, 1, 1,
c     .              xzn, xec, yzn, yec, 
c     .              us,  vs,  cj1, cj2)
      call ConvCoef(nx, ny, 4, 1,         !---> Upwinding
     .              xzu, xeu, yzu, yeu,
     .              us,  vs,  cj1, cj2)
      call ConvCoef(nx, ny, 1, 0,
     .              xzn, xec, yzn, yec,
     .              us,  vs,  c1s, c2s)
      call ConvCoef(nx, ny, 1, 0,
     .              xzn, xec, yzn, yec,
     .              un,  vn,  c1n, c2n)


C**** Needs revision
C---: In porous regions, divide convective terms by porosity
C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

          if(nRegType(ireg,jreg).eq.RM_POROUS) then  !---: If a porous region

            dLocPoros = dPRporos(ireg,jreg)

            do j=jS+1, jN
               do i=iW, iE
                 cj1(i,j) = cj1(i,j)/dLocPoros
                 cj2(i,j) = cj2(i,j)/dLocPoros
                 c1s(i,j) = c1s(i,j)/dLocPoros
                 c2s(i,j) = c2s(i,j)/dLocPoros
                 c1n(i,j) = c1n(i,j)/dLocPoros
                 c2n(i,j) = c2n(i,j)/dLocPoros
              end do
            end do

          end if   !---: If not a porous region, just go to the next

        end do   !---: End J-region sweep
      end do   !---: End I-region sweep


      call DConvU(nx, ny, c1s, c2s, us, cnvs)
      call DConvU(nx, ny, c1n, c2n, un, cnvn)
      call DDiffU(nx, ny, rac, rbc, rbn, rgn, us, difs)
      call DDiffU(nx, ny, rac, rbc, rbn, rgn, un, difn)

      call PorosCoef(nx, ny, 1, 1,
     .               nReg, nRegBrd, nRegType,
     .               dPRporos, dPRporc1, dPRporc2,
     .               us, vs, cpj)
      call PorosCoef(nx, ny, 1, 0,
     .               nReg, nRegBrd, nRegType,
     .               dPRporos, dPRporc1, dPRporc2,
     .               us, vs, cps)
      call PorosCoef(nx, ny, 1, 0,
     .               nReg, nRegBrd, nRegType,
     .               dPRporos, dPRporc1, dPRporc2,
     .               un, vn, cpn)


C---:                    +------------------+
C---:                    | First Split Step |
C---:                    +------------------+

      do j=2,ny
        do i=1,nx

          ind = (j-2)*nx+i
          rkj = dk2*dju(i,j)

C-------> Coefficients of the left-hand side tri-diagonal matrix

CC-------> No upwinding:
C          a(1,ind) = rkj*(-cj1(i,j)-re1*rac(i,j))
C          a(2,ind) = rkj*(cj1(i+1,j)-cj1(i,j)
C     .      +re1*(rac(i+1,j)+rac(i,j))) + dk2*cpj(i,j) + dOne
C          a(3,ind) = rkj*(cj1(i+1,j)-re1*rac(i+1,j))

C-------> Upwinding:
          if(cj1(i,j).ge.dZero) then
            a(1,ind) = rkj*(-cj1(i-1,j)-re1*rac(i,j))
            a(2,ind) = dOne + rkj*(cj1(i,j)
     .        +re1*(rac(i+1,j)+rac(i,j))) + dk2*cpj(i,j)
            a(3,ind) = rkj*(-re1*rac(i+1,j))
          else
            a(1,ind) = rkj*(-re1*rac(i,j))
            a(2,ind) = dOne +  rkj*(-cj1(i,j)
     .        +re1*(rac(i+1,j)+rac(i,j))) + dk2*cpj(i,j)
            a(3,ind) = rkj*(cj1(i+1,j)-re1*rac(i+1,j))
            end if

C-------> Coefficients of the right-hand side vector

          b(ind) = un(i,j) - us(i,j) + rkj*(-cnvs(i,j) - cnvn(i,j))
     .           + rkj*re1*(difs(i,j)+difn(i,j))
     .           - dk2*(cps(i,j)*us(i,j)+cpn(i,j)*un(i,j))

        end do
      end do


C---: Solve tridiagonal system via LU factorization.
C---: The returning vector is loaded onto RHS vector b.
      call AltTridLU(nx*(ny-1), a, b)


C--->                     +-------------------+
C--->                     | Second Split Step |
C--->                     +-------------------+

      do j=2,ny
        do i=1,nx

          ind = (j-2)*nx+i
          rkj = dk2*dju(i,j)

C-------> Coefficients of the left-hand side tri-diagonal matrix

CC-------> No upwinding:
C          a(1,ind) = rkj*(-cj2(i,j-1)-re1*rgn(i,j-1))
C          a(2,ind) = rkj*(cj2(i,j)-cj2(i,j-1)
C     .      +re1*(rgn(i,j)+rgn(i,j-1))) + dk2*cpj(i,j) + dOne
C          a(3,ind) = rkj*(cj2(i,j)-re1*rgn(i,j))

C-------> Upwinding:
          if(cj2(i,j).ge.dZero) then
            a(1,ind) = rkj*(-cj2(i,j-1)-re1*rgn(i,j-1))
            a(2,ind) = dOne + rkj*(cj2(i,j)
     .        +re1*(rgn(i,j)+rgn(i,j-1))) + dk2*cpj(i,j)
            a(3,ind) = rkj*(-re1*rgn(i,j))
          else
            a(1,ind) = rkj*(-re1*rgn(i,j-1))
            a(2,ind) = dOne + rkj*(-cj2(i,j)
     .        +re1*(rgn(i,j)+rgn(i,j-1))) + dk2*cpj(i,j)
            a(3,ind) = rkj*(cj2(i,j+1)-re1*rgn(i,j))
          end if

C-------> Note:  The coefficients of the right-hand side vector
C-------> on the second split step correspond to those of the
C-------> resulting vector from the first split step

        end do
      end do



C**** Account for blockage in tridiagonal matrix.
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

C---------: Interior of blockage
            do j=jS+1,jN
              do i=iW,iE
                ind = (j-2)*nx+i
                a(1,ind) = dZero
                a(2,ind) = dOne
                a(3,ind) = dZero
                b(ind)   = dZero
              end do
            end do

          end if   !---: If not a blockage

C-------: W E S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,WEST))
          case(BM_INTERN, BM_OUTLT1, BM_OUTLT2)  !---> Internal, outlets

          case(BM_WALL1, BM_WALL2, BM_INLET)   !---> Walls, inlet
            do j=jS+1, jN
              ind = (j-2)*nx+iW
              a(1,ind) = dZero
              a(2,ind) = dOne
              a(3,ind) = dZero
              b(ind)   = dZero
            end do

          case default
            print*, 'Wrong nBdTypeW flag in region ',ireg,',',jreg
            stop
          end select

C-------: E A S T   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,EAST))
          case(BM_INTERN, BM_OUTLT1, BM_OUTLT2)  !---> Internal, outlets

          case(BM_WALL1, BM_WALL2, BM_INLET)   !---> Walls, inlet
            do j=jS+1, jN
              ind = (j-2)*nx+iE
              a(1,ind) = dZero
              a(2,ind) = dOne
              a(3,ind) = dZero
              b(ind)   = dZero
            end do

          case default
            print*, 'Wrong nBdTypeE flag in region ',ireg,',',jreg
            stop
          end select

        end do   !---: End J-region sweep
      end do   !---: End I-region sweep


C---: Solve tridiagonal system via LU factorization.
C---: The returning vector is loaded onto RHS vector b.
      call AltTridLU(nx*(ny-1), a, b)


C---: Load results back into returning array
      do j=2,ny
        do i=1,nx
          ind = (j-2)*nx+i
          dus(i,j) = b(ind)
        end do
      end do


      return
      end

*----------------------------------------------------------------------*
*                              YMomentum                               *
*----------------------------------------------------------------------*

      subroutine YMomentum(nx, ny,
     .                     nReg, nRegBrd,
     .                     nRegType, nMomBdTp,
     .                     dk, re, fr,
     .                     dPRporos, dPRporc1, dPRporc2,
     .                     ran, rbn, rbc, rgc, djv,
     .                     xen, yen, xzc, yzc,
     .                     xev, yev, xzv, yzv,
     .                     d, dn, 
     .                     us, vs, un, vn,
     .                     dvs)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn    ! Maximum number of total points
      parameter ( mn  = mnx*mny )

C---: Subroutine arguments
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

      REAL    dk, re, fr

C---: Porous region dimensionless coefficients and local porosity
      REAL    dPRporos(mgri,mgrj), dPRporc1(mgri,mgrj), 
     .        dPRporc2(mgri,mgrj) 

      REAL    d(0:mnx,0:mny),  dn(0:mnx,0:mny),
     .        us(0:mnx,0:mny), vs(0:mnx,0:mny),
     .        un(0:mnx,0:mny), vn(0:mnx,0:mny),
     .        dvs(0:mnx,0:mny)

      REAL    ran(0:mnx,0:mny), rbn(0:mnx,0:mny),
     .        rbc(0:mnx,0:mny), rgc(0:mnx,0:mny),
     .        djv(0:mnx,0:mny)

      REAL    xen(0:mnx,0:mny), yen(0:mnx,0:mny),
     .        xzc(0:mnx,0:mny), yzc(0:mnx,0:mny),
     .        xev(0:mnx,0:mny), yev(0:mnx,0:mny),
     .        xzv(0:mnx,0:mny), yzv(0:mnx,0:mny)


C---: Local variables
      INTEGER i, j, ind, ireg, jreg, iW, iE, jS, jN
      REAL    re1, dk2, fr1, rkj, dLocPoros, buoy

      REAL    cj1(0:mnx,0:mny),  cj2(0:mnx,0:mny), 
     .        c1s(0:mnx,0:mny),  c2s(0:mnx,0:mny),  
     .        c1n(0:mnx,0:mny),  c2n(0:mnx,0:mny)
      REAL    cps(0:mnx,0:mny),  cpn(0:mnx,0:mny),
     .        cpj(0:mnx,0:mny)
      REAL    cnvs(0:mnx,0:mny), difs(0:mnx,0:mny),
     .        cnvn(0:mnx,0:mny), difn(0:mnx,0:mny)

      REAL    a(3,mn), b(mn)

C---: Double precision constants
      REAL   dZero, dOne, dTwo, dFour, dHalf
      parameter ( dZero = 0.00d+00 ,
     .            dOne  = 1.00d+00 ,
     .            dTwo  = 2.00d+00 ,
     .            dFour = 4.00d+00 ,
     .            dHalf = 0.50d+00 )


      re1 = dOne/re
      dk2 = dk*dHalf
      fr1 = dOne/fr

c      call ConvCoef(nx, ny, 2, 1,
c     .              xzc, xen, yzc, yen,
c     .              us,  vs,  cj1, cj2)
      call ConvCoef(nx, ny, 5, 1,         !---> Upwinding
     .              xzv, xev, yzv, yev,
     .              us,  vs,  cj1, cj2)
      call ConvCoef(nx, ny, 2, 0,
     .              xzc, xen, yzc, yen,
     .              us,  vs,  c1s, c2s)
      call ConvCoef(nx, ny, 2, 0,
     .              xzc, xen, yzc, yen,
     .              un,  vn,  c1n, c2n)


C**** Needs revision
C---: In porous regions, divide convective terms by porosity
C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

          if(nRegType(ireg,jreg).eq.RM_POROUS) then  !---: If a porous region

            dLocPoros = dPRporos(ireg,jreg)

            do j=jS, jN
               do i=iW+1, iE
                 cj1(i,j) = cj1(i,j)/dLocPoros
                 cj2(i,j) = cj2(i,j)/dLocPoros
                 c1s(i,j) = c1s(i,j)/dLocPoros
                 c2s(i,j) = c2s(i,j)/dLocPoros
                 c1n(i,j) = c1n(i,j)/dLocPoros
                 c2n(i,j) = c2n(i,j)/dLocPoros
              end do
            end do

          end if   !---: If not a porous region, just go to the next

        end do   !---: End J-region sweep
      end do   !---: End I-region sweep


      call DConvV(nx, ny, c1s,c2s,vs,cnvs)
      call DConvV(nx, ny, c1n,c2n,vn,cnvn)
      call DDiffV(nx, ny, ran, rbc, rbn, rgc, vs, difs)
      call DDiffV(nx, ny, ran, rbc, rbn, rgc, vn, difn)

      call PorosCoef(nx, ny, 2, 1,
     .               nReg, nRegBrd, nRegType,
     .               dPRporos, dPRporc1, dPRporc2,
     .               us, vs, cpj)
      call PorosCoef(nx, ny, 2, 0,
     .               nReg, nRegBrd, nRegType,
     .               dPRporos, dPRporc1, dPRporc2,
     .               us, vs, cps)
      call PorosCoef(nx, ny, 2, 0,
     .               nReg, nRegBrd, nRegType,
     .               dPRporos, dPRporc1, dPRporc2,
     .               un, vn, cpn)


C---:                    +------------------+
C---:                    | First Split Step |
C---:                    +------------------+

      do j=1,ny
        do i=2,nx
          ind = (j-1)*(nx-1)+i-1
          rkj = dk2*djv(i,j)

C-------> Coefficients of the left-hand side tri-diagonal matrix

CC-------> No upwinding:
C          a(1,ind) = rkj*(-cj1(i-1,j)-re1*ran(i-1,j))
C          a(2,ind) = rkj*(cj1(i,j)-cj1(i-1,j)
C     .      +re1*(ran(i,j)+ran(i-1,j))) + dk2*cpj(i,j) + dOne
C          a(3,ind) = rkj*(cj1(i,j)-re1*ran(i,j))

C-------> Upwinding:
          if(cj1(i,j).ge.dZero) then
            a(1,ind) = rkj*(-cj1(i-1,j)-re1*ran(i-1,j))
            a(2,ind) = rkj*(cj1(i,j)
     .        +re1*(ran(i,j)+ran(i-1,j))) + dk2*cpj(i,j) + dOne
            a(3,ind) = rkj*(-re1*ran(i,j))
          else
            a(1,ind) = rkj*(-re1*ran(i-1,j))
            a(2,ind) = rkj*(-cj1(i,j)
     .        +re1*(ran(i,j)+ran(i-1,j))) + dk2*cpj(i,j) + dOne
            a(3,ind) = rkj*(cj1(i+1,j)-re1*ran(i,j))
          end if

C-------> Coefficients of the righ-hand side vector

          buoy = dk*(d(i,j+1)+d(i,j)+dn(i,j+1)+dn(i,j))/(dFour*fr)

          b(ind) = vn(i,j) - vs(i,j) + rkj*(-cnvs(i,j) - cnvn(i,j))
     .           + rkj*re1*(difs(i,j) + difn(i,j))
     .           - dk2*(cps(i,j)*vs(i,j)+cpn(i,j)*vn(i,j))
     .           - buoy

        end do
      end do


C---: Solve tridiagonal system via LU factorization.
C---: The returning vector is loaded onto RHS vector b.
      call AltTridLU((nx-1)*ny, a, b)


C---:                    +-------------------+
C---:                    | Second Split Step |
C---:                    +-------------------+

      do j=1,ny
        do i=2,nx
          ind = (j-1)*(nx-1)+i-1
          rkj = dk2*djv(i,j)

C-------> Coefficients of the left-hand side tri-diagonal matrix

CC-------> No upwinding:
C          a(1,ind) = rkj*(-cj2(i,j)-re1*rgc(i,j))
C          a(2,ind) = rkj*(cj2(i,j+1)-cj2(i,j)
C     .      +re1*(rgc(i,j+1)+rgc(i,j))) + dk2*cpj(i,j) + dOne
C          a(3,ind) = rkj*(cj2(i,j+1)-re1*rgc(i,j+1))

C-------> Upwinding:
          if(cj2(i,j).ge.dZero) then
            a(1,ind) = rkj*(-cj2(i,j-1)-re1*rgc(i,j))
            a(2,ind) = rkj*(cj2(i,j)
     .        +re1*(rgc(i,j+1)+rgc(i,j))) + dk2*cpj(i,j) + dOne
            a(3,ind) = rkj*(-re1*rgc(i,j+1))
          else
            a(1,ind) = rkj*(-re1*rgc(i,j))
            a(2,ind) = rkj*(-cj2(i,j)
     .        +re1*(rgc(i,j+1)+rgc(i,j))) + dk2*cpj(i,j) + dOne
            a(3,ind) = rkj*(cj2(i,j+1)-re1*rgc(i,j+1))
          end if

C-------> Note:  The coefficients of the right-hand side vector
C-------> on the second split step correspond to those of the
C-------> resulting vector from the first split step

        end do
      end do



C**** Account for blockage in tridiagonal matrix.
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
              do i=iW+1,iE
                ind = (j-1)*(nx-1)+i-1
                a(1,ind) = dZero
                a(2,ind) = dOne
                a(3,ind) = dZero
                b(ind)   = dZero
              end do
            end do

          end if   !---: If not a blockage

C-------: S O U T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,SOUTH))
          case(BM_INTERN, BM_OUTLT1, BM_OUTLT2)  !---> Internal, outlets 

          case(BM_WALL1, BM_WALL2, BM_INLET)   !---> Walls, inlet
            do i=iW+1, iE
              ind = (jS-1)*(nx-1)+i-1
              a(1,ind) = dZero
              a(2,ind) = dOne
              a(3,ind) = dZero
              b(ind)   = dZero
            end do

          case default
            print*, 'Wrong nBdTypeS flag in region ',ireg,',',jreg
            stop
          end select

C-------: N O R T H   Boundary :-------]
          select case(nMomBdTp(ireg,jreg,NORTH))
          case(BM_INTERN, BM_OUTLT1, BM_OUTLT2)  !---> Internal, outlets 

          case(BM_WALL1, BM_WALL2, BM_INLET)   !---> Walls, inlet
            do i=iW+1, iE
              ind = (jN-1)*(nx-1)+i-1
              a(1,ind) = dZero
              a(2,ind) = dOne
              a(3,ind) = dZero
              b(ind)   = dZero
            end do

          case default
            print*, 'Wrong nBdTypeN flag in region ',ireg,',',jreg
            stop
          end select

        end do   !---: End J-region sweep
      end do   !---: End I-region sweep


C---: Solve tridiagonal system via LU factorization.
C---: The returning vector is loaded onto RHS vector b.
      call AltTridLU((nx-1)*ny, a, b)

C---: Load results back into returning array
      do j=1,ny
        do i=2,nx
          ind = (j-1)*(nx-1)+i-1
          dvs(i,j) = b(ind)
        end do
      end do


      return
      end

*----------------------------------------------------------------------*
*                              ConvCoef                                *
*----------------------------------------------------------------------*
*
* ncomp -> Component (variable) flag:
*   1: x-momentum equation
*   2: y-momentum equation
*   3: thermal energy equation
*   4: x-momentum equation, upwinding
*   5: y-momentum equation, upwinding
*   6: thermal energy equation, upwinding
*
* njacob -> Jacobian flag (LHS of auxiliary momentum equations)
*
* Metric coefficients to be passed depending on component flag:
*   nloc     xzi   xet   yzi   yet
*  ------   ----- ----- ----- -----
*     1      xzn   xec   yzn   yec
*     2      xzc   xen   yzc   yen
*     3      xzv   xeu   yzv   yeu
*     4      xzu   xeu   yzu   yeu
*     5      xzv   xev   yzv   yev
*     6      xzc   xec   yzc   yec
*
      subroutine ConvCoef(nx, ny, ncomp, njacob,
     .                    xzi, xet, yzi, yet,
     .                    u, v, cc1, cc2)

      implicit none

C---: Configurable parameters
      include "config.f"

C---: Arguments
      INTEGER nx, ny
      INTEGER ncomp, njacob

      REAL   xzi(0:mnx,0:mny),  xet(0:mnx,0:mny),
     .       yzi(0:mnx,0:mny),  yet(0:mnx,0:mny)

      REAL   u(0:mnx,0:mny),   v(0:mnx,0:mny),
     .       cc1(0:mnx,0:mny), cc2(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j
      REAL    djac

C---: Double precision constants
      REAL   dZero, dOne, dTwo, dFour, dHalf
      parameter ( dZero = 0.00d+00 ,
     .            dOne  = 1.00d+00 ,
     .            dTwo  = 2.00d+00 ,
     .            dFour = 4.00d+00 ,
     .            dHalf = 0.50d+00 ) 

      djac = dOne                   !---> Convective coefficient
      if(njacob.eq.1) djac = dTwo   !---> LHS jacobian coefficient

      select case(ncomp)

C---> Convective coefficients for x-momentum equation
      case(1)
        do j=1, ny+1
          do i=1, nx+1
            cc1(i,j) = (djac*yet(i,j)*(u(i,j)+u(i-1,j))
     .                -xet(i,j)*(v(i,j)+v(i,j-1)))*dHalf
          end do
        end do
        do j=1,ny
          do i=1,nx
            cc2(i,j) = (xzi(i,j)*(v(i+1,j)+v(i,j))
     .                 -djac*yzi(i,j)*(u(i,j+1)+u(i,j)))*dHalf
          end do
        end do

C---> Convective coefficients for y-momentum equation
      case(2)
        do j=1,ny
          do i=1,nx
            cc1(i,j) = (yet(i,j)*(u(i,j+1)+u(i,j))
     .                 -djac*xet(i,j)*(v(i+1,j)+v(i,j)))*dHalf
          end do
        end do
        do j=1,ny+1
          do i=1,nx+1
            cc2(i,j) = (djac*xzi(i,j)*(v(i,j)+v(i,j-1))
     .                 -yzi(i,j)*(u(i,j)+u(i-1,j)))*dHalf
          end do
        end do

C---> Convective coefficients for thermal energy equation
      case(3)
        do j=1,ny
          do i=1,nx
            cc1(i,j)=(yet(i,j)*u(i,j)-xet(i,j)*(v(i+1,j)+v(i,j)
     .        +v(i+1,j-1)+v(i,j-1))/dFour)*dHalf
            cc2(i,j)=(xzi(i,j)*v(i,j)-yzi(i,j)*(u(i,j+1)+u(i-1,j+1)
     .        +u(i,j)+u(i-1,j))/dFour)*dHalf
          end do
        end do

C---> Upwinding LHS x-momentum equation
      case(4)
        do j=1,ny
          do i=0,nx
            cc1(i,j)=djac*yet(i,j)*u(i,j)-xet(i,j)*(v(i+1,j)+v(i,j)
     .              + v(i+1,j-1)+v(i,j-1))/dFour
            cc2(i,j)=xzi(i,j)*(v(i+1,j)+v(i,j)+v(i+1,j-1)
     .              + v(i,j-1))/dFour - djac*yzi(i,j)*u(i,j)
          end do
        end do

C---> Upwinding LHS y-momentum equation
      case(5)
        do j=0,ny
          do i=1,nx
            cc1(i,j)=yet(i,j)*(u(i,j+1)+u(i-1,j+1)+u(i,j)
     .              +u(i-1,j))/dFour - djac*xet(i,j)*v(i,j)
            cc2(i,j)=djac*xzi(i,j)*v(i,j)-yzi(i,j)*(u(i,j+1)
     .              +u(i-1,j+1)+u(i,j)+u(i-1,j))/dFour
          end do
        end do

C---> Upwinding LHS thermal energy equation
      case(6)
        do j=1,ny+1
          do i=1,nx+1
          cc1(i,j) = (yet(i,j)*(u(i,j)+u(i-1,j))
     .              -xet(i,j)*(v(i,j)+v(i,j-1)))*dHalf
          cc2(i,j) = (xzi(i,j)*(v(i,j)+v(i,j-1))
     .               -yzi(i,j)*(u(i,j)+u(i-1,j)))*dHalf
          end do
        end do

      case default
        print*, 'Error: Wrong component flag passed to ConvCoef: ', 
     .    ncomp 
        stop
      end select

      return
      end

*----------------------------------------------------------------------*
*                               DConvU                                 *
*----------------------------------------------------------------------*

      subroutine DConvU(nx, ny, c1, c2, u, c)

      implicit none

      include "config.f"

      INTEGER nx, ny

      REAL    u(0:mnx,0:mny),  c(0:mnx,0:mny),
     .        c1(0:mnx,0:mny), c2(0:mnx,0:mny)

      INTEGER i, j

      do j=2, ny 
        do i=1, nx
          c(i,j) = -c2(i,j-1)*u(i,j-1)-c1(i,j)*u(i-1,j)
     .      +(c1(i+1,j)-c1(i,j)+c2(i,j)-c2(i,j-1))*u(i,j)
     .      +c1(i+1,j)*u(i+1,j)+c2(i,j)*u(i,j+1)
        end do
      end do

      return
      end

*----------------------------------------------------------------------*
*                               DDiffU                                 *
*----------------------------------------------------------------------*

      subroutine DDiffU(nx, ny, ac, bc, bn, gn, u, d)

      implicit none

      include "config.f"

      INTEGER nx, ny
      REAL    ac(0:mnx,0:mny), bc(0:mnx,0:mny), 
     .        bn(0:mnx,0:mny), gn(0:mnx,0:mny)
      REAL    u(0:mnx,0:mny)
      REAL    d(0:mnx,0:mny)

      INTEGER i, j
      REAL    s1, s2

      do j=2, ny
        do i=1, nx
          s1 =ac(i+1,j)*(u(i+1,j)-u(i,j))
     .       -ac(i,j)*(u(i,j)-u(i-1,j))
     .       +bc(i+1,j)*(u(i+1,j+1)+u(i,j+1)-u(i+1,j-1)-u(i,j-1))
     .       -bc(i,j)*(u(i,j+1)+u(i-1,j+1)-u(i,j-1)-u(i-1,j-1))
          s2 =bn(i,j)*(u(i+1,j+1)+u(i+1,j)-u(i-1,j+1)-u(i-1,j))
     .       -bn(i,j-1)*(u(i+1,j)+u(i+1,j-1)-u(i-1,j)-u(i-1,j-1))
     .       +gn(i,j)*(u(i,j+1)-u(i,j))
     .       -gn(i,j-1)*(u(i,j)-u(i,j-1))
         d(i,j) = s1 + s2
        end do
      end do

      return
      end

*----------------------------------------------------------------------*
*                               DConvV                                 *
*----------------------------------------------------------------------*

      subroutine DConvV(nx, ny, c1, c2, v, c)

      implicit none

      include "config.f"

      INTEGER nx, ny

      REAL    v(0:mnx,0:mny),  c(0:mnx,0:mny),
     .        c1(0:mnx,0:mny), c2(0:mnx,0:mny)

      INTEGER i, j

      do j=1, ny
        do i=2, nx
          c(i,j) = -c2(i,j)*v(i,j-1)-c1(i-1,j)*v(i-1,j)
     .      +(c1(i,j)-c1(i-1,j)+c2(i,j+1)-c2(i,j))*v(i,j)
     .      +c1(i,j)*v(i+1,j)+c2(i,j+1)*v(i,j+1)
        end do
      end do

      return
      end

*----------------------------------------------------------------------*
*                               DDiffV                                 *
*----------------------------------------------------------------------*

      subroutine DDiffV(nx, ny, an, bc, bn, gc, v, d)

      implicit none

      include "config.f"

      INTEGER nx, ny
      REAL    an(0:mnx,0:mny), bc(0:mnx,0:mny), 
     .        bn(0:mnx,0:mny), gc(0:mnx,0:mny)
      REAL    v(0:mnx,0:mny)
      REAL    d(0:mnx,0:mny)

      INTEGER i, j
      REAL    s1, s2

      do j=1, ny
        do i=2, nx
          s1 =an(i,j)*(v(i+1,j)-v(i,j))
     .       -an(i-1,j)*(v(i,j)-v(i-1,j))
     .       +bn(i,j)*(v(i+1,j+1)+v(i,j+1)-v(i+1,j-1)-v(i,j-1)) 
     .       -bn(i-1,j)*(v(i,j+1)+v(i-1,j+1)-v(i,j-1) -v(i-1,j-1))
          s2 =bc(i,j+1)*(v(i+1,j+1)+v(i+1,j)-v(i-1,j+1)-v(i-1,j)) 
     .       -bc(i,j)*(v(i+1,j)+v(i+1,j-1)-v(i-1,j)-v(i-1,j-1)) 
     .       +gc(i,j+1)*(v(i,j+1)-v(i,j))
     .       -gc(i,j)*(v(i,j)-v(i,j-1))
         d(i,j) = s1 + s2
        end do
      end do

      return
      end

*----------------------------------------------------------------------*
*                              PorosCoef                               *
*----------------------------------------------------------------------*

      subroutine PorosCoef(nx, ny, ncomp, njacob,
     .                     nReg, nRegBrd, nRegType,
     .                     dPRporos, dPRporc1, dPRporc2,
     .                     u, v, cp)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, ncomp, njacob

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Porous region dimensionless coefficients and local porosity
      REAL    dPRporos(mgri,mgrj), dPRporc1(mgri,mgrj), 
     .        dPRporc2(mgri,mgrj) 

      REAL u(0:mnx,0:mny), v(0:mnx,0:mny), cp(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, ireg, jreg, iW, iE, jS, jN
      REAL    porc1, porc2, unorm, unrm1

C---: Double precision constants
      REAL   dZero, dOne, dTwo, dFour, dHalf
      parameter ( dZero = 0.00d+00 ,
     .            dOne  = 1.00d+00 ,
     .            dTwo  = 2.00d+00 ,
     .            dFour = 4.00d+00 ,
     .            dHalf = 0.50d+00 )


C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

          if(nRegType(ireg,jreg).ne.RM_POROUS) then

            do j=jS, jN
              do i=iW, iE
                cp(i,j) = dZero
              end do
            end do

!---------: If not a porous region try the next one
            cycle

          end if


C***: Temporary
C***: The following is done just to keep the same notation as in
C***: previous versions of this subroutine (for debugging purposes)

          porc1 = dPRporc1(ireg,jreg)
          porc2 = dPRporc2(ireg,jreg)

          select case(ncomp)
          case(1)   !---> X-momentum equation

            do j=jS+1, jN
              do i=iW, iE

                unorm = dsqrt(u(i,j)**2 + (v(i,j) + v(i+1,j)
     .            +v(i,j-1) + v(i+1,j-1)/dFour)**2)
                unrm1 = dZero
                if(unorm.gt.1.d-8) unrm1 = (u(i,j)**2)/unorm
                if(njacob.eq.1) unorm = unrm1 + unorm
                cp(i,j) = porc1 + porc2*unorm

              end do
            end do

          case(2)   !---> Y-momentum equation

            do j=jS, jN
              do i=iW+1, iE
                unorm = dsqrt(((u(i-1,j+1) + u(i,j+1) + u(i-1,j)
     .            +u(i,j))/dFour)**2 + v(i,j)**2)
                unrm1 = dZero
                if(unorm.gt.1.d-8) unrm1 = (v(i,j)**2)/unorm
                if(njacob.eq.1) unorm = unrm1 + unorm
                cp(i,j) = porc1 + porc2*unorm
              end do
            end do

          case default
            print*, 'Error: Wrong ncomp flag passed to PorosCoef: ', 
     .        ncomp 
            stop
          end select

        end do   !---> End horizontal sweep
      end do   !---> End vertical sweep

      return
      end

*----------------------------------------------------------------------*
*                               TridLU                                 *
*----------------------------------------------------------------------*
*
* Tridiagonal LU decomposition (Thomas algorithm)
*
      subroutine TridLU (nmax, n, a, b, x)

      implicit none

C---: Arguments
      INTEGER nmax, n
      REAL a(nmax,3), b(nmax), x(nmax)

C---: Local variables
      INTEGER i

C---: Construction of L and U from elements of a
      a(1,3) = a(1,3)/a(1,2)
      do i=2, n-1
        a(i,2) = a(i,2) - (a(i,1)*a(i-1,3))
        a(i,3) = a(i,3)/a(i,2)
      end do
      a(n,2) = a(n,2) - (a(n,1)*a(n-1,3))
			
C---: Forward substitution (Solve Ly=b)
      b(1) = b(1)/a(1,2)
      do i=2, n
        b(i) = (b(i) - a(i,1)*b(i-1))/a(i,2)
      end do

C---: Backward substitution (Solve Ux=y)
      x(n) = b(n)
      do i=(n-1), 1, -1
        x(i) = b(i) - (a(i,3)*x(i+1))
      end do
			
      return
      end


*----------------------------------------------------------------------*
*                              AltTridLU                               *
*----------------------------------------------------------------------*
*
* Alternate Tridiagonal LU decomposition.
*
c      subroutine AltTridLU(n, a, b)
c      implicit none
cC---: Arguments
c      INTEGER n
c      REAL    a(3,n), b(n)
cC---: Local variables
c      INTEGER i
cC---: Construction of L and U from elements of a
c      a(3,1) = a(3,1)/a(2,2)
c      do i=2, n-1
c        a(2,i) = a(2,i) - (a(1,i)*a(3,i-1))
c        a(3,i) = a(3,i)/a(2,i)
c      end do
c      a(2,n) = a(2,n) - (a(1,n)*a(3,n-1))
cC---: Forward substitution (Solve Ly=b)
c      b(1) = b(1)/a(2,1)
c      do i=2, n
c        b(i) = (b(i) - a(1,i)*b(i-1))/a(2,i)
c      end do
cC---: Backward substitution (Solve Ux=y)
c      do i=(n-1), 1, -1
c        b(i) = b(i) - (a(3,i)*b(i+1))
c      end do
c      return
c      end

*----------------------------------------------------------------------*
*                              AltTridLU                               *
*----------------------------------------------------------------------*
*
* Alternate Tridiagonal LU decomposition.
*
      subroutine AltTridLU(n, a, b)

      implicit none

C---: Arguments
      INTEGER n
      REAL    a(3,n), b(n)

C---: Local variables
      INTEGER i

C---: Construction of L and U from elements of a
      a(3,1) = a(3,1)/a(2,2)
      b(1) = b(1)/a(2,1)

      do i=2, n-1
        a(2,i) = a(2,i) - (a(1,i)*a(3,i-1))
        a(3,i) = a(3,i)/a(2,i)

C---: Forward substitution (Solve Ly=b)
        b(i) = (b(i) - a(1,i)*b(i-1))/a(2,i)
      end do

      a(2,n) = a(2,n) - (a(1,n)*a(3,n-1))
      b(n) = (b(n) - a(1,n)*b(n-1))/a(2,n)

C---: Backward substitution (Solve Ux=y)
      do i=(n-1), 1, -1
        b(i) = b(i) - (a(3,i)*b(i+1))
      end do
			
      return
      end

*-----------------------|---|---|---V---|---|---|----------------------*
