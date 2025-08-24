*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                              thermal.F                               *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* ThermEnergy - Thermal energy equation (solve for temperature)
* EqState     - Ideal gas equation of state
*
*----------------------------------------------------------------------*

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"

*----------------------------------------------------------------------*
*                             ThermEnergy                              *
*----------------------------------------------------------------------*

      subroutine ThermEnergy(nx, ny,
     .                       nReg, nRegBrd,
     .                       nTRgType, nTemBdTp,
     .                       dk, pe,
     .                       dTRgVal, dHGSTval, dBCVal,
     .                       rau, rbu, rbv, rgv, djc,
     .                       xeu, yeu, xzv, yzv,
     .                       xec, yec, xzc, yzc,
     .                       un, vn, u, v, tn, t)

      implicit none

C---: Configurable parameters
      include "config.f"

C---> Non-configurable parameters
      INTEGER mn    ! Maximum number of total points
      parameter ( mn  = mnx*mny )

C---: Arguments
      INTEGER nx, ny

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

      REAL    dk, pe

      REAL    u(0:mnx,0:mny),  v(0:mnx,0:mny),  t(0:mnx,0:mny),
     .        un(0:mnx,0:mny), vn(0:mnx,0:mny), tn(0:mnx,0:mny)

      REAL   rau(0:mnx,0:mny), rbu(0:mnx,0:mny),
     .       rbv(0:mnx,0:mny), rgv(0:mnx,0:mny),
     .       djc(0:mnx,0:mny)

      REAL    xeu(0:mnx,0:mny), yeu(0:mnx,0:mny),
     .        xzv(0:mnx,0:mny), yzv(0:mnx,0:mny),
     .        xec(0:mnx,0:mny), yec(0:mnx,0:mny),
     .        xzc(0:mnx,0:mny), yzc(0:mnx,0:mny)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Temperature Boundary type
      INTEGER nTemBdTp(mgri,mgrj,*)

C---: Temperature region values
      REAL    dTRgVal(mgri,mgrj)

C---: Dimensionless heat generation source term values
      REAL    dHGSTval(mgri,mgrj)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL dBCVal(mgri,mgrj,4,*)


C---: Local variables
      INTEGER i, j, ind, ireg, jreg, iW, iE, jS, jN
      REAL    pe1, rkj, c, d
      REAL   cu1(0:mnx,0:mny), cun(0:mnx,0:mny),
     .       cv1(0:mnx,0:mny), cvn(0:mnx,0:mny), 
     .       s(0:mnx,0:mny)

      REAL   a(3,mn), b(mn)

C---: Commonly used double precision constants
      REAL   dZero, dOne, dTwo, dHalf
      parameter ( dZero = 0.00d+00 ,
     .            dOne  = 1.00d+00 ,
     .            dTwo  = 2.00d+00 ,
     .            dHalf = 0.50d+00 ) 

      pe1 = dOne/pe

C---: Apply boundary conditions
      call TempBoundCond(nx, ny,
     .                   nReg, nRegBrd,
     .                   nTRgType, nTemBdTp,
     .                   dTRgVal, dBCVal,
     .                   t)

C---: Convective Coefficients
c      call ConvCoef(nx, ny, 3, 0,         !---> No upwinding
c     .              xzv, xeu, yzv, yeu,
c     .              u,   v,   cu1, cv1)
      call ConvCoef(nx, ny, 6, 0,          !---> Upwinding
     .             xzc, xec, yzc, yec,
     .             u,   v,   cu1, cv1)
      call ConvCoef(nx, ny, 3, 0,
     .             xzv, xeu, yzv, yeu,
     .             un,   vn, cun, cvn)


C---: Heat generation term
      do jreg=1, nReg(_J_)   !--- For all regions in J direction
        do ireg=1, nReg(_I_)   !--- For all regions in I direction

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

          if(nTRgType(ireg,jreg).eq.RT_HEATGN) then

C----------: Assign heat generation value to all points in reg.
             do j=jS+1, jN
               do i=iW+1,iE
                 s(i,j) = dHGSTval(ireg,jreg)
               end do
             end do
          else   !---: If heat geration not specified
             do j=jS+1, jN
               do i=iW+1,iE
                 s(i,j) = dZero
               end do
             end do

          end if

        end do   !---: End J-region sweep
      end do   !---: End I-region sweep


C-----------------------
C---: First Split Step :
C-----------------------
      do j=2, ny
        do i=2, nx
          ind = (j-2)*(nx-1)+i-1
          rkj = dk * djc(i,j) * dHalf

C--------> Coefficients of the Left-hand side matrix

CC--------: No upwinding
C          a(1,ind) = rkj*(-cu1(i-1,j)-pe1*rau(i-1,j))
C          a(2,ind) = dOne + rkj*(cu1(i,j)-cu1(i-1,j)
C     .             + pe1*(rau(i,j)+rau(i-1,j)))
C          a(3,ind) = rkj*(cu1(i,j)-pe1*rau(i,j))

C-------: Upwinding
          if(cu1(i,j).ge.dZero) then
            a(1,ind) = rkj*(-cu1(i-1,j)-pe1*rau(i-1,j))
            a(2,ind) = dOne + rkj*(cu1(i,j)
     .               + pe1*(rau(i,j)+rau(i-1,j)))
            a(3,ind) = rkj*(-pe1*rau(i,j))
          else 
            a(1,ind) = rkj*(-pe1*rau(i-1,j))
            a(2,ind) = dOne + rkj*(-cu1(i,j)
     .               + pe1*(rau(i,j)+rau(i-1,j)))
            a(3,ind) = rkj*(cu1(i+1,j)-pe1*rau(i,j))
          end if

C-------> Coefficients of the right-hand side vector
          c = -cvn(i,j-1)*tn(i,j-1)-cun(i-1,j)*tn(i-1,j)
     .      +(cun(i,j)-cun(i-1,j)+cvn(i,j)-cvn(i,j-1))*tn(i,j)
     .      +cun(i,j)*tn(i+1,j)+cvn(i,j)*tn(i,j+1)
          d = rau(i,j)*(tn(i+1,j)-tn(i,j))
     .      -rau(i-1,j)*(tn(i,j)-tn(i-1,j))
     .      +rbu(i,j)*(tn(i+1,j+1)+tn(i,j+1)-tn(i+1,j-1)-tn(i,j-1))
     .      -rbu(i-1,j)*(tn(i,j+1)+tn(i-1,j+1)-tn(i,j-1)-tn(i-1,j-1))
     .      +rbv(i,j)*(tn(i+1,j+1)+tn(i+1,j)-tn(i-1,j+1)-tn(i-1,j))
     .      -rbv(i,j-1)*(tn(i+1,j)+tn(i+1,j-1)-tn(i-1,j)-tn(i-1,j-1))
     .      +rgv(i,j)*(tn(i,j+1)-tn(i,j))
     .      -rgv(i,j-1)*(tn(i,j)-tn(i,j-1))

          b(ind) = dk*s(i,j)*dHalf + dTwo*rkj*(-c + pe1*d)

        end do
      end do


C---: Solve the tridiagonal system by LU factorization
      call AltTridLU((nx-1)*(ny-1),a,b)

C------------------------
C---: Second Split Step :
C------------------------
      do j=2, ny
        do i=2, nx
          ind = (j-2)*(nx-1)+i-1
          rkj = dk*djc(i,j)*dHalf

C--------> Coefficients of the Left-hand side matrix

CC--------: No upwinding
C          a(1,ind) = rkj*(-cv1(i,j-1)-pe1*rgv(i,j-1))
C          a(2,ind) = dOne + rkj*(cv1(i,j)-cv1(i,j-1)
C     .             +pe1*(rgv(i,j)+rgv(i,j-1)))
C          a(3,ind) = rkj*(cv1(i,j)-pe1*rgv(i,j))

C-------: Upwinding
          if(cv1(i,j).ge.dZero) then
            a(1,ind) = rkj*(-cv1(i,j-1)-pe1*rgv(i,j-1))
            a(2,ind) = dOne + rkj*(cv1(i,j)
     .               + pe1*(rgv(i,j)+rgv(i,j-1)))
            a(3,ind) = rkj*(-pe1*rgv(i,j))
          else
            a(1,ind) = rkj*(-pe1*rgv(i,j-1))
            a(2,ind) = dOne + rkj*(-cv1(i,j)
     .               + pe1*(rgv(i,j)+rgv(i,j-1)))
            a(3,ind) = rkj*(cv1(i,j+1)-pe1*rgv(i,j))
          end if

C-------: Note that the Right-hand vector coefficients are those
C-------: of the returning vector from previous step

        end do
      end do


C---: Interior of fixed-temperature regions
      do jreg=1, nReg(_J_)   !--- For all regions in J direction
        do ireg=1, nReg(_I_)   !--- For all regions in I direction

          if(nTRgType(ireg,jreg).eq.RT_TEMPER) then

C---------: Extract boundary indeces
            iW = nRegBrd(ireg,jreg,WEST)
            iE = nRegBrd(ireg,jreg,EAST)
            jS = nRegBrd(ireg,jreg,SOUTH)
            jN = nRegBrd(ireg,jreg,NORTH)

            do j=jS+1, jN
              do i=iW+1, iE
                ind = (j-2)*(nx-1)+i-1
                a(1,ind) = dZero
                a(2,ind) = dOne
                a(3,ind) = dZero
                b(ind)   = dZero

              end do
            end do

          end if

        end do   !---: End J-region sweep
      end do   !---: End I-region sweep


C---: Solve the tridiagonal system by LU factorization
      call AltTridLU((nx-1)*(ny-1),a,b)


C---: Update temperature array
      do j=2, ny
        do i=2, nx
          ind = (j-2)*(nx-1)+i-1
          t(i,j) = t(i,j) + b(ind)
        end do
      end do


      return
      end

*----------------------------------------------------------------------*
*                               EqState                                *
*----------------------------------------------------------------------*

      subroutine EqState(nx, ny,
     .                   uref, densref, tmax, tref, rconst, 
     .                   p, t, den)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny
      REAL    uref, densref, tmax, tref, rconst
      REAL    p(0:mnx,0:mny), t(0:mnx,0:mny), den(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j
      REAL    pref, c1, c2

C---: Double precision constants
      REAL    dZero, dOne
      parameter ( dZero = 0.00d+00 ,
     .            dOne  = 1.00d+00 )


      pref = densref * rconst * tref

c***note: needs revision
c      do j=1,ny
c        do i=1,nx

      do j=2, ny
        do i=2, nx

          c1 = p(i,j) * densref * uref**2 + pref
          c2 = densref*rconst*(t(i,j)*(tmax-tref)+tref)

          den(i,j) = c1/c2 - dOne

          if(dabs(den(i,j)).lt.1.d-10) den(i,j) = dZero

        end do
      end do

      return
      end

*-----------------------|---|---|---V---|---|---|----------------------*
