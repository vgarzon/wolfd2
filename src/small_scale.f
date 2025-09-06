*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                             smallScale.F                             *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* SmallScale  - Small scale components of solution vector
* (datanh)    - Hyperbolic arctangent
*
*----------------------------------------------------------------------*

/*
  Only declare the ATD subroutine if the preprocessor flag _SMALLSCL_
  has been declared 
*/

#ifdef _SMALLSCL_

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"

*----------------------------------------------------------------------*
*                              SmallScale                              *
*----------------------------------------------------------------------*

      subroutine SmallScale(nx, ny, initflg, nthermen,
     .                      lCartesGrid,
     .                      nReg, nRegBrd,
     .                      nRegType, nTRgType, nMomBdTp, nTemBdTp,
     .                      nPpeSolver, msorit,
     .                      dlref, uref,tref, tmax,
     .                      dka, re, pe,
     .                      sortol, sorrel,
     .                      fp,
     .                      cu0, TsCoef, HsCoef, TemCoef,
     .                      bnumc, rmax, rlc,
     .                      dTRgVal, dBCVal,
     .                      rau, rbu, rbv, rgv,
     .                      dju, djv, djc,
     .                      xeu, yeu, xzv, yzv,
     .                      xzu, yzu, xev, yev,
     .                      xec, yec, xzc, yzc,
     .                      u1, v1, t1, 
     .                      uss, vss, pss, tss)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, initflg, nthermen

C---: Cartesian grid flag
      LOGICAL lCartesGrid

C---: Number of regions in each direction
      INTEGER nReg(*)

C---: Location of region borders
      INTEGER nRegBrd(mgri,mgrj,*)

C---: Momentum region types
      INTEGER nRegType(mgri,mgrj)

C---: Region types: Thermal energy
      INTEGER nTRgType(mgri,mgrj)

C---: Boundary type (Ireg, Jreg, Boundary)
      INTEGER nMomBdTp(mgri,mgrj,*)

C---: Momentum Boundary type
      INTEGER nTemBdTp(mgri,mgrj,*)

      INTEGER nPpeSolver, msorit

      REAL    dlref, uref,tref, tmax,
     .        dka, re, pe, sortol, sorrel

      REAL    cu0, TsCoef, HsCoef, TemCoef,
     .        bnumc, rmax, rlc

      REAL   u1(0:mnx,0:mny),  v1(0:mnx,0:mny),  t1(0:mnx,0:mny)
      REAL   uss(0:mnx,0:mny), vss(0:mnx,0:mny),
     .       pss(0:mnx,0:mny), tss(0:mnx,0:mny)


C---: Small-scale Schuman filter parameters
      REAL   fp(*)

C---: Temperature region values
      REAL    dTRgVal(mgri,mgrj)

C---: Boundary condition values (Ireg, Jreg, Boundary, Variable)
      REAL    dBCVal(mgri,mgrj,4,*)


C---: Metric information
      REAL   rau(0:mnx,0:mny), rbu(0:mnx,0:mny),
     .       rbv(0:mnx,0:mny), rgv(0:mnx,0:mny),
     .       dju(0:mnx,0:mny), djv(0:mnx,0:mny), 
     .       djc(0:mnx,0:mny),
     .       xeu(0:mnx,0:mny), yeu(0:mnx,0:mny), 
     .       xzv(0:mnx,0:mny), yzv(0:mnx,0:mny),
     .       xzu(0:mnx,0:mny), yzu(0:mnx,0:mny), 
     .       xev(0:mnx,0:mny), yev(0:mnx,0:mny),
     .       xec(0:mnx,0:mny), yec(0:mnx,0:mny), 
     .       xzc(0:mnx,0:mny), yzc(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j, k, l, ireg, jreg, iW, iE, jS, jN
      INTEGER nmap, nSorConv
      REAL    umpsd, vmpsd, tmpsd, rehmin, pehmin,
     .        hs, dk, pr, rnu, dkappa,
     .        ump, vmp, tmp,
     .        tArea, cuU, cuV, cuT,
     .        hxy, uz, ue, vz, ve,
     .        uxsq, uysq, vxsq, vysq,
     .        delu2n, reh, ts, bnum, rmap,
     .        grduf1, grduf2, grdvf1, grdvf2,
     .        grduf, grdvf, grd,
     .        s1, s2, s3, rnrmjs,
     .        zeta1, zeta2, zeta3, 
     .        alf11, alf12, alf21, alf22, alf13, alf23,
     .        alf31, alf32, alf33, 
     .        um1, um2, um3, vm1, vm2, vm3, tm1, tm2, tm3,
     .        um, vm, tm, au, av, at, uscon, vscon,
     .        tz, te, txsq, tysq, delt2n, peh,
     .        grdtf1, grdtf2, grdtf

C---: Local arrays
      REAL  umap(0:mnx,0:mny,3), vmap(0:mnx,0:mny,3),
     .      tmap(0:mnx,0:mny,3)

      REAL  uc(0:mnx,0:mny), vc(0:mnx,0:mny)
      REAL  uf(0:mnx,0:mny), vf(0:mnx,0:mny), tf(0:mnx,0:mny)
      REAL  xz(0:mnx,0:mny), xe(0:mnx,0:mny), yz(0:mnx,0:mny),
     .      ye(0:mnx,0:mny), rj(0:mnx,0:mny)
      REAL  ul(0:mnx,0:mny), vl(0:mnx,0:mny), tl(0:mnx,0:mny)

C---: Cell areas
      REAL  area(0:mnx,0:mny)

C---: Double precision constants
      REAL   dZero, dOne, dTwo, dThree, dFour, dSix, dHalf
      parameter ( dZero  = 0.00d+00 ,
     .            dOne   = 1.00d+00 ,
     .            dTwo   = 2.00d+00 ,
     .            dThree = 3.00d+00 ,
     .            dFour  = 4.00d+00 ,
     .            dSix   = 6.00d+00 ,
     .            dHalf  = 0.50d+00 )

C---: Logistic map constants
      REAL  dAr, dAm, rc,  piosr2
      parameter ( dAr = 4.82842712474d+00 ,
     .            dAm = 1.47839783948d+00 ,
     .            rc  = 0.20710678119d+00 )

cC---: Save local map iterates
c      save umap, vmap, tmap
      save

C---: Note:  Local variable names vs. variable names in main
C---:   cu0   <- ssCu0
C---:   bnumc <- ssBnCrit
C---:   rmax  <- ssRMpMax
C---:   rlc   <- ssRMpExp

C---: Pi/sqrt(2)
      piosr2 = dacos(-dOne)/sqrt(dTwo)

C---: Map seeds (note: Should read from file?)
      umpsd = 0.92d0
      vmpsd = 0.31d0
      tmpsd = 0.50d0

C---: Minimum reh and peh at which to start iterating the maps
      rehmin = 3.d0
      pehmin = 3.d0

C---: Small-Scale reference length and time step size (???)
      hs     = dlref*HsCoef
      dk     = dka*dlref/uref

C---: Prandtl number
      pr     = pe/re

C---: Kinematic viscosity (nu) and thermal diffusivity (kappa)
      rnu    = uref*dlref/re
      dkappa = rnu/pr

C---: Initialize maps and cell areas
      if(initflg.le.0) then
        ump = umpsd
        vmp = vmpsd
        tmp = tmpsd

        do i=0, nx+1
          do j=0, ny+1

            do l=1, 3
              ump = rc * dAr * ump * ( dOne - dAm * dabs(ump) )
              vmp = rc * dAr * vmp * ( dOne - dAm * dabs(vmp) )
              tmp = rc * dAr * tmp * ( dOne - dAm * dabs(tmp) )
              umap(i,j,l) = ump
              vmap(i,j,l) = vmp
              tmap(i,j,l) = tmp
            end do

            uss(i,j) = dZero
            vss(i,j) = dZero
            tss(i,j) = dZero

          end do
        end do

C-----: Initialize total dimensionless grid areas
        tArea = dZero

C-----: Calculate cell areas for each group of finite volumes
        do i=1,nx
          do j=1,ny

C---------: Dimensionless cell areas, assuming Delta_xi = Delta_eta = 1
            area(i,j) = dOne/djc(i,j)

C---------: Accumulate cell areas
            tArea = tArea + area(i,j)

          end do
        end do

cC-----: Print total dimensional areas
c        write(STDOUT,'(a)') '* Total grid areas:'
c        write(STDOUT,'(4x,a,e13.6)') 'Area: ', tArea * dlref**2
c        write(STDOUT,*)

        if(initflg.lt.0) return

      end if

C------------------------------| U, V, T |------------------------------

C---: III. Calculate small-scale temperature fluctuations at center of 
C---: [p] C.V. [grid point i-1/2,j-1/2]

      do i=1,nx
        do j=1,ny
C-------: Unscale metric information
          xz(i,j) = dlref*xzc(i,j)
          xe(i,j) = dlref*xec(i,j)
          yz(i,j) = dlref*yzc(i,j)
          ye(i,j) = dlref*yec(i,j)
          rj(i,j) = djc(i,j)/(dlref*dlref)

C-------: Average vertical component of velocity on grid point i-1/2,j-1/2
C-------: and unscale large-scale velocity components
          vl(i,j) = uref*(v1(i,j)+v1(i,j-1))*dHalf
          ul(i,j) = uref*(u1(i,j)+u1(i-1,j))*dHalf

C-------: Unscale temperature
          tl(i,j) = (tmax-tref)*t1(i,j)+tref

C-------: Calculate contravariant velocity components 
          uc(i,j) = ye(i,j)*ul(i,j)-xe(i,j)*vl(i,j)
          vc(i,j) = xz(i,j)*vl(i,j)-yz(i,j)*ul(i,j)

C-------: Initialize arrays for filtering
          uf(i,j) = uc(i,j)
          vf(i,j) = vc(i,j)
          tf(i,j) = tl(i,j)

        end do
      end do

C---: Enforce velocity boundary conditions (???)

C---: Enforce temperature boundary conditions
      call TempBoundCond(nx, ny,
     .                   nReg, nRegBrd,
     .                   nTRgType, nTemBdTp,
     .                   dTRgVal, dBCVal,
     .                   tf)


C---: Perform high-pass filter
      call Filter(nx, ny, _U_,
     .            nReg, nRegBrd,
     .            nRegType, nMomBdTp, nTRgType,
     .            fp(_U_), uf)

      call Filter(nx, ny, _V_,
     .            nReg, nRegBrd,
     .            nRegType, nMomBdTp, nTRgType,
     .            fp(_V_), vf)

      call Filter(nx, ny, _T_,
     .            nReg, nRegBrd,
     .            nRegType, nMomBdTp, nTRgType,
     .            fp(_T_), tf)


C---: Keep high-frequency part
      do i=1,nx
        do j=1,ny
          uf(i,j) = uc(i,j)-uf(i,j)
          vf(i,j) = vc(i,j)-vf(i,j)
          tf(i,j) = tl(i,j)-tf(i,j)
        end do
      end do


C---: Scan regions: Vertical sweep
      do jreg=1, nReg(_J_)      
C-----: Scan regions: Horizontal sweep
        do ireg=1, nReg(_I_)

C-------: Extract boundary indeces
          iW = nRegBrd(ireg,jreg,WEST)
          iE = nRegBrd(ireg,jreg,EAST)
          jS = nRegBrd(ireg,jreg,SOUTH)
          jN = nRegBrd(ireg,jreg,NORTH)

      do i=iW+1, iE         !---> Note: Start unindented loops
        do j=jS+1, jN


         if(nRegType(ireg,jreg).eq.RM_BLOCKG .or.
     .       nTRgType(ireg,jreg).eq.RT_TEMPER    ) cycle


C-------: Calculate sub-grid reference length (unscaled)
          hxy=dsqrt(xz(i,j)**2+yz(i,j)**2+xe(i,j)**2+ye(i,j)**2)

C-------: Calculate 2-norm of del u (large-scale, unfiltered):
          uz = (ul(i+1,j)-ul(i-1,j))*dHalf
          ue = (ul(i,j+1)-ul(i,j-1))*dHalf
          vz = (vl(i+1,j)-vl(i-1,j))*dHalf
          ve = (vl(i,j+1)-vl(i,j-1))*dHalf
          uxsq = (rj(i,j)*(ye(i,j)*uz-yz(i,j)*ue))**2
          uysq = (rj(i,j)*(xz(i,j)*ue-xe(i,j)*uz))**2
          vxsq = (rj(i,j)*(ye(i,j)*vz-yz(i,j)*ve))**2
          vysq = (rj(i,j)*(xz(i,j)*ve-xe(i,j)*vz))**2
          delu2n = dsqrt(uxsq + uysq + vxsq + vysq)

C-------: Calculate 2-norm of grad T (large-scale, unfiltered):
          tz = (tl(i+1,j)-tl(i-1,j))*dHalf
          te = (tl(i,j+1)-tl(i,j-1))*dHalf
          txsq = (rj(i,j)*(ye(i,j)*tz-yz(i,j)*te))**2
          tysq = (rj(i,j)*(xz(i,j)*te-xe(i,j)*tz))**2
          delt2n = dsqrt(txsq+tysq)

C-------: Calculate sub-grid Reynolds number
          reh = delu2n*hxy**2/rnu
          peh = pr*reh

          if(peh.gt.pehmin) then

C---------: Calculate cuU, cuV, cuT:
C*** Needs revision
            cuU = cu0 * TemCoef
            cuV = cu0
            cuT = cu0

C---------: Check for small cu, and set number of iterates accordingly
            if (cu0.gt.1.e-10) then

C-----------: Calculate local small-scale time scale:
c              ts = TsCoef*piosr2*fm*reh**(dOne/dThree)/(cu0*delu2n)

C*** Note:
              ts = TsCoef*piosr2*reh**(dOne/dThree)/(cu0*delu2n)

C-----------: Calculate local number of map iterates:
              nmap = dOne + dk/ts
            else
C-----------: Don't iterate the maps if cu0 is less than tolerance
              nmap = 0
            end if

C*** Needs revision
            if(nmap.gt.50) nmap = 50

C---------: Calculate flow bifurcation parameter:
            bnum = sqrt(15.d0)*(hs**2*delu2n/rnu)**(dOne/dSix)

C---------: Calculate chaotic map bifurcation parameter:
            rmap = rmax*dtanh((bnum/bnumc)**rlc*datanh(rc/rmax))

C---------: Calculate anisotropy corrections:
C---------: a. Physical anisotropy of contravariant velocity components
            uz = (uf(i+1,j)-uf(i-1,j))*dHalf
            ue = (uf(i,j+1)-uf(i,j-1))*dHalf
            vz = (vf(i+1,j)-vf(i-1,j))*dHalf
            ve = (vf(i,j+1)-vf(i,j-1))*dHalf
            tz = (tf(i+1,j)-tf(i-1,j))*dHalf
            te = (tf(i,j+1)-tf(i,j-1))*dHalf
            grduf1 = rj(i,j)*(ye(i,j)*uz-yz(i,j)*ue)
            grduf2 = rj(i,j)*(xz(i,j)*ue-xe(i,j)*uz)
            grdvf1 = rj(i,j)*(ye(i,j)*vz-yz(i,j)*ve)
            grdvf2 = rj(i,j)*(xz(i,j)*ve-xe(i,j)*vz)
            grdtf1 = rj(i,j)*(ye(i,j)*tz-yz(i,j)*te)
            grdtf2 = rj(i,j)*(xz(i,j)*te-xe(i,j)*tz)
            grduf = dsqrt(grduf1*grduf1 + grduf2*grduf2)
            grdvf = dsqrt(grdvf1*grdvf1 + grdvf2*grdvf2)
            grdtf = dsqrt(grdtf1*grdtf1 + grdtf2*grdtf2)
            grd   = dsqrt(grduf*grduf + grdvf*grdvf + grdtf*grdtf)
            s1 = grduf/grd
            s2 = grdvf/grd
            s3 = grdtf/grd

C---------: b. Account for grid spacing anisotropy
c            rnrmjs = dsqrt((xzu(i,j)*s1)**2+(yeu(i,j)*s2)**2)
            rnrmjs = dsqrt((xz(i,j)*s1+xe(i,j)*s2)**2
     .               +(yz(i,j)*s1+ye(i,j)*s2)**2)
            zeta1 = dsqrt(dTwo)*s1/rnrmjs
            zeta2 = dsqrt(dTwo)*s2/rnrmjs
            zeta3 = dOne

            if(dabs(grduf).gt.dZero) then
              alf11 = grduf1/grduf
              alf12 = grduf2/grduf
              alf13 = grdtf1/grduf
            else
              alf11 = dZero
              alf12 = dZero
              alf13 = dZero
            end if      
            if(dabs(grdvf).gt.dZero) then
              alf21 = grdvf1/grdvf
              alf22 = grdvf2/grdvf
              alf23 = grdtf2/grdvf
            else
              alf21 = dZero
              alf22 = dZero
              alf23 = dZero
            end if
            if(dabs(grdtf).gt.dZero) then
              alf31 = grduf1/grdtf
              alf32 = grduf2/grdtf
              alf33 = dsqrt(grdtf1**2 + grdtf2**2)/grdtf
            else
              alf31 = dZero
              alf32 = dZero
              alf33 = dZero
            end if      

C---------: Iterate chaotic map:
            um1 = umap(i,j,1)
            um2 = umap(i,j,2)
            um3 = umap(i,j,3)
            vm1 = vmap(i,j,1)
            vm2 = vmap(i,j,2)
            vm3 = vmap(i,j,3)
            tm1 = tmap(i,j,1)
            tm2 = tmap(i,j,2)
            tm3 = tmap(i,j,3)

            do k=1,nmap
              um1 = rmap * dAr * um1 * (dOne - dAm * dabs(um1) )
              um2 = rmap * dAr * um2 * (dOne - dAm * dabs(um2) )
              um3 = rmap * dAr * um3 * (dOne - dAm * dabs(um3) )
              vm1 = rmap * dAr * vm1 * (dOne - dAm * dabs(vm1) )
              vm2 = rmap * dAr * vm2 * (dOne - dAm * dabs(vm2) )
              vm3 = rmap * dAr * vm3 * (dOne - dAm * dabs(vm3) )
              tm1 = rmap * dAr * tm1 * (dOne - dAm * dabs(tm1) )
              tm2 = rmap * dAr * tm2 * (dOne - dAm * dabs(tm2) )
              tm3 = rmap * dAr * tm3 * (dOne - dAm * dabs(tm3) )
            end do

C---------: Save current map iterate to initialize calculations 
C---------: at next large-scale time step
            umap(i,j,1) = um1
            umap(i,j,2) = um2
            umap(i,j,3) = um3
            vmap(i,j,1) = vm1
            vmap(i,j,2) = vm2
            vmap(i,j,3) = vm3
            tmap(i,j,1) = tm1
            tmap(i,j,2) = tm2
            tmap(i,j,3) = tm3

C---------: Construct map combinations (contravariant)
            um = alf11*um1 + alf12*um2
            vm = alf21*um1 + alf22*um2
            tm = alf31*tm1 + alf32*tm2 + alf33*tm3

C---------: Amplitude factors
            au = cuU*reh**(dOne/dSix)*dsqrt(rnu*delu2n)
            av = cuV*reh**(dOne/dSix)*dsqrt(rnu*delu2n)
            at = (dThree*cuT**dFour*peh/pr)**(dOne/dSix)*dsqrt(dkappa)
            at = at*delt2n/dsqrt(delu2n)*TemCoef

C---------: Scale amplitude factor wrt. area ratio and hxy
            au = au * dsqrt(area(i,j) / tArea) * hxy**(dOne/dThree)
            av = av * dsqrt(area(i,j) / tArea) * hxy**(dOne/dThree)
            at = at * dsqrt(area(i,j) / tArea) * hxy**(dOne/dThree)

C---------: Small scale temperature in physical coordinates (scaled)
            uscon = av*zeta1*um
            vscon = av*zeta2*vm

C---------: Small scale velocity in physical coordinates (scaled)
            uss(i,j) = rj(i,j)*(xz(i,j)*uscon+ye(i,j)*vscon)/uref
            vss(i,j) = rj(i,j)*(ye(i,j)*vscon+yz(i,j)*uscon)/uref
            tss(i,j) = (at*tm*zeta3)/(tmax-tref)

C---------: Discard very small numbers
            if(dabs(uss(i,j)).lt.1.d-14) uss(i,j) = dZero
            if(dabs(vss(i,j)).lt.1.d-14) vss(i,j) = dZero
            if(dabs(tss(i,j)).lt.1.d-14) tss(i,j) = dZero      

          end if

        end do   !---> End unindented loop
      end do

        end do   !---> End regions sweep 
      end do

C---: Apply boundary conditions to small scale variables
      call SmlSclBC(nx, ny,
     .              nReg, nRegBrd, 
     .              nRegType, nMomBdTp, nTRgType, nTemBdTp,
     .              dBCVal,
     .              uss, vss, pss, tss)


C---------------------------------| P |---------------------------------

C---: Solve small-scale ppe (i.e. enforce mass conservation)

      do j=0,ny
        do i=0,nx
          pss(i,j) = dZero
        end do
      end do

C---: Apply boundary conditions to small scale variables
      call SmlSclBC(nx, ny,
     .              nReg, nRegBrd, 
     .              nRegType, nMomBdTp, nTRgType, nTemBdTp,
     .              dBCVal,
     .              uss, vss, pss, tss)

C---: Small-scale Pressure-Poisson equation
      call Ppe(nx, ny,
     .         nReg, nRegBrd,
     .         nRegType, 
     .         lCartesGrid,
     .         nPpeSolver, msorit, nSorConv,
     .         dka, sortol, sorrel,
     .         rau, rbu, rbv, rgv,
     .         xeu, yeu, xzv, yzv,
     .         uss, vss, pss)

C---: Apply boundary conditions to small scale variables
      call SmlSclBC(nx, ny,
     .              nReg, nRegBrd,
     .              nRegType, nMomBdTp, nTRgType, nTemBdTp,
     .              dBCVal,
     .              uss, vss, pss, tss)

C---: Project small-scale velocity onto a divergence free-space
C---: to enforce mass conservation
      call Project(nx, ny, 
     .             nReg, nRegBrd,
     .             nRegType, nMomBdTp,
     .             dka,
     .             dju, djv,
     .             yeu, xzv, yzu, xev,
     .             pss, uss, vss)


      return
      end

#endif


#ifdef _DATANH_

/*
  Define the double precision hyperbolic arctangent function if
  not available in the standard math library.
*/

*----------------------------------------------------------------------*
*                                datanh                                *
*----------------------------------------------------------------------*
*      
      real*8 function datanh(x)

      implicit none

      real*8 x

      datanh = dlog((1.00d0 + x)/(1.00d0 - x)) / 2.00d0

      return
      end
*
*----------------------------------------------------------------------*

#endif _DATANH_


*-----------------------|---|---|---V---|---|---|----------------------*
