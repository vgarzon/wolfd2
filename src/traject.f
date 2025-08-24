*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                              traject.F                               *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* PartParam    - Load individual particle parameters and Eom coeff.
* Traject      - Driver function to integrate position of all particles
* HeunTrap     - Heun/Trapezoidal integration
* FwdEuler     - Forward Euler integrations
* TrajFunc     - Function of Heun/Trapezoidal algorithm
* TrajJac      - Jacobian of Heun/Trapezoidal algorithm
* FindPos      - Find nearest indeces to a given position
* BiLinInterp  - Bilinear interpolation of grid function to given pos.
* Gauss        - Gaussian elimination
*
*----------------------------------------------------------------------*

/*
  Only compile the contents of this file when _TRAJECT_ is defined
*/

#ifdef _TRAJECT_

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"

*----------------------------------------------------------------------*
*                             InitTraject                              *
*----------------------------------------------------------------------*

      subroutine InitTraject(partFile,
     .                       ntr, nts,
     .                       iTrajFreq, nTSubStp,
     .                       nTOutBnd,
     .                       nTrMethod, nTrCdEq, mTrHTmit,
     .                       dlref, densref, re,
     .                       dTrHTtol, dTrHTdel,
     .                       xp, yp, up, vp,
     .                       cpartx, cparty, repc)

      implicit none

      include "config.f"

C---: Arguments
      character*(*) partFile

      INTEGER ntr, nts, iTrajFreq, nTSubStp,
     .        nTrMethod, nTrCdEq, mTrHTmit

      INTEGER nTOutBnd(*)

      REAL    dlref, densref, re, dTrHTtol, dTrHTdel

      REAL    xp(*), yp(*), up(*), vp(*),
     .        cpartx(*), cparty(*), repc(*)

C---: Local variables
      character*(mfnmlgth) trajFile
      INTEGER i, j, l, ind
      INTEGER ntrx, ntry
      REAL    tinitx1, tinity1, tinitxinc, tinityinc,
     .        partdiam, partdens, dfx, dfy, c1

C---: Double precision constants
      REAL   dZero, dThree, dFour
      parameter ( dZero  = 0.00d+00 ,
     .            dThree = 3.00d+00 ,
     .            dFour  = 4.00d+00 )


C---: Read particle information from file [particles.dat]
      open(14,file=partFile,status='old')
         read(14,'(16x,i6)') ntrx
         read(14,'(16x,i6)') ntry
         read(14,'(16x,i6)') iTrajFreq
         read(14,'(16x,i6)') nTrMethod
         read(14,'(16x,i6)') nTrCdEq
         read(14,'(16x,i6)') nTSubStp
         read(14,'(16x,i6)') mTrHTmit
         read(14,'(16x,i6)') dTrHTtol
         read(14,'(16x,i6)') dTrHTdel
         read(14,'(16x,i6)') tinitx1
         read(14,'(16x,i6)') tinity1
         read(14,'(16x,i6)') tinitxinc
         read(14,'(16x,i6)') tinityinc
         read(14,'(16x,i6)') partdiam
         read(14,'(16x,i6)') partdens
         read(14,'(16x,i6)') dfx
         read(14,'(16x,i6)') dfy
      close(14)

C---: Number of trajectories 
      ntr = ntrx*ntry

C---: Check trajectory array bounds
      if(ntr.gt.mntr) then
        print*, 'Requested number of trajectories is too large',
     .  '(mntr=',mntr,').'
        print*, 'Please increase size of parameter mntr and recompile.'
        stop
      end if

C---: Initial position distribution
      do j=1, ntry
        do i=1, ntrx
          ind=(j-1)*ntrx+i
          xp(ind) = tinitx1 + dble(i-1)*tinitxinc
          yp(ind) = tinity1 + dble(j-1)*tinityinc
        end do
      end do

C---: Assume particles start from rest
      do l=1, ntr
        up(l) = dZero
        vp(l) = dZero
      end do

C---: Initialize nTOutBnd: Assume all particles are inside domain
      do l=1, ntr
        nTOutBnd(l) = 0
      end do

C---: Compute particle EoM and Rep coefficients
      do l=1, ntr
        c1 = dlref/(partdiam*partdens)
        cpartx(l) = dfx*c1*dThree/dFour
        cparty(l) = dfy*c1*dThree/dFour
C-----: Particle Reynolds number coefficient
        repc(l) = re*partdiam/dlref
      end do

C---: Open file to save trajectories
      trajfile='traject.out'
      open(13,file=trajFile,status='unknown')

C---: Initialize trajectory file with ntr and nts
      write(13,*) ntr, (nts/itrajfreq)

      return
      end


*----------------------------------------------------------------------*
*                               Traject                                *
*----------------------------------------------------------------------*

      subroutine Traject(nx, ny,
     .                   ntr, ntsubstp,
     .                   nTrMethod, nTrCdEq, mTrHTmit,
     .                   nTOutBnd,
     .                   dkflow,
     .                   densref, fr,
     .                   dTrHTtol, dTrHTdel,
     .                   cpartx, cparty, repc,
     .                   x, y, u, v, un, vn,
     .                   dens, densn,
     .                   xp, yp, up, vp)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny, ntr, ntsubstp,
     .        nTrMethod, nTrCdEq, mTrHTmit

      INTEGER nTOutBnd(*)

      REAL    dkflow, densref, fr, dTrHTtol, dTrHTdel
      REAL    cpartx(*), cparty(*), repc(*)

      REAL    x(0:mnx,0:mny),    y(0:mnx,0:mny),
     .        u(0:mnx,0:mny),    un(0:mnx,0:mny),
     .        v(0:mnx,0:mny),    vn(0:mnx,0:mny),
     .        dens(0:mnx,0:mny), densn(0:mnx,0:mny)

      REAL    xp(*), yp(*), up(*), vp(*)

C---: External procedures
      INTEGER iFindPos
      REAL    BiLinInterp

C---: Local variables
      INTEGER neqns
      parameter (neqns=4)

      INTEGER k, l, ip, jp
      REAL    dk, rep, uf1, vf1, ufn, vfn, df1, dfn,
     .        dStokes, cd, cpx1, cpy1, cpxn, cpyn

C---: Double precision constants
      REAL   dZero, dOne, dTwo, dThree, dSix
      parameter ( dZero  = 0.00d+00 ,
     .            dOne   = 1.00d+00 ,
     .            dTwo   = 2.00d+00 ,
     .            dThree = 3.00d+00 ,
     .            dSix   = 6.00d+00 )

C---: Working array
      REAL    w(neqns)

C---: Split the fluid-phase time step if requested
      if(ntsubstp.ne.1) then
         dk = dkflow/dble(ntsubstp)
      else
         dk = dkflow
      end if

C---: Start substeps
      do k=1, ntsubstp

C-----: Traverse all particles
        do l=1, ntr

C-------: If out of bounds do not integrate further: Cycle
          if(nTOutBnd(l).gt.0) cycle

C-------: Find position in grid or signal out-of-domain
          nTOutBnd(l) = iFindPos(nx, ny, ip, jp, xp(l), yp(l), x, y)

C-------: Find local velocity and density via bilinear interpolation
          uf1 = BiLinInterp(ip, jp, xp(l), yp(l), x, y, u)
          vf1 = BiLinInterp(ip, jp, xp(l), yp(l), x, y, v)
          ufn = BiLinInterp(ip, jp, xp(l), yp(l), x, y, un)
          vfn = BiLinInterp(ip, jp, xp(l), yp(l), x, y, vn)
          df1 = BiLinInterp(ip, jp, xp(l), yp(l), x, y, dens)
          dfn = BiLinInterp(ip, jp, xp(l), yp(l), x, y, densn)

C-------: Dimensional density of the fluid
          df1 = densref*(df1 + dOne)
          dfn = densref*(dfn + dOne)

C-------: Particle Reynolds number
          rep = repc(l)*dsqrt((uf1-up(l))**2 + (vf1-vp(l))**2)

C-------: Compute drag coefficient
          dstokes = 24.00d0/rep   !---> Stokes drag

C-------: Select drag coefficient formula
          select case(nTrCdEq)
          case (1)   !---> Stokes
            cd = dstokes

          case (2)   !---> Chein
            cd = dstokes*(dOne + rep**(dTwo/dThree)/dSix)

          case (3)   !---> White
            cd = 0.40d0 + dstokes + dSix/(dOne + dsqrt(rep))

          case (4)   !---> Tilly
            cd = dstokes*(dOne + 0.1970d0 * rep**(0.63d0)
     .         + (0.26d-3)*rep**(1.38d0))

          case default
            print*, 'Error: Wrong nTrCdEq flag passed to Traject'
            stop
          end select

C-------: Coefficients of EoM
          cpx1 = cd*cpartx(l)*df1
          cpxn = cd*cpartx(l)*dfn
          cpy1 = cd*cparty(l)*df1
          cpyn = cd*cparty(l)*dfn

C-------: Initialize working array
          w(1) = xp(l)
          w(2) = up(l)
          w(3) = yp(l)
          w(4) = vp(l)

C-------: Select integration method
          select case(nTrMethod)
          case (1)
             call HeunTrap(mTrHTmit, k, dTrHTtol, dTrHTdel, fr,
     .                     uf1, vf1, ufn, vfn,
     .                     cpx1, cpxn, cpy1, cpyn, w)
          case (2)
             call FwdEuler(dk, fr, ufn, vfn, cpxn, cpyn, w)

          case default
             print*, 'Error: Wrong nTrMethod flag passed to Traject'
             stop
          end select

C-------: Find position in grid or signal out-of-domain
          nTOutBnd(l) = iFindPos(nx, ny, ip, jp, xp(l), yp(l), x, y)

C-------: If out of bounds do not update xp, yp, etc.
          if(nTOutBnd(l).gt.0) cycle

C-------: Return updated values to working array
          xp(l) = w(1)
          up(l) = w(2)
          yp(l) = w(3)
          vp(l) = w(4)

        end do

      end do

      return
      end

*----------------------------------------------------------------------*
*                              FwdEuler                                *
*----------------------------------------------------------------------*

      subroutine FwdEuler(h, fr, uf, vf, cpx, cpy, u)

      implicit none

      include "config.f"

      REAL    h, fr, uf, vf, cpx, cpy
      REAL    u(*)

      u(2) = u(2) + h*(cpx*dabs(uf-u(2))*(uf-u(2)))
      u(1) = u(1) + h*u(2)
      u(4) = u(4) + h*(cpy*dabs(vf-u(4))*(vf-u(4))-(1.0d0/fr))
      u(3) = u(3) + h*u(4)

      return
      end

*----------------------------------------------------------------------*
*                              HeunTrap                                *
*----------------------------------------------------------------------*

      subroutine HeunTrap(maxit, h, toler, delta, fr,
     .                    uf1, vf1, ufn, vfn, 
     .                    cpx1, cpxn, cpy1, cpyn, u)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER neqns
      parameter (neqns=4)

      INTEGER maxit
      REAL    h, toler, delta, fr
      REAL    uf1, vf1, ufn, vfn, 
     .        cpx1, cpxn, cpy1, cpyn
      REAL    u(neqns)

C---: External procecures
      REAL    TrajFunc

C---: Local variables
      REAL    g(neqns), udel(neqns), f(neqns), us(neqns) 
      REAL    dj(neqns,neqns)

      INTEGER i, j, m
      REAL    h2, dkrondel, dumax, du

C---: Constants
      REAL    dZero, dOne, dTwo
      parameter ( dZero = 0.00d+00 ,
     .            dOne  = 1.00d+00 ,
     .            dTwo  = 2.00d+00 )

      LOGICAL conv


      h2 = h/dTwo
      conv = .false.

C---: Begin Newton iterations
      do m=1, maxit
        if (m.eq.1) then
C-------: Evaluate u*(i,n+1) for use in Heun's method
          do i=1, neqns
            g(i) = TrajFunc(i, fr, ufn, vfn, cpxn, cpyn, u)
            us(i) = u(i) + h * g(i)
          end do
          if (maxit.eq.1) then   !---> Heun's method
            do i=1, neqns
              u(i) = us(i)
            end do
            conv = .true.
            exit
          end if
C-------: Calculate initial guess for trapezoidal rule
C-------: from Heun's method
          do i=1, neqns
            us(i) = u(i) + h2*(g(i) + 
     .        TrajFunc(i,fr, uf1, vf1, cpx1, cpy1, us))
          end do
        end if
C-----: Load jacobian for Newton iterations 
        call TrajJac(uf1, vf1, cpx1, cpy1, us, dj)
        do i=1, neqns
          f(i) = us(i) - h2*TrajFunc(i,fr, uf1, vf1, cpx1, cpy1, us) 
     .         - (u(i) + h2*g(i))
          do j=1, neqns
            dkrondel = dZero
            if(i.eq.j) dkrondel = dOne
            dj(i,j) = dkrondel - h2*dj(i,j)
          end do
        end do
        do i=1, neqns
          f(i) = -f(i)
        end do
C-----: Gaussian elimination
        call Gauss(neqns, neqns, dj, f, udel)

        dumax = dZero
        do i=1, neqns
          du = dabs(udel(i))
          if(du.gt.dumax) dumax = du
          us(i) = us(i) + delta * udel(i)
        end do
        if (dumax.lt.toler) then
          conv = .true.
          do i=1, neqns
            u(i) = us(i)
          end do
          exit
        end if
      end do   !---> End Newton iterations
      if (.not. conv)
     .  print*, 'Newton iterations failed to converge after ', maxit

      return
      end

*----------------------------------------------------------------------*
*                              TrajFunc                                *
*----------------------------------------------------------------------*
      
      REAL function TrajFunc(i, fr, uf, vf, cpx, cpy, w)

      implicit none

      include "config.f"

      INTEGER i
      REAL    fr, uf, vf, cpx, cpy
      REAL    w(*)

C---: Local variables
      REAL    temp

      select case (i)
        case (1)
          temp = w(2)

        case (2)
          temp = cpx*dabs(uf - w(2))*(uf - w(2))

        case (3)
          temp = w(4)

        case (4)
          temp = cpy*dabs(vf - w(4))*(vf - w(4)) - 1.00d0/fr
      end select

      TrajFunc = temp

      return
      end
      
*----------------------------------------------------------------------*
*                               TrajJac                                *
*----------------------------------------------------------------------*

      subroutine TrajJac(uf, vf, cpx, cpy, w, dj)
            
      implicit none

      include "config.f"

      INTEGER neqns
      parameter (neqns=4)

      REAL    uf, vf, cpx, cpy
      REAL    w(neqns), dj(neqns,neqns)

      INTEGER i, j

      do j=1, neqns
        do i=1, neqns
          dj(i,j) = 0.00d0
        end do
      end do

      dj(1,2) = 1.00d0
      dj(2,2) = cpx*((uf - w(2)) - dabs(uf - w(2)))
      dj(3,4) = 1.00d0
      dj(4,4) = cpy*((vf - w(4)) - dabs(vf - w(4)))

      return
      end

*----------------------------------------------------------------------*
*                             iFindPos                                 *
*----------------------------------------------------------------------*

      INTEGER function iFindPos(nx, ny, ip, jp, xp, yp, x, y)

      implicit none

      include "config.f"

      INTEGER nx, ny
      INTEGER ip, jp
      REAL    xp, yp
      REAL    x(0:mnx,0:mny), y(0:mnx,0:mny)

      LOGICAL foundxpos, foundypos

C---: Local variables
      INTEGER i, j
      REAL    xrel, yrel
     

C---: Default status
      iFindPos = 0

      foundxpos = .false.
      foundypos = .false.

      ip = -1
      jp = -1
      do j=1, ny
        do i=1, nx
          if(x(i,j).gt.xp.and.y(i,j).gt.yp) then
            ip = i
            jp = j
            xrel = x(i,j) - xp
            yrel = y(i,j) - yp
            foundxpos = .true.
            foundypos = .true.
            exit
          end if
        end do
        if(foundypos) exit
      end do

C---: Check to see if still inside domain
      if(ip.le.1.and.xrel.ge.0.00d0)  iFindPos = 1  !-> Left bndy
      if(jp.le.1.and.yrel.ge.0.00d0)  iFindPos = 3  !-> Bottom bndy
      if(ip.lt.0.or.jp.lt.0)      iFindPos = 2  !-> Right or Top bnds

cC---: Traverse entire grid until nearest point is found
c      dmindist = dabs(x(nx,ny)-x(1,1))
c      dminxrel = 0.00d0
c      dminyrel = 0.00d0
c      ip = 0
c      jp = 0
c      do j=1,ny
c         do i=1,nx
c            xrel = x(i,j) - xp
c            yrel = y(i,j) - yp
c            distrel = dsqrt(xrel**d2 + yrel**d2)
c            if(distrel.lt.dmindist) then
c               dmindist = distrel
c               dminxrel = xrel
c               dminyrel = yrel
c               ip = i
c               jp = j
c            end if
c         end do
c      end do
c      iFindPos = 0
cC---: Check to see if still inside domain
c      if(ip.le.1.and.dminxrel.ge.0.00d0)  iFindPos = 1   !-> Left bndy
c      if(ip.ge.nx.and.dminxrel.le.0.00d0) iFindPos = 2   !-> Right bndy
c      if(jp.le.1.and.dminyrel.ge.0.00d0)  iFindPos = 3   !-> Bottom bndy
c      if(jp.ge.ny.and.dminyrel.le.0.00d0) iFindPos = 4   !-> Top bndy
cC---: If still inside, correct ip and jp in necessary
c      if(iFindPos.le.0) then
c         if(xrel.le.0.00d0) ip = ip + 1
c         if(yrel.le.0.00d0) jp = jp + 1
c      end if

      return
      end

*----------------------------------------------------------------------*
*                             BiLinInterp                              *
*----------------------------------------------------------------------*

      REAL function BiLinInterp(i, j, xs, ys, x, y, u)

      implicit none

      include "config.f"

      INTEGER i, j
      REAL    xs, ys
      REAL    x(0:mnx,0:mny), y(0:mnx,0:mny), u(0:mnx,0:mny)

C---: Local variables
      REAL    x1, x2, x3, x4, y1, y2, y3, y4,
     .        f1, f2, f3, f4, fa, fb, fc, fd,
     .        ya, xb, yc, xd, fsx, fsy

c      x11 = x(i-1,j-1)
c      x12 = x(i-1,j)
c      x21 = x(i,j-1)
c      x22 = x(i,j)
c      x1 = (x11 + x12)/d2
c      x2 = (x21 + x22)/d2
c      y11 = y(i-1,j-1)
c      y12 = y(i-1,j)
c      y21 = y(i,j-1)
c      y22 = y(i,j)
c      y1 = (y11 + y21)/d2
c      y2 = (y12 + y22)/d2
c      f11 = u(i-1,j-1)
c      f12 = u(i-1,j)
c      f21 = u(i,j-1)
c      f22 = u(i,j)
c      BiLinInterp = f11 + (f21 - f11) * (xs - x1) / (x2 - x1) +
c     .    (f12 - f11) * (ys - y1) / (y2 - y1) +
c     .    (f22 - f21 - f12 + f11) * (xs - x1) * (ys - y1) /
c     .    ((x2 - x1)*(y2 - y1))


      x1 = x(i-1,j-1)
      x2 = x(i,j-1)
      x3 = x(i,j)
      x4 = x(i-1,j)
      y1 = y(i-1,j-1)
      y2 = y(i,j-1)
      y3 = y(i,j)
      y4 = y(i-1,j)
      f1 = u(i-1,j-1)
      f2 = u(i,j-1)
      f3 = u(i,j)
      f4 = u(i-1,j)

      fa = f1 + (xs-x1) * (f2-f1) / (x2-x1)
      fb = f2 + (ys-y2) * (f3-f2) / (y3-y2)
      fc = f4 + (xs-x4) * (f3-f4) / (x3-x4)
      fd = f1 + (ys-y1) * (f4-f1) / (y4-y1)

      ya = y1 + (xs-x1) * (y2-y1) / (x2-x1)
      xb = x2 + (ys-y2) * (x3-x2) / (y3-y2)
      yc = y4 + (xs-x4) * (y3-y4) / (x3-x4)
      xd = x1 + (ys-y1) * (x4-x1) / (y4-y1)

      fsx = fd + (xs-xd) * (fb-fd) / (xb-xd)
      fsy = fa + (ys-ya) * (fc-fa) / (yc-ya)

c      fs = (fsx + fsy)/2.00d0

      BiLinInterp = (fsx + fsy)/2.00d0

      return
      end

*----------------------------------------------------------------------*
*                               Gauss                                  *
*----------------------------------------------------------------------*

      subroutine Gauss (nmax, n, a, b, x)

      implicit none

      INTEGER nmax, n
      REAL    a(nmax,nmax), b(nmax), x(nmax)

C---: Local variables
      INTEGER i, j, k, imax
      REAL    da, amax, atemp, btemp
      REAL    dm(nmax, nmax)

      do k=1, n-1   !---> Row interchange
        amax = dabs(a(k,k))   ! Locate largest element (in absolute
        imax = k              ! value) in column beneath pivot element
        do i = k+1, n
          da = dabs(a(i,k))
          if(da.gt.amax) then
            amax = da
            imax = i
          endif
        end do

        if (imax.ne.k) then   ! Interchange rows to place largest
          do j = k, n         ! element of current column in pivot
            atemp = a(k,j)    ! position
            a(k,j) = a(imax,j)
            a(imax,j) = atemp
          end do
          btemp = b(k)
          b(k) = b(imax)
          b(imax) = btemp
        endif

        do i = k+1, n   !---> Forward elimination
          dm(i,k) = a(i,k)/a(k,k)
          b(i) = b(i) - dm(i,k)*b(k)
          do j=k+1, n
            a(i,j) = a(i,j) - dm(i,k)*a(k,j)
          end do
        end do
      end do      

      x(n) = b(n)/a(n,n)
      do i = n-1, 1, -1   !---> Backward Substitution
        x(i) = 0.00d0
          do j=i+1, n
          x(i) = x(i) + a(i,j)*x(j)
        end do
        x(i) = (b(i) - x(i))/a(i,i)
      end do

      return
      end

#endif

*-----------------------|---|---|---V---|---|---|----------------------*
