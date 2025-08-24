*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                                grid.f                                *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* Grid          - Driver function
* ReadGrid      - Read grid information from file using Plot3D format
* MirrorPts     - Calculate mirror points in original grid
* FullGrid      - Obtain intermediate (averaged) coordinate points and
*                 produce full grid
* Metric        - Compute metric information, including inverse of
*                 transformation jacobian
* CheckGridSize - Check that requested grid does not maximum parameters
*
*----------------------------------------------------------------------*

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"


*----------------------------------------------------------------------*
*                                Grid                                  *
*----------------------------------------------------------------------*

      subroutine Grid(gridFile, 
     .                nx, ny, 
     .                dlref,
     .                gx, gy,
     .                rau, rbu, rbv, rgv,
     .                ran, rbn, rgn,
     .                rac, rbc, rgc,
     .                dju, djv, djc, djn,
     .                xen, yen, xzn, yzn,
     .                xec, yec, xzc, yzc,
     .                xeu, yeu, xzv, yzv,
     .                xzu, yzu, xev, yev)

      implicit none

C---: Configurable parameters
      include "config.f"

C---: Arguments
      INTEGER nx, ny

      REAL    dlref
      character*(*) gridFile

C---: Original (natural location) grid points
      REAL    gx(0:mnx,0:mny),   gy(0:mnx,0:mny)

C---: Metric tensor components and other metric coefficients
      REAL   ran(0:mnx,0:mny), rbn(0:mnx,0:mny), rgn(0:mnx,0:mny),
     .       rac(0:mnx,0:mny), rbc(0:mnx,0:mny), rgc(0:mnx,0:mny),
     .       dju(0:mnx,0:mny), djv(0:mnx,0:mny),
     .       djc(0:mnx,0:mny), djn(0:mnx,0:mny),
     .       rau(0:mnx,0:mny), rbu(0:mnx,0:mny),
     .       rbv(0:mnx,0:mny), rgv(0:mnx,0:mny),
     .       xen(0:mnx,0:mny), yen(0:mnx,0:mny),
     .       xzn(0:mnx,0:mny), yzn(0:mnx,0:mny),
     .       xec(0:mnx,0:mny), yec(0:mnx,0:mny),
     .       xzc(0:mnx,0:mny), yzc(0:mnx,0:mny),
     .       xeu(0:mnx,0:mny), yeu(0:mnx,0:mny),
     .       xzv(0:mnx,0:mny), yzv(0:mnx,0:mny),
     .       xzu(0:mnx,0:mny), yzu(0:mnx,0:mny),
     .       xev(0:mnx,0:mny), yev(0:mnx,0:mny)

C---: Local variables:
C---: Interpolated staggered grids (centers of U, V and P cells
C---: respectively)
      REAL    xu(0:mnx,0:mny), yu(0:mnx,0:mny),
     .        xv(0:mnx,0:mny), yv(0:mnx,0:mny),
     .        xc(0:mnx,0:mny), yc(0:mnx,0:mny)

C---: Local indeces
      INTEGER i, j

C---: Read grid (into 2D array from 2D or 3D file)
      call ReadGrid(gridFile, nx, ny, gx, gy)

      write (STDOUT,*)
      write (STDOUT,'(a,2(i5,a),i8,a)') '* Grid Size: ', 
     .  nx,' * ',ny,' = ', (nx*ny), ' points'
      write (STDOUT,*)

      write (STDLOG,*)
      write (STDLOG,'(a,2(i5,a),i8,a)') '* Grid Size: ', 
     .  nx,' * ',ny,' = ', (nx*ny), ' points'
      write (STDLOG,*)

C---: Scale grid with respect to the reference length
      do j=1,ny
        do i=1,nx
          gx(i,j) = gx(i,j)/dlref
          gy(i,j) = gy(i,j)/dlref
        end do
      end do

C---: Calculate mirror points
      call MirrorPts(nx, ny, gx, gy)

C---: Interpolate staggered grid points
      call FullGrid(nx, ny, gx, gy, xu, yu, xv, yv, xc, yc)

C---: Compute metric information (transformation parameters)
      call Metric(nx, ny,
     .            gx, gy,             ! Original grid points
     .            xu,  yu,            !
     .            xv,  yv,            ! Staggered grids
     .            xc,  yc,            !
     .            rau, rbu, rbv, rgv,
     .            ran, rbn, rgn,
     .            rac, rbc, rgc,
     .            dju, djv, djc, djn,
     .            xen, yen, xzn, yzn,
     .            xec, yec, xzc, yzc,
     .            xeu, yeu, xzv, yzv,
     .            xzu, yzu, xev, yev)

      return
      end


*----------------------------------------------------------------------*
*                              ReadGrid                                *
*----------------------------------------------------------------------*
* 
* ndim -> Grid file dimension:
*   2: 2-dimensional (2D)
*   3: 3-dimensional (3D)
*

C---: Note that the following compiler directive is required for jnum, 
C---: dnum, etc..
C---: This implies that the function as it stands now might only work
C---: with HP fortran.

#ifdef _HP_

$HP9000_800 INTRINSICS

#endif _HP_

      subroutine ReadGrid(gFile, nx, ny, x, y)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny
      character*(*) gFile
      REAL    x(0:mnx,0:mny), y(0:mnx,0:mny)

C---: External procedures
      LOGICAL lCheckFile, lStIsNumb
      INTEGER mTrimLeft, nStrLen, mExtractW

#ifndef _HP_

      INTEGER jnum

#endif _HP_

C---: Local variables
      character*(mlinelgt) lineStr, word

      INTEGER indeces(3)
      INTEGER ncw, nwords, nerror, nstat, lcount
      INTEGER i, j, nz, k

C---: Check that grid file is readable
      if(.not.lCheckFile(gFile,nError)) then
        print*, 'Error code :', nError
        stop
      end if

C---: Find out the grid dimension (2D or 3D)
C---: Open the grid file
      open(11,file=gFile,status='old')

C---: Read the first line in the file
      read(11,'(a)',iostat=nstat) lineStr

      if(mTrimLeft(lineStr).lt.1) then
        print*, 'Syntax error in first line of grid file [',
     .          gFile(1:nStrLen(gFile)),']'
        stop
      end if

C---: Start parsing words
      nwords = 0
      do

C-----: If string is empty, stop parsing words
        if(mTrimLeft(lineStr).lt.1) exit

        ncw = mExtractW(lineStr,word)
        if(ncw.lt.1)  exit 
        nwords = nwords + 1   !---> Increase word counter

        if (nwords.gt.3) exit 

        if(lStIsNumb(word)) then
          indeces(nwords) = jnum(word(1:ncw))
        else
          call SyntaxError(lcount,nwords)
        end if

      end do   !---> End parsing words


C---: default indeces
      nx = indeces(1)
      ny = indeces(2)

C---: Read grid information from file
      select case(nwords)
      case(2)   !---> Expect 2D file
        call CheckGridSize(nx,ny)
        read(11,*) ((x(i,j), i=1,nx), j=1,ny),
     .             ((y(i,j), i=1,nx), j=1,ny)

      case(3)   !---> Expect 3D file, but only read first k-plane

        nz = indeces(3)

        write(STDLOG,*)
        write(STDLOG,'(a)') '- Note: The grid file is in 3D format.'
        write(STDLOG,'(a)') '        Using only the first k-plane.'

        call CheckGridSize(nx,ny)
        read(11,*) (((x(i,j), i=1,nx), j=1,ny), k=1,nz),
     .             (((y(i,j), i=1,nx), j=1,ny), k=1,nz)

      case default
        print*, 'Syntax error in first line of grid file [',
     .          gFile(1:nStrLen(gFile)),']'
        stop
      end select

C---: Close grid file
      close(11)

      return
      end

*----------------------------------------------------------------------*
*                             MirrorPts                                *
*----------------------------------------------------------------------*

      subroutine MirrorPts(nx, ny, x, y)

      implicit none

      include "config.f"

C---: Arguments
      INTEGER nx, ny
      REAL    x(0:mnx,0:mny), y(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j

C---: Double precision constants
      REAL    dTwo
      parameter ( dTwo = 2.00d+00 )

C---: Calculate mirror points
      do i=1,nx   !---> Bottom and top boundaries
         x(i,0)    = dTwo*x(i,1)  - x(i,2)
         x(i,ny+1) = dTwo*x(i,ny) - x(i,ny-1)
         y(i,0)    = dTwo*y(i,1)  - y(i,2)
         y(i,ny+1) = dTwo*y(i,ny) - y(i,ny-1)
      end do
      do j=1,ny   !---> Left and right boundaries
         x(0,j)    = dTwo*x(1,j)  - x(2,j)
         x(nx+1,j) = dTwo*x(nx,j) - x(nx-1,j)
         y(0,j)    = dTwo*y(1,j)  - y(2,j)
         y(nx+1,j) = dTwo*y(nx,j) - y(nx-1,j)
      end do

C---: Lower left corner
      x(0,0) = dTwo*x(1,1) - x(2,2)
      y(0,0) = dTwo*y(1,1) - y(2,2)

C---: Upper left corner
      x(0,ny+1) = dTwo*x(1,ny) - x(2,ny-1)
      y(0,ny+1) = dTwo*y(1,ny) - y(2,ny-1)

C---: Lower right corner
      x(nx+1,0) = dTwo*x(nx,1) - x(nx-1,2)
      y(nx+1,0) = dTwo*y(nx,1) - y(nx-1,2)

C---: Upper right corner
      x(nx+1,ny+1) = dTwo*x(nx,ny) - x(nx-1,ny-1)
      y(nx+1,ny+1) = dTwo*y(nx,ny) - y(nx-1,ny-1)

      return
      end

*----------------------------------------------------------------------*
*                              FullGrid                                *
*----------------------------------------------------------------------*

      subroutine FullGrid(nx, ny, x, y, xu, yu, xv, yv, xc, yc)

      implicit none

C---: Configurable parameters
      include "config.f"

C---: Arguments
      INTEGER  nx, ny

C---: Original (natural location) grid points
      REAL      x(0:mnx,0:mny),  y(0:mnx,0:mny)

C---: Interpolated grids (centers of U, V and P cells respectively)
      REAL      xu(0:mnx,0:mny), yu(0:mnx,0:mny),
     .          xv(0:mnx,0:mny), yv(0:mnx,0:mny),
     .          xc(0:mnx,0:mny), yc(0:mnx,0:mny)

C---: Local variables
      INTEGER i, j

C---: Double precision constants
      REAL    dHalf
      parameter ( dHalf = 0.50d+00 )

C---: Average points in x-direction (center of U cell)
      do j=1, ny+1
        do i=0, nx+1
          xu(i,j) = dHalf*(x(i,j) + x(i,j-1))
          yu(i,j) = dHalf*(y(i,j) + y(i,j-1))
        end do
      end do

C---: Average points in y-direction (center of V cell)
      do j=0, ny+1
        do i=1, nx+1
          xv(i,j) = dHalf*(x(i,j) + x(i-1,j))
          yv(i,j) = dHalf*(y(i,j) + y(i-1,j))
        end do
      end do

C---: Average points at cell center (P or natural cell)
      do j=1, ny+1
        do i=1, nx+1
          xc(i,j) = dHalf*(xu(i,j) + xu(i-1,j))
          yc(i,j) = dHalf*(yv(i,j) + yv(i,j-1))
        end do
      end do

      return
      end

*----------------------------------------------------------------------*
*                               Metric                                 *
*----------------------------------------------------------------------*

      subroutine Metric(nx, ny,
     .                  xn,  yn,               ! Original grid points
     .                  xu,  yu,               !
     .                  xv,  yv,               ! Staggered grids
     .                  xc,  yc,               !
     .                  rau, rbu, rbv, rgv,
     .                  ran, rbn, rgn,
     .                  rac, rbc, rgc,
     .                  dju, djv, djc, djn,
     .                  xen, yen, xzn, yzn,
     .                  xec, yec, xzc, yzc,
     .                  xeu, yeu, xzv, yzv,
     .                  xzu, yzu, xev, yev)

      implicit none

C---: Read configurable parameters
      include "config.f"

C---: Arguments
      INTEGER nx, ny

C---: Orignal grid points
      REAL      xn(0:mnx,0:mny), yn(0:mnx,0:mny)

C---: Interpolated grids (centers of U, V and P cells respectively)
      REAL      xu(0:mnx,0:mny), yu(0:mnx,0:mny),
     .          xv(0:mnx,0:mny), yv(0:mnx,0:mny),
     .          xc(0:mnx,0:mny), yc(0:mnx,0:mny)

C---: Metric information
      REAL   ran(0:mnx,0:mny), rbn(0:mnx,0:mny), rgn(0:mnx,0:mny),
     .       rac(0:mnx,0:mny), rbc(0:mnx,0:mny), rgc(0:mnx,0:mny),
     .       rau(0:mnx,0:mny), rbu(0:mnx,0:mny),
     .       rbv(0:mnx,0:mny), rgv(0:mnx,0:mny),
     .       xen(0:mnx,0:mny), yen(0:mnx,0:mny),
     .       xzn(0:mnx,0:mny), yzn(0:mnx,0:mny),
     .       xec(0:mnx,0:mny), yec(0:mnx,0:mny),
     .       xzc(0:mnx,0:mny), yzc(0:mnx,0:mny),
     .       xeu(0:mnx,0:mny), yeu(0:mnx,0:mny),
     .       xzv(0:mnx,0:mny), yzv(0:mnx,0:mny),
     .       xzu(0:mnx,0:mny), yzu(0:mnx,0:mny),
     .       xev(0:mnx,0:mny), yev(0:mnx,0:mny),
     .       dju(0:mnx,0:mny), djv(0:mnx,0:mny),
     .       djc(0:mnx,0:mny), djn(0:mnx,0:mny)

C---: Local scalar variables
      INTEGER i, j
      REAL    g11, g12, g22

C---: Double precision constants
      REAL     dOne, dQrtr
      parameter ( dOne  = 1.00d+00 ,
     .            dQrtr = 0.25d+00 )


C---: Transformation parameters:  At original grid points (N)
      do j=1, ny
        do i=1, nx

C-------: Basic metric coefficients
          xzn(i,j) = xv(i+1,j) - xv(i,j)
          xen(i,j) = xu(i,j+1) - xu(i,j)
          yzn(i,j) = yv(i+1,j) - yv(i,j)
          yen(i,j) = yu(i,j+1) - yu(i,j)

C-------: Determinant of transformation Jacobian matrix
          djn(i,j) = dOne/(xzn(i,j)*yen(i,j)-xen(i,j)*yzn(i,j))

C-------: Metric tensor coefficients g11, g12, g22
          g11 = xzn(i,j)*xzn(i,j) + yzn(i,j)*yzn(i,j)
          g12 = xzn(i,j)*xen(i,j) + yzn(i,j)*yen(i,j)
          g22 = xen(i,j)*xen(i,j) + yen(i,j)*yen(i,j)

C-------: Additional metric coefficients
          ran(i,j) =  djn(i,j)*g22
          rbn(i,j) = -djn(i,j)*g12*dQrtr
          rgn(i,j) =  djn(i,j)*g11

        end do
      end do

C---: Transformation parameters:  Center of U cells (U)
      do j=1, ny
        do i=1, nx

C-------: Basic metric coefficients
          xzu(i,j) = xc(i+1,j) - xc(i,j)
          xeu(i,j) = xn(i,j)   - xn(i,j-1)
          yzu(i,j) = yc(i+1,j) - yc(i,j)
          yeu(i,j) = yn(i,j)   - yn(i,j-1)

C-------: Determinant of transformation Jacobian matrix
          dju(i,j) = dOne/(xzu(i,j)*yeu(i,j)-xeu(i,j)*yzu(i,j))

C-------: Metric tensor coefficients g11, g12, g22
          g11 = xzu(i,j)*xzu(i,j) + yzu(i,j)*yzu(i,j)
          g12 = xzu(i,j)*xeu(i,j) + yzu(i,j)*yeu(i,j)
          g22 = xeu(i,j)*xeu(i,j) + yeu(i,j)*yeu(i,j)

C-------: Additional metric coefficients
          rau(i,j) =  dju(i,j)*g22
          rbu(i,j) = -dju(i,j)*g12*dQrtr

C Not required at present:
C          rgu(i,j) =  dju(i,j)*g11
C

        end do
      end do

C---: Transformation parameters:  Center of V cells (V)
      do j=1, ny
        do i=1, nx

C-------: Basic metric coefficients
          xzv(i,j) = xn(i,j)   - xn(i-1,j)
          xev(i,j) = xc(i,j+1) - xc(i,j)
          yzv(i,j) = yn(i,j)   - yn(i-1,j)
          yev(i,j) = yc(i,j+1) - yc(i,j)

C-------: Determinant of transformation Jacobian matrix
          djv(i,j) = dOne/(xzv(i,j)*yev(i,j)-xev(i,j)*yzv(i,j))

C-------: Metric tensor coefficients g11, g12, g22
          g11 = xzv(i,j)*xzv(i,j) + yzv(i,j)*yzv(i,j)
          g12 = xzv(i,j)*xev(i,j) + yzv(i,j)*yev(i,j)
          g22 = xev(i,j)*xev(i,j) + yev(i,j)*yev(i,j)

C-------: Additional metric coefficients

C Not required at present:
C          rav(i,j) =  djv(i,j)*g22
C
          rbv(i,j) = -djv(i,j)*g12*dQrtr
          rgv(i,j) =  djv(i,j)*g11

        end do
      end do

C---: Transformation parameters:  Center of P cells (C)
      do j=1, ny
        do i=1, nx

C-------: Basic metric coefficients
          xzc(i,j) = xu(i,j) - xu(i-1,j)
          xec(i,j) = xv(i,j) - xv(i,j-1)
          yzc(i,j) = yu(i,j) - yu(i-1,j)
          yec(i,j) = yv(i,j) - yv(i,j-1)

C-------: Determinant of transformation Jacobian matrix
          djc(i,j) = dOne/(xzc(i,j)*yec(i,j)-xec(i,j)*yzc(i,j))

C-------: Metric tensor coefficients g11, g12, g22
          g11 = xzc(i,j)*xzc(i,j) + yzc(i,j)*yzc(i,j)
          g12 = xzc(i,j)*xec(i,j) + yzc(i,j)*yec(i,j)
          g22 = xec(i,j)*xec(i,j) + yec(i,j)*yec(i,j)

C-------: Additional metric coefficients
          rac(i,j) =  djc(i,j)*g22
          rbc(i,j) = -djc(i,j)*g12*dQrtr
          rgc(i,j) =  djc(i,j)*g11

        end do
      end do

      return
      end

*----------------------------------------------------------------------*
*                           CheckGridSize                              *
*----------------------------------------------------------------------*

      subroutine CheckGridSize(nx, ny)

      implicit none

C---: Read configurable parameters
      include "config.f"

C---: Arguments
      INTEGER nx, ny

      if(((nx+1).gt.mnx) .or. ((ny+1).gt.mny)) then

        write(STDOUT,*)
        write(STDOUT,'(1x,a)')
     .    '-> Error:  The requested grid file is too large.'
        write(STDOUT,'(1x,2a)')
     .    '-> Please increase the parameters mnx and mny in ',
     .    'the file config.f '
        write(STDOUT,'(1x,a,i4,a,i4,a)') 
     .    '-> to at least ', (nx+1), ' and ', (ny+1),
     .    ' respectively and recompile.'
        write(STDOUT,*)        

        write(STDLOG,*)
        write(STDLOG,'(1x,a)')
     .    '-> Error:  The requested grid file is too large.'
        write(STDLOG,'(1x,2a)')
     .    '-> Please increase the parameters mnx and mny in ',
     .    'the file config.f '
        write(STDLOG,'(1x,a,i4,a,i4,a)') 
     .    '-> to at least ', (nx+1), ' and ', (ny+1),
     .    ' respectively and recompile.'
        write(STDLOG,*)        

        stop

      end if

      return
      end

*-----------------------|---|---|---V---|---|---|----------------------*
