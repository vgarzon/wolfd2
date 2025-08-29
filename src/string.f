*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
*                                                                      *
*                               string.F                               *
*                                                                      *
*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*

*----------------------------------------------------------------------*
*
* string.F    -
* mTrimLeft   -
* mExtractW   -
* mToLowerC   -
* mToUpperC   -
* lStIsNumb   -
* lStIsReal   -
* lStIsInt    -
* nStrLen     -
* mInitStr    -
* SyntaxError -
* (inum)      -
* (jnum)      -
* (rnum)      -
* (dnum)      -
* PrintTitle  -
* PrintDiff   -
* LogInputPar -
*
*----------------------------------------------------------------------*

/* 
  Preprocessor header file(s) 
*/

#include "wolfd2.h"


*----------------------------------------------------------------------*
*                              mTrimLeft                               *
*----------------------------------------------------------------------*
*
* Revised 98.07.24
*
      INTEGER function mTrimLeft(str)

      implicit none
      character*(*) str
      INTEGER i, l, loc

C---: Length of string
      l = len(str)

C---: Search for first non-blank character
      loc = 0
      do i=1, l
        if(ichar(str(i:i)).ge.33.and.ichar(str(i:i)).le.126) then
          loc = i
          exit
        end if
      end do

C---: If empty string (blank or nullc)
      if(loc.eq.0) then
        mTrimLeft = 0
        return
      end if

C---: Shift non-blank characters to the left
      str(1:l-loc) = str(loc:l)

      mTrimLeft = loc

      end


*----------------------------------------------------------------------*
*                              mExtractW                               *
*----------------------------------------------------------------------*
*
* Revised 98.07.24
* 
      INTEGER function mExtractW(str, wrd)

      implicit none
      character*(*) str, wrd
      INTEGER mTrimLeft
      INTEGER i, l, loc, ntrim

C---: Default value
      mExtractW = 0

C---: Trim any blank padding on the left
      ntrim = mTrimLeft(str)

C---: If string is empty, return default value
      if (ntrim.lt.1) return

C---: Length of string
      l = len(str)

C---: Search for first blank character
      loc = 0
      do i=1, l
        if(ichar(str(i:i)).lt.33.or.ichar(str(i:i)).gt.126) then
          loc = i
          exit
        end if
      end do

C---: If empty string (blank or nullc) return default value
      if(loc.eq.0) return

C---: Initialize word
      do i=1,len(wrd)
        wrd(i:i) = char(0)
      end do

C---: Copy non-blank characters to 'word'
      wrd(1:loc-1) = str(1:loc-1)

C---: Remove word from string
      str(1:l-loc-1) = str(loc:l)

C---: Return length of word
      mExtractW = loc-1

      return
      end

*----------------------------------------------------------------------*
*                              mToLowerC                               *
*----------------------------------------------------------------------*
*
* Revised 98.07.24
* 
      subroutine mToLowerC(str)

      implicit none
      character*(*) str
      INTEGER i, nbuf, l

C---: Length of string
      l = len(str)

C---: Distance between upper and lower case characters
      nbuf = ichar('a') - ichar('A')

C---: Convert to lower case
      do i=1, l
        if(ichar(str(i:i)).ge.ichar('A').and.
     .     ichar(str(i:i)).le.ichar('Z'))
     .    str(i:i) = char(ichar(str(i:i)) + nbuf)
      end do

      return
      end

*----------------------------------------------------------------------*
*                              mToUpperC                               *
*----------------------------------------------------------------------*
*
* Revised 98.07.24
* 
      subroutine mToUpperC(str)

      implicit none
      character*(*) str
      INTEGER i, nbuf, l

C---: Length of string
      l = len(str)

C---: Distance between upper and lower case characters
      nbuf = ichar('A') - ichar('a')

C---: Convert to upper case
      do i=1, l
        if(ichar(str(i:i)).ge.ichar('a').and.
     .     ichar(str(i:i)).le.ichar('z'))
     .    str(i:i) = char(ichar(str(i:i)) + nbuf)
      end do

      return
      end


*----------------------------------------------------------------------*
*                              lStIsNumb                               *
*----------------------------------------------------------------------*
*
* Revised 98.07.24
* 
      LOGICAL function lStIsNumb(s)

      implicit none
      character*(*) s
      character*1 c
      INTEGER nStrLen
      INTEGER i, l

C---: Length of string
      l = nStrLen(s)   !---> Note: Only check until blank is found

      lStIsnumb = .true.

C---: Check if possibly a number.
C---: Check the first character in s
      c = s(1:1)
      if((c.lt.'0'.or.c.gt.'9').and.c.ne.'+'.and.c.ne.'-'.and.
     .    c.ne.'.') then
        lStIsnumb = .false.
        return
      end if

C---: Now check the rest of string for unwanted characters
      do i=2, l
        c = s(i:i)
        if((c.lt.'0'.or.c.gt.'9').and.
     .      c.ne.'+'.and.c.ne.'-'.and.c.ne.'e'.and.c.ne.'E'.and.
     .      c.ne.'d'.and.c.ne.'D'.and.c.ne.'.'.and.c.ne.',')
     .  lStIsnumb = .false.
      end do

      return
      end


*----------------------------------------------------------------------*
*                              lStIsReal                               *
*----------------------------------------------------------------------*
*
* Revised 98.08.20
* 
      LOGICAL function lStIsReal(s)

      implicit none
      character*(*) s
      character*1 c
      INTEGER nStrLen
      INTEGER i, l

C---: Length of string
      l = nStrLen(s)   !---> Note: Only check until blank is found

      lStIsReal = .true.

C---: Check if possibly a number.
C---: Check the first character in s
      c = s(1:1)
      if((c.lt.'0'.or.c.gt.'9').and.c.ne.'+'.and.c.ne.'-'.and.
     .    c.ne.'.') then
        lStIsReal = .false.
        return
      end if

C---: Now check the rest of string for unwanted characters
      do i=2, l
        c = s(i:i)
        if((c.lt.'0'.or.c.gt.'9').and.
     .      c.ne.'+'.and.c.ne.'-'.and.c.ne.'e'.and.c.ne.'E'.and.
     .      c.ne.'d'.and.c.ne.'D'.and.c.ne.'.'.and.c.ne.',')
     .  lStIsReal = .false.
      end do

      return
      end


*----------------------------------------------------------------------*
*                              lStIsInt                               *
*----------------------------------------------------------------------*
*
* Revised 98.08.20
* 
      LOGICAL function lStIsInt(s)

      implicit none
      character*(*) s
      character*1 c
      INTEGER nStrLen
      INTEGER i, l

C---: Length of string
      l = nStrLen(s)   !---> Note: Only check until blank is found

      lStIsInt = .true.

C---: Check if possibly a number (first character in string)
      c = s(1:1)
      if((c.lt.'0'.or.c.gt.'9').and.c.ne.'+'.and.c.ne.'-') then
        lStIsInt = .false.
        return
      end if

C---: Now check the rest of string for unwanted characters
      do i=2, l
        c = s(i:i)
        if(c.lt.'0'.or.c.gt.'9') then
          lStIsInt = .false.
          exit
        end if
      end do

      return
      end


*----------------------------------------------------------------------*
*                               nStrLen                                *
*----------------------------------------------------------------------*
*
* Revised 98.07.24
* 
      INTEGER function nStrLen(str)

      implicit none
      character*(*) str
      INTEGER i

C---: Start with last character and find the first nonblank
      do i=len(str), 1, -1
        if(ichar(str(i:i)).ge.33.and.ichar(str(i:i)).le.126) then
          nstrlen=i
          return
        end if
      end do

C---: All characters are blanks
      nstrlen = 0

      return
      end


*----------------------------------------------------------------------*
*                               mInitStr                               *
*----------------------------------------------------------------------*
*
* Revised 98.07.27
*
      subroutine mInitStr(str)

      implicit none
      character*(*) str
      INTEGER i

C---: Fill string with NULL characters
      do i=1,len(str)
        str(i:i) = char(0)
      end do

      return
      end


*----------------------------------------------------------------------*
*                             SyntaxError                              *
*----------------------------------------------------------------------*

      subroutine SyntaxError(nline, nword)

      implicit none
      INTEGER nline, nword

      print '(a,i4,a,i4)',
     .  '* Syntax Error in line ', nline, ', word ', nword
      stop

      end


#ifndef _HP_

*----------------------------------------------------------------------*
*                                 inum                                 *
*----------------------------------------------------------------------*

      integer*2 function inum(cNum)

      implicit none
      character*(*) cNum
      read(cNum,*) inum

      return
      end

*----------------------------------------------------------------------*
*                                 jnum                                 *
*----------------------------------------------------------------------*

      INTEGER function jnum(cNum)

      implicit none
      character*(*) cNum
      read(cNum,*) jnum

      return
      end

*----------------------------------------------------------------------*
*                                 rnum                                 *
*----------------------------------------------------------------------*

      real*4 function rnum(cNum)

      implicit none
      character*(*) cNum
      read(cNum,*) rnum

      return
      end

*----------------------------------------------------------------------*
*                                 dnum                                 *
*----------------------------------------------------------------------*

      REAL function dnum(cNum)

      implicit none
      character*(*) cNum
      read(cNum,*) dnum

      return
      end


#endif _HP_


*----------------------------------------------------------------------*
*                              PrintTitle                              *
*----------------------------------------------------------------------*

      subroutine PrintTitle(nUnit)

      implicit none

      character*39 tl(13)
      INTEGER nUnit
      INTEGER i

      data tl(1) /'------------------------------ #   ####'/
      data tl(2) /'@     @        @  @@@@        #  #    #'/
      data tl(3) /'@     @        @ @           #       # '/
      data tl(4) /'@  @  @  @@@@  @ @@@    #####  #####  .'/
      data tl(5) /'@  @  @ @    @ @ @    #    # #       ..'/
      data tl(6) /'@  @  @ @    @ @ @   #    # #       ...'/
      data tl(7) /' @@ @@   @@@@  @ @   ##### ######  ....'/
      data (tl(8)(i:i),i=1,39) / 39 *'-' /
      data tl(9)  /'                 v0.3i'/
      data tl(10) /'Computational Fluid Dynamics Laboratory'/
      data tl(11) /' Department  of Mechanical Engineering'/
      data tl(12) /'        University  of Kentucky'/
      data (tl(13)(i:i),i=1,39) / 39 *'-' /

      do i=1, 13
        write(nUnit,'(20x,a)') tl(i)
      end do

      return
      end


*----------------------------------------------------------------------*
*                              PrintDiff                               *
*----------------------------------------------------------------------*

      subroutine PrintDiff(nPrDfFreq, lPrtBanner, nBannFreq,
     .                     nthermen, nsmallscl,
     .                     k, ks, nQLIter, nMEIter, nSorIter,
     .                     time, dif)

      implicit none

      INTEGER nPrDfFreq, lPrtBanner, nBannFreq, 
     .        nthermen, nsmallscl,
     .        k, ks, nQLIter, nMEIter, nSorIter
      REAL    time
      REAL    dif(4)

C---: Local variables
      character*70 cHLine1
      character*78 cHLine2

C---: Convergence indicators
      character*1 cQLC

      INTEGER nBanCnt, m, i

      data cHLine1 / '----------------------------------------' /
      data cHLine2 / '----------------------------------------' /

      data nBanCnt / 0 /

      save    nBanCnt

C---: Default status of lQLConv indicator
      cQLC = ' '

C---: If have reached frequency,
      if(mod((k-ks),nPrDfFreq).eq.0) then

        if(lPrtBanner.eq.1) then  !---> If banner is requested
          if(nBanCnt.eq.0.or.nBanCnt.ge.nBannFreq) then

            if(nthermen.eq.1) then

              write(STDOUT,'(1x,a)') cHLine2
              write(STDOUT,'(1x,a6,a13,a5,a7,a10,4(a12))')
     .          'Tstp','Time   ','QlIt','PpeIt',
     .          '|p1-pn|','|u1-un|','|v1-vn|', '|t1-tn|'
              write(STDOUT,'(1x,a)') cHLine2

              write(STDLOG,'(1x,a)') cHLine2
              write(STDLOG,'(1x,a6,a13,a5,a7,a10,4(a12))')
     .          'Tstp','Time   ','QlIt','PpeIt',
     .          '|p1-pn|','|u1-un|','|v1-vn|', '|t1-tn|'
              write(STDLOG,'(1x,a)') cHLine2

            else

              write(STDOUT,'(1x,a)') cHLine1
              write(STDOUT,'(1x,a6,a13,a5,a7,a10,3(a13))')
     .          'Tstp','Time   ','QlIt','PpeIt',
     .          '|p1-pn|','|u1-un|','|v1-vn|'
              write(STDOUT,'(1x,a)') cHLine1

              write(STDLOG,'(1x,a)') cHLine1
              write(STDLOG,'(1x,a6,a13,a5,a7,a10,3(a13))')
     .          'Tstp','Time   ','QlIt','PpeIt',
     .          '|p1-pn|','|u1-un|','|v1-vn|'
              write(STDLOG,'(1x,a)') cHLine1

            end if

            nBanCnt = 1
          else
C---------: Increment banner counter
            nBanCnt = nBanCnt + 1
          end if
        end if

C-----: If QL iterations did not converge, give indication
        if(nQLiter.lt.0) cQLC = '*'

        if(nthermen.eq.1) then

          write(STDOUT,'(1x,i6,e13.6,i4,a1,i6,4(e12.5))')
     .      k,time,nQLIter,cQLC,nSorIter,(dif(m),m=1,4)

          write(STDLOG,'(1x,i6,e13.6,i4,a1,i6,4(e12.5))')
     .      k,time,nQLIter,cQLC,nSorIter,(dif(m),m=1,4)

        else

          write(STDOUT,'(1x,i6,e13.6,i4,a1,i6,4(e13.6))')
     .      k,time,nQLIter,cQLC,nSorIter,(dif(m),m=1,3)

          write(STDLOG,'(1x,i6,e13.6,i4,a1,i6,4(e13.6))')
     .      k,time,nQLIter,cQLC,nSorIter,(dif(m),m=1,3)

        end if
  
      end if

      return
      end


*----------------------------------------------------------------------*
*                            LogInputPar                               *
*----------------------------------------------------------------------*

      subroutine LogInputPar(nUnit, 
     .               inputFile, bndryFile, gridFile,
     .               nthermen, neqstate, nsmallscl, ntraject,
     .               nts, mqiter, nmeiter, nfiltu, nfiltv, nfiltt, 
     .               nPpeSolver, msorit,
     .               dk, qtol, dmeittol, fpu, fpv, fpt, 
     .               sortol, sorrel,
     .               ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .               dlref, uref, tref, tmax,
     .               re, pe, fr)

      implicit none

      INTEGER nUnit
      character*(*) inputFile, bndryFile, gridFile
      INTEGER nthermen, neqstate, nsmallscl, ntraject,
     .        nts, mqiter, nmeiter, nfiltu, nfiltv, nfiltt, 
     .        nPpeSolver, msorit
      REAL    dk, qtol, dmeittol, fpu, fpv, fpt, 
     .        sortol, sorrel,
     .        ssCu0, ssTsCoef, ssHsCoef, ssTemCoef,
     .        dlref, uref, tref, tmax,
     .        re, pe, fr

      INTEGER nStrLen


      if (nUnit.ne.STDOUT) then
        write (nUnit,*)
        write (nUnit,'(4a)') '* Input parameters file:    [',
     .        inputFile(1:nStrLen(inputFile)),']'
        write (nUnit,'(4a)') '* Boundary conditions file: [',
     .        bndryFile(1:nStrLen(bndryFile)),']'
        write (nUnit,'(4a)') '* Grid information file:    [',
     .        gridFile(1:nStrLen(gridFile)),']'
      end if

      write (nUnit,*)
      write (nUnit, '(a)') '* Control flags'
      if(nthermen.eq.1) then
        write (nUnit,'(4x,a)') '+ Thermal-Energy'
      else
        write (nUnit,'(4x,a)') '+ Cold Flow'
      end if
      if(nsmallscl.eq.1) then 
        write (nUnit,'(4x,a)') '+ Turbulence Modeling: (ATD)'
      else
        write (nUnit,'(4x,a)') '+ Laminar Flow'
      end if

      if(neqstate.eq.1)  write (nUnit,'(3x,a)') '+ Buoyancy'
      if(ntraject.eq.1)  write (nUnit,'(3x,a)') '+ Trajectories'

      write (nUnit,*)
      write (nUnit,'(a)')       '* Integration parameters'
      write (nUnit,'(4x,a,i7)')    'Number of time steps: ', nts
      write (nUnit,'(4x,a,e12.6)') 'Time-step size:       ',
     .  dk*dlref/uref

      write (nUnit,*)
      write (nUnit,'(a)')       '* Scaling parameters'
      write (nUnit,'(4x,a,e12.6)')   'Reynolds number:      ', re
      if(nthermen.eq.1)
     .  write (nUnit,'(4x,a,e12.6)') 'Peclet number:        ', pe
      if(neqstate.eq.1)
     .  write (nUnit,'(4x,a,e12.6)') 'Froude number:        ', fr


      return
      end


*-----------------------|---|---|---V---|---|---|----------------------*

