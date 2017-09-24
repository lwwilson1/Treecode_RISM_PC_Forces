!
!   Author:  Henry A. Boateng  (boateng@umich.edu)
!   Department of Mathematics
!   University of Michigan, Ann Arbor
!
!   Copyright (c) 2013. The Regents of the University of Michigan.
!   All Rights Reserved.
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      PROGRAM TREEDRIVER
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! runtime parameters

      INTEGER :: numparsS,numparsT,order,maxparnodeS
      REAL(KIND=r8) :: theta 

      INTEGER :: ftype, pot_type
      REAL(KIND=r8) :: kappa, eta, eps, T

! arrays for coordinates and charges and energy of target particles

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xS,yS,zS,qS  !source particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: denergy,tenergy !exact energy and energy via treecode

      REAL(KIND=r8),DIMENSION(6) :: xyzminmax
      INTEGER,DIMENSION(3) :: xyzdim

! variables for potential energy computation

      REAL(KIND=r8) :: tpeng,dpeng

! variables needed for f90 DATE_AND_TIME intrinsic

      CHARACTER (LEN=25) :: sampin1,sampin3,sampout
      REAL(KIND=r8)      :: timedirect,timetree

! variables for error calculations

      REAL(KIND=r8) :: inferr,relinferr
      REAL(KIND=r8) :: n2err,reln2err

! local variables

      INTEGER :: i,err,stat
      CHARACTER (LEN=10) :: c1, c2, c3, c4, c5, c6, c7
      REAL(KIND=r8) :: a1, a2, a3, a4

! EXECUTABLE STATEMENTS
      WRITE(6,*) 'Enter the name of input file 1 (xyzq target data)'
      READ(5,*) sampin1

      WRITE(6,*) 'Enter the type of input file 1 (0 is pqr, 1 is flat)'
      READ(5,*) ftype
      
      WRITE(6,*) 'Enter the name of input file 2 (source grid correlation data)'
      READ(5,*) sampin3

      WRITE(6,*) 'Enter direct forces results'
      READ(5,*) sampin3
      
      !WRITE(6,*) 'Enter the name of output file'
      !READ(5,*) sampout
      sampout = 'out.txt'
      

      WRITE(6,*) 'Enter particle number for targets'
      READ(5,*) numparsT

      WRITE(6,*) 'Enter THETA : ' ! The multipole acceptability criterion
      READ(5,*) theta         

      WRITE(6,*) 'Enter ORDER : ' ! The order of the approximation 
      READ(5,*) order

      WRITE(6,*) 'Enter MAXPARNODE : '  ! maximum number of particles in a leaf
      READ(5,*) maxparnodeS ! maxparnodeS is for leaves of target tree


      WRITE(6,*) 'Enter voxvol value'
      READ(5,*) voxvol 
      
      
      WRITE(6,*) 'Enter xyzminmax source grid limits'
      READ(5,*) xyzminmax(1),xyzminmax(2),xyzminmax(3)
      READ(5,*) xyzminmax(4),xyzminmax(5),xyzminmax(6)

      WRITE(6,*) 'Enter xyzdim source grid dimensions'
      READ(5,*) xyzdim(1),xyzdim(2),xyzdim(3)

      numparsS=xyzdim(1)*xyzdim(2)*xyzdim(3)


      ALLOCATE(xT(numparsT),yT(numparsT),zT(numparsT),qT(numparsT),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for xT, yT, zT and qT! '
          STOP
      END IF

      ALLOCATE(guv(numparsS),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for xT, yT, zT and qT! '
          STOP
      END IF

      ALLOCATE(tforces(3,numparsT),dforces(3,numparsT),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating tenergy or denergy! '
          STOP
      END IF

! Read in coordinates and charges for the target particles
      OPEN(unit=82,file=sampin1,status='old',action='read')
 
      WRITE(6,*) ' ' 
      WRITE(6,*) "Reading in targets..."

      IF (ftype == 0) THEN
          i = 1
          DO
              READ(82, 100, iostat=stat) c1, c2, c3, c4, c5, a1, a2, a3, a4, c6, c7
              IF (stat /= 0) EXIT
              IF (c1(1:4) == 'ATOM') THEN
                  WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                        char(13), " Reading in target ", i, " of ", numparsT
                  xT(i) = a1
                  yT(i) = a2
                  zT(i) = a3
                  qT(i) = a4
                  i = i + 1
              END IF
          END DO

      ELSE
          DO i=1,numparsT
              WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                    char(13), " Reading in target ", i, " of ", numparsT
              READ(82,*) xT(i),yT(i),zT(i),qT(i)
!             READ(82,'(F15.10,3F16.10)') xS(i),yS(i),zS(i),qS(i)
          END DO
      END IF
     
      CLOSE(82)


! Read in the values for the correlation at each source grid point       
      OPEN(unit=83,file=sampin2,status='old',action='read')

      WRITE(6,*) ' ' 
      WRITE(6,*) "Reading in source grid point correlation..."

      DO i=1,numparsS
          READ(83,*) guv(i)
          IF ( MOD(i, 100000) == 0) THEN
             WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                    char(13), " Reading in direct ", i, " of ", numparsS
          END IF
      END DO

      CLOSE (83)


! Read in the values for the exact forces at each target       
      OPEN(unit=84,file=sampin3,status='old',action='read')

      WRITE(6,*) ' ' 
      WRITE(6,*) "Reading in direct energy results..."

      READ(84,13) timedirect
      DO i=1,numparsT
          READ(84,*) dforces(1,i), dforces(2,i), dforces(3,i)
          IF ( MOD(i, 100000) == 0) THEN
             WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                    char(13), " Reading in direct ", i, " of ", numparsT
          END IF
      END DO

      CLOSE (84)


      OPEN(unit=85,file=sampout,status='replace',action='write')

! Calling main subroutine to approximate the energy

      CALL TREECODE(xT, yT, zT, qT, &                !target particle info
                    guv, xyzminmax, xyzdim, &        !source grid info
                    numparsS, numparsT, &            !number of sources, targets
                    order, theta, maxparnodeS, &     !tree info
                    voxvol, &                        !time, parameter
                    tforces, timetree)               !output


      WRITE(6,*) ' '
      WRITE(6,'("                      Direct time (s): ",ES15.8)') timedirect
      WRITE(6,'("                        Tree time (s): ",ES15.8)') timetree
      WRITE(6,*) ' '


! compute energy errors

      inferr = MAXVAL(ABS(dforces-tforces))
      relinferr = inferr / MAXVAL(ABS(dforces))
      n2err = SQRT(DOT_PRODUCT(dforces-tforces,dforces-tforces))
      reln2err = n2err / SQRT(DOT_PRODUCT(dforces,dforces))

! output errors to standard out

      WRITE(6,*) ' '
      WRITE(6,'(" Absolute inf norm error in potential: ",ES15.8)') inferr
      WRITE(6,'(" Relative inf norm error in potential: ",ES15.8)') relinferr
      WRITE(6,*) ' ' 
      WRITE(6,'("   Absolute 2 norm error in potential: ",ES15.8)') n2err
      WRITE(6,'("   Relative 2 norm error in potential: ",ES15.8)') reln2err
      WRITE(6,*) ' '
         
      write(85,15)numparsS,numparsT,maxparnodeS,order,theta
      write(85,16)relinferr,reln2err,timedirect,timetree

      CLOSE(unit=85)

 13   FORMAT(E24.16)
 15   FORMAT(I8,2X,I8,2X,I4,2X,I3,2X,F12.8)
 16   FORMAT(E24.16,2X,E24.16,2X,F24.16,2X,F24.16) 
 17   FORMAT(14X,E24.16)
100   FORMAT(A7, A5, A5, A7, A6, F8.3, F8.3, F8.3, F8.4, A8, A9) 

      END PROGRAM TREEDRIVER

!!!!!!!!!!!!!!

