
      PROGRAM TREEDRIVER
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! runtime parameters

      INTEGER :: numparsS, numparsT, ftype, numsolv
      REAL(KIND=r8) :: voxvol

! arrays for coordinates, charge, force calculations 

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: xT,yT,zT,qT !target particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: dforces
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: qvdens, qvtemp
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: guv
      REAL(KIND=r8),DIMENSION(6) :: xyzminmax
      INTEGER,DIMENSION(3) :: xyzdim

! variables needed for output files and cpu time

      CHARACTER (LEN=25) :: sampin1,sampin2,sampout
      REAL(KIND=r8)      :: timedirect,timebeg,timeend

! local variables

      INTEGER :: i,j,err,stat
      CHARACTER (LEN=10) :: c1, c2, c3, c4, c5, c6, c7
      REAL(KIND=r8) :: a1, a2, a3, a4

! EXECUTABLE STATEMENTS
      WRITE(6,*) 'Enter the name of input file 1 (xyzq target data)'
      READ(5,*) sampin1
      WRITE(6,*) 'Enter the type of input file 1 (0 is pqr, 1 is flat)'
      READ(5,*) ftype
      WRITE(6,*) 'Enter the name of input file 2 (source grid correlation data)'
      READ(5,*) sampin2
      WRITE(6,*) 'Enter the name of output file'
      READ(5,*) sampout

      WRITE(6,*) 'Enter particle number for targets'
      READ(5,*) numparsT
      WRITE(6,*) 'Enter voxvol value'
      READ(5,*) voxvol
      WRITE(6,*) 'Enter number of solvent types'
      READ(5,*) numsolv

      ALLOCATE(qvdens(numsolv),qvtemp(numsolv),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for qvdens,qvtemp! '
          STOP
      END IF
      WRITE(6,*) 'Enter qvdens, qvtemp for each solvent atom type'
      DO i=1,numsolv
          READ(5,*) qvdens(i), qvtemp(i)
      END DO

      WRITE(6,*) 'Enter xyzminmax source grid limits'
      READ(5,*) xyzminmax(1),xyzminmax(2),xyzminmax(3)
      READ(5,*) xyzminmax(4),xyzminmax(5),xyzminmax(6)
      WRITE(6,*) 'Enter xyzdim source grid dimensions'
      READ(5,*) xyzdim(1),xyzdim(2),xyzdim(3)

      numparsS=xyzdim(1)*xyzdim(2)*xyzdim(3)

      ALLOCATE(xT(numparsT),yT(numparsT),zT(numparsT),qT(numparsT),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for xS, yS, zS and qS! '
          STOP
      END IF

      ALLOCATE(guv(numparsS,numsolv),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for guv! '
          STOP
      END IF

      ALLOCATE(dforces(3,numparsT),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating tenergy or denergy! '
          STOP
      END IF

      OPEN(unit=82,file=sampin1,status='old',action='read')
 
      print *, "Reading in sources..."
      IF (ftype == 0) THEN
          i = 1
          DO
              READ(82, 100, iostat=stat) c1, c2, c3, c4, c5, a1, a2, a3, a4, c6, c7
              IF (stat /= 0) EXIT
              IF (c1(1:4) == 'ATOM') THEN
                  WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  & 
                        char(13), " Reading in source ", i, " of ", numparsT
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
                    char(13), " Reading in source ", i, " of ", numparsT
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
          READ(83,*) guv(i,:)
          IF ( MOD(i, 100000) == 0) THEN
             WRITE(6,'(A1,A,I12,A,I12)',ADVANCE='NO')  &
                    char(13), " Reading in guv ", i, " of ", numparsS
          END IF
      END DO

      CLOSE (83)


      OPEN(unit=85,file=sampout,status='replace',action='write')

      print *
      print *, "Calculating potentials by direct summation..."

      CALL CPU_TIME(timebeg)

      CALL DIRECT_ENG(xT, yT, zT, qT, guv, qvtemp, qvdens, &
                      xyzminmax, xyzdim, &
                      numparsS, numparsT, numsolv, voxvol, dforces)

      CALL CPU_TIME(timeend)
      timedirect = timeend - timebeg

      print *, "Writing direct results to file..."

      WRITE(85,13)timedirect
      DO j=1,numparsT
         WRITE(85,*) dforces(1,j), dforces(2,j), dforces(3,j)
      END DO

      CLOSE(unit=85)

      print *
      print *,   "Total time (s): ", timedirect
      print *

 13   FORMAT(E24.16)
100   FORMAT(A7, A5, A5, A7, A6, F8.3, F8.3, F8.3, F8.4, A8, A9) 

      END PROGRAM TREEDRIVER

!!!!!!!!!!!!!!

      SUBROUTINE DIRECT_ENG(xT, yT, zT, qT, guv, qvtemp, qvdens, &
                            xyzminmax, xyzdim, &
                            numparsS, numparsT, numsolv, voxvol, dforces)
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

      INTEGER,INTENT(IN) :: numparsS, numparsT, numsolv
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(IN) :: xT, yT, zT, qT
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzminmax
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim
      REAL(KIND=r8),DIMENSION(numparsS,numsolv),INTENT(IN) :: guv
      REAL(KIND=r8),DIMENSION(numsolv),INTENT(IN) :: qvdens
      REAL(KIND=r8),DIMENSION(numsolv),INTENT(IN) :: qvtemp
      REAL(KIND=r8),INTENT(IN) :: voxvol

      REAL(KIND=r8),DIMENSION(3,numparsT),INTENT(OUT) :: dforces

! local variables 
   
      INTEGER :: i,j,k,kk,jj,nn,yzdim
      REAL(KIND=r8) :: tx,ty,tz,xi,yi,zi,dist
      REAL(KIND=r8),DIMENSION(3) :: xyz_dd

      xyz_dd(1) = (xyzminmax(2)-xyzminmax(1)) / (xyzdim(1)-1)
      xyz_dd(2) = (xyzminmax(4)-xyzminmax(3)) / (xyzdim(2)-1)
      xyz_dd(3) = (xyzminmax(6)-xyzminmax(5)) / (xyzdim(3)-1)

      xi = xyzminmax(1)-xyz_dd(1)
      yzdim = xyzdim(2)*xyzdim(3)
      dforces = 0.0_r8

      DO i=0,xyzdim(1)-1
          xi = xi + xyz_dd(1)
          yi = xyzminmax(3)-xyz_dd(2)
          DO j=0,xyzdim(2)-1
              yi = yi + xyz_dd(2)
              zi = xyzminmax(5)-xyz_dd(3)
              DO k=0,xyzdim(3)-1
                  zi = zi + xyz_dd(3)
                  nn = (i*yzdim) + (j*xyzdim(3)) + k + 1

                  DO jj=1,numparsT
                      tx=xi-xT(jj)
                      ty=yi-yT(jj)
                      tz=zi-zT(jj)

                      dist = 1.0_r8 / (tx*tx + ty*ty + tz*tz)**(1.5_r8)
                      DO kk=1,numsolv
                          dforces(1,jj) = dforces(1,jj) + tx * dist &
                                        * guv(nn,kk) * qvdens(kk) * qvtemp(kk)
                          dforces(2,jj) = dforces(2,jj) + ty * dist &
                                        * guv(nn,kk) * qvdens(kk) * qvtemp(kk)
                          dforces(3,jj) = dforces(3,jj) + tz * dist &
                                        * guv(nn,kk) * qvdens(kk) * qvtemp(kk)
                      END DO
                  END DO

              END DO
          END DO
      END DO

      DO jj=1,numparsT
          dforces(1,jj) = qT(jj) * dforces(1,jj)
          dforces(2,jj) = qT(jj) * dforces(2,jj)
          dforces(3,jj) = qT(jj) * dforces(3,jj)
      END DO

      dforces = dforces * voxvol
      
      END SUBROUTINE DIRECT_ENG
