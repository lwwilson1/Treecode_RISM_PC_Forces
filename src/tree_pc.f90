!
!   Authors:  Leighton W. Wilson (lwwilson@umich.edu)
!             Henry A. Boateng  (boateng@umich.edu)
!
!   Department of Mathematics
!   University of Michigan, Ann Arbor
!
!   Copyright (c) 2013-2017. The Regents of the University of Michigan.
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
      MODULE treecode_procedures
      IMPLICIT NONE

! r8 is 8-byte (double precision) real

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      REAL(KIND=r8),PARAMETER :: pi = 4_r8*ATAN(1.D0)

! global variables for taylor expansions

      INTEGER :: torder, torderlim, torderflat
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: cf,cf1,cf2,cf3
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:,:) :: a, b
      
! global variables used when computing potential/force

      REAL(KIND=r8),DIMENSION(3) :: tarpos
      REAL(KIND=r8) :: thetasq,tarposq

! global variables for postition and charge storage

      REAL(KIND=r8),DIMENSION(3) :: xyz_ddglob
      INTEGER,DIMENSION(3) :: xyz_dimglob

! global variables for tree level tracking

      INTEGER :: maxlevel, minlevel

! node pointer and node type declarations

      TYPE tnode_pointer
           TYPE(tnode), POINTER :: p_to_tnode
      END TYPE tnode_pointer
      TYPE tnode
           INTEGER          :: numpar
           REAL(KIND=r8),DIMENSION(3) :: xyz_min, xyz_max, xyz_mid

           REAL(KIND=r8)    :: radius,sqradius,aspect
           INTEGER          :: level,num_children,exist_ms
           REAL(KIND=r8),DIMENSION(:),POINTER :: ms
           TYPE(tnode_pointer), DIMENSION(8) :: child

           INTEGER,DIMENSION(3) :: xyz_dim, xyz_lowindex, xyz_highindex
      END TYPE tnode

      CONTAINS


!!!!!!!!!!!!!!!


      SUBROUTINE SETUP(xyzminmax,xyzdim,xyzind,order,theta)
      IMPLICIT NONE
!
! SETUP allocates and initializes arrays needed for the Taylor expansion.
! Also, global variables are set and the Cartesian coordinates of
! the smallest box containing the particles is determined.
!
      INTEGER,INTENT(IN) :: order
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzminmax
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim
      INTEGER,DIMENSION(6),INTENT(INOUT) :: xyzind
      REAL(KIND=r8),INTENT(IN) :: theta

! local variables

      INTEGER :: err,i
      REAL(KIND=r8) :: t1

! global integers and reals:  TORDER, TORDERLIM and THETASQ

      torder=order
      torderlim=torder+1
      torderflat=torderlim*(torderlim+1)*(torderlim+2)/6
      thetasq=theta*theta

      xyz_dimglob = xyzdim

      xyz_ddglob(1) = (xyzminmax(2)-xyzminmax(1)) / (xyz_dimglob(1)-1)
      xyz_ddglob(2) = (xyzminmax(4)-xyzminmax(3)) / (xyz_dimglob(2)-1)
      xyz_ddglob(3) = (xyzminmax(6)-xyzminmax(5)) / (xyz_dimglob(3)-1)

! allocate global Taylor expansion variables

      ALLOCATE(cf(0:torder), cf1(torderlim), cf2(torderlim), cf3(torderlim), &
               a(-2:torderlim, -2:torderlim, -2:torderlim), &
               b(-2:torderlim, -2:torderlim, -2:torderlim), STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating Taylor variables! '
         STOP
      END IF

! initialize arrays for Taylor sums and coeffs
      DO i = 0, torder
         cf(i) = REAL(i,KIND=r8) + 1.0_r8
      END DO

      DO i = 1, torderlim
         t1 = 1.0_r8 / REAL(i,KIND=r8)
         cf1(i) = t1;
         cf2(i) = 1.0_r8 - 0.5_r8 * t1
         cf3(i) = 1.0_r8 - t1
      END DO

! find bounds of Cartesian box enclosing the particles

      xyzind(1) = 0
      xyzind(2) = xyz_dimglob(1)-1
      xyzind(3) = 0
      xyzind(4) = xyz_dimglob(2)-1
      xyzind(5) = 0
      xyzind(6) = xyz_dimglob(3)-1

      RETURN
      END SUBROUTINE SETUP


!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE CREATE_TREE_N0(p,maxparnode,xyzmm,xyzdim,xyzind,level)
      IMPLICIT NONE
!
! CREATE_TREE_N0 recursively creates the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box. The division of a cluster terminates
! when the number of particles in a cluster are is less or equal to maxparnode

      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: level,maxparnode
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzmm
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim
      INTEGER,DIMENSION(6),INTENT(IN) :: xyzind

! local variables
      REAL(KIND=r8), DIMENSION(3) :: xyz_len
      REAL(KIND=r8) :: lmax,t2
      INTEGER :: i, err, loclev, numposchild

      REAL(KIND=r8), DIMENSION(6,8) :: xyzmms
      INTEGER, DIMENSION(3,8) :: xyzdims
      INTEGER, DIMENSION(6,8) :: xyzinds

      REAL(KIND=r8), DIMENSION(6) ::  lxyzmm
      INTEGER, DIMENSION(3) :: lxyzdim
      INTEGER, DIMENSION(6) :: lxyzind

      xyzmms = 0.0_r8
      xyzdims = 0
      xyzinds = 0
      lxyzmm = 0.0_r8
      lxyzdim = 0
      lxyzind = 0
     
! allocate pointer 

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF

! set node fields: number of particles, exist_ms
! and xyz bounds 

      p%numpar=product(xyzdim)
      p%exist_ms=0

      p%xyz_min=xyzmm(1:5:2)
      p%xyz_max=xyzmm(2:6:2)

      p%xyz_dim=xyzdim

      p%xyz_lowindex=xyzind(1:5:2)
      p%xyz_highindex=xyzind(2:6:2)

! compute aspect ratio

      xyz_len=p%xyz_max-p%xyz_min

      lmax=MAXVAL(xyz_len)
      t2=MINVAL(xyz_len)

      IF (t2 .NE. 0.0_r8) THEN
         p%aspect=lmax/t2
      ELSE
         p%aspect=0.0_r8
      END IF

! midpoint coordinates , RADIUS and SQRADIUS 

      p%xyz_mid=(p%xyz_max+p%xyz_min)/2.0_r8
      p%sqradius=SUM(xyz_len**2)/4.0_r8
      p%radius=SQRT(p%sqradius)

! set particle limits, tree level of node, and nullify children pointers

      p%level=level
      p%num_children=0

      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO

      IF (p%numpar .GT. maxparnode) THEN
!
! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
!
         xyzmms(:,1)=xyzmm
         xyzdims(:,1)=xyzdim
         xyzinds(:,1)=xyzind

         CALL PARTITION_8(xyzmms,xyzdims,xyzinds,xyz_len,lmax,numposchild)
!
! create children if indicated and store info in parent
!
         loclev=level+1

         DO i=1,numposchild
            IF (((xyzinds(1,i) .LE. xyzinds(2,i)) .AND. &
                 (xyzinds(3,i) .LE. xyzinds(4,i))) .AND. &
                 (xyzinds(5,i) .LE. xyzinds(6,i))) THEN

               p%num_children=p%num_children+1

               lxyzmm=xyzmms(:,i)
               lxyzdim=xyzdims(:,i)
               lxyzind=xyzinds(:,i)

               CALL CREATE_TREE_N0(p%child(p%num_children)%p_to_tnode, &
                                   maxparnode,lxyzmm,lxyzdim,lxyzind,loclev)
            END IF
         END DO

      END IF   

      END SUBROUTINE CREATE_TREE_N0


!!!!!!!!!!!!!!!


      SUBROUTINE PARTITION_8(xyzmms,xyzdims,xyzinds,xyz_len,lmax,numposchild)

      IMPLICIT NONE
!
! PARTITION_8 determines the particle indices of the eight sub boxes
! containing the particles after the box defined by particles I_BEG
! to I_END is divided by its midpoints in each coordinate direction.
! The determination of the indices is accomplished by the subroutine
! PARTITION. A box is divided in a coordinate direction as long as the
! resulting aspect ratio is not too large. This avoids the creation of
! "narrow" boxes in which Talyor expansions may become inefficient.
! On exit the INTEGER array IND (dimension 8 x 2) contains
! the indice limits of each new box (node) and NUMPOSCHILD the number 
! of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
! that box J is empty.
!
      REAL(KIND=r8),DIMENSION(6,8),INTENT(INOUT) :: xyzmms
      INTEGER,DIMENSION(3,8),INTENT(INOUT) :: xyzdims
      INTEGER,DIMENSION(6,8),INTENT(INOUT) :: xyzinds

      REAL(KIND=r8),DIMENSION(3),INTENT(IN) :: xyz_len
      REAL(KIND=r8),INTENT(IN) :: lmax
      INTEGER,INTENT(INOUT) :: numposchild

! local variables

      INTEGER :: i
      REAL(KIND=r8) :: critlen
      INTEGER :: xdim,ydim,zdim,xn,yn,zn
      INTEGER :: xlowind,xhighind,ylowind,yhighind,zlowind,zhighind
      REAL(KIND=r8) :: xlowmid,xhighmid,ylowmid,yhighmid,zlowmid,zhighmid


      numposchild=1
      critlen=lmax/sqrt(2.0_r8)

      xdim=xyzdims(1,1)
      ydim=xyzdims(2,1)
      zdim=xyzdims(3,1)

      xn=xdim/2
      yn=ydim/2
      zn=zdim/2

      IF (xyz_len(1) .GE. critlen) THEN

         xlowmid=xyzmms(1,1)+(xn-1)*xyz_ddglob(1)
         xhighmid=xyzmms(2,1)-(xdim-xn-1)*xyz_ddglob(1)

         xlowind=xyzinds(1,1)+(xn-1)
         xhighind=xyzinds(2,1)-(xdim-xn-1)

         xyzmms(:,2)=xyzmms(:,1)
         xyzinds(:,2)=xyzinds(:,1)

         xyzmms(2,1)=xlowmid
         xyzmms(1,2)=xhighmid

         xyzinds(2,1)=xlowind
         xyzinds(1,2)=xhighind

         numposchild=2*numposchild

      END IF 

      IF (xyz_len(2) .GE. critlen) THEN

         ylowmid=xyzmms(3,1)+(yn-1)*xyz_ddglob(2)
         yhighmid=xyzmms(4,1)-(ydim-yn-1)*xyz_ddglob(2)

         ylowind=xyzinds(3,1)+(yn-1)
         yhighind=xyzinds(4,1)-(ydim-yn-1)

         DO i=1,numposchild
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzinds(:,numposchild+i)=xyzinds(:,i)

            xyzmms(4,i)=ylowmid
            xyzmms(3,numposchild+i)=yhighmid

            xyzinds(4,i)=ylowind
            xyzinds(3,numposchild+i)=yhighind
         END DO

         numposchild=2*numposchild

      END IF

      IF (xyz_len(3) .GE. critlen) THEN

         zlowmid=xyzmms(5,1)+(zn-1)*xyz_ddglob(3)
         zhighmid=xyzmms(6,1)-(zdim-zn-1)*xyz_ddglob(3)

         zlowind=xyzinds(5,1)+(zn-1)
         zhighind=xyzinds(6,1)-(zdim-zn-1)

         DO i=1,numposchild

            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzinds(:,numposchild+i)=xyzinds(:,i)

            xyzmms(6,i)=zlowmid
            xyzmms(5,numposchild+i)=zhighmid

            xyzinds(6,i)=zlowind
            xyzinds(5,numposchild+i)=zhighind

         END DO

         numposchild=2*numposchild

      END IF

      xyzdims(1,:)=xyzinds(2,:)-xyzinds(1,:)+1
      xyzdims(2,:)=xyzinds(4,:)-xyzinds(3,:)+1
      xyzdims(3,:)=xyzinds(6,:)-xyzinds(5,:)+1

      RETURN 
      END SUBROUTINE PARTITION_8


!!!!!!!!!!!!!!


      SUBROUTINE PC_TREECODE_FORCES(p, guv, &              !source info
                                    xT, yT, zT, qT,  &     !target info
                                    numparsS, numparsT, &  !number sources, targets
                                    forcesT)               !output
      IMPLICIT NONE
!
! CP_TREECODE is the driver routine which calls COMPUTE_CP1 for each
! source particle, setting the global variable TARPOS before the call. After
! the calls to COMPUTE_CP1, CP_TREECODE calls COMPUTE_CP2 to compute
! the energy using power series. 
! P is the root node of the target tree, xT,yT,zT are the coordinates of
! the target particles while xS,yS,zS,qS are the coordinates and charges of the
! source particles. The energy at target 'i', is stored in EnP(i).
!
 
      INTEGER,INTENT(IN) :: numparsS,numparsT
      TYPE(tnode),POINTER :: p  
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(IN) :: xT,yT,zT,qT
      REAL(KIND=r8),DIMENSION(numparsS),INTENT(IN) :: guv
      REAL(KIND=r8),DIMENSION(3,numparsT),INTENT(INOUT) :: forcesT
 
! local variables

      INTEGER :: i
      REAL(KIND=r8),DIMENSION(3) :: force(3)

      forcesT=0.0_r8


      DO i=1,numparsT
          force=0.0_r8
          tarpos(1)=xT(i)
          tarpos(2)=yT(i)
          tarpos(3)=zT(i)
          tarposq=qT(i)

          CALL COMPUTE_PC(p,guv,numparsS,force)
          !DO j=1,p%num_children
          !    CALL COMPUTE_PC(p%child(j)%p_to_tnode,guv,numparsS,forcetemp)
          !    force = force + forcetemp
          !END DO

          forcesT(1,i) = -tarposq * force(1)
          forcesT(2,i) = -tarposq * force(2)
          forcesT(3,i) = -tarposq * force(3)

      END DO

      RETURN
      END SUBROUTINE PC_TREECODE_FORCES


!!!!!!!!!!!!!!!!!!!!!!!!


      RECURSIVE SUBROUTINE COMPUTE_PC(p,guv,numpars,force)
      IMPLICIT NONE

! COMPUTE_CP1 is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target 
! cluster are updated. If the MAC is not satisfied then the algorithm 
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: numpars
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: guv
      REAL(KIND=r8),DIMENSION(3),INTENT(OUT) :: force(3)

! local variables

      REAL(KIND=r8),DIMENSION(3) :: xyz_t, forcetemp
      REAL(KIND=r8) :: distsq
      INTEGER :: i, err, kk, k1, k2, k3

! determine DISTSQ for MAC test

      xyz_t=tarpos-p%xyz_mid
      distsq=SUM(xyz_t**2)

! intialize potential energy and force

      force = 0.0_r8

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .NE. 0.0_r8)) THEN

         IF (p%exist_ms .EQ. 0) THEN
             ALLOCATE(p%ms(torderflat),STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error allocating node moments! '
                STOP
             END IF

             p%ms=0.0_r8
             p%exist_ms=1
             CALL COMP_MS(p, guv, numpars)
         END IF

         a=0.0_r8
         CALL COMP_TCOEFF(xyz_t(1), xyz_t(2), xyz_t(3))

         kk = 0
         DO k3 = 0, torder
            DO k2 = 0, torder-k3
               DO k1 = 0, torder-k3-k2
                  kk = kk + 1
                  force(1) = force(1) + cf(k1) * a(k1+1,k2,k3) * p%ms(kk)
                  force(2) = force(2) + cf(k2) * a(k1,k2+1,k3) * p%ms(kk)
                  force(3) = force(3) + cf(k3) * a(k1,k2,k3+1) * p%ms(kk)
               END DO
            END DO
         END DO

      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
!
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT(p, guv, numpars, force)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_PC(p%child(i)%p_to_tnode, guv, &
                               numpars, forcetemp)
               force = force + forcetemp
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_PC


!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_TCOEFF(dx,dy,dz)
      IMPLICIT NONE
!
! COMP_TCOEFF computes the Taylor coefficients of the potential
! using a recurrence formula.  The center of the expansion is the
! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.
!
      REAL(KIND=r8),INTENT(IN) :: dx,dy,dz

! local variables

      REAL(KIND=r8) :: ddx,ddy,ddz,fac,sqfac
      INTEGER :: i,j,k,i1,i2,j1,j2,k1,k2

! setup variables

      ddx=2.0_r8*dx
      ddy=2.0_r8*dy
      ddz=2.0_r8*dz
      fac=1.0_r8/(dx*dx+dy*dy+dz*dz)
      sqfac=SQRT(fac)

! 0th coeff or function val 

      a(0,0,0)=-sqfac

! 2 indices are 0

      a(1,0,0)=fac*dx*a(0,0,0)
      a(0,1,0)=fac*dy*a(0,0,0)
      a(0,0,1)=fac*dz*a(0,0,0)
  
      DO i=2,torderlim
         i1=i-1; i2=i-2
         a(i,0,0)=fac*(ddx*cf2(i)*a(i1,0,0)-cf3(i)*a(i2,0,0))
         a(0,i,0)=fac*(ddy*cf2(i)*a(0,i1,0)-cf3(i)*a(0,i2,0))
         a(0,0,i)=fac*(ddz*cf2(i)*a(0,0,i1)-cf3(i)*a(0,0,i2))
      END DO

! 1 index 0, 1 index 1, other >=1

      a(1,1,0)=fac*(dx*a(0,1,0)+ddy*a(1,0,0))
      a(1,0,1)=fac*(dx*a(0,0,1)+ddz*a(1,0,0))
      a(0,1,1)=fac*(dy*a(0,0,1)+ddz*a(0,1,0))

      DO i=2,torderlim-1
         i1=i-1; i2=i-2
         a(1,0,i)=fac*(dx*a(0,0,i)+ddz*a(1,0,i1)-a(1,0,i2))
         a(0,1,i)=fac*(dy*a(0,0,i)+ddz*a(0,1,i1)-a(0,1,i2))
         a(0,i,1)=fac*(dz*a(0,i,0)+ddy*a(0,i1,1)-a(0,i2,1))
         a(1,i,0)=fac*(dx*a(0,i,0)+ddy*a(1,i1,0)-a(1,i2,0))
         a(i,1,0)=fac*(dy*a(i,0,0)+ddx*a(i1,1,0)-a(i2,1,0))
         a(i,0,1)=fac*(dz*a(i,0,0)+ddx*a(i1,0,1)-a(i2,0,1))
      END DO
  
! 1 index 0, others >= 2
      
      DO i=2,torderlim-2 
         i1=i-1; i2=i-2
         DO j=2,torderlim-i
            j1=j-1; j2=j-2
            a(i,j,0)=fac*(ddx*cf2(i)*a(i1,j,0)+ddy*a(i,j1,0) &
                            -cf3(i)*a(i2,j,0)-a(i,j2,0))
            a(i,0,j)=fac*(ddx*cf2(i)*a(i1,0,j)+ddz*a(i,0,j1) &
                            -cf3(i)*a(i2,0,j)-a(i,0,j2))
            a(0,i,j)=fac*(ddy*cf2(i)*a(0,i1,j)+ddz*a(0,i,j1) &
                            -cf3(i)*a(0,i2,j)-a(0,i,j2))
         END DO
      END DO

 
! 2 indices 1, other >= 1

      a(1,1,1)=fac*(dx*a(0,1,1)+ddy*a(1,0,1)+ddz*a(1,1,0))

      DO i=2,torderlim-2
         i1=i-1; i2=i-2
         a(1,1,i)=fac*(dx*a(0,1,i)+ddy*a(1,0,i)+ddz*a(1,1,i1) &
                        -a(1,1,i2))
         a(1,i,1)=fac*(dx*a(0,i,1)+ddy*a(1,i1,1)+ddz*a(1,i,0) &
                        -a(1,i2,1))
         a(i,1,1)=fac*(dy*a(i,0,1)+ddx*a(i1,1,1)+ddz*a(i,1,0) &
                        -a(i2,1,1))
      END DO

! 1 index 1, others >=2
      DO i=2,torderlim-3
         i1=i-1; i2=i-2 
         DO j=2,torderlim-i
            j1=j-1; j2=j-2
            a(1,i,j)=fac*(dx*a(0,i,j)+ddy*a(1,i1,j)+ddz*a(1,i,j1) &
                           -a(1,i2,j)-a(1,i,j2))
            a(i,1,j)=fac*(dy*a(i,0,j)+ddx*a(i1,1,j)+ddz*a(i,1,j1) &
                           -a(i2,1,j)-a(i,1,j2))
            a(i,j,1)=fac*(dz*a(i,j,0)+ddx*a(i1,j,1)+ddy*a(i,j1,1) &
                           -a(i2,j,1)-a(i,j2,1))
         END DO
      END DO

! all indices >=2
      DO k=2,torderlim-4
         k1=k-1; k2=k-2
         DO j=2,torderlim-2-k
            j1=j-1; j2=j-2
            DO i=2,torderlim-k-j
               a(i,j,k)=fac*(ddx*cf2(i)*a(i-1,j,k)+ddy*a(i,j1,k) &
                           +ddz*a(i,j,k1)-cf3(i)*a(i-2,j,k) &
                           -a(i,j2,k)-a(i,j,k2)) 
            END DO
         END DO
      END DO

      RETURN

      END SUBROUTINE COMP_TCOEFF


!!!!!!!!!!!!!!!


      SUBROUTINE COMP_MS(ap,guv,numpars)
      IMPLICIT NONE
!
! COMP_CMS computes the moments for node P needed in the Taylor approximation
!
      INTEGER,INTENT(IN) :: numpars
      TYPE(tnode),POINTER :: ap
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: guv

      REAL(KIND=r8) :: tx,ty,tz,dx,dy,dz,xl,yl,zl,xm,ym,zm
      INTEGER :: xlind, ylind, zlind, xhind, yhind, zhind, yzhind
      INTEGER :: i,nn,j,k,k1,k2,k3,kk


      xl=ap%xyz_min(1)
      yl=ap%xyz_min(2)
      zl=ap%xyz_min(3)

      xm=ap%xyz_mid(1)
      ym=ap%xyz_mid(2)
      zm=ap%xyz_mid(3)

      xlind=ap%xyz_lowindex(1)
      ylind=ap%xyz_lowindex(2)
      zlind=ap%xyz_lowindex(3)

      xhind=ap%xyz_highindex(1)
      yhind=ap%xyz_highindex(2)
      zhind=ap%xyz_highindex(3)

      yzhind=xyz_dimglob(2)*xyz_dimglob(3)
      ap%ms=0.0_r8

      DO i=xlind,xhind
         DO j=ylind,yhind
            DO k=zlind,zhind

               nn=(i*yzhind)+(j*xyz_dimglob(3))+k+1

               dx=xl+(i-xlind)*xyz_ddglob(1)-xm
               dy=yl+(j-ylind)*xyz_ddglob(2)-ym
               dz=zl+(k-zlind)*xyz_ddglob(3)-zm

               kk = 0
               tz = 1.0_r8

               DO k3 = 0, torder
                  ty = 1.0_r8
                  DO k2 = 0, torder-k3
                     tx = 1.0_r8
                     DO k1 = 0, torder-k3-k2
                        kk = kk + 1
                        ap%ms(kk) = ap%ms(kk) + guv(nn)*tx*ty*tz
                        tx = dx * tx
                     END DO
                     ty = dy * ty
                  END DO
                  tz = tz * dz
               END DO

            END DO
         END DO
      END DO
         
      RETURN
      END SUBROUTINE COMP_MS


!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE COMP_DIRECT(p, guv, numpars, force)

      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the targets
! in the current cluster due to the 
! current source  (determined by the global variable TARPOS). 
!
      INTEGER,INTENT(IN) :: numpars
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: guv
      REAL(KIND=r8),DIMENSION(3),INTENT(OUT) :: force(3)

! local variables

      INTEGER :: i,j,k,nn
      REAL(KIND=r8) :: tx,ty,tz,xl,yl,zl,dist
      INTEGER :: xlind,ylind,zlind,xhind,yhind,zhind,yzhind

      xl=p%xyz_min(1)
      yl=p%xyz_min(2)
      zl=p%xyz_min(3)

      xlind=p%xyz_lowindex(1)
      ylind=p%xyz_lowindex(2)
      zlind=p%xyz_lowindex(3)

      xhind=p%xyz_highindex(1)
      yhind=p%xyz_highindex(2)
      zhind=p%xyz_highindex(3)

      force = 0.0_r8
      yzhind=xyz_dimglob(2)*xyz_dimglob(3)

      DO i=xlind,xhind
         tx=xl+(i-xlind)*xyz_ddglob(1)-tarpos(1)
         DO j=ylind,yhind
            ty=yl+(j-ylind)*xyz_ddglob(2)-tarpos(2)
            DO k=zlind,zhind
               tz=zl+(k-zlind)*xyz_ddglob(3)-tarpos(3)
               nn=(i*yzhind)+(j*xyz_dimglob(3))+k+1

               !print *, tx,ty,tz, nn, guv(nn)
               !print *, i,j,k, nn, tx+tarpos(1),ty+tarpos(2),tz+tarpos(3)

               dist = 1.0_r8 / (tx*tx + ty*ty + tz*tz)**(1.5_r8)

               force(1) = force(1) + guv(nn) * tx * dist
               force(2) = force(2) + guv(nn) * ty * dist
               force(3) = force(3) + guv(nn) * tz * dist

            END DO
         END DO
      END DO   

      RETURN
      END SUBROUTINE COMP_DIRECT


!!!!!!!!!!!!!!!!!


      SUBROUTINE CLEANUP(p)
      IMPLICIT NONE
!
! CLEANUP deallocates allocated global variables and then
! calls recursive routine REMOVE_NODE to delete the tree.
!
      TYPE(tnode),POINTER :: p      

! local variables
  
      INTEGER :: err

      DEALLOCATE(cf,cf1,cf2,cf3,a,b, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating Taylor variables! '
         STOP
      END IF

      CALL REMOVE_NODE(p)
      DEALLOCATE(p, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating root node! '
         STOP
      END IF 
      NULLIFY(p)         

      RETURN
      END SUBROUTINE CLEANUP


!!!!!!!!!!!


      RECURSIVE SUBROUTINE REMOVE_NODE(p)
      IMPLICIT NONE
!
! REMOVE_NODE recursively removes each node from the
! tree and deallocates its memory for MS array if it
! exits.
!
      TYPE(tnode),POINTER :: p 

! local variables

      INTEGER :: i,err

      IF (p%exist_ms .EQ. 1) THEN
         DEALLOCATE(p%ms,STAT=err)
         IF (err .NE. 0) THEN
            WRITE(6,*) 'Error deallocating node MS! '
            STOP
         END IF               
      END IF

      IF (p%num_children .GT. 0) THEN
          DO i=1,p%num_children
            CALL REMOVE_NODE(p%child(i)%p_to_tnode)
            DEALLOCATE(p%child(i)%p_to_tnode,STAT=err)
            IF (err .NE. 0) THEN
               WRITE(6,*) 'Error deallocating node child! '
               STOP
            END IF                           
          END DO
      END IF 

      RETURN                
      END SUBROUTINE REMOVE_NODE      

      END MODULE treecode_procedures

!!!!!!!!!!!!!!!

      SUBROUTINE TREECODE_FOR(xT, yT, zT, qT, &
                      guv, xyzminmax, xyzdim, &
                      numparsS, numparsT, &
                      order,theta, maxparnodeS, &
                      voxvol, &
                      tEn, timetree)

      USE treecode_procedures
      IMPLICIT NONE

!====================================================================
!                                                                   
! xT,yT,zT,qT    :: x,y,z coordinates and charges of targets
! xyzdim         :: source grid dimensions
! xyzminmax      :: min, max source grid limits
! numparS        :: number of sources
! numparT        :: number of targets
! tEn            :: array of dimension (3,numparsT) for storing forces
!                   at each target
! maxparnodeS    :: maximum number of particles in a leaf
! timetree       :: The total time for the treecode computation
!=====================================================================

      INTEGER,INTENT(IN) :: numparsS,numparsT,order,maxparnodeS
      REAL(KIND=r8),DIMENSION(numparsT),INTENT(IN) :: xT,yT,zT,qT

      REAL(KIND=r8),DIMENSION(numparsS),INTENT(IN) :: guv 
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzminmax
      INTEGER,DIMENSION(3),INTENT(IN) :: xyzdim

      REAL(KIND=r8),INTENT(IN) :: theta
      REAL(KIND=r8),INTENT(IN) :: voxvol

      REAL(KIND=r8),DIMENSION(3,numparsT),INTENT(OUT) :: tEn
      REAL(KIND=r8),INTENT(OUT) :: timetree

! local variables

      TYPE(tnode),POINTER :: trootS
      INTEGER :: level
      INTEGER,DIMENSION(6) :: xyzind

! variables needed for f90 DATE_AND_TIME intrinsic

      REAL(KIND=r8) :: totaltime, timebeg, timeend


! Call SETUP to allocate arrays for Taylor expansions
! and setup global variables. 

      CALL SETUP(xyzminmax,xyzdim,xyzind,order,theta)

! nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(trootS)

! set global variables to track tree levels during construction

      level=0
      minlevel=50000
      maxlevel=0

      WRITE(6,*) ' '
      WRITE(6,*) 'Creating tree...'

      CALL CPU_TIME(timebeg)
      CALL CREATE_TREE_N0(trootS,maxparnodeS,xyzminmax,xyzdim,xyzind,level)
      CALL CPU_TIME(timeend)

      timetree = timeend - timebeg

! print tree information to stdout 

         WRITE(6,*) ' '
         WRITE(6,*) 'Tree created. '
         WRITE(6,*) 'Tree parameters: '
         WRITE(6,*) ' '
         WRITE(6,*) '         numpar: ',trootS%numpar
         WRITE(6,*) '          x_mid: ',trootS%xyz_mid(1)
         WRITE(6,*) '          y_mid: ',trootS%xyz_mid(2)
         WRITE(6,*) '          z_mid: ',trootS%xyz_mid(3)
         WRITE(6,*) '         radius: ',trootS%radius   
         WRITE(6,*) '         torder: ',torder
         WRITE(6,*) '          theta: ',theta
         WRITE(6,*) '     maxparnode: ',maxparnodeS
         WRITE(6,*) ' '
 
      CALL CPU_TIME(timebeg)

!Call driver routine for cluster-particle
      CALL PC_TREECODE_FORCES(trootS, guv, &         !source info
                              xT, yT, zT, qT,  &     !target info
                              numparsS, numparsT, &  !number sources, targets
                              tEn)                   !output

      tEn = tEn * voxvol

      CALL CPU_TIME(timeend)
      timetree = timetree + timeend - timebeg

         WRITE(6,*) ' '
         WRITE(6,*) '   Finished calculation. '
         WRITE(6,*) '    Tree timing results: '
         WRITE(6,*) ' '
         WRITE(6,*) ' Tree creation time (s): ', timetree-totaltime
         WRITE(6,*) '  Treecode run time (s): ', totaltime
         WRITE(6,*) 'Treecode total time (s): ', timetree

! Call CLEANUP to deallocate global variables and tree structure.

         WRITE(6,*) ' '
         WRITE(6,*) 'Deallocating tree structure...'
         WRITE(6,*) ' '

      CALL CLEANUP(trootS)

      END SUBROUTINE TREECODE_FOR
