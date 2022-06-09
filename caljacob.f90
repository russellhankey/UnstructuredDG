SUBROUTINE calcjacob
      use setup
      IMPLICIT NONE
      INCLUDE 'MESH2D.INC'

      INTEGER :: is,js,K,ic,ifp,jfp,I,J,m
      DOUBLE PRECISION :: psi,eta,xpsi,xeta,ypsi,yeta
      INTEGER,DIMENSION(4) :: IC2V
	  
	  double precision, dimension(2) :: tecnode(2)
	  tecnode(1)=-1/6
	  tecnode(2)=1/6

      DO ic=1,NCELL
          DO m=1,4
          IC2V(m)=IVCELL(ic,m)
          END DO
		 

      DO K=1,4
      DO J=1,3
      IF(K.EQ.1 .AND. J.EQ.1) THEN
              psi=gp(1)
              eta=-0.5
      ELSE IF (K.EQ.1 .AND. J.EQ.2) THEN
             psi=gp(2)
             eta=-0.5
     ELSE IF(K.EQ.1 .AND. J.EQ.3) THEN
             psi=gp(3) 
              eta=-0.5
      END IF
     
      IF(K.EQ.2 .AND. J.EQ.1) THEN
              psi=0.5
              eta=gp(1)
      ELSE IF (K.EQ.2 .AND. J.EQ.2) THEN
             eta=gp(2)
             psi=0.5
     ELSE IF(K.EQ.2 .AND. J.EQ.3) THEN
             eta=gp(3)
             psi=0.5
      END IF

      IF(K.EQ.3 .AND. J.EQ.1) THEN
              psi=gp(3)
              eta=0.5
      ELSE IF (K.EQ.3 .AND. J.EQ.2) THEN
             psi=gp(2)
             eta=0.5
     ELSE IF(K.EQ.3 .AND. J.EQ.3) THEN
             psi=gp(1)
              eta=0.5
      END IF

      IF(K.EQ.4 .AND. J.EQ.1) THEN
              psi=-0.5
              eta=gp(3)
      ELSE IF (K.EQ.4 .AND. J.EQ.2) THEN
             eta=gp(2)
             psi=-0.5
     ELSE IF(K.EQ.4 .AND. J.EQ.3) THEN
             eta=gp(1)
             psi=-0.5
      END IF
      
      xpsi=0.0
      xeta=0.0
      ypsi=0.0
      yeta=0.0
      xpsi=(eta-0.5)*XV(IC2V(1))+(0.5-eta)*XV(IC2V(2)) &
              + (0.5+eta)*XV(IC2V(3))-(eta+0.5)*XV(IC2V(4))
      xeta=(psi-0.5)*XV(IC2V(1))-(0.5+psi)*XV(IC2V(2)) &
              + (0.5+psi)*XV(IC2V(3))+(0.5-psi)*XV(IC2V(4))
      ypsi=(eta-0.5)*YV(IC2V(1))+(0.5-eta)*YV(IC2V(2)) &
              + (0.5+eta)*YV(IC2V(3))-(eta+0.5)*YV(IC2V(4))
      yeta=(psi-0.5)*YV(IC2V(1))-(0.5+psi)*YV(IC2V(2)) &
              + (0.5+psi)*YV(IC2V(3))+(0.5-psi)*YV(IC2V(4))
			  
			  
         djacobf(1,J,K,ic)=xpsi
	 djacobf(2,J,K,ic)=xeta
	 djacobf(3,J,K,ic)=ypsi
	 djacobf(4,J,K,ic)=yeta

      jacobf(J,K,ic)=xpsi*yeta-xeta*ypsi
	  
      END DO
      END DO
	  
!	  DO K=1,4
!      DO J=1,2
!      IF(K.EQ.1 .AND. J.EQ.1) THEN
!              psi=-1/6
!              eta=-0.5
!      ELSE IF (K.EQ.1 .AND. J.EQ.2) THEN
!             psi=1/6
!             eta=-0.5
!
!      END IF
!     
!      IF(K.EQ.2 .AND. J.EQ.1) THEN
!              psi=0.5
!              eta=-1/6
!      ELSE IF (K.EQ.2 .AND. J.EQ.2) THEN
!             eta=1/6
!             psi=0.5
!     
!      END IF
!
!      IF(K.EQ.3 .AND. J.EQ.1) THEN
!              psi=1/6
!              eta=0.5
!      ELSE IF (K.EQ.3 .AND. J.EQ.2) THEN
!             psi=-1/6
!             eta=0.5
!     
!      END IF
!
!      IF(K.EQ.4 .AND. J.EQ.1) THEN
!              psi=-0.5
!              eta=1./.6
!      ELSE IF (K.EQ.4 .AND. J.EQ.2) THEN
!             eta=-1./6.
!             psi=-0.5
!     
!      END IF
!      
!      xpsi=0.0
!      xeta=0.0
!      ypsi=0.0
!      yeta=0.0
!      xpsi=(eta-0.5)*XV(IC2V(1))+(0.5-eta)*XV(IC2V(2)) &
!              + (0.5+eta)*XV(IC2V(3))-(eta+0.5)*XV(IC2V(4))
!      xeta=(psi-0.5)*XV(IC2V(1))-(0.5+psi)*XV(IC2V(2)) &
!              + (0.5+psi)*XV(IC2V(3))+(0.5-psi)*XV(IC2V(4))
!      ypsi=(eta-0.5)*YV(IC2V(1))+(0.5-eta)*YV(IC2V(2)) &
!              + (0.5+eta)*YV(IC2V(3))-(eta+0.5)*YV(IC2V(4))
!      yeta=(psi-0.5)*YV(IC2V(1))-(0.5+psi)*YV(IC2V(2)) &
!              + (0.5+psi)*YV(IC2V(3))+(0.5-psi)*YV(IC2V(4))
!			  
!			  
!     !djacobf(1,J,K,ic)=xpsi
!	 !djacobf(2,J,K,ic)=xeta
!	 !djacobf(3,J,K,ic)=ypsi
!	 !djacobf(4,J,K,ic)=yeta
!
!      jacobtecint(J,K,ic)=xpsi*yeta-xeta*ypsi
!	 
!	  
!      END DO
!      END DO
      
	  
      DO I=1,3
        DO J=1,3
           psi=gp(I)
           eta=gp(J)
           xpsi=0.0
      xeta=0.0
      ypsi=0.0
      yeta=0.0
      xpsi=(eta-0.5)*XV(IC2V(1))+(0.5-eta)*XV(IC2V(2)) &
              + (0.5+eta)*XV(IC2V(3))-(eta+0.5)*XV(IC2V(4))
      xeta=(psi-0.5)*XV(IC2V(1))-(0.5+psi)*XV(IC2V(2)) &
              + (0.5+psi)*XV(IC2V(3))+(0.5-psi)*XV(IC2V(4))
      ypsi=(eta-0.5)*YV(IC2V(1))+(0.5-eta)*YV(IC2V(2)) &
              + (0.5+eta)*YV(IC2V(3))-(eta+0.5)*YV(IC2V(4))
      yeta=(psi-0.5)*YV(IC2V(1))-(0.5+psi)*YV(IC2V(2)) &
              + (0.5+psi)*YV(IC2V(3))+(0.5-psi)*YV(IC2V(4))
      
	  djacobv(1,I,J,ic)=xpsi
	  djacobv(2,I,J,ic)=xeta
	  djacobv(3,I,J,ic)=ypsi
	  djacobv(4,I,J,ic)=yeta
      jacobv(I,J,ic)=xpsi*yeta-xeta*ypsi
      END DO
      END DO
	  
	  
	  DO I=1,4
	     if (I.eq.1) then
		     psi=-0.5
			 eta=-0.5
		 else if (I.eq.2) then
		    psi=0.5
			eta=-0.5
		else if (I.eq.3) then
		    psi=0.5
			eta=0.5
		else if (I.eq.4) then
		    psi=-0.5
			eta=0.5  
		end if 
		xpsi=(eta-0.5)*XV(IC2V(1))+(0.5-eta)*XV(IC2V(2)) &
              + (0.5+eta)*XV(IC2V(3))-(eta+0.5)*XV(IC2V(4))
        xeta=(psi-0.5)*XV(IC2V(1))-(0.5+psi)*XV(IC2V(2)) &
              + (0.5+psi)*XV(IC2V(3))+(0.5-psi)*XV(IC2V(4))
        ypsi=(eta-0.5)*YV(IC2V(1))+(0.5-eta)*YV(IC2V(2)) &
              + (0.5+eta)*YV(IC2V(3))-(eta+0.5)*YV(IC2V(4))
        yeta=(psi-0.5)*YV(IC2V(1))-(0.5+psi)*YV(IC2V(2)) &
              + (0.5+psi)*YV(IC2V(3))+(0.5-psi)*YV(IC2V(4))
			  
        djacobvertex(1,I,ic)=xpsi
	    djacobvertex(2,I,ic)=xeta
	    djacobvertex(3,I,ic)=ypsi
	    djacobvertex(4,I,ic)=yeta
        jacobvertex(I,ic)=xpsi*yeta-xeta*ypsi

      END DO
	  
!	   DO I=1,2
!        DO J=1,2
!           psi=tecnode(I)
!           eta=tecnode(J)
!           xpsi=0.0
!      xeta=0.0
!      ypsi=0.0
!      yeta=0.0
!      xpsi=(eta-0.5)*XV(IC2V(1))+(0.5-eta)*XV(IC2V(2)) &
!              + (0.5+eta)*XV(IC2V(3))-(eta+0.5)*XV(IC2V(4))
!      xeta=(psi-0.5)*XV(IC2V(1))-(0.5+psi)*XV(IC2V(2)) &
!              + (0.5+psi)*XV(IC2V(3))+(0.5-psi)*XV(IC2V(4))
!      ypsi=(eta-0.5)*YV(IC2V(1))+(0.5-eta)*YV(IC2V(2)) &
!              + (0.5+eta)*YV(IC2V(3))-(eta+0.5)*YV(IC2V(4))
!      yeta=(psi-0.5)*YV(IC2V(1))-(0.5+psi)*YV(IC2V(2)) &
!              + (0.5+psi)*YV(IC2V(3))+(0.5-psi)*YV(IC2V(4))
!      
!	  !djacobv(1,I,J,ic)=xpsi
!	  !djacobv(2,I,J,ic)=xeta
!	  !djacobv(3,I,J,ic)=ypsi
!	  !djacobv(4,I,J,ic)=yeta
!      jacobtecv(I,J,ic)=xpsi*yeta-xeta*ypsi
!
!      END DO
!      END DO
	  
END DO

END 

