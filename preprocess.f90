SUBROUTINE preprocess
      use setup
      IMPLICIT NONE
      INCLUDE 'MESH2D.INC'

      INTEGER :: I1,IC1,IC2,IF2,K,NF,NF2,II,K2
      INTEGER :: VF1,VF2,IB
	  double precision:: XFIN1,XFIN2,YFIN1,YFIN2,XFOU1,XFOU2,YFOU1,YFOU2, eps1,eps2,eps3,eps4, YFCYCL1,YFCYCL2,YFCYCL3,YFCYCL4,YF1,YF2
	double precision :: XFCYCL1,XFCYCL2,XFCYCL3,XFCYCL4,XF1,XF2
	  allocate(BOUNFACEINDEX(NFACES))
	  
	  DO K=1,NFACES
	      BOUNFACEINDEX(K)=0
      END DO
	  

      NINLET=0
      NOUTLET=0
      NINTER=0
      NSYMM=0
      NWALL=0
      NWALC=0
      NWAC1=0
      NWAC2=0
	  NCYCL=0
      
      OPEN(110,FILE=NAMEBND2D)
      DO 110 II=1,NBOUN
      READ(110,*,END=1101) IB,IF2,IVBOUN(IB,1),IVBOUN(IB,2),IBOUNTYPE(IB)
	  
      IF (IBOUNTYPE(II).EQ.'INLE') THEN
              NINLET=NINLET+1
			  IBFINL(NINLET)= IF2
			  FACEMARKER(IF2)='INLE'

      END IF
      IF (IBOUNTYPE(II).EQ.'OUTL') THEN
              NOUTLET=NOUTLET+1
			  IBFOUT(NOUTLET) = IF2 
			  FACEMARKER(IF2)='OUTL'

      END IF
      IF (IBOUNTYPE(II).EQ.'SYMM') THEN
              NSYMM=NSYMM+1
			   IBFSYMM(NSYMM) = IF2
			   FACEMARKER(IF2)='SYMM'

      END IF
      IF (IBOUNTYPE(II).EQ.'WALL') THEN
              NWALL=NWALL+1
			  IBFWAL(NWALL) = IF2
			  FACEMARKER(IF2)='WALL'

      END IF
	  IF (IBOUNTYPE(II).EQ.'CYCL') THEN
              NCYCL=NCYCL+1
			  IBFCYCL(NCYCL) = IF2
			  FACEMARKER(IF2)='CYCL'

      END IF
     
110 CONTINUE

1101 CONTINUE



      CLOSE(110)
	  
	  
	    NINTER=0
      DO NF=1,NFACES
      IC2=IF2C(NF,2)
      IF(IC2 .le. NCELL) THEN
              NINTER=NINTER+1
              IBFINTER(NINTER)=NF
			  
      END IF
      END DO
	  
	  WRITE(*,*)'ninlet =',NINLET
      WRITE(*,*)'noutlet=',NOUTLET
      WRITE(*,*)'nsymm  =',NSYMM
      WRITE(*,*)'nwall  =',NWALL
      WRITE(*,*)'ncycl  =',NCYCL
      WRITE(*,*)'ninter =',NINTER
	  
	  NPAIR = 0
    !ALLOCATE(CYC2FPAIR(NCYCL/2,2))
	DO NF=1,NCYCL
	    K=IBFCYCL(NF)
		YFIN1=YV(IF2V(K,1))
		YFIN2=YV(IF2V(K,2))
			  DO NF2=1,NCYCL
			  IF (YFIN1 .NE. YFIN2) THEN
			     IF(NF2 .GT. NF) THEN
				     K2=IBFCYCL(NF2)
					 YFOU1=YV(IF2V(K2,1))
		             YFOU2=YV(IF2V(K2,2))
				     eps1=abs(YFIN1-YFOU1)
					 eps2=abs(YFIN1-YFOU2)
					 eps3=abs(YFIN2-YFOU1)
					 eps4=abs(YFIN2-YFOU2)
					 IF ((eps1 .lt. tolCYC .AND. eps4 .lt. tolCYC) .or.  (eps2 .lt. tolCYC .AND. eps3 .lt. tolCYC)) THEN
					       NPAIR=NPAIR+1
						   CYC2FPAIR(NPAIR,1)=K
						   CYC2FPAIR(NPAIR,2)=K2
				     END IF
				END IF
			END IF
			END DO
   END DO
   
   DO NF=1,NCYCL
	    K=IBFCYCL(NF)
		XFIN1=XV(IF2V(K,1))
		XFIN2=XV(IF2V(K,2))
			  DO NF2=1,NCYCL
			  IF (XFIN1 .NE. XFIN2) THEN
			     IF(NF2 .GT. NF) THEN
				     K2=IBFCYCL(NF2)
					 XFOU1=XV(IF2V(K2,1))
		             XFOU2=XV(IF2V(K2,2))
				     eps1=abs(XFIN1-XFOU1)
					 eps2=abs(XFIN1-XFOU2)
					 eps3=abs(XFIN2-XFOU1)
					 eps4=abs(XFIN2-XFOU2)
					 IF ((eps1 .lt. tolCYC .AND. eps4 .lt. tolCYC) .or.  (eps2 .lt. tolCYC .AND. eps3 .lt. tolCYC)) THEN
					       NPAIR=NPAIR+1
						   CYC2FPAIR(NPAIR,1)=K
						   CYC2FPAIR(NPAIR,2)=K2
				     END IF
				END IF
			END IF
			END DO
   END DO
					 
		write(*,*) 'npair =', NPAIR
	  
	 ! do nf =1,ncycl
	 !     YFCYCL1=YV(IF2V(IBFCYCL(nf),1))
	!	  YFCYCL2=YV(IF2V(IBFCYCL(nf),2))
	!	 ! YF1=0.5*(YFCYCL1+YFCYCL2)
	!	  do nf2=1,ncycl
	!	     YFCYCL3=YV(IF2V(IBFCYCL(nf2),1))
	!	     YFCYCL4=YV(IF2V(IBFCYCL(nf2),2))
	!	 ! YF2=0.5*(YFCYCL3+YFCYCL4)
	!	  eps1=abs(YFCYCL1-YFCYCL3)
	!	  eps2=abs(YFCYCL1-YFCYCL4)
	!	  eps3=abs(YFCYCL2-YFCYCL3)
	!	  eps4=abs(YFCYCL2-YFCYCL4)
	!	  if (YFCYCL1 .NE. YFCYCL2 .AND. IBFCYCL(nf).NE.IBFCYCL(nf2)) then
	!	       if ((eps1.LT. 1.e-3 .AND. eps4.LT. 1.e-3).OR.(eps2.LT. 1.e-3 .AND. eps3.LT. 1.e-3))  then
	!	          FPAIRCYCL(nf)=IBFCYCL(nf2)
	!			  IF2C(IBFCYCL(nf),2) = IF2C(FPAIRCYCL(nf),1)
	!			  IF2C(IBFCYCL(nf),4)=IF2C(FPAIRCYCL(nf),3)
	!	!		  write (*,*)'boundary face=', IBFCYCL(nf), 'boundary cell=', IF2C(IBFCYCL(nf),1), &
	!	!'boundary face pair', FPAIRCYCL(nf), 'boundary pair cell=', IF2C(IBFCYCL(nf),2)
	!	       end if
	!	  end if
	!	 end do
	 !end do
	 
	 
	 !do nf =1,ncycl
	 !     XFCYCL1=XV(IF2V(IBFCYCL(nf),1))
	!	  XFCYCL2=XV(IF2V(IBFCYCL(nf),2))
	!	 ! YF1=0.5*(YFCYCL1+YFCYCL2)
	!	  do nf2=1,ncycl
	!	     XFCYCL3=XV(IF2V(IBFCYCL(nf2),1))
	!	     XFCYCL4=XV(IF2V(IBFCYCL(nf2),2))
	!	 ! YF2=0.5*(YFCYCL3+YFCYCL4)
	!	  eps1=abs(XFCYCL1-XFCYCL3)
	!	  eps2=abs(XFCYCL1-XFCYCL4)
	!	  eps3=abs(XFCYCL2-XFCYCL3)
	!	  eps4=abs(XFCYCL2-XFCYCL4)
	!	  if (XFCYCL1 .NE. XFCYCL2 .AND. IBFCYCL(nf).NE.IBFCYCL(nf2)) then
	!	       if ((eps1.LT. 1.e-3 .AND. eps4.LT. 1.e-3).OR.(eps2.LT. 1.e-3 .AND. eps3.LT. 1.e-3))  then
	!	          FPAIRCYCL(nf)=IBFCYCL(nf2)
	!			  IF2C(IBFCYCL(nf),2) = IF2C(FPAIRCYCL(nf),1)
	!			  IF2C(IBFCYCL(nf),4)=IF2C(FPAIRCYCL(nf),3)
	!	!		  write (*,*)'boundary face=', IBFCYCL(nf), 'boundary cell=', IF2C(IBFCYCL(nf),1), &
	!	!'boundary face pair', FPAIRCYCL(nf), 'boundary pair cell=', IF2C(IBFCYCL(nf),2)
	!	       end if
	!	  end if
	!	 end do
	 !end do
	  
	  
	! do nf=1,ncycl
	!    !IF2C(IBFCYCL(nf),2) = IF2C(FPAIRCYCL(nf),1)
	!	!IF2C(IBFCYCL(nf),4)=IF2C(FPAIRCYCL(nf),3)
	!	!write (*,*) 'boundary face=', IBFCYCL(nf)
	!	write (*,*)'boundary face=', IBFCYCL(nf), 'boundary cell=', IF2C(IBFCYCL(nf),1), &
	!	'boundary face pair', FPAIRCYCL(nf), 'boundary pair cell=', IF2C(IBFCYCL(nf),2)
	!end do
	!  
	! 
	! do nf=1,ncycl
	!   write(*,*) 'boundary cell', IF2C(IBFCYCL(nf),1),'pair cell', IF2C(IBFCYCL(nf),2)
	! end do

	
    
   !write(*,*) '  PREPROCESS NINTER=', NINTER

	  
	  END SUBROUTINE preprocess










