SUBROUTINE READ_INPUT
      use setup
      IMPLICIT NONE
      INCLUDE 'MESH2D.INC'
      OPEN(30,FILE='QUAD.INP')
	      READ(30,*)
		  READ(30,*)
		  READ(30,*)  NAMECEL
		  READ(30,*)  NAMEVRT
		  READ(30,*)  NAMEBND2D
		  READ(30,*)
		  READ(30,*)  N
		  READ(30,*)
		  READ(30,*)  NRE
		  READ(30,*)
		  READ(30,*) rinf
		  READ(30,*) pinf
		  READ(30,*) uinf
		  READ(30,*) vinf
		  READ(30,*)
		  READ(30,*) MAXITER
          READ(30,*)	
          READ(30,*) CFL	
          READ(30,*)	
          READ(30,*) NDT
		  READ(30,*)
		  READ(30,*)	restart
          READ(30,*) NAMERESTART	
          READ(30,*) restart_ord	
	  CLOSE(30)
	  RETURN
      END SUBROUTINE READ_INPUT



      SUBROUTINE READ_VRT_DATA
!----Read coordinates of vertices (XV,YV,ZV)
      use setup
      IMPLICIT NONE
      INTEGER :: I,IV
      INCLUDE 'MESH2D.INC'
!
      OPEN(10,FILE=NAMEVRT)
      DO 10 I=1,MAXVRT
      READ(10,*,END=1001) IV,XV(IV),YV(IV)
      IF (MOD(IV,10).EQ.0) WRITE(*,*) 'reading vertex=',IV
10    CONTINUE
1000  FORMAT(I9,6X,2(D20.12,1X))
1001  CONTINUE
      CLOSE(10)
      NVERT=IV
      WRITE(*,*) 'Maximum value of vertex number=',NVERT
!
      RETURN
      END


	  
	  SUBROUTINE READ_CELL_DATA
      use setup
      IMPLICIT NONE
      INCLUDE 'MESH2D.INC'
      INTEGER :: I,IC,K,ICKEY
	  

!
!----Read fluid cell definitions
!----Make sure that ONLY FLUID CELLS ARE STORED

      OPEN(20,FILE=NAMECEL)
      DO 20 I=1,MAXCEL
! IVCELL for vertex numbers of one cell
      READ(20,*,END=2001) IC,(IVCELL(I,K),K=1,4), &
      ICTYPE(I),ICKEY
      ICLMAP(I)=IC
! Currently linear mapping (1 to 1)
      IF (MOD(I,100).EQ.0) WRITE(*,*) 'reading cell=',I
20    CONTINUE
2001  CONTINUE
      CLOSE(20)
2000  FORMAT(I9,6X,5(I9,1X),I4)
      NCELL=I-1
      WRITE(*,*) 'Number of cells=',NCELL
	     
!
      RETURN
      END


      !SUBROUTINE READ_CELL_DATA
      !        IMPLICIT NONE
      !        INCLUDE 'MESH2D.INC'
      !        INTEGER :: I,IC,K,ICKEY

      !        OPEN(200,FILE=NAMECEL)
      !        DO 20 I=1,MAXCEL
      !            READ(20,*,END=2001) IC,(IVCELL(I,K),K=1,4),&
      !                    ICTYPE(I),ICKEY
      !            ICLMAP(I)=IC
      !            IF (MOD(I,1000).EQ.0) WRITE(*,*)'reading cell=',I
      !        20 CONTINUE
      !        2001 CONTINUE
      !        CLOSE(200)
      !        NCELL=I-1
      !        WRITE(*,*) 'Number of cells=',NCELL
      ! RETURN
      ! END

       
