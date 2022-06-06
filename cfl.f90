SUBROUTINE calc_cfl
use setup
implicit none
include 'MESH2D.INC'
     double precision :: rho,u,v,ww,Bx,By,Bz,rhoE,p
	 double precision :: temp1,temp2,temp3,temp
	 double precision :: rlambdamax,NORM_U,NORM_B
	 integer :: ic
	 
	 rlambdamax=0
	 do ic=1,NCELL
	    rho=u0(1,ic)
		u=u0(2,ic)/rho
		v=u0(3,ic)/rho
		ww=u0(4,ic)/rho
		Bx=u0(6,ic)
		By=u0(7,ic)
		Bz=u0(8,ic)
		rhoE=u0(5,ic)
		NORM_U = u**2 + v**2 + ww**2
		NORM_B = Bx**2 + By**2 + Bz**2
		p=(gama-1)*(rhoE-0.5*rho*NORM_U-0.5*NORM_B)

!		write(*,*)
!		write(*,*)
!		write(*,*)
!		write(*,*)"iteration:",iter
!		write(*,*)"rho:",rho
!		write(*,*)"u:",u
!		write(*,*)"v:",v
!		write(*,*)"w:",ww
!		write(*,*)"Bx:",Bx
!		write(*,*)"By:",By
!		write(*,*)"Bz:",Bz
!		write(*,*)"rhoE:",rhoE
!		write(*,*)"p:",p
!		write(*,*)
!		write(*,*)
!		write(*,*)

		temp1=abs(u)+sqrt(gama*p/rho)
		temp2=abs(v)+sqrt(gama*p/rho)
		temp3=abs(ww)+sqrt(gama*p/rho)
		temp=sqrt(temp1*temp1 + Bx*Bx + temp2*temp2 + By*By + temp3*temp3 + Bz*Bz)
!		write(*,*)
!		write(*,*)
!		write(*,*)
!		write(*,*)"temp1:",temp1
!		write(*,*)"temp2:",temp2
!		write(*,*)"temp:",temp
!		write(*,*)
!		write(*,*)
!		write(*,*)
		if (temp .gt. rlambdamax) rlambdamax=temp
	 end do
	 
!	 write(*,*)"rlambdamax:",rlambdamax

	 dt=minl*CFL/rlambdamax

!	 write(*,*)"dt:",dt
	 end
	 
	  SUBROUTINE calc_minl
	     use setup
		 implicit none
		 include 'MESH2D.INC'
	     integer :: ic
		 DOUBLE PRECISION :: ax,ay,bx,by
	     DOUBLE PRECISION :: vecprod
		 
		 do ic=1,NCELL
		     ax=XV(IVCELL(ic,3))-XV(IVCELL(ic,1))
			 ay=YV(IVCELL(ic,3))-YV(IVCELL(ic,1))
			 bx=XV(IVCELL(ic,4))-XV(IVCELL(ic,2))
			 by=YV(IVCELL(ic,4))-YV(IVCELL(ic,2))
			 
			 vecprod=ax*by-ay*bx
			 if (vecprod < 0) write(*,*) 'something wrong my friend'
			 area(ic)=0.5*vecprod
	     end do
			 
		
	     minl=sqrt(area(1))
		 do ic=1,NCELL
		    if(minl.gt.sqrt(area(ic))) minl=sqrt(area(ic))
		 end do
		  
	 END SUBROUTINE calc_minl
	  
		 

