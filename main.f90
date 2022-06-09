PROGRAM main

      use setup
      IMPLICIT NONE
      INCLUDE 'MESH2D.INC'

      INTEGER :: i,j,k,ic,checkstat,flag
	  double precision :: t,t2,residnorm
	  double precision :: proctime,timestart,timestop
	  

! initialization	  
	  CALL READ_INPUT
	   CALL READ_CELL_DATA
	  CALL READ_VRT_DATA
	  !CALL READ_CELL_DATA
	  CALL CONNECTIVITY
	  CALL preprocess
	  
	  CALL initialize
	   
	  !CALL calcjacob
	  
	  CALL calc_minl
	 
! opening files and writing headers
	  open(11,file="resid.out")
	  open(9001,file="record.out")
           
	write(11,*) "    ITER        RESIDNORM         TIME"


	write(9001,*) "For some reason, Fortran won't run without a residnorm NAN unless"
	write(9001,*) "we report the dt for each iteration, so this file records that."
	write(9001,*) "Even though this file is intended to prevent a residnorm NAN,"
	write(9001,*) "it can still be used to record other useful measurements as well."
	write(9001,*)
	write(9001,*) "       iter      residnorm             dt                         t"

! specify final time and initialize local variables
	  ft=600
	  t=0
	  kcount=0
	  checkstat = 0
	  iter = 0
	  flag = 0
	  call CPU_TIME(timestart)


!		save original values if the t is decreased later in the program	  
!
!     MHD
!	  do ic=1,NCELL
!		do k=1,8
!		   o0(k,ic)=u0(k,ic)
!		   ox(k,ic)=ux(k,ic)
!		   oy(k,ic)=uy(k,ic)
!		   oxx(k,ic)=uxx(k,ic)
!		   oxy(k,ic)=uxy(k,ic)
!		   oyy(k,ic)=uyy(k,ic)
!		end do
!	  end do

!	  Euler
	  do ic=1,NCELL
		do k=1,4
		   o0(k,ic)=u0(k,ic)
		   ox(k,ic)=ux(k,ic)
		   oy(k,ic)=uy(k,ic)
		   oxx(k,ic)=uxx(k,ic)
		   oxy(k,ic)=uxy(k,ic)
		   oyy(k,ic)=uyy(k,ic)
		end do
	  end do


! loop of iterations
	  do while(t .lt. ft)
	       iter=iter+1
		   if(flag.eq.0) call calc_cfl
!dt = ndt
		if (t+dt .gt. ft) dt=ft-t
		t=t+dt
		!kcount=kcount+1

		write(9001,*)iter,residnorm,dt,t

!		MHD
! save previous iteration's result
!		do ic=1,NCELL
!		   do k=1,8
!		      w10(k,ic)=u0(k,ic)
!			  w1x(k,ic)=ux(k,ic)
!			  w1y(k,ic)=uy(k,ic)
!			  w1xx(k,ic)=uxx(k,ic)
!			  w1xy(k,ic)=uxy(k,ic)
!			  w1yy(k,ic)=uyy(k,ic)
!		   end do
!		 end do
!
!
!
!! Runge-Kutta step 1
!		!CALL Detector
!	    CALL calright
!		do ic=1,NCELL
!		   do k=1,8
!		     u0(k,ic)=w10(k,ic)+dt*right0(k,ic)
!			 ux(k,ic)=w1x(k,ic)+dt*rightx(k,ic)
!			 uy(k,ic)=w1y(k,ic)+dt*righty(k,ic)
!			 uxx(k,ic)=w1xx(k,ic)+dt*rightxx(k,ic)
!			 uxy(k,ic)=w1xy(k,ic)+dt*rightxy(k,ic)
!			 uyy(k,ic)=w1yy(k,ic)+dt*rightyy(k,ic)
!		   end do
!		end do
!		   
!
!
!! Runge-Kutta step 2
!		! CALL Detector  
!		 CALL calright
!
!		 do ic=1,NCELL
!		   do k=1,8
!		     t2=u0(k,ic)+dt*right0(k,ic)
!		     u0(k,ic)=0.75*w10(k,ic)+0.25*t2
!			 
!			 t2=ux(k,ic)+dt*rightx(k,ic)
!			 ux(k,ic)=0.75*w1x(k,ic)+0.25*t2
!			 
!			 t2=uy(k,ic)+dt*righty(k,ic)
!			 uy(k,ic)=0.75*w1y(k,ic)+0.25*t2
!			 
!			 t2=uxx(k,ic)+dt*rightxx(k,ic)
!			 uxx(k,ic)=0.75*w1xx(k,ic)+0.25*t2
!			 
!			 t2=uxy(k,ic)+dt*rightxy(k,ic)
!			 uxy(k,ic)=0.75*w1xy(k,ic)+0.25*t2
!			 
!			 t2=uyy(k,ic)+dt*rightyy(k,ic)
!			 uyy(k,ic)=0.75*w1yy(k,ic)+0.25*t2
!		   end do
!		end do
!		
!
!
!! Runge-Kutta step 3
!	!	CALL Detector 
!		CALL calright
!		 do ic=1,NCELL
!		   do k=1,8
!		     t2=u0(k,ic)+dt*right0(k,ic)
!		     u0(k,ic)=1./3.*w10(k,ic)+2./3.*t2
!			 
!			 t2=ux(k,ic)+dt*rightx(k,ic)
!			 ux(k,ic)=1./3.*w1x(k,ic)+2./3.*t2
!			 
!			 t2=uy(k,ic)+dt*righty(k,ic)
!			 uy(k,ic)=1./3.*w1y(k,ic)+2./3.*t2
!			 
!			 t2=uxx(k,ic)+dt*rightxx(k,ic)
!			 uxx(k,ic)=1./3.*w1xx(k,ic)+2./3.*t2
!			 
!			 t2=uxy(k,ic)+dt*rightxy(k,ic)
!			 uxy(k,ic)=1./3.*w1xy(k,ic)+2./3.*t2
!			 
!			 t2=uyy(k,ic)+dt*rightyy(k,ic)
!			 uyy(k,ic)=1./3.*w1yy(k,ic)+2./3.*t2
!		   end do
!		end do
!
!
!		  
!! calculate residnorm
!		residnorm =0
!		do ic=1,NCELL
!		   do k=1,5
!		      residnorm=residnorm +abs( right0(k,ic))
!		   end do
!		end do
!		residnorm=residnorm/(NCELL*5)

!		Euler
! save previous iteration's result
		do ic=1,NCELL
			do k=1,4
			   w10(k,ic)=u0(k,ic)
			   w1x(k,ic)=ux(k,ic)
			   w1y(k,ic)=uy(k,ic)
			   w1xx(k,ic)=uxx(k,ic)
			   w1xy(k,ic)=uxy(k,ic)
			   w1yy(k,ic)=uyy(k,ic)
			end do
		  end do
 
 
 
 ! Runge-Kutta step 1
		  RK = 1
		 !CALL Detector
		 CALL calright
		 do ic=1,NCELL
			do k=1,4
			  u0(k,ic)=w10(k,ic)+dt*right0(k,ic)
			  ux(k,ic)=w1x(k,ic)+dt*rightx(k,ic)
			  uy(k,ic)=w1y(k,ic)+dt*righty(k,ic)
			  uxx(k,ic)=w1xx(k,ic)+dt*rightxx(k,ic)
			  uxy(k,ic)=w1xy(k,ic)+dt*rightxy(k,ic)
			  uyy(k,ic)=w1yy(k,ic)+dt*rightyy(k,ic)
			end do
		 end do
 
 
 ! Runge-Kutta step 2
		 RK = 2
		 ! CALL Detector  
		  CALL calright
 
		  do ic=1,NCELL
			do k=1,4
			  t2=u0(k,ic)+dt*right0(k,ic)
			  u0(k,ic)=0.75*w10(k,ic)+0.25*t2
			  
			  t2=ux(k,ic)+dt*rightx(k,ic)
			  ux(k,ic)=0.75*w1x(k,ic)+0.25*t2
			  
			  t2=uy(k,ic)+dt*righty(k,ic)
			  uy(k,ic)=0.75*w1y(k,ic)+0.25*t2
			  
			  t2=uxx(k,ic)+dt*rightxx(k,ic)
			  uxx(k,ic)=0.75*w1xx(k,ic)+0.25*t2
			  
			  t2=uxy(k,ic)+dt*rightxy(k,ic)
			  uxy(k,ic)=0.75*w1xy(k,ic)+0.25*t2
			  
			  t2=uyy(k,ic)+dt*rightyy(k,ic)
			  uyy(k,ic)=0.75*w1yy(k,ic)+0.25*t2
			end do
		 end do
		 
 
 
 ! Runge-Kutta step 3
		 RK = 3
	 !	CALL Detector 
		 CALL calright
		  do ic=1,NCELL
			do k=1,4
			  t2=u0(k,ic)+dt*right0(k,ic)
			  u0(k,ic)=1./3.*w10(k,ic)+2./3.*t2
			  
			  t2=ux(k,ic)+dt*rightx(k,ic)
			  ux(k,ic)=1./3.*w1x(k,ic)+2./3.*t2
			  
			  t2=uy(k,ic)+dt*righty(k,ic)
			  uy(k,ic)=1./3.*w1y(k,ic)+2./3.*t2
			  
			  t2=uxx(k,ic)+dt*rightxx(k,ic)
			  uxx(k,ic)=1./3.*w1xx(k,ic)+2./3.*t2
			  
			  t2=uxy(k,ic)+dt*rightxy(k,ic)
			  uxy(k,ic)=1./3.*w1xy(k,ic)+2./3.*t2
			  
			  t2=uyy(k,ic)+dt*rightyy(k,ic)
			  uyy(k,ic)=1./3.*w1yy(k,ic)+2./3.*t2
			end do
		 end do
 
 
		   
 ! calculate residnorm
		 residnorm =0
		 do ic=1,NCELL
			do k=1,4
			   residnorm=residnorm + abs(right0(k,ic))
			end do
		 end do
		 residnorm=residnorm/(NCELL*4)


!		if(flag.eq.0 .and. residnorm.ne.residnorm) then
!			t = 0
!			dt = ndt
!			iter = 0
!			write(*,*)"dt switched to defined dt:",ndt
!			flag = 1
!
!			do ic=1,NCELL
!				do k=1,8
!				   u0(k,ic)=o0(k,ic)
!				   ux(k,ic)=ox(k,ic)
!				   uy(k,ic)=oy(k,ic)
!				   uxx(k,ic)=oxx(k,ic)
!				   uxy(k,ic)=oxy(k,ic)
!				   uyy(k,ic)=oyy(k,ic)
!				end do
!			  end do
!
!			cycle
!		end if

! stop the program if something goes wrong

		if(abs(residnorm) .le. 1.e-16) then!checkstat=1
			write(*,*)"|Residnorm| less than 1E-16    (ie: Nothing seems to be changing)"
			exit
	    else if(iter .gt. MAXITER) then!checkstat=1
			write(*,*)"Maximum allowed iterations"
			exit
		else if(residnorm .ne. residnorm) then
			write(*,*)"Residnorm NAN"
			write(*,*)"Residnorm:",residnorm
			exit
		end if
		

! calculate the time since the program was started
		call CPU_TIME(timestop)
		proctime=timestop-timestart

! update on progress every 20 iterations
		if(mod(iter,20)==0) then
			write(11,111) iter, residnorm, proctime
            write(*,*) 'At iteration ',iter,'residual norm is ',residnorm
		end if
     
! plot every 100 iterations
		 if (mod(iter,100)==0) then 
            call tecplotter(NRE)
		 !call tecplotnorefinement
         end if
		end do


! calculate time program took to run
		call CPU_TIME(timestop)
		proctime=timestop-timestart

! write info after the program finishes
				 WRITE(*,*)
	     write(*,*) " COMPLETED ", iter, " ITERATIONS"
		 write(*,*) "t = ",t,"after running for",proctime,"seconds"

	close(11)
	close(9001)
111     format(I9,1X,3(G16.9,1X))
		
		END PROGRAM main
		
		
		   
		   
		   
				 
		
	  
	  

      
