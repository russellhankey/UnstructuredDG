SUBROUTINE initialize
      use setup
	  implicit none
	  INCLUDE 'MESH2D.INC'

      INTEGER :: I,K,ic,igp,jgp,iv
      DOUBLE PRECISION :: rho,u,v,ww,p,rhoE,Bx,By,Bz,x,y
	  DOUBLE PRECISION :: NORM_B,pr,NORM_U
	  double precision, dimension(:,:)	::   xxs(2,4)
	  double precision,dimension(2) :: xx
	  double precision :: xr, yr
	  double precision :: psi,eta
	  double precision :: ur,vr,V2,theta
	  double precision :: tt

!	  MHD
!	  double precision :: Qvi(8)

!     Euler
	  double precision :: Qvi(4)
	  
	  
	  !double precision	:: rinf,pinf,uinf,vinf
	  !rinf=1.0
	  !pinf=7.9365
	  !uinf=1.0
	  !vinf=0.0
	  
	  theta=ATAN(0.5)
	  
	  
!	  MHD
!	  ALLOCATE (u0(1:8,1:NCELL))
!      ALLOCATE (ux(1:8,1:NCELL))
!      ALLOCATE (uy(1:8,1:NCELL))
!      ALLOCATE (uxx(1:8,1:NCELL))
!      ALLOCATE (uxy(1:8,1:NCELL))
!      ALLOCATE (uyy(1:8,1:NCELL))
!
!      ALLOCATE (w10(1:8,1:NCELL))
!      ALLOCATE (w1x(1:8,1:NCELL))
!      ALLOCATE (w1y(1:8,1:NCELL))
!      ALLOCATE (w1xx(1:8,1:NCELL))
!      ALLOCATE (w1xy(1:8,1:NCELL))
!      ALLOCATE (w1yy(1:8,1:NCELL))
!
!      ALLOCATE (right0(1:8,1:NCELL))
!      ALLOCATE (rightx(1:8,1:NCELL))
!      ALLOCATE (righty(1:8,1:NCELL))
!      ALLOCATE (rightxx(1:8,1:NCELL))
!      ALLOCATE (rightxy(1:8,1:NCELL))
!      ALLOCATE (rightyy(1:8,1:NCELL))
!
!	  allocate (o0(1:8,1:NCELL))
!	  allocate (ox(1:8,1:NCELL))
!	  allocate (oy(1:8,1:NCELL))
!	  allocate (oxx(1:8,1:NCELL))
!	  allocate (oxy(1:8,1:NCELL))
!	  allocate (oyy(1:8,1:NCELL))
!	  
!	  allocate(gp(1:3))
!	  allocate(w(1:3))
!	  allocate(jacobvertex(1:4,1:NCELL))
!	  allocate(djacobvertex(1:4,1:4,1:NCELL))
!	  allocate(djacobf(1:4,1:3,1:4,1:NCELL))
!	  allocate(jacobf(1:3,1:4,1:NCELL))
!	  allocate(djacobv(1:4,1:3,1:3,1:NCELL))
!	  allocate(jacobv(1:3,1:3,1:NCELL))
!	  allocate(f_edge(1:8,1:3,1:4,1:NCELL))
!	  allocate(area(1:NCELL))
!	  allocate(rmass(6))
!	  
!!	  allocate(jacobtecint(1:2,1:4,1:NCELL))
!!	  allocate(jacobtecv(1:2,1:2,1:NCELL))
!	  
!	  allocate(Qv(1:8,1:3,1:3,1:NCELL))
!	  allocate(Qub(1:8,1:3,1:2,1:NCELL))
!	  allocate(Qlr(1:8,1:3,1:2,1:NCELL))
!	  allocate(Qa(1:8,1:3,1:4,1:NCELL))
!	  
!	  allocate(iface2fp(N,4))
!	  allocate(jface2fp(N,4))

!	  Euler
	  ALLOCATE (u0(1:4,1:NCELL))
      ALLOCATE (ux(1:4,1:NCELL))
      ALLOCATE (uy(1:4,1:NCELL))
      ALLOCATE (uxx(1:4,1:NCELL))
      ALLOCATE (uxy(1:4,1:NCELL))
      ALLOCATE (uyy(1:4,1:NCELL))

      ALLOCATE (w10(1:4,1:NCELL))
      ALLOCATE (w1x(1:4,1:NCELL))
      ALLOCATE (w1y(1:4,1:NCELL))
      ALLOCATE (w1xx(1:4,1:NCELL))
      ALLOCATE (w1xy(1:4,1:NCELL))
      ALLOCATE (w1yy(1:4,1:NCELL))

      ALLOCATE (right0(1:4,1:NCELL))
      ALLOCATE (rightx(1:4,1:NCELL))
      ALLOCATE (righty(1:4,1:NCELL))
      ALLOCATE (rightxx(1:4,1:NCELL))
      ALLOCATE (rightxy(1:4,1:NCELL))
      ALLOCATE (rightyy(1:4,1:NCELL))

	  allocate (o0(1:4,1:NCELL))
	  allocate (ox(1:4,1:NCELL))
	  allocate (oy(1:4,1:NCELL))
	  allocate (oxx(1:4,1:NCELL))
	  allocate (oxy(1:4,1:NCELL))
	  allocate (oyy(1:4,1:NCELL))

	  allocate(gp(1:3))
	  allocate(w(1:3))
	  allocate(jacobvertex(1:4,1:NCELL))
	  allocate(djacobvertex(1:4,1:4,1:NCELL))
	  allocate(djacobf(1:4,1:3,1:4,1:NCELL))
	  allocate(jacobf(1:3,1:4,1:NCELL))
	  allocate(djacobv(1:4,1:3,1:3,1:NCELL))
	  allocate(jacobv(1:3,1:3,1:NCELL))
	  allocate(f_edge(1:4,1:3,1:4,1:NCELL))
	  allocate(area(1:NCELL))
	  allocate(rmass(6))
	  
!	  allocate(jacobtecint(1:2,1:4,1:NCELL))
!	  allocate(jacobtecv(1:2,1:2,1:NCELL))
	  
	  allocate(Qv(1:4,1:3,1:3,1:NCELL))
	  allocate(Qub(1:4,1:3,1:2,1:NCELL))
	  allocate(Qlr(1:4,1:3,1:2,1:NCELL))
	  allocate(Qa(1:4,1:3,1:4,1:NCELL))
	  
	  allocate(iface2fp(N,4))
	  allocate(jface2fp(N,4))
	  
	  
	  w(1)=5./18.
	  w(2)=8./18.
	  w(3)=5./18.
	  gp(1)=-sqrt(15.)/10.
      gp(2)=0.0
      gp(3)=sqrt(15.)/10.
     

!	  MHD
!      DO I=1,NCELL
!	  DO K=1,8
!      u0(K,I)=0
!      u0(K,I)=0
!      u0(K,I)=0
!      u0(K,I)=0
!      ux(K,I)=0
!      uy(K,I)=0
!      uxx(K,I)=0
!      uxy(K,I)=0
!      uyy(K,I)=0
!      END DO
!      END DO
!	  
!	 ! CALL output
!	  !call tecplotnorefinement
!	  !call tecplotter(N)
!	  rmass(1)=1
!	  rmass(2)=12
!	  rmass(3)=12
!	  rmass(4)=180
!	  rmass(5)=144
!	  rmass(6)=180
!	  
!	   CALL calcjacob 
!	   
!	   do ic=1,NCELL
!	      do iv=1,4
!		     xxs(1,iv)=XV(IVCELL(ic,iv))
!			 xxs(2,iv)=YV(IVCELL(ic,iv))
!		  end do
!		  
!		  do igp=1,3
!		     do jgp=1,3
!			    CALL xyCoor_atGps(4,xxs,gp(igp),gp(jgp),xx)
!				x = xx(1)
!				
!            y = xx(2)
!            rho = gama**2
!            u = 0.01!-sin(y)
!            v = 0!sin(x)
!            ww = 0.d0
!            pr = gama
!            Bx = 0!-sin(y)
!            By = 0!sin(x*2)
!            Bz = 0.d0
!            NORM_U = u**2+v**2+ww**2
!            NORM_B = Bx**2+By**2+Bz**2
!            Qvi(1) = rho
!            Qvi(2) = rho*u
!            Qvi(3) = rho*v
!            Qvi(4) = rho*ww
!            Qvi(5) = pr/(gama-1) + 0.5d0*rho*NORM_U + 0.5d0*NORM_B
!            Qvi(6) = Bx
!            Qvi(7) = By
!            Qvi(8) = Bz
!			Qv(1:8,igp,jgp,ic)=Qvi(1:8)
!		    end do
!		end do
!		do igp=1,3
!		    do jgp=1,3
!			   do k=1,8
!			   u0(k,ic)=u0(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*rmass(1)
!			   ux(k,ic)=ux(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*gp(igp)*rmass(2)
!			   uy(k,ic)=uy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*gp(jgp)*rmass(3)
!			   uxx(k,ic)=uxx(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(igp)*gp(igp)-1./12.)*rmass(4)
!			   uxy(k,ic)=uxy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(igp)*gp(jgp))*rmass(5)
!			   uyy(k,ic)=uyy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(jgp)*gp(jgp)-1./12.)*rmass(6)
!			   end do
!			 end do
!		end do
		
!		Euler
      DO I=1,NCELL
		DO K=1,4
		u0(K,I)=0
		ux(K,I)=0
		uy(K,I)=0
		uxx(K,I)=0
		uxy(K,I)=0
		uyy(K,I)=0
		END DO
		END DO
		
	   ! CALL output
		!call tecplotnorefinement
		!call tecplotter(N)
		rmass(1)=1
		rmass(2)=12
		rmass(3)=12
		rmass(4)=180
		rmass(5)=144
		rmass(6)=180
		
		 CALL calcjacob 
		 
		 do ic=1,NCELL
			do iv=1,4
			   xxs(1,iv)=XV(IVCELL(ic,iv))
			   xxs(2,iv)=YV(IVCELL(ic,iv))
			end do
			
			do igp=1,3
			   do jgp=1,3
				  CALL xyCoor_atGps(4,xxs,gp(igp),gp(jgp),xx)
				  x = xx(1)
				  
			  y = xx(2)
			  rho = gama**2
			  u = 0.01!-sin(y)
			  v = 0!sin(x)
			  pr = gama

			  NORM_U = u**2+v**2

			  Qvi(1) = rho
			  Qvi(2) = rho*u
			  Qvi(3) = rho*v
			  Qvi(4) = pr/(gama-1) + 0.5d0*rho*NORM_U
			  Qv(1:4,igp,jgp,ic)=Qvi(1:4)
			  end do
		  end do
		  do igp=1,3
			  do jgp=1,3
				 do k=1,4
				 u0(k,ic)=u0(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*rmass(1)
				 ux(k,ic)=ux(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*gp(igp)*rmass(2)
				 uy(k,ic)=uy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*gp(jgp)*rmass(3)
				 uxx(k,ic)=uxx(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(igp)*gp(igp)-1./12.)*rmass(4)
				 uxy(k,ic)=uxy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(igp)*gp(jgp))*rmass(5)
				 uyy(k,ic)=uyy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(jgp)*gp(jgp)-1./12.)*rmass(6)
				 end do
			   end do
		  end do

		end do
		!call tecplotnorefinement
		call tecplotter(N)
		! ! for face 4
	! do jgp=1,3
	!   CALL xyCoor_atGps(4,xxs,-0.5,gp(N-jgp+1),xx)
	!          x=xx(1)
	!			y=xx(2)
	!			rho=gama**2
	!			u=-sin(y)
	!			v=sin(x)
	!			w=0.d0
	!			pr=gama
	!			Bx=-sin(y)
	!			By=sin(x*2)
	!			Bz=0.d0
	!			NORM_U=u**2+v**2+w**2
	!			NORM_B=Bx**2+By**2+Bz**2
	!			Qvi(1)=rho
	!			Qvi(2)=rho*u
	!			Qvi(3)=rho*v
	!			Qvi(4)=rho*w
	!			Qvi(5)=pr/(gama-1)+0.5*rho*NORM_U+0.5*NORM_B
	!			Qvi(6)=Bx
	!			Qvi(7)=By
	!			Qvi(8)=Bz
	!			Qa(1:8,jgp,1,ic)=Qvi(1:8)
	!end do
	
	!!for face 2
	!do jgp=1,3
	!   CALL xyCoor_atGps(4,xxs,0.5,gp(jgp),xx)
	!   x=xx(1)
	!			y=xx(2)
	!			rho=gama**2
	!			u=-sin(y)
	!			v=sin(x)
	!			w=0.d0
	!			pr=gama
	!			Bx=-sin(y)
	!			By=sin(x*2)
	!			Bz=0.d0
	!			NORM_U=u**2+v**2+w**2
	!			NORM_B=Bx**2+By**2+Bz**2
	!			Qvi(1)=rho
	!			Qvi(2)=rho*u
	!			Qvi(3)=rho*v
	!			Qvi(4)=rho*w
	!			Qvi(5)=pr/(gama-1)+0.5*rho*NORM_U+0.5*NORM_B
	!			Qvi(6)=Bx
	!			Qvi(7)=By
	!			Qvi(8)=Bz
	!			Qa(1:8,jgp,2,ic)=Qvi(1:8)
	!end do
	
	!!for face 3
	!do igp=1,3
	!   CALL xyCoor_atGps(4,xxs,gp(N-igp+1),0.5,xx)
	!   x=xx(1)
	!			y=xx(2)
	!			rho=gama**2
	!			u=-sin(y)
	!			v=sin(x)
	!			w=0.d0
	!			pr=gama
	!			Bx=-sin(y)
	!			By=sin(x*2)
	!			Bz=0.d0
	!			NORM_U=u**2+v**2+w**2
	!			NORM_B=Bx**2+By**2+Bz**2
	!			Qvi(1)=rho
	!			Qvi(2)=rho*u
	!			Qvi(3)=rho*v
	!			Qvi(4)=rho*w
	!			Qvi(5)=pr/(gama-1)+0.5*rho*NORM_U+0.5*NORM_B
	!			Qvi(6)=Bx
	!			Qvi(7)=By
	!			Qvi(8)=Bz
	!			Qa(1:8,igp,3,ic)=Qvi(1:8)
	!end do
	
	!! !for face 1
	!  do igp=1,3
	!   CALL xyCoor_atGps(4,xxs,gp(igp),-0.5,xx)
	!   x=xx(1)
	!			y=xx(2)
	!			rho=gama**2
	!			u=-sin(y)
	!			v=sin(x)
	!			w=0.d0
	!			pr=gama
	!			Bx=-sin(y)
	!			By=sin(x*2)
	!			Bz=0.d0
	!			NORM_U=u**2+v**2+w**2
	!			NORM_B=Bx**2+By**2+Bz**2
	!			Qvi(1)=rho
	!			Qvi(2)=rho*u
	!			Qvi(3)=rho*v
	!			Qvi(4)=eho*w
	!			Qvi(5)=pr/(gama-1)+0.5*rho*NORM_U+0.5*NORM_B
	!			Qvi(6)=Bx
	!			Qvi(7)=By
	!			Qvi(8)=Bz
	!			Qa(1:8,igp,4,ic)=Qvi(1:8)
	!end do
	!end do
	
		
		
		
	  ! the following code is for euler vortex case
	!  do ic=1, NCELL
	!     do iv=1,4
	!	   xxs(1,iv)=XV(IVCELL(ic,iv))
	!	   xxs(2,iv)=YV(IVCELL(ic,iv))
	!	 end do
	!	 
	!	do igp =1,3
	!	   do jgp=1,3
	!	      CALL xyCoor_atGps(4,xxs,gp(igp),gp(jgp),xx)
	!		  xr=xx(1)-5
	!		  yr=xx(2)-5
	!		  ur = COS(theta)-yr*EXP((1-xr**2-yr**2)/2)
	!		   vr = SIN(theta)+yr*EXP((1-xr**2-yr**2)/2)
	!		   Qvi(1) = rinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(1/(gama-1))
	!		   !write(*,*) Qvi(1)
	!		   pr = pinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(gama/(gama-1))
	!		   V2 = ur**2+vr**2
	!		   Qvi(2)=Qvi(1)*ur
	!		   Qvi(3)=Qvi(1)*vr
	!		   Qvi(4)=pr/(gama-1)+0.5*Qvi(1)*V2
	!		   Qv(1:4,igp,jgp,ic)=Qvi(1:4)
	!	 end do
	! end do
	! !write(*,*) Qv(1,1,1,ic)
	! 
	!     do igp=1,3
	!	    do jgp=1,3
	!		   do k=1,4
	!		   u0(k,ic)=u0(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*rmass(1)
	!		   ux(k,ic)=ux(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*gp(igp)*rmass(2)
	!		   uy(k,ic)=uy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*gp(jgp)*rmass(3)
	!		   uxx(k,ic)=uxx(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(igp)*gp(igp)-1./12.)*rmass(4)
	!		   uxy(k,ic)=uxy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(igp)*gp(jgp))*rmass(5)
	!		   uyy(k,ic)=uyy(k,ic) + w(igp)*w(jgp)*Qv(k,igp,jgp,ic)*(gp(jgp)*gp(jgp)-1./12.)*rmass(6)
	!		   end do
	!		 end do
	!	end do
	!	
	! ! for face 4
	! do jgp=1,3
	!   CALL xyCoor_atGps(4,xxs,-0.5,gp(N-jgp+1),xx)
	!   xr=xx(1)-5
	!		  yr=xx(2)-5
	!		  ur = COS(theta)-yr*EXP((1-xr**2-yr**2)/2)
	!		   vr = SIN(theta)+yr*EXP((1-xr**2-yr**2)/2)
	!		   Qvi(1) = rinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(1/(gama-1))
	!		   pr = pinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(gama/(gama-1))
	!		   V2 = ur**2+vr**2
	!		   Qvi(2)=Qvi(1)*ur
	!		   Qvi(3)=Qvi(1)*vr
	!		   Qvi(4)=pr/(gama-1)+0.5*Qvi(1)*V2
	!		   Qa(1:4,jgp,4,ic)=Qvi(1:4)
	!end do
	
	!!for face 2
	!do jgp=1,3
	!   CALL xyCoor_atGps(4,xxs,0.5,gp(jgp),xx)
	!   xr=xx(1)-5
	!		  yr=xx(2)-5
	!		  ur = COS(theta)-yr*EXP((1-xr**2-yr**2)/2)
	!		   vr = SIN(theta)+yr*EXP((1-xr**2-yr**2)/2)
	!		   Qvi(1) = rinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(1/(gama-1))
	!		   pr = pinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(gama/(gama-1))
	!		   V2 = ur**2+vr**2
	!		   Qvi(2)=Qvi(1)*ur
	!		   Qvi(3)=Qvi(1)*vr
	!		   Qvi(4)=pr/(gama-1)+0.5*Qvi(1)*V2
	!		   Qa(1:4,jgp,2,ic)=Qvi(1:4)
	!end do
	
	!!for face 3
	!do igp=1,3
	!   CALL xyCoor_atGps(4,xxs,gp(N-igp+1),0.5,xx)
	!   xr=xx(1)-5
	!		  yr=xx(2)-5
	!		  ur = COS(theta)-yr*EXP((1-xr**2-yr**2)/2)
	!		   vr = SIN(theta)+yr*EXP((1-xr**2-yr**2)/2)
	!		   Qvi(1) = rinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(1/(gama-1))
	!		   pr = pinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(gama/(gama-1))
	!		   V2 = ur**2+vr**2
	!		   Qvi(2)=Qvi(1)*ur
	!		   Qvi(3)=Qvi(1)*vr
	!		   Qvi(4)=pr/(gama-1)+0.5*Qvi(1)*V2
	!		   Qa(1:4,igp,3,ic)=Qvi(1:4)
	!end do
	! !for face 1
	!  do igp=1,3
	!   CALL xyCoor_atGps(4,xxs,gp(igp),-0.5,xx)
	!   xr=xx(1)-5
	!		  yr=xx(2)-5
	!		  ur = COS(theta)-yr*EXP((1-xr**2-yr**2)/2)
	!		   vr = SIN(theta)+yr*EXP((1-xr**2-yr**2)/2)
	!		   Qvi(1) = rinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(1/(gama-1))
	!		   pr = pinf*(1-(gama-1)*0.3**2/2*EXP(1-xr**2-yr**2))**(gama/(gama-1))
	!		   V2 = ur**2+vr**2
	!		   Qvi(2)=Qvi(1)*ur
	!		   Qvi(3)=Qvi(1)*vr
	!		   Qvi(4)=pr/(gama-1)+0.5*Qvi(1)*V2
	!		   Qa(1:4,igp,1,ic)=Qvi(1:4)
	!end do
	
	!end do

    !call tecplotter(N)
	!call tecplotnorefinement
	 
	
      END SUBROUTINE initialize
	  
	  
 FUNCTION calc_Q(q0,q1,q2,q3,q4,q5,psi,eta)

	      IMPLICIT NONE
	      DOUBLE PRECISION,intent(in) :: q0,q1,q2,q3,q4,q5,psi,eta
		  DOUBLE PRECISION :: calc_Q

	  
	      calc_Q=q0*1+q1*psi+q2*eta+q3*(psi*psi-1./12.)+q4*psi*eta+q5*(eta*eta-1./12.)
		  
		  			 
	 END FUNCTION calc_Q

	 SUBROUTINE MapBaseFunc(func,xi,eta)
	 IMPLICIT NONE
	 DOUBLE PRECISION :: func(4),xi,eta
	       
	      func(1)=(0.5-xi)*(0.5-eta)
		  func(2)=(0.5+xi)*(0.5-eta)
		  func(3)=(0.5+xi)*(0.5+eta)
		  func(4)=(0.5-xi)*(0.5+eta)
		  RETURN
     END
	 
	 SUBROUTINE face2fpconnec
	 
	 use setup
	 IMPLICIT NONE
	 INCLUDE 'MESH2D.INC'
	 
	 integer ::nfp,iface
	 
	 iface=1
	 do nfp=1,N
	    iface2fp(nfp,iface)=nfp
		jface2fp(nfp,iface)=1
     end do
	 
	 iface=2
	 do nfp=1,N
	    iface2fp(nfp,iface)=N+1
		jface2fp(nfp,iface)=nfp
     end do
	 
	 iface=3
	 do nfp=1,N
	    iface2fp(nfp,iface)=N-nfp+1
		jface2fp(nfp,iface)=N+1
     end do
	 
	 iface=4
	 do nfp=1,N
	    iface2fp(nfp,iface)=1
		jface2fp(nfp,iface)=N-nfp+1
     end do
	 
	 END SUBROUTINE face2fpconnec
	    


      

