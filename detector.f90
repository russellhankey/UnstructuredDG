SUBROUTINE Detector

    use setup
	IMPLICIT NONE
	INCLUDE 'MESH2D.INC'
	
	double precision :: mm
	integer :: ic,igp,jgp,flag01,flag02,iface
	integer :: nfp
	double precision :: psi,eta,Qbar,h,deltarho,rho
	mm=0.00002
	
	do 111 ic=1,NCELL
	    Qbar=0
		do igp=1,3
		   do jgp=1,3
		      psi=gp(igp)
			  eta=gp(jgp)
			  Qbar= Qbar + calc_Q(u0(0,ic),ux(0,ic),uy(0,ic),uxx(0,ic),uxy(0,ic),uyy(0,ic),psi,eta)
		   end do
		 end do
		 Qbar=Qbar/9
		 flag01=1
		 flag02=0
		 h=sqrt(area(ic))
		 
		 do 112 iface=1,4
		    if (flag02 .eq.1) GOTO 111
			do nfp=1,3
			   if(flag01 .eq.1) GOTO 112
			      if (iface.eq.1) then
				     psi=gp(nfp)
					 eta=-0.5
					 rho=calc_Q(u0(0,ic),ux(0,ic),uy(0,ic),uxx(0,ic),uxy(0,ic),uyy(0,ic),psi,eta)
					 deltarho=abs(rho-Qbar)
					 if(deltarho .gt. mm*h*h) then
					     flag01=1
						 CALL limiter(ic)
					end if
				end if
				
			      if (iface.eq.2) then
				     psi=0.5
					 eta=gp(nfp)
					 rho=calc_Q(u0(0,ic),ux(0,ic),uy(0,ic),uxx(0,ic),uxy(0,ic),uyy(0,ic),psi,eta)
					 deltarho=abs(rho-Qbar)
					 if(deltarho .gt. mm*h*h) then
					     flag01=1
						 CALL limiter(ic)
					end if
				end if	 
				
				if (iface.eq.3) then
				     psi=gp(nfp)
					 eta=0.5
					 rho=calc_Q(u0(0,ic),ux(0,ic),uy(0,ic),uxx(0,ic),uxy(0,ic),uyy(0,ic),psi,eta)
					 deltarho=abs(rho-Qbar)
					 if(deltarho .gt. mm*h*h) then
					     flag01=1
						 CALL limiter(ic)
					end if
				end if	
				
				if (iface.eq.4) then
				     psi=-0.5
					 eta=gp(nfp)
					 rho=calc_Q(u0(0,ic),ux(0,ic),uy(0,ic),uxx(0,ic),uxy(0,ic),uyy(0,ic),psi,eta)
					 deltarho=abs(rho-Qbar)
					 if(deltarho .gt. mm*h*h) then
					     flag01=1
						 CALL limiter(ic)
					end if
				end if
				
				if (flag01 .eq.1) flag02=1
		 end do
112 continue	
		
	
111 continue

		 	 contains
	  DOUBLE PRECISION FUNCTION calc_Q(q0,q1,q2,q3,q4,q5,psi2,eta2)

	      IMPLICIT NONE
	      DOUBLE PRECISION,intent(in) :: q0,q1,q2,q3,q4,q5,psi2,eta2

	  
	      calc_Q=q0*1+q1*psi2+q2*eta2+q3*(psi2*psi2-1./12.)+q4*psi2*eta2+q5*(eta2*eta2-1./12.)
		  
		  			 
	 END 
		 
END SUBROUTINE Detector

    SUBROUTINE limiter(ic)
	use setup
	implicit none
	   include 'MESH2D.INC'

        integer :: ic,flag,i,icc,k,iface,igp
		double precision,dimension(4) :: MMK,m,Qbar,Qb,Qq
		double precision,dimension(3,4,4) :: a
		integer,dimension(4) :: flag2,mina
		double precision :: psi,eta,Q,u,v,p,calc_Q
		
		
		do k=1,4
		  flag2(k)=0
		  mina(k)=100000
		end do
		
		
		flag=0
		do i=1,4
		    icc=c2c(i,ic)
			if(icc .ne.0) then
			     Qb(1)=u0(1,icc)
				 Qb(2)=u0(2,icc)
				 Qb(3)=u0(3,icc)
				 Qb(4)=u0(4,icc)
				 
				 if (flag .eq.0) then
				      do k=1,4
					     MMK(k)=Qb(k)
						 m(k)=Qb(k)
						 flag = flag+1
					  end do
				 end if
				 
				 do k=1,4
				    if(MMK(k) .lt. Qb(k)) MMK(k)= Qb(k)
					if(m(k) .gt. Qb(k)) m(k)= Qb(k)
				end do
				
			end if
		end do
		
		Qbar(1)=u0(1,ic)
		Qbar(2)=u0(2,ic)
		Qbar(3)=u0(3,ic)
		Qbar(4)=u0(4,ic)
		
		do k=1,4
		   do iface =1,4
		      do igp=1,3
			     a(igp,iface,k)=1
			   end do
			end do
	    end do
		
		do k=1,4
		   do iface=1,4
		      do igp=1,3
			     if (iface .eq. 1) then
				     psi=gp(igp)
					 eta=-0.5
				 else if(iface .eq.2) then
				     psi=0.5
					 eta=gp(igp)
				 else if(iface .eq.3) then
				     eta=0.5
					 psi=gp(igp)
				else
				      psi=-0.5
					 eta=gp(igp)
				end if
				
				Q=calc_Q(u0(k,ic),ux(k,ic),uy(k,ic),uxx(k,ic),uxy(k,ic),uyy(k,ic),psi,eta)
				if(Q>MMK(k)) then
				  a(igp,iface,k)=max(0.0, (MMK(k)-Qbar(k))/(Q-Qbar(k)))
				  flag2(k)=flag2(k)+1
				end if
				
				if(Q<m(k)) then
				  a(igp,iface,k)=max(0.0, (m(k)-Qbar(k))/(Q-Qbar(k)))
				  flag2(k)=flag2(k)+1
				end if
				
				if(mina(k).gt.a(igp,iface,k)) mina(k)=a(igp,iface,k)
			end do
	end do
	end do
	
	do k=1,4
	   if(flag2(k) .ne.0) then
	      u0(k,ic)=u0(k,ic)
		  ux(k,ic)=mina(k)*ux(k,ic)
		  uy(k,ic)=mina(k)*uy(k,ic)
		  uxx(k,ic)=0
		  uxy(k,ic)=0
		  uyy(k,ic)=0
		end if
	end do
	
	do iface=1,4
	   do igp=1,3
	      if (iface .eq. 1) then
				     psi=gp(igp)
					 eta=-0.5
				 else if(iface .eq.2) then
				     psi=0.5
					 eta=gp(igp)
				 else if(iface .eq.3) then
				     eta=0.5
					 psi=gp(igp)
				else
				      psi=-0.5
					 eta=gp(igp)
				end if
				
			do k=1,4
			   Qq(k)=calc_Q(u0(k,ic),ux(k,ic),uy(k,ic),uxx(k,ic),uxy(k,ic),uyy(k,ic),psi,eta)
			end do
			
			u=Qq(2)/Qq(1)
			v=Qq(3)/Qq(1)
			p=(gama-1)*(Qq(4)-0.5*Qq(1)*(u*u+v*v))
			
			if(p .lt.0) then
			    do k=1,4
				    u0(k,ic)=u0(k,ic)
					ux(k,ic)=0
					uy(k,ic)=0
					uxx(k,ic)=0
					uxy(k,ic)=0
					uyy(k,ic)=0
				end do
			end if
	end do
	end do
	
	!contains
	!  DOUBLE PRECISION FUNCTION calc_Q(q0,q1,q2,q3,q4,q5,psi,eta)

	!      IMPLICIT NONE
	!      DOUBLE PRECISION,intent(in) :: q0,q1,q2,q3,q4,q5,psi,eta

	!  
	!      calc_Q=q0*1+q1*psi+q2*eta+q3*(psi*psi-1./12.)+q4*psi*eta+q5*(eta*eta-1./12.)
	!	  
	!	  			 
	! END
	
	END SUBROUTINE limiter
			
			
	      
				
		
		
		      
		
		
		
		
			
		   
		
		
		 

