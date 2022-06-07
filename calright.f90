SUBROUTINE calright
      use setup
      IMPLICIT NONE
      INCLUDE 'MESH2D.INC'

      DOUBLE PRECISION,DIMENSION(8,3,3) :: FC,GC
	  INTEGER :: ic,k,m
      DOUBLE PRECISION,DIMENSION(8) :: f1_int,f1x_int,f1y_int
      DOUBLE PRECISION,DIMENSION(8) :: g1_int,g1x_int,g1y_int
      DOUBLE PRECISION,DIMENSION(0:5,8) :: rint
      DOUBLE PRECISION,DIMENSION(0:5,8) :: fluxr,fluxl,fluxu,fluxb,flux
	  DOUBLE PRECISION :: w_1,w0,w1
	  DOUBLE PRECISION :: a1,ax,ay,axx,axy,ayy
	  !DOUBLE PRECISION,DIMENSION(3) :: w
	  DOUBLE PRECISION :: psi,eta
      INTEGER :: igp,jgp


	  
	  !w(1)=5.0/18.
	  !w(2)=8./18.
	  !w(3)=5.0/18.

      CALL FLUXedge 

	  CALL BCflux
	  do ic=1,NCELL
	      CALL calc_FG(ic,FC,GC)
		  
		  do k=1,8
		      f1_int(k)=Rint_area(FC(k,1,1),FC(k,1,3),FC(k,3,1), FC(k,3,3), &
								  FC(k,2,1),FC(k,2,3),FC(k,1,2),FC(k,3,2),FC(k,2,2))
			  f1x_int(k)=Rint_area(FC(k,1,1)*gp(1),FC(k,1,3)*gp(1),FC(k,3,1)*gp(3),FC(k,3,3)*gp(3), &
								   FC(k,2,1)*gp(2),FC(k,2,3)*gp(2),FC(k,1,2)*gp(1),FC(k,3,2)*gp(3),FC(k,2,2)*gp(2))
			  f1y_int(k)=Rint_area(FC(k,1,1)*gp(1),FC(k,1,3)*gp(3),FC(k,3,1)*gp(1),FC(k,3,3)*gp(3), &
								  FC(k,2,1)*gp(1),FC(k,2,3)*gp(3),FC(k,1,2)*gp(2),FC(k,3,2)*gp(2),FC(k,2,2)*gp(2))
			  
			  
			  g1_int(k)=Rint_area(GC(k,1,1),GC(k,1,3),GC(k,3,1), GC(k,3,3), &
								  GC(k,2,1),GC(k,2,3),GC(k,1,2),GC(k,3,2),GC(k,2,2))
			  g1x_int(k)=Rint_area(GC(k,1,1)*gp(1),GC(k,1,3)*gp(1),GC(k,3,1)*gp(3),GC(k,3,3)*gp(3), &
								   GC(k,2,1)*gp(2),GC(k,2,3)*gp(2),GC(k,1,2)*gp(1),GC(k,3,2)*gp(3),GC(k,2,2)*gp(2))
			  g1y_int(k)=Rint_area(GC(k,1,1)*gp(1),GC(k,1,3)*gp(3),GC(k,3,1)*gp(1),GC(k,3,3)*gp(3), &
								  GC(k,2,1)*gp(1),GC(k,2,3)*gp(3),GC(k,1,2)*gp(2),GC(k,3,2)*gp(2),GC(k,2,2)*gp(2))
			  
		  end do
		  
		  do k=1,8
		     rint(0,k)=0.0
		     rint(1,k)=f1_int(k)
			 rint(2,k)=g1_int(k)
			 rint(3,k)=2*f1x_int(k)
			 rint(4,k)=f1y_int(k)+g1x_int(k)
			 rint(5,k)=2*g1y_int(k)
		  end do
		  
		  w_1=gp(1)*gp(1)-1./12.
		  w0=gp(2)*gp(2)-1./12.
		  w1=gp(3)*gp(3)-1./12.
		  
		  do k=1,8
		    fluxr(0,k)=Rint_edge(f_edge(k,1,2,ic),f_edge(k,2,2,ic),f_edge(k,3,2,ic))
			fluxr(1,k)=0.5*Rint_edge(f_edge(k,1,2,ic),f_edge(k,2,2,ic),f_edge(k,3,2,ic))
			fluxr(2,k)=Rint_edge(f_edge(k,1,2,ic)*gp(1),f_edge(k,2,2,ic)*gp(2),f_edge(k,3,2,ic)*gp(3))
			fluxr(3,k)=1./6.*Rint_edge(f_edge(k,1,2,ic),f_edge(k,2,2,ic),f_edge(k,3,2,ic))
		    fluxr(4,k)=0.5*Rint_edge(f_edge(k,1,2,ic)*gp(1),f_edge(k,2,2,ic)*gp(2),f_edge(k,3,2,ic)*gp(3))
			
			fluxr(5,k)=Rint_edge(f_edge(k,1,2,ic)*w_1,f_edge(k,2,2,ic)*w0,f_edge(k,3,2,ic)*w1)
		 end do
		 
		 do k=1,8
		    fluxl(0,k)=Rint_edge(f_edge(k,1,4,ic),f_edge(k,2,4,ic),f_edge(k,3,4,ic))
			fluxl(1,k)=-0.5*Rint_edge(f_edge(k,1,4,ic),f_edge(k,2,4,ic),f_edge(k,3,4,ic))
			fluxl(2,k)=Rint_edge(f_edge(k,1,4,ic)*gp(3),f_edge(k,2,4,ic)*gp(2),f_edge(k,3,4,ic)*gp(1))
			fluxl(3,k)=1./6.*Rint_edge(f_edge(k,1,4,ic),f_edge(k,2,4,ic),f_edge(k,3,4,ic))
		    fluxl(4,k)=-0.5*Rint_edge(f_edge(k,1,4,ic)*gp(3),f_edge(k,2,4,ic)*gp(2),f_edge(k,3,4,ic)*gp(1))
			
			fluxl(5,k)=Rint_edge(f_edge(k,1,2,ic)*w_1,f_edge(k,2,2,ic)*w0,f_edge(k,3,2,ic)*w1)
		 end do
		 
		 do k=1,8
		    fluxu(0,k)=Rint_edge(f_edge(k,1,3,ic),f_edge(k,2,3,ic),f_edge(k,3,3,ic))
			fluxu(1,k)=Rint_edge(f_edge(k,1,3,ic)*gp(3),f_edge(k,2,3,ic)*gp(2),f_edge(k,3,3,ic)*gp(1))
			fluxu(2,k)=0.5*Rint_edge(f_edge(k,1,3,ic),f_edge(k,2,3,ic),f_edge(k,3,3,ic))
			fluxu(3,k)=Rint_edge(f_edge(k,1,3,ic)*w_1,f_edge(k,2,3,ic)*w0,f_edge(k,3,3,ic)*w1)
		    fluxu(4,k)=0.5*Rint_edge(f_edge(k,1,3,ic)*gp(3),f_edge(k,2,3,ic)*gp(2),f_edge(k,3,3,ic)*gp(1))
			
			fluxu(5,k)=1./6.*Rint_edge(f_edge(k,1,3,ic),f_edge(k,2,3,ic),f_edge(k,3,3,ic))
		 end do
		 
		 do k=1,8
		    fluxb(0,k)=Rint_edge(f_edge(k,1,1,ic),f_edge(k,2,1,ic),f_edge(k,3,1,ic))
			fluxb(1,k)=Rint_edge(f_edge(k,1,1,ic)*gp(1),f_edge(k,2,1,ic)*gp(2),f_edge(k,3,1,ic)*gp(3))
			fluxb(2,k)=-0.5*Rint_edge(f_edge(k,1,1,ic),f_edge(k,2,1,ic),f_edge(k,3,1,ic))
			fluxb(3,k)=Rint_edge(f_edge(k,1,1,ic)*w_1,f_edge(k,2,1,ic)*w0,f_edge(k,3,1,ic)*w1)
		    fluxb(4,k)=-0.5*Rint_edge(f_edge(k,1,1,ic)*gp(1),f_edge(k,2,1,ic)*gp(2),f_edge(k,3,1,ic)*gp(3))
			
			fluxb(5,k)=1./6.*Rint_edge(f_edge(k,1,1,ic),f_edge(k,2,1,ic),f_edge(k,3,1,ic))
		 end do
		 
		 do m=0,5
		    do k=1,8
			  flux(m,k)=-fluxl(m,k)+fluxr(m,k)+fluxu(m,k)-fluxb(m,k)
			end do
		end do
		
		do igp=1,3
		  do jgp=1,3
		     psi=gp(igp)
			 eta=gp(jgp)
			 a1 = a1+jacobv(igp,jgp,ic)*w(igp)*w(jgp)*rmass(1)
			 ax =ax+jacobv(igp,jgp,ic)*w(igp)*w(jgp)*psi*psi*rmass(2)
			 ay=ay+jacobv(igp,jgp,ic)*w(igp)*w(jgp)*eta*eta*rmass(3)
			 axx =axx+jacobv(igp,jgp,ic)*w(igp)*w(jgp)*(psi*psi-1.0/12)*(psi*psi-1.0/12)*rmass(4)
			 axy=axy+jacobv(igp,jgp,ic)*w(igp)*w(jgp)*psi*eta*psi*eta*rmass(5)
			 ayy =ayy+jacobv(igp,jgp,ic)*w(igp)*w(jgp)*(eta*eta-1.0/12)*(eta*eta-1.0/12)*rmass(6)
		 end do
		end do
		
		do k=1,8
		  right0(k,ic)=(rint(0,k)-flux(0,k))/a1
		  rightx(k,ic)=(rint(1,k)-flux(1,k))/ax
		  righty(k,ic)=(rint(2,k)-flux(2,k))/ay
		  rightxx(k,ic)=(rint(3,k)-flux(3,k))/axx
		  rightxy(k,ic)=(rint(4,k)-flux(4,k))/axy
		  rightyy(k,ic)=(rint(5,k)-flux(5,k))/ayy

!		  if(right0 .ne. right0) then
!			write(*,*)"right0 nan"
!			call abort
!		else if(rightx .ne. rightx) then
!			write(*,*)"rightx nan"
!			call abort
!		else if(righty .ne. righty) then
!			write(*,*)"righty nan"
!			call abort
!		else if(rightxx .ne. rightxx) then
!			write(*,*)"rightxx nan"
!			call abort
!		else if(rightxy .ne. rightxy) then
!			write(*,*)"rightxy nan"
!			call abort
!		else if(rightyy .ne. rightyy) then
!			write(*,*)"rightyy nan"
!			call abort
!		end if

		end do
		end do

		contains
	  DOUBLE PRECISION FUNCTION Rint_area(f_1_1,f_11,f1_1,f11, f0_1,f01,f_10,f10, f00)

	      IMPLICIT NONE
	      DOUBLE PRECISION,intent(in) :: f_1_1,f_11,f1_1,f11,f0_1,f00,f01,f10,f_10

	  
	      Rint_area= 25./324.*(f_1_1+f_11+f1_1+f11)+40./324.*(f0_1+f01+f_10+f10)+64./324.*f00
		  END
		  
	       DOUBLE PRECISION FUNCTION Rint_edge(fv_1,fv0,fv1)

	        IMPLICIT NONE
	        DOUBLE PRECISION,intent(in) :: fv_1,fv0,fv1
		  
	         Rint_edge=5.0/18.*fv_1+8./18.*fv0+5.0/18.*fv1
		  END

		END SUBROUTINE calright
		
		
		
		 
		
		
				  
		  
	  
	  



