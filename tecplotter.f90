SUBROUTINE output
   use setup
   IMPLICIT NONE 
   INCLUDE 'MESH2D.INC'
   integer :: ic,i,j,k,vert,count,iv
   double precision :: rho,u,v,rhoE,p
   double precision :: psi,eta,x,y,calc_Q
   double precision :: temprho,tempu,tempv,temprhoE
   double precision :: mach,c,inte
   
   character(80) :: filename=' '
   character*16  ::tecfile
   
   write(filename,'(a,i6.6,a)')'testec',iter,'.dat'
   tecfile=trim(filename)
   open(unit=1,file=tecfile,form='formatted')
   
   write(1,*) 'TITLE = "DG 2D QUADMESH"'
   write(1,*) 'VARIABLES = "X","Y","rho","u","v","p","Mach","Temp"'
   write(1,*) 'ZONE T="REFINED MESHES"  N=',NVERT, "E=",NCELL,"DATAPACKING=POINT"
   write(1, *) 'ZONETYPE=FEQuadrilateral'
   
   do iv=1,NVERT
      x=XV(iv)
	  y=YV(iv)
	  
	  do k=ICVSTA(iv),ICVSTA(iv+1)-1
	      ic=ICVERT(k)
		  if (iv .eq. IVCELL(ic,1)) then
		    psi=-0.5
		    eta=-0.5
		    vert=1
		  else if (iv .eq. IVCELL(ic,2))   then
		     psi=0.5
			 eta=-0.5
			 vert=2
		   else if  (iv .eq. IVCELL(ic,3))  then
		     psi=0.5
			 eta=0.5
			 vert=3
		   else
		     psi=-0.5
			 eta=0.5
			 vert=4
		   end if
		   
		   temprho =calc_Q(u0(1,ic),ux(1,ic),uy(1,ic),uxx(1,ic), &
						   uxy(1,ic),uyy(1,ic),psi,eta)
		   rho = rho+temprho
		   
		   tempu =calc_Q(u0(2,ic),ux(2,ic),uy(2,ic),uxx(2,ic), &
						   uxy(2,ic),uyy(2,ic),psi,eta)/temprho
		   u = u+tempu
		   
		   tempv =calc_Q(u0(3,ic),ux(3,ic),uy(3,ic),uxx(3,ic), &
						   uxy(3,ic),uyy(3,ic),psi,eta)/temprho
		   v = v+tempv
		   
		   temprhoE = calc_Q(u0(4,ic),ux(4,ic),uy(4,ic),uxx(4,ic), &
						   uxy(4,ic),uyy(4,ic),psi,eta)
		   rhoE =rhoE+temprhoE
		   
		   count =count+1
		   
		end do
		  
		  rho=rho/count
		  u=u/count
		  v=v/count
		  rhoE=rhoE/count
		  p=(gama-1)*(rhoE-0.5*rho*(u*u+v*v))
		  c=sqrt(gama*p/rho)
		  mach=sqrt(u**2+v**2)/c
		  inte=p/(gama-1)/rho
		  
		  write(1,1000) XV(iv),YV(iv), rho,u,v,p,mach,inte
	end do
		  do ic=1,NCELL
		      write(1,1001) (IVCELL(ic,i),i=1,4)
		  end do
1000       format(8(G16.8,1X))
1001      format (4(I9,1X))   
           close(1)
		   
		   
	!	   contains
	 ! DOUBLE PRECISION FUNCTION calc_Q(q0,q1,q2,q3,q4,q5,psi,eta)

	 !     IMPLICIT NONE
	 !     DOUBLE PRECISION,intent(in) :: q0,q1,q2,q3,q4,q5,psi,eta

	 ! 
	 !     calc_Q=q0*1+q1*psi+q2*eta+q3*(psi*psi-1./12.)+q4*psi*eta+q5*(eta*eta-1./12.)
	!	  
	!	  			 
	 !END 
		   
		   end 
		   
		   
			 
		  
		  
	  


