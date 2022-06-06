

      SUBROUTINE FLUXedge
      use setup
      IMPLICIT NONE
      INCLUDE 'MESH2D.INC'

      INTEGER :: icleft,icright,sign_l,NF,faml,famr,k,sign_r
      INTEGER :: ifacelc,ifacerc,nfp,nfpr,i,ic,iface
	  INTEGER :: flagl,flagr

      DOUBLE PRECISION,DIMENSION(8) :: Fnl,Fnr
      DOUBLE PRECISION,DIMENSION(8) :: Ql,Qr
      DOUBLE PRECISION,DIMENSION(2) :: normfl
	  

      DOUBLE PRECISION :: psi,eta,calc_Q
	  
	  
	  
      DO NF=1,NINTER
	     
	     iface=IBFINTER(NF)

		 
         icleft=IF2C(iface,1)

		 
         icright=IF2C(iface,2)

         ifacelc=IF2C(iface,3)
         ifacerc=IF2C(iface,4)
		 
		 faml=mod(ifacelc,2)+1
		 famr=mod(ifacerc,2)+1

         DO nfp=1,3

         ic=icleft

         sign_l=1
         if (ifacelc.eq.1) then
                 psi=gp(nfp)
                 eta=-0.5
                 do k=1,8
                    Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft), &
                            uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
                  end do
                  sign_l=-1
                  normfl(1)=-djacobf(3,nfp,ifacelc,icleft)
                  normfl(2)=djacobf(1,nfp,ifacelc,icleft)
          end if

          if (ifacelc.eq.2) then
                 psi=0.5
                 eta=gp(nfp)
                 do k=1,8
                    Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),&
                            uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
                  end do
                  
                  normfl(1)=djacobf(4,nfp,ifacelc,icleft)
                  normfl(2)=-djacobf(2,nfp,ifacelc,icleft)
          end if

          if (ifacelc.eq.3) then
                 psi=gp(N-nfp+1)
                 eta=0.5
                 do k=1,4
                    Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),&
                            uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
                  end do
                  
                  normfl(1)=-djacobf(3,nfp,ifacelc,icleft)
                  normfl(2)=djacobf(1,nfp,ifacelc,icleft)
          end if

          if (ifacelc.eq.4) then
                 psi=-0.5
                 eta=gp(N-nfp+1)
                 do k=1,8
                    Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),&
                            uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
                  end do
                  sign_l=-1

                  normfl(1)=djacobf(4,nfp,ifacelc,icleft)
                  normfl(2)=-djacobf(2,nfp,ifacelc,icleft)
          end if
		  
		  !flagl=0
		  !flagr=0
		  !if(ifacelc .eq. 1. .or. ifacelc .eq. 2.) then
		  !         flagl=0
	      !else
		  !         flagl=1
	      !end if
		  
		  !if(ifacerc .eq. 1. .or. ifacerc .eq. 2.) then
		  !         flagl=1
	      !else
		  !         flagl=0
	      !end if
		  
		  !if (flagl .eq. flagr) then
		  !     nfpr=nfp
		  !else
		  !     nfpr=N+1-nfp
		  !end if
		  nfpr=N+1-nfp
		  sign_r=1
		  if (ifacerc.eq.1) then
                 psi=gp(nfpr)
                 eta=-0.5
                 do k=1,8
                    Qr(k)=calc_Q(u0(k,icright),ux(k,icright),uy(k,icright), &
                            uxx(k,icright),uxy(k,icright),uyy(k,icright),psi,eta)
                  end do

          end if

          if (ifacerc.eq.2) then
                 psi=0.5
                 eta=gp(nfpr)
                 do k=1,8
                    Qr(k)=calc_Q(u0(k,icright),ux(k,icright),uy(k,icright),&
                            uxx(k,icright),uxy(k,icright),uyy(k,icright),psi,eta)
                  end do
				  sign_r=-1

          end if

          if (ifacerc.eq.3) then
                 psi=gp(N-nfpr+1)
                 eta=0.5
                 do k=1,8
                    Qr(k)=calc_Q(u0(k,icright),ux(k,icright),uy(k,icright),&
                            uxx(k,icright),uxy(k,icright),uyy(k,icright),psi,eta)
                  end do
				  sign_r=-1

          end if

          if (ifacerc.eq.4) then
                 psi=-0.5
                 eta=gp(N-nfpr+1)
                 do k=1,8
                    Qr(k)=calc_Q(u0(k,icright),ux(k,icright),uy(k,icright),&
                            uxx(k,icright),uxy(k,icright),uyy(k,icright),psi,eta)
                  end do
          end if
		  

		  CALL getrusanovflux(Ql,Qr,Fnl,Fnr,normfl,sign_l,sign_r)

		  
		  do i=1,8
		      f_edge(i,nfp,ifacelc,icleft)=Fnl(i)
			  f_edge(i,nfpr,ifacerc,icright)=Fnr(i)
		  end do
		  		  
		  END DO

		  END DO
		  
		  
	!	  contains
	 ! DOUBLE PRECISION FUNCTION calc_Q(q0,q1,q2,q3,q4,q5,psi,eta)

	 !     IMPLICIT NONE
	 !     DOUBLE PRECISION,intent(in) :: q0,q1,q2,q3,q4,q5,psi,eta

	 ! 
	 !     calc_Q=q0*1+q1*psi+q2*eta+q3*(psi*psi-1./12.)+q4*psi*eta+q5*(eta*eta-1./12.)
	!	  
	!	 
     !END 
	 
       
		  END SUBROUTINE FLUXedge
		  
		  
		  SUBROUTINE calc_FG(ic,FC,GC)
		  use setup
		  implicit none
		  include 'MESH2D.INC'
		  
		  double precision,dimension(8,3,3), intent(out) :: FC,GC
		  double precision :: rho,u,v,VnF,VnG,rhoE,p
		  double precision,dimension(8,3,3) :: Q
		  double precision :: psi,eta,calc_Q
		  integer :: igp,jgp,i,j,k
		  integer, intent(in) :: ic
		  double precision :: s11,s12,s21,s22
		  double precision :: ww,Bx,By,Bz,NORM_U,NORM_B
		  double precision :: Fi(8),Gi(8)
		  
		  do igp=1,3
		     do jgp=1,3
			     psi=gp(igp)
				 eta=gp(jgp)
			   do k=1,8 !Q gauss points volume integral
			       Q(k,igp,jgp)=calc_Q(u0(k,ic),ux(k,ic),uy(k,ic),uxx(k,ic),uxy(k,ic),uyy(k,ic),psi,eta)
			   end do
			   
			   rho=Q(1,igp,jgp)
			   u=Q(2,igp,jgp)/rho
			   v=Q(3,igp,jgp)/rho
			   ww=Q(4,igp,jgp)/rho
			   rhoE=Q(5,igp,jgp)
			   Bx=Q(6,igp,jgp)
			   By=Q(7,igp,jgp)
			   Bz=Q(8,igp,jgp)
			   NORM_U=u**2+v**2+ww**2
			   NORM_B=Bx**2+By**2*Bz**2
			   p=(gama-1)*(rhoE-0.5*rho*NORM_U-0.5*NORM_B)
			   
			   Fi(1)=rho*u
			   Fi(2)=rho*u*u+p+(0.5d0*NORM_B-Bx**2)
			   Fi(3)=rho*u*v-Bx*By
			   Fi(4)=rho*u*ww-Bx*Bz
			   Fi(5)=u*(rhoE+p+0.5d0*NORM_B)-(u*Bx+v*By+ww*Bz)*Bx
			   Fi(6)=0.d0
			   Fi(7)=u*By-v*Bx
			   Fi(8)=u*Bz-ww*Bx
			   
			   Gi(1)=rho*v
			   Gi(2)=rho*u*v-Bx*By
			   Gi(3)=rho*v*v+p+(0.5d0*NORM_B-By**2)
			   Gi(4)=rho*v*ww-By*Bz
			   Gi(5)=v*(rhoE+p+0.5d0*NORM_B)-(u*Bx+v*By+ww*Bz)*By
			   Gi(6)=-(u*By-v*Bx)
			   Gi(7)=0.d0
			   Gi(8)=v*Bz-ww*By
			   
			   s11=djacobv(4,igp,jgp,ic)
			   s12=-djacobv(2,igp,jgp,ic)
			   s21=-djacobv(3,igp,jgp,ic)
			   s22=djacobv(1,igp,jgp,ic)
			   
			   VnF=s11*u+s12*v
			   VnG=s21*u+s22*v
			   
			   FC(1,igp,jgp)=s11*Fi(1)+s12*Gi(1)
			   FC(2,igp,jgp)=s11*Fi(2)+s12*Gi(2)
			   FC(3,igp,jgp)=s11*Fi(3)+s12*Gi(3)
			   FC(4,igp,jgp)=s11*Fi(4)+s12*Gi(4)
			   FC(5,igp,jgp)=s11*Fi(5)+s12*Gi(5)
			   FC(6,igp,jgp)=s11*Fi(6)+s12*Gi(6)
			   FC(7,igp,jgp)=s11*Fi(7)+s12*Gi(7)
			   FC(8,igp,jgp)=s11*Fi(8)+s12*Gi(8)
			   
			   GC(1,igp,jgp)=s21*Fi(1)+s22*Gi(1)
			   GC(2,igp,jgp)=s21*Fi(2)+s22*Gi(2)
			   GC(3,igp,jgp)=s21*Fi(3)+s22*Gi(3)
			   GC(4,igp,jgp)=s21*Fi(4)+s22*Gi(4)
			   GC(5,igp,jgp)=s21*Fi(5)+s22*Gi(5)
			   GC(6,igp,jgp)=s21*Fi(6)+s22*Gi(6)
			   GC(7,igp,jgp)=s21*Fi(7)+s22*Gi(7)
			   GC(8,igp,jgp)=s21*Fi(8)+s22*Gi(8)
			   
			 end do
		end do
		
	!	contains
	 ! DOUBLE PRECISION FUNCTION calc_Q(q0,q1,q2,q3,q4,q5,psi,eta)

	 !     IMPLICIT NONE
	 !     DOUBLE PRECISION,intent(in) :: q0,q1,q2,q3,q4,q5,psi,eta

	 ! 
	 !     calc_Q=q0*1+q1*psi+q2*eta+q3*(psi*psi-1./12.)+q4*psi*eta+q5*(eta*eta-1./12.)
	!	  
	!	
     !END 
	 
		
   END  SUBROUTINE calc_FG
   
   
			   
		  
		  
		  
		  
		  




