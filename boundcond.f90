SUBROUTINE BCflux
      use setup
      IMPLICIT NONE
      INCLUDE 'MESH2D.INC'
	  integer :: icleft,icright,ifacelc,ifacerc,nfp,nfpr,sign_l,sign_r,ifacelcB,icleftB,iface,ipair,k1,k2
	  double precision :: psi,eta,calc_Q
	  double precision,dimension(8) :: Ql,Qr
	  double precision,dimension(2) :: normfl
	  double precision,dimension(8) :: Fnl,Fnr
	  integer :: k,nf,i
	  integer :: flagl,flagr
	  
	  sign_r=1
	  do nf=1,NINLET
	      iface=IBFINL(nf)
	      
	      icleft=IF2C(iface,1)
		  icright=IF2C(iface,2)
		  
		  ifacelc=IF2C(iface,3)
		  ifacerc=IF2C(iface,4)
		  
		  do nfp=1,3
		     sign_l=1
			 if(ifacelc .eq. 1.) then
			      psi=gp(nfp)
				  eta=-0.5
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  sign_l=-1
				  normfl(1)=-djacobf(3,nfp,ifacelc,icleft)
				  normfl(2)=djacobf(1,nfp,ifacelc,icleft)
			 end if  
			
			 if(ifacelc .eq. 2.) then
			      psi=0.5
				  eta=gp(nfp)
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  normfl(1)=djacobf(4,nfp,ifacelc,icleft)
				  normfl(2)=-djacobf(2,nfp,ifacelc,icleft)
			 end if  
			 
			 if(ifacelc .eq. 3.) then
			      psi=gp(N-nfp+1)
				  eta=0.5
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  normfl(1)=-djacobf(3,nfp,ifacelc,icleft)
				  normfl(2)=djacobf(1,nfp,ifacelc,icleft)
			 end if 
			 
			 if(ifacelc .eq. 4.) then
			      psi=-0.5
				  eta=gp(N-nfp+1)
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  sign_l=-1
				  normfl(1)=djacobf(4,nfp,ifacelc,icleft)
				  normfl(2)=-djacobf(2,nfp,ifacelc,icleft)
			 end if
			 
			 Qr(1:8)=Ql(1:8)
			 !Qr(1)=rinf
			 !Qr(2)=uinf*rinf
			 !Qr(3)=rinf*vinf
			 !Qr(4)=pinf/(gama-1)+0.5*(Qr(2)**2+Qr(3)**2)/Qr(1)
			 CALL getrusanovflux(Ql,Qr,Fnl,Fnr,normfl,sign_l,sign_r)
			 
			 do k=1,8
			    f_edge(k,nfp,ifacelc,icleft)=Fnl(k)
		     end do
			 
			 end do
			 end do
			 
			 do nf=1,NOUTLET
			     
			     iface=IBFOUT(nf)
	              icleft=IF2C(iface,1)
		          icright=IF2C(iface,2)
		  
		          ifacelc=IF2C(iface,3)
		          ifacerc=IF2C(iface,4)
		  
		          do nfp=1,3
		              sign_l=1
			          if(ifacelc .eq. 1.) then
			         psi=gp(nfp)
				  eta=-0.5
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  sign_l=-1
				  normfl(1)=-djacobf(3,nfp,ifacelc,icleft)
				  normfl(2)=djacobf(1,nfp,ifacelc,icleft)
			 end if  
			
			 if(ifacelc .eq. 2.) then
			      psi=0.5
				  eta=gp(nfp)
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  normfl(1)=djacobf(4,nfp,ifacelc,icleft)
				  normfl(2)=-djacobf(2,nfp,ifacelc,icleft)
			 end if  
			 
			 if(ifacelc .eq. 3.) then
			      psi=gp(N-nfp+1)
				  eta=0.5
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  normfl(1)=-djacobf(3,nfp,ifacelc,icleft)
				  normfl(2)=djacobf(1,nfp,ifacelc,icleft)
			 end if 
			 
			 if(ifacelc .eq. 4.) then
			      psi=-0.5
				  eta=gp(N-nfp+1)
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  sign_l=-1
				  normfl(1)=djacobf(4,nfp,ifacelc,icleft)
				  normfl(2)=-djacobf(2,nfp,ifacelc,icleft)
			 end if
			 
			 Qr=Ql
			
			 CALL getrusanovflux(Ql,Qr,Fnl,Fnr,normfl,sign_l,sign_r)
			 
			 do k=1,8
			    f_edge(k,nfp,ifacelc,icleft)=Fnl(k)
		     end do
			 
			 end do
			 end do
		
			 
			 do ipair=1,NPAIR
			    k1=CYC2FPAIR(ipair,1)
				k2=CYC2FPAIR(ipair,2)
				
				icleft=IF2C(k1,1)
		        icright=IF2C(k2,1)
		  
		        ifacelc=IF2C(k1,3)
		        ifacerc=IF2C(k2,3)
				
				do nfp=1,3
		           sign_l=1
			      if(ifacelc .eq. 1.) then
			          psi=gp(nfp)
				      eta=-0.5
				      do k=1,8
				         Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			          end do
				     sign_l=-1
				     normfl(1)=-djacobf(3,nfp,ifacelc,icleft)
				     normfl(2)=djacobf(1,nfp,ifacelc,icleft)
			    end if  
			
			 if(ifacelc .eq. 2.) then
			      psi=0.5
				  eta=gp(nfp)
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  normfl(1)=djacobf(4,nfp,ifacelc,icleft)
				  normfl(2)=-djacobf(2,nfp,ifacelc,icleft)
			 end if  
			 
			 if(ifacelc .eq. 3.) then
			      psi=gp(N-nfp+1)
				  eta=0.5
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  normfl(1)=-djacobf(3,nfp,ifacelc,icleft)
				  normfl(2)=djacobf(1,nfp,ifacelc,icleft)
			 end if 
			 
			 if(ifacelc .eq. 4.) then
			      psi=-0.5
				  eta=gp(N-nfp+1)
				  do k=1,8
				     Ql(k)=calc_Q(u0(k,icleft),ux(k,icleft),uy(k,icleft),uxx(k,icleft),uxy(k,icleft),uyy(k,icleft),psi,eta)
			      end do
				  sign_l=-1
				  normfl(1)=djacobf(4,nfp,ifacelc,icleft)
				  normfl(2)=-djacobf(2,nfp,ifacelc,icleft)
			 end if
			 
		!	 flagl=0
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
		  
		end do
		end do
			    
			 !contains
	 ! FUNCTION calc_Q(q0,q1,q2,q3,q4,q5,psi,eta)


	 !     IMPLICIT NONE
	 !     DOUBLE PRECISION,intent(in) :: q0,q1,q2,q3,q4,q5,psi,eta
	!	  DOUBLE PRECISION :: calc_Q

	 ! 
	 !     calc_Q=q0*1+q1*psi+q2*eta+q3*(psi*psi-1./12.)+q4*psi*eta+q5*(eta*eta-1./12.)
	!	  
	!	  			 
	 !END FUNCTION calc_Q
			 
			 end subroutine BCflux
			 
			 
			 
			
			 
			 
			 
			 
			 
			 
			 
		  
		  
	  

