SUBROUTINE getrusanovflux(Ql,Qr,Fnl,Fnr,normfl,sign_l,sign_r)

      use setup
      IMPLICIT NONE
	  double precision,dimension(8),intent(in) :: Ql,Qr
	  double precision,dimension(8),intent(out) :: Fnl,Fnr
	  double precision,dimension(2),intent(in) :: normfl
	  integer,intent(in):: sign_l,sign_r
	  double precision :: rhol,rmxl,rmyl,rhoEl,pl,eigv
	  double precision :: rhor,rmxr,rmyr,rhoEr,pr
	  double precision :: sonicl,sonicr,lambdal,lambdar
	  double precision :: ul,ur,vl,vr,lambda,magn
	  double precision :: Vn_l,Vn_r,unit_x,unit_y
	  double precision :: rhoi,c_a,c_d,c_n,c_tmp,c_f
	  double precision :: Bn_l,Bn_r,NORM_Bl,NORM_Br,Bx,By,Bz,wl,wr, &
	                       Bxl,Byl,Bzl,Bxr,Byr,Bzr
	  integer :: k
	  
	  magn=sqrt(normfl(1)**2+normfl(2)**2)
	  unit_x=normfl(1)/magn*sign_l
	  unit_y=normfl(2)/magn*sign_l
	  
	  rhoi=1.d0/Ql(1)
	  ul=rhoi*Ql(2)
	  vl=rhoi*Ql(3)
	  wl=rhoi*Ql(4)
	  Bxl=Ql(6)
	  Byl=Ql(7)
	  Bzl=Ql(8)
	  
	  NORM_Bl=Bxl**2+Byl**2+Bzl**2
	  pl=(gama-1.d0)*(Ql(5)-0.5*Ql(1)*(ul**2+vl**2+wl**2))
	  pl=pl+0.5*NORM_Bl
	  Vn_l=ul*unit_x+vl*unit_y
	  Bn_l=Bxl*unit_x+Byl*unit_y
	  
	  rhoi=1.d0/Qr(1)
	  ur=rhoi*Qr(2)
	  vr=rhoi*Qr(3)
	  wr=rhoi*Qr(4)
	  Bxr=Qr(6)
	  Byr=Qr(7)
	  Bzr=Qr(8)
	  
	  NORM_Br=Bxr**2+Byr**2+Bzr**2
	  pr=(gama-1.d0)*(Qr(5)-0.5*Qr(1)*(ur**2+vr**2+wr**2))
	  pr=pr+0.5*NORM_Br
	  Vn_r=ur*unit_x+vr*unit_y
	  Bn_r=Bxr*unit_x+Byr*unit_y
	  
	  
	  Fnl(1)=Ql(1)*Vn_l
	  Fnl(2)=Ql(2)*Vn_l + pl*unit_x - Bxl*Bn_l
	  Fnl(3)=Ql(3)*Vn_l + pl*unit_y - Byl*Bn_l
	  Fnl(4)=Ql(4)*Vn_l -Bzl*Bn_l
	  Fnl(5) = (Ql(5)+pl)*Vn_l - Bn_l*(ul*Bxl+vl*Byl+wl*Bzl)
	  Fnl(6) = (vl*Bxl-ul*Byl)*unit_y
      Fnl(7) = (ul*Byl-vl*Bxl)*unit_x
      Fnl(8) = Bzl*Vn_l - wl*Bn_l 
	  
	  Fnr(1)=Qr(1)*Vn_r
	  Fnr(2)=Qr(2)*Vn_r + pr*unit_x - Bxr*Bn_r
	  Fnr(3)=Qr(3)*Vn_r + pr*unit_y - Byr*Bn_l
	  Fnr(4)=Qr(4)*Vn_r -Bzr*Bn_r
	  Fnr(5) = (Qr(5)+pr)*Vn_r - Bn_r*(ur*Bxr+vr*Byr+wr*Bzr)
	  Fnr(6) = (vr*Bxr-ur*Byr)*unit_y
      Fnr(7) = (ur*Byr-vr*Bxr)*unit_x
      Fnr(8) = Bzr*Vn_r - wr*Bn_r 
	  
	  pl = pl -  0.5d0*NORM_Bl
      pr = pr -  0.5d0*NORM_Br
      Bx = 0.5d0*(Bxl+Bxr)
      By = 0.5d0*(Byl+Byr)
      Bz = 0.5d0*(Bzl+Bzr)
	  
	  c_a  = gama*(pl+pr)/(Ql(1)+Qr(1))
      c_d  = (Bx**2+By**2+Bz**2)/(0.5d0*(Ql(1)+Qr(1)))    
      c_n  = (Bx*unit_x+By*unit_y)**2/(0.5d0*(Ql(1)+Qr(1)))
      c_tmp= dsqrt((c_a+c_d)**2-4.0d0*c_a*c_n)
    
      c_f  = dsqrt(0.5d0*(c_a+c_d+c_tmp))
    
     eigv = 0.5d0*dabs(Vn_l+Vn_r) + c_f
	  
	  
	  do k=1,8
	      Fnl(k)=0.5*(Fnl(k)+Fnr(k)-eigv*(Qr(k)-Ql(k)))*magn*sign_l
		  Fnr(k)=sign_l*sign_r*Fnl(k)
	  end do

	  
	  END SUBROUTINE getrusanovflux
	  
	  
	  
		  
	  
