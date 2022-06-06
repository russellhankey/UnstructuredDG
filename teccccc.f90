SUBROUTINE tecplotnorefinement
       

use setup
IMPLICIT NONE
INCLUDE 'MESH2D.INC'

integer :: i,j,is,k
integer :: IC2V(4),icleft,icright,ifacelc,ifacerc,iface,icell,gnid
double precision :: psi,eta,calc_Q,Qleft(8),xxf(2,4),gam1,v2,xx2(2),NORM_U,NORM_B
double precision,allocatable :: plotX(:,:),plotQ(:,:),mach(:),intene(:)
integer,allocatable :: gnode_cellcount(:)
character(80) :: filename=''
character*16 :: tecfile

gam1=gama-1

allocate(plotX(2,NVERT))
allocate(plotQ(8,NVERT))
allocate(mach(NVERT))
allocate(intene(NVERT))
allocate(gnode_cellcount(NVERT))

plotQ(:,:)=0
gnode_cellcount(:)=0

do icell =1,NCELL
   do k=1,4
      IC2V(k)=IVCELL(icell,k)
   end do
   do k=1,4
      gnid=IC2V(k)
      select case(k)
      case(1)
         psi=-0.5
         eta=-0.5
		 do i=1,8
          plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
                               uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
       end do
      case(2)
          psi=0.5
          eta=-0.5
		  do i=1,8
          plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
                               uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
       end do
       case(3)
          psi=0.5
          eta=0.5
		  do i=1,8
          plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
                               uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
       end do
       case(4)
          psi=-0.5
          eta=0.5
		  do i=1,8
          plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
                               uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
       end do
       end select
       
       !do i=1,4
       !   plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
       !                        uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)/jacobvertex(I,ic)
       !end do
       
       if(gnode_cellcount(gnid)==0) then
            plotX(1,gnid)=XV(gnid)
            plotX(2,gnid)=YV(gnid)
        end if
        gnode_cellcount(gnid)=gnode_cellcount(gnid)+1
     end do
end do

do i=1,NVERT
   plotQ(1:8,i)=plotQ(1:8,i)/real(gnode_cellcount(i))
   plotQ(2,i)=plotQ(2,i)/plotQ(1,i)
   plotQ(3,i)=plotQ(3,i)/plotQ(1,i)
   plotQ(4,i) = plotQ(4,i)/plotQ(1,i)
	  NORM_U=plotQ(2,i)**2+plotQ(3,i)**2+plotQ(4,i)**2
	  NORM_B=plotQ(6,i)**2+plotQ(7,i)**2+plotQ(8,i)**2
	  plotQ(5,i) = (gama-1)*(plotQ(5,i)-0.5*plotQ(1,i)*NORM_U-0.5*NORM_B)
   !v2=plotQ(2,i)**2+plotQ(3,i)**2
   !plotQ(4,i)=gam1*(plotQ(4,i)-0.5*plotQ(1,i)*v2)
   !mach(i)=sqrt(v2)/sqrt(gama*plotQ(4,i)/plotQ(1,i))
   !intene(i)=plotQ(4,i)/(gam1*plotQ(1,i))
end do

write(filename,'(a,i6.6,a)')'testec',iter,'.dat'
tecfile=trim(filename)

open(unit=1,file=tecfile,form='formatted')
  write(1,*)'TITLE="DG2G QUAD MESH"'
  write(1,*) 'VARIABLES= "X", "Y", "ROU", "U", "V", "P", "Mach", "Temp"'
  write(1,*) 'ZONE T="REFINED MESHES"    N=',NVERT,"E=", NCELL, "DATAPACKING=POINT"
  write(1, *) 'ZONETYPE=FEQuadrilateral'
  
 do i=1,NVERT
    write(1,1000)(plotX(j,i),j=1,2),(plotQ(is,i),is=1,8)
 end do
 
 do i=1,NCELL
     write(1,1001) (IVCELL(i,j),j=1,4)
 end do
 
1000  format (10(G16.8,1X))
1001  format (4(I9,1X))

    close(1)
    deallocate(plotQ,plotX,mach,intene)
    
 

END

