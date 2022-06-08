SUBROUTINE tecplotter(Npart)

use setup
IMPLICIT NONE
INCLUDE 'MESH2D.INC'

integer,intent(in) :: Npart
integer :: plotnodes,plotcells,Nfacenodes
integer :: k,is,js,newnode,newnode_r,nd,nodeid,nface,ic,jc
integer :: i,j,inode,jnode
integer :: IC2V(4),gnid,icleft,icright,ifacelc,ifacerc,iface,icell,gcid
integer,allocatable :: gnode_cellcount(:),g_index_node(:,:)
integer,allocatable :: ivnewcells(:,:),newfacenodes(:,:)
double precision,allocatable :: plotX(:,:),plotQ(:,:),mach(:),intene(:)
character(80) :: filename=''
character*16 :: tecfile

! MHD
!double precision :: psi,eta,calc_Q,Qleft(8),Qright(8),xxf(2,4),gam1,v2,xx(2),NORM_U,NORM_B

! Euler
double precision :: psi,eta,calc_Q,Qleft(4),Qright(4),xxf(2,4),gam1,v2,xx(2),NORM_U

plotnodes=0
plotcells=0
Nfacenodes=Npart-1
gam1=gama-1

allocate(plotX(2,MAXVRT))
allocate(mach(MAXVRT))
allocate(intene(MAXVRT))
allocate(gnode_cellcount(MAXVRT))
allocate(newfacenodes(Nfacenodes,NFACES))
allocate(ivnewcells(4,NCELL*Npart*Npart))
allocate(g_index_node(Npart+1,Npart+1))

! MHD
!allocate(plotQ(8,MAXVRT))

! Euler
allocate(plotQ(4,MAXVRT))

plotQ(:,:)=0
gnode_cellcount(:)=0

do icell=1,NCELL

   do k=1,4
       IC2V(k)=IVCELL(icell,k)
   end do
   do k=1,4
      gnid=IC2V(k)
	  select case(k)
	  case(1)
	      psi=-0.5
		  eta=-0.5
	  case(2)
	      psi=0.5
		  eta=-0.5
	  case(3)
	      psi=0.5
		  eta=0.5	  
	  case(4)
	      psi=-0.5
		  eta=0.5	  
	  end select

!	  MHD
!	  do i=1,8
!          plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
!                               uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
!       end do

!	  Euler	  
	  do i=1,4
	  plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell), &
                            uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
	  end do
	  
      if (gnode_cellcount(gnid)==0) then
	       plotnodes = plotnodes+1
		   plotX(1,gnid)=XV(gnid)
		   plotX(2,gnid)=YV(gnid)
	   end if
	   gnode_cellcount(gnid)=gnode_cellcount(gnid)+1
	 end do
	 
end do

do iface=1,NFACES
   icleft  = IF2C(iface,1)
   icright = IF2C(iface,2)
   ifacelc = IF2C(iface,3)
   ifacerc = IF2C(iface,4)
   
   do newnode=1,Nfacenodes
       plotnodes = plotnodes + 1
	   gnid = plotnodes
	   newfacenodes(newnode,iface) = gnid
	   
	   Qleft = 0.0
	   Qright = 0.0
		
		select case(ifacelc)
			case(1) 
				psi = real(newnode)/real(Npart)-0.5
				eta = -0.5
				
			case(2)
				psi = 0.5
				eta = real(newnode)/real(Npart)-0.5
				
			case(3)
				psi = real(Npart-newnode)/real(Npart)-0.5		
				eta = 0.5
				
			case(4)
				psi = -0.5
				eta = real(Npart-newnode)/real(Npart)-0.5		
				
			end select
			
!	  		MHD
!	  		do i=1,8
!    		      plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
!                               uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
!      		end do

!			  Euler	  
			do i=1,4
				plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell), &
									  uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
			end do

			do k=1,4
			   nodeid=IVCELL(icleft,k)
			   xxf(1,k)=XV(nodeid)
			   xxf(2,k)=YV(nodeid)
			end do
			
			CALL xyCoor_atGps(4,xxf,psi,eta,xx)
			plotX(1,gnid)=xx(1)
			plotX(2,gnid)=xx(2)
			
			if(trim(FACEMARKER(iface)) == 'INTE') then
				! Find the right cell coord and Q
				newnode_r = Npart - newnode
				select case(ifacerc)
				case(1) 
					psi = real(newnode_r)/real(Npart)-0.5
					eta = -0.5
				case(2)
					psi = 0.5
					eta = real(newnode_r)/real(Npart)-0.5
				case(3)
					psi = real(Npart-newnode_r)/real(Npart)-0.5
					eta = 0.5
				case(4)
					psi = -0.5
					eta = real(Npart-newnode_r)/real(Npart)-0.5	
				end select
						
				do is=1,N
				do js=1,N
!	  		MHD
!	  		do i=1,8
!    		      plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
!                               uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
!      		end do

!			  Euler	  
			do i=1,4
				plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell), &
									  uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
			end do
				end do
				end do
				
			else
			Qright=Qleft
			end if
			
			plotQ(:,gnid)=Qleft(:)+Qright(:)
			gnode_cellcount(gnid)=2
    end do
end do

do icell=1,NCELL
    do k=1,4
	   nodeid=IVCELL(icell,k)
	   xxf(1,k)=XV(nodeid)
	   xxf(2,k)=YV(nodeid)
	 end do
	 
	 do inode=1,Nfacenodes
	    do jnode=1,Nfacenodes
		    plotnodes=plotnodes+1
			gnid=plotnodes
			psi=real(inode)/real(Npart)-0.5
			eta=real(jnode)/real(Npart)-0.5
			
			plotQ(1:4,gnid)=0.0
			
!	  		MHD
!	  		do i=1,8
!    		      plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell),&
!                               uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
!      		end do

!			  Euler	  
			do i=1,4
				plotQ(i,gnid)=plotQ(i,gnid)+calc_Q(u0(i,icell),ux(i,icell),uy(i,icell), &
									  uxx(i,icell),uxy(i,icell),uyy(i,icell),psi,eta)
			end do
			gnode_cellcount(gnid)=1
			g_index_node(inode+1,jnode+1)=gnid
			CALL xyCoor_atGps(4,xxf,psi,eta,xx)
			plotX(1,gnid)=xx(1)
			plotX(2,gnid)=xx(2)
		end do
	end do
	 
	 g_index_node(1,1) = IVCELL(icell,1)
	 g_index_node(Npart+1,1) = IVCELL(icell,2)
	 g_index_node(Npart+1,Npart+1) = IVCELL(icell,3)
	 g_index_node(1,Npart+1) = IVCELL(icell,4)
	 
	 nface=1
		iface = IC2F(icell,nface)
		icleft  = IF2C(iface,1)
		if(icell==icleft) then
			do newnode=1,Nfacenodes
				g_index_node(newnode+1,1) = newfacenodes(newnode,iface)
			end do
		else
			do newnode=1,Nfacenodes
				g_index_node(Npart-newnode+1,1) = newfacenodes(newnode,iface)
			end do
		end if	
		
			
		nface=2
		iface = IC2F(icell,nface)
		icleft  = IF2C(iface,1)
		
		if(icell==icleft) then
			do newnode=1,Nfacenodes
				g_index_node(Npart+1,newnode+1) = newfacenodes(newnode,iface)
			end do
		else
			do newnode=1,Nfacenodes
				g_index_node(Npart+1,Npart-newnode+1) = newfacenodes(newnode,iface)
			end do
		end if	
		
		nface=3
		iface = IC2F(icell,nface)
		icleft  = IF2C(iface,1)
		
		if(icell==icleft) then
			do newnode=1,Nfacenodes
				g_index_node(Npart-newnode+1,Npart+1) = newfacenodes(newnode,iface)
			end do
		else
			do newnode=1,Nfacenodes
				g_index_node(newnode+1,Npart+1) = newfacenodes(newnode,iface)
			end do
		end if		
		
																						
		nface=4
		iface = IC2F(icell,nface)
		icleft  = IF2C(iface,1)
		
		if(icell==icleft) then
			do newnode=1,Nfacenodes
				g_index_node(1,Npart-newnode+1) = newfacenodes(newnode,iface)
			end do
		else
			do newnode=1,Nfacenodes
				g_index_node(1,newnode+1) = newfacenodes(newnode,iface)
			end do
		end if	
		
	
		do ic=1,Npart
		do jc=1,Npart
			plotcells = plotcells + 1
			gcid = plotcells
			
			ivnewcells(1,gcid) = g_index_node(ic,jc)
			ivnewcells(2,gcid) = g_index_node(ic+1,jc)
			ivnewcells(3,gcid) = g_index_node(ic+1,jc+1)
			ivnewcells(4,gcid) = g_index_node(ic,jc+1)
		end do
		end do
		
end do

! MHD
!  do nd=1,plotnodes
!      plotQ(:,nd)=plotQ(:,nd)/real(gnode_cellcount(nd))
!	  plotQ(2,nd) = plotQ(2,nd)/plotQ(1,nd)
!	  plotQ(3,nd) = plotQ(3,nd)/plotQ(1,nd)
!	  plotQ(4,nd) = plotQ(4,nd)/plotQ(1,nd)
!	  NORM_U=plotQ(2,nd)**2+plotQ(3,nd)**2+plotQ(4,nd)**2
!	  NORM_B=plotQ(6,nd)**2+plotQ(7,nd)**2+plotQ(8,nd)**2
!	  plotQ(5,nd) = (gama-1)*(plotQ(5,nd)-0.5*plotQ(1,nd)*NORM_U-0.5*NORM_B)
!!	  v2=plotQ(2,nd)**2+plotQ(3,nd)**2
!!	  plotQ(4,nd) = gam1*(plotQ(4,nd)-0.5*plotQ(1,nd)*v2)
!!	  mach(nd) = sqrt(plotQ(2,nd)**2+plotQ(3,nd)**2) / &
!!        		   sqrt(gama*plotQ(4,nd)/plotQ(1,nd))
!!	 intene(nd) = plotQ(4,nd)/(gam1*plotQ(1,nd))
!   end do

! Euler
	do nd=1,plotnodes
		plotQ(:,nd)=plotQ(:,nd)/real(gnode_cellcount(nd))
		plotQ(2,nd) = plotQ(2,nd)/plotQ(1,nd)
		plotQ(3,nd) = plotQ(3,nd)/plotQ(1,nd)
		v2=plotQ(2,nd)**2+plotQ(3,nd)**2
		plotQ(4,nd) = gam1*(plotQ(4,nd)-0.5*plotQ(1,nd)*v2)
		mach(nd) = sqrt(plotQ(2,nd)**2+plotQ(3,nd)**2) / &
					 sqrt(gama*plotQ(4,nd)/plotQ(1,nd))
	    intene(nd) = plotQ(4,nd)/(gam1*plotQ(1,nd))
	end do
   
   write(filename,'(a,i6.6,a)')'testec',iter,'.dat'
   tecfile=trim(filename)

!   MHD   
!   open (unit=1,file =tecfile, form='formatted')
!    write(1, *) 'TITLE = "DG2D QUAD MESH"'
!    write(1, *) 'VARIABLES = "X", "Y", "ROU", "U", "V", "W", "P", "Bx", "By", "Bz"'
!    write(1, *) 'ZONE T="REFINED MESHES"    N=', plotnodes, "E=", plotcells, "DATAPACKING=POINT"
!    write(1, *) 'ZONETYPE=FEQuadrilateral'
!	
!	do i=1,plotnodes
!       	write(1,1000) (plotX(j,i),j=1,2),(plotQ(is,i),is=1,8)
!    end do
!
!    do i = 1,plotcells
!    	write(1,1001) (ivnewcells(j,i),j=1,4)
!    end do
!
!1000    format (10(G16.8,1X))
!1001    format (4(I9,1X))
!
!    close(1)

!	Euler
   open (unit=1,file =tecfile, form='formatted')
   write(1, *) 'TITLE = "DG2D QUAD MESH"'
   write(1, *) 'VARIABLES = "X", "Y", "ROU", "U", "V", "P", "Mach", "Temp"'
   write(1, *) 'ZONE T="REFINED MESHES"    N=', plotnodes, "E=", plotcells, "DATAPACKING=POINT"
   write(1, *) 'ZONETYPE=FEQuadrilateral'
   
   do i=1,plotnodes
		  write(1,1000) (plotX(j,i),j=1,2),(plotQ(is,i),is=1,4),mach(i),intene(i)
   end do

   do i = 1,plotcells
	   write(1,1001) (ivnewcells(j,i),j=1,4)
   end do

1000    format (8(G16.8,1X))
1001    format (4(I9,1X))

   close(1)
	deallocate(plotQ,plotX,mach)
	
	END SUBROUTINE tecplotter
			
   
