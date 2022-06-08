MODULE setup
      IMPLICIT NONE
	  !character(64) :: NAMECEL,NAMEVRT,NAMEBND2D
	  !integer, allocatable :: IVCELL(:,:),IVBOUN(:,:)
	  !integer :: NCELL,NVERT,NBOUN,NFACES,NCELLTOT
	  !character(8),allocatable ::IBOUNTYPE(:)
	  !double precision,allocatable ::XV(:),YV(:)
	  
	  !integer :: NINLET,NOUTLET,NSYMM,NWALL,NCYCL
	  !INTEGER,ALLOCATABLE :: IBFINL(:),IBFOUT(:),IBFSYMM(:),IBFCYCL(:),IBFWAL(:)
	  !integer,allocatable :: BOUNFACEINDEX(:)
	  !integer,allocatable :: IVCELL(:,:),IC2F(:,:),IF2V(:,:),IF2C(:,:)
	  !integer,allocatable :: ICVERT(:),ICVSTA(:)
	  !integer,dimension(4,2) :: IVFACE
	  
      double precision,dimension(:,:,:,:),allocatable :: F1,F2,resid,SRE
	  double precision	:: rinf,pinf,uinf,vinf

      double precision,dimension(:,:,:,:),allocatable :: etiny,imres,imresold
      double precision :: gama,minl,ft,dt,ndt
      double precision,dimension(:,:),allocatable :: u0,ux,uy
      double precision,dimension(:,:),allocatable :: uxx,uxy,uyy
      double precision,dimension (:,:),allocatable :: w10,w1x,w1y
      double precision,dimension(:,:),allocatable :: w1xx,w1xy,w1yy
      double precision,dimension(:,:),allocatable :: right0,rightx,righty
      double precision,dimension(:,:),allocatable :: rightxx,rightxy,rightyy
      double precision,dimension(:,:),allocatable :: o0,ox,oy
      double precision,dimension(:,:),allocatable :: oxx,oxy,oyy
      double precision,dimension(:,:,:,:),allocatable :: fv_vedge,fv_hedge
      double precision,dimension(:,:,:,:),allocatable :: djacobf,djacobv
	  double precision,dimension(:,:,:),allocatable :: djacobvertex,jacobf,jacobv,jacobtecint,jacobtecv
	  double precision,dimension(:,:),allocatable :: jacobvertex
	  double precision,dimension(:,:,:,:),allocatable :: f_edge
      double precision :: umesh,CFL
	  double precision,dimension(:),allocatable :: area
	  double precision,dimension(:,:,:,:),allocatable :: Qv  ! for volume integral
	  double precision,dimension(:,:,:,:),allocatable :: Qub !for area integral upper and bottom edge
	  double precision,dimension(:,:,:,:),allocatable :: Qlr !for area integral left and right edge
	  double precision,dimension(:,:,:,:),allocatable :: Qa
	  double precision,dimension(:,:),allocatable:: iface2fp,jface2fp

      double precision,dimension (:),allocatable :: rmass,gp,xmod,w
      integer,dimension (:),allocatable :: jmap

      double precision :: r13,r23,haf,r12,r06,dx,dy
	  double precision,parameter :: tolCYC=1d-5
      integer :: nx,ny,np,ndof,RK
      integer :: kcount,MAXITER
      parameter(r13=1./3.,r23=2./3.)
      parameter(r12=1./12.,r06=1./6.,haf=0.5)
      parameter(np=3)
      parameter(gama=5.d0/3)
      END MODULE setup

