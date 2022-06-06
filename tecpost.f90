SUBROUTINE xyCoor_atGps(ns,xxs,xi,eta,xx)

IMPLICIT NONE
integer,intent(in) :: ns
integer :: i,nb
double precision :: xxs(2,ns),xx(2),xi,eta
double precision,pointer,dimension(:) :: func
allocate(func(ns))
if(ns==4) call MapBaseFunc(func,xi,eta)
!if(ns==8) call MapBaseFunc8(func,xi,eta)
!if(ns==12) call MapBaseFunc(func,xi,eta)

do i=1,2
  xx(i)=0.0
  do nb=1,ns
     xx(i)=xx(i)+xxs(i,nb)*func(nb)
  end do
 end do
 deallocate(func)
 return 
 END SUBROUTINE xyCoor_atGps
