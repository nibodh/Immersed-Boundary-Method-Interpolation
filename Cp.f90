! ! !!!_______________PM________________!!!
program airfoil

IMPLICIT NONE
integer,parameter :: mytype = 8
integer, parameter :: nx=481,ny=288,nz=14, nclx=2, ncly=0, nclz=0
!  integer, parameter :: nx1=500 , ny1=500
 real(mytype), parameter:: dt=0.0005
real(mytype) :: dx, dy,theta,xo,yo,m, xmax, xmin, ymax, ymin, d2min,x1,x2, y1, y2,p1,p2,p3,p4
real(mytype) :: p11,p21,p12,p22,p, delta
integer :: i,j,out,k,l, counter, count, q,r,s, no, case
! real(4), dimension(nx1,ny1) :: ux1, uy1, x1, y1
integer, parameter :: n=159, dims=2
real(mytype), dimension(n) :: x,y
real(mytype), dimension(n,nz) :: pre, pavg
real(mytype), dimension(n) :: preFinal
real(mytype), dimension(n,dims) :: Po ! POINTS ARRAY	
real(mytype), dimension(nx,ny,nz)::pp
real(mytype), parameter::  xlx=20 , yly=12, zlz=1
real(mytype), dimension(nx,ny):: xgrid,ygrid

integer::im1,ip1
real(mytype)::dist,dx1,dy1

integer :: kAvg

if (nclx==0) dx=xlx/nx 
if (nclx==1 .or. nclx==2) dx=xlx/(nx-1.) 
if (ncly==0) dy=yly/ny !(ny-1.)
if (ncly==1.or.ncly==2) dy=yly/(ny-1.)
delta = 0.00005 + (dx**2 + dy**2)**0.5

OPEN(10,FILE='Cylinder-test.dat',FORM='FORMATTED',&
        RECL=4, STATUS='OLD')
  
     DO i= 1,n
           READ(10,*) Po(i,1), Po(i,2)						! reading points file
     ENDDO
    
  CLOSE(10)
  
   OPEN(10,FILE='pmean-2.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO k=1,nz
     DO j=ny,1, -1
        DO i=1,nx
           READ(10,REC=COUNT) pp(i,j,k)						! reading pressure values
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(10)
  
!  print *, '1,1', pp(1,1,19)/4000
!print *, '2,1', pp(2,1,19)/4000
!print *, '1,2', pp(1,2,19)/4000
!print *, '241,337', pp(481,337,19)/4000
!print *, '1,337', pp(1,337,19)/4000
!print *, '481,1', pp(481,1,1)/4000

 !    DO j=165,173
  !      DO i=104,114
!	    print *, i, j, pp(i,j,19)/4000	
 !       ENDDO
  !   ENDDO
  
  
!   do i=1,400
!   
!   im1=i-1
!   ip1=i+1
!   if(i==1)im1=399
!   if(i==400)ip1=2
!   
!   dx1=Po(ip1,1)-Po(im1,1)
!   dy1=Po(im1,2)-Po(ip1,2)
! 
!   dist = dsqrt(dx1**2+dy1**2)
! 
!   dx1 = dx1/dist
!   dy1 = dy1/dist
! 
!   x(i) = Po(i,1) + dx1 * delta
!   y(i) = Po(i,2) + dy1 * delta
!   
!   end do
  
!**************************************************************************

!  delta = 3.5e-2

  do i=1,n

  im1=i-1
  ip1=i+1
  if(i==1)im1= n-1
  if(i==n)ip1= 2

  !if(i==200)ip1 = i
  !if(i==201)im1 = i



  dx1=Po(ip1,1)-Po(im1,1)
  dy1=Po(ip1,2)-Po(im1,2)

  dist = dsqrt(dx1**2+dy1**2)

  dx1 = dx1/dist
  dy1 = dy1/dist

  x(i) = Po(i,1) - dy1 * delta
  y(i) = Po(i,2) + dx1 * delta

  
  end do

OPEN(UNIT=3,FILE='points_boundary.dat', FORM = 'FORmatted', status = 'replace')
do i=1,n
  WRITE(3,*) x(i), y(i)!-pp(2,ny-4,nz/2))
enddo

CLOSE(3)
!-------------------------------------------------------------added from interp_pre.f90----------------------------
! 
d2min=dx*dx


xmax=3.0
xmin=2.000036
ymax=1.060009
ymin=0.939991

! defining actual coordinates
s=(xmax+8*dx)/dx+1
r=(xmin-8*dx)/dx+1
!p=(ymin-8*dy)/dy+1
q=(ymax+8*dy)/dy+1

 !print *, r, s, p,q


! STOP

do j=1,ny
do i=1,nx

xgrid(i,j)=(i-1)*dx
ygrid(i,j)=(j-1)*dy

enddo
enddo
OPEN(UNIT=3,FILE='interpolant_points.dat', FORM = 'FORmatted', status = 'replace')

do kAvg=1,nz
print *,'k loop:',kAvg,'out of',nz
do i=1,n !no
count=1
!   do l=p,q
!   do k=r,s
  do l= 5*ny/12, 7*ny/12
  do k= 5*nx/24, 7*nx/24
	x1 = xgrid(k,l)
	y1 = ygrid(k,l)
	x2 = xgrid(k+1,l)
	y2 = ygrid(k,l+1)
	p11 = pp(k,l,kAvg)
	p21 = pp(k+1,l,kAvg)
	p12 = pp(k,l+1,kAvg)
	p22 = pp(k+1,l+1,kAvg)

	call mesh_cell(x1,y1,x2,y2,x(i),y(i),p11,p21,p12,p22,p,case,out)

!      d2=(x(i)-xgrid(k,l))*(x(i)-xgrid(k,l))+(y(i)-ygrid(k,l))*(y(i)-ygrid(k,l))
!  !  print *, d2
!!  if (sqrt((x(k,l)-xc)*(x(k,l)-xc)+(x(k,l)-xc)*(x(k,l)-xc))>rad) then
!      if (d2<=d2min .and. count==1) then
!!     x1=x(k,l)
!!     y1=y(k,l)
!      call polygon_test(xgrid(k,l),ygrid(k,l),out)
!      if (out==1) then
!            d_1=d2
!      p1=pp(k,l,kAvg)
!      count=2
!      
!! 	  print *, x(i),',',y(i),',',xgrid(k,l),',',ygrid(k,l)
!      endif
      
!    elseif (d2<=d2min .and. count==2) then
!!     x2=x(k,l)
!! 	    y2=y(k,l)
!	call polygon_test(xgrid(k,l),ygrid(k,l),out)
!	if (out==1) then
!	  d_2=d2
!	  p2=pp(k,l,kAvg)
!	  count=3
!      WRITE(3,*) xgrid(k,l),ygrid(k,l)
!! 	  print *, x(i),',',y(i),',',xgrid(k,l),',',ygrid(k,l)
!	endif

!    elseif (d2<=d2min .and. count==3) then
!! 	    x3=x(k,l)
!! 	    y3=y(k,l)
!	  call polygon_test(xgrid(k,l),ygrid(k,l),out)
!	   if (out==1) then
!	      d_3=d2
!	        p3=pp(k,l,kAvg)
!	      count=4
!      WRITE(3,*) xgrid(k,l),ygrid(k,l)
!	      ! 	      print *, x(i),',',y(i),',',xgrid(k,l),',',ygrid(k,l)
!	    endif
   
!    elseif (d2<=d2min .and. count==4) then
!!     x4=x(k,l)
!! 	    y4=y(k,l)
!	    call polygon_test(xgrid(k,l),ygrid(k,l),out)
!	    if (out==1) then
!	      d_4=d2
!	      p4=pp(k,l,kAvg)
!	      count=5
!      WRITE(3,*) xgrid(k,l),ygrid(k,l)
!! 	      print *, x(i),',',y(i),',',xgrid(k,l),',',ygrid(k,l)
!	    endif

! !    print *, xl(i,j),yl(i,j), x1, y1,x2,y2,x3,y3, x4, y4
!    endif
    
 !    pre(i,j)=(p1+p2+p3+p4)/4/dt				! basic averaging
  !   print *, pre(i,j)
    if (out == 0) then
    pre(i,kAvg) = (1/dt)*p
	if (case ==1 ) then
	WRITE(3,*) x1 ,y1
	elseif (case ==2 ) then
	WRITE(3,*) x1, y1
	WRITE(3,*) x2, y1
	elseif (case ==3 ) then
	WRITE(3,*) x1, y1
	WRITE(3,*) x1, y2
	elseif (case ==4 ) then
	WRITE(3,*) x1, y1
	WRITE(3,*) x2, y1
	WRITE(3,*) x1, y2
	WRITE(3,*) x2, y2
	endif	
    exit
    end if
   ! pre(i,kAvg)=(1/dt)*((p1/d_1)+(p2/d_2)+(p3/d_3)+(p4/d_4))/((1./d_1)+(1./d_2)+(1./d_3)+(1./d_4))
   !    pre(i,j)=1./dt*(p1/(d_1*d_1)+p2/(d_2*d_2)+p3/(d_3*d_3)+p4/(d_4*d_4))/(1./(d_1*d_1)+1./(d_2*d_2)+1./(d_3*d_3)+1./(d_4*d_4))		! IDW p=2
!    endif
    enddo
    enddo
enddo
enddo!kAvg
CLOSE(3)
! i=1
! d2min=dx*dx
! count=1
!   do k=1,nx/2
!   do l=1,ny
!     d2=(xl(i,j)-x(k,l))*(xl(i,j)-x(k,l))+(yl(i,j)-y(k,l))*(yl(i,j)-y(k,l))
!   !  print *, d2
!   if (sqrt((x(k,l)-xc)*(x(k,l)-xc)+(y(k,l)-yc)*(y(k,l)-yc))>rad) then
!     if (d2<=d2min .and. count==1) then
! !     x1=x(k,l)
! !     y1=y(k,l)
!     d_1=d2
!     p1=pp(k,l,nz/2)
!     count=2
!     elseif (d2<=d2min .and. count==2) then
! !     x2=x(k,l)
! !     y2=y(k,l)
!     d_2=d2
!     p2=pp(k,l,nz/2)
!     count=3
! !     elseif (d2<=d2min .and. count==3) then
! ! !     x3=x(k,l)
! ! !     y3=y(k,l)
! !     d_3=d2
! !     p3=pp(k,l,nz/2)
! !     count=4
! !     elseif (d2<=d2min .and. count==4) then
! ! !     x4=x(k,l)
! ! !     y4=y(k,l)
! !     d_4=d2
! !     p4=pp(k,l,nz/2)
! !     count=5
! !  !    print *, xl(i,j),yl(i,j), x1, y1,x2,y2,x3,y3, x4, y4
!     endif
!     
!  !    pre(i,j)=(p1+p2+p3+p4)/4/dt				! basic averaging
!   !   print *, pre(i,j)
!     pre(i,j)=1./dt*(p1/d_1+p2/d_2)/(1./d_1+1./d_2)
!     endif
!     enddo
!     enddo
! enddo
! print *, no
! do i=1,no
! print *, pre(i)
! enddo
print *, 'dx', dx, 'dy', dy, 'nx', nx,'ny', ny, 'X', xgrid(1, ny), 'to', xgrid(nx,ny), 'Y', ygrid(nx,1), 'to', ygrid (nx, ny)
print *, 'geometry center', xgrid(121,145), ygrid(121, 145)
print *, 'geometry horizontal', xgrid(109,145), xgrid(133,145)
print *, 'geometry vertical', ygrid(121, 133), ygrid(121, 157)
print *, 'pressure', pp(121, 157, 7), pp(121, 133, 7), pp(113, 156, 7), pp(113, 134, 7), pp(129, 156, 7), pp(129, 134, 7)
print *, 'pressure reflection', pp(121, 158, 7), pp(121, 130, 7), pp(113, 156, 7), pp(113, 132, 7), pp(129, 156, 7), pp(129, 132, 7)
print *, 'pressure reflection', pp(121, 158, 7), pp(121, 131, 7), pp(113, 156, 7), pp(113, 133, 7), pp(129, 156, 7), pp(129, 133, 7)

OPEN(UNIT=3,FILE='Cp.dat', FORM = 'FORmatted', status = 'replace')

preFinal = 0.d0
do kAvg=1,nz
do i=1,n!no
  preFinal(i) = preFinal(i)+pre(i,kAvg) 
enddo  
enddo
preFinal = preFinal/nz

do i=1,n!no
	WRITE(3,*) x(i), 2*preFinal(i)/0.33/4000!-pp(2,ny-4,nz/2))
enddo

CLOSE(3)

end program airfoil


subroutine mesh_cell(x1,y1,x2,y2,x,y,p11,p21,p12,p22,p,case, out)
IMPLICIT NONE
integer,parameter :: mytype = 8
integer, parameter :: dims=2,n=100		! n- no of points, dims- no of dimensions
integer :: crossings, COUNT, i, j, onbound=0, case,out!,n
real(mytype) :: x,y,x1,y1,x2,y2,p11,p21,p12,p22,p
out = 0

if (x == x1 .and. y == y1 ) then !point same as grid point
case =1
p = p11

elseif (y == y1 .and. x > x1 .and. x < x2 ) then! linear x interpolation
case =2
p = (p11*(x2-x) + p21*(x-x1))/(x2-x1)

elseif (x == x1 .and. y > y1 .and. y < y2 ) then ! linear y interpolation
case =3
p = (p11*(y2-x) + p12*(y-y1))/(y2-y1)

elseif (x > x1 .and. x < x2 .and. y > y1 .and. y < y2) then !point perfectly between 4 points
case =4
p = (p11*(x2-x)*(y2-y) + p21*(x-x1)*(y2-y) + p12*(x2-x)*(y-y1) + p22*(x-x1)*(y-y1) )/((x2-x1)*(y2-y1))

Else
out = 1

endif

end subroutine mesh_cell




  

subroutine polygon_test(x,y,out)

! use BSI_points
IMPLICIT NONE
integer,parameter :: mytype = 8
integer, parameter :: dims=2,n=100		! n- no of points, dims- no of dimensions
integer :: crossings, COUNT, i, j, onbound=0, out!,n
real(mytype) :: t,x,y, yc
real(mytype), dimension(n,dims) :: P 		! POINTS ARRAY	

 OPEN(10,FILE='Cylinder.dat',FORM='FORMATTED',&
        RECL=4, STATUS='OLD')
  !COUNT = 1
     DO i=1,n
           READ(10,*) P(i,1), P(i,2)						! reading points file
         !  COUNT = COUNT + 1
!          print *, P(i,1), P(i,2)
     ENDDO
    
  CLOSE(10)
  
 crossings=0!; out=10        
do i=1,n

    if (i<n) then
	
	if (P(i,1)<x .and. x<P(i+1,1) .or. P(i,1)>x .and. x>P(i+1,1)) then		! points with x,y coor different from those of the vertex points
	  t= (x-P(i+1,1))/(P(i,1)-P(i+1,1))
	  yc= t*P(i,2)+(1-t)*P(i+1,2)
    
	    if (y==yc) then
	      onbound=1 !return (on boundary)
	      print *, 'onbound'
	    else if (y>yc) then
	      crossings=crossings+1
! 	      print *, '1~'
	    else 
	      crossings=crossings
	    endif
	endif
    
	if (P(i,1)==x .and. P(i,2) <=y) then				! for vertex points or on the vertical line
	    if (P(i,2)==y) then					! vertex point
		 onbound=1!return (on boundary)
		 print *, 'onbound'
	    endif
	    if (P(i+1,1)==x) then
		if ((P(i,2) <= y .and. y <= P(i+1,2)) .or. (P(i,2) >= y .and. y >= P(i+1,2))) then 		! on the vertical line joining the vertices
		  onbound=1 !return (on boundary)
		  print *, 'onbound'
		endif
	    else if (P(i+1,1) >x) then 									! the point just doesnt graze the surface but cuts the curve
	      crossings=crossings+1
! 	      print *, '2~'
	    else
	      crossings=crossings
	    endif												! the same case but for the right part of the curve
	    if (i==1) then 	
	      if (P(n,1)>x) then										! assuming the points are in the order of their position on the curve 
		crossings=crossings+1
! 	    print *, '3~'
	      endif
	    else
	      if (P(i-1,1)>x) then										
		crossings=crossings+1
! 		print *, '4~'
	      endif
	    endif
	endif       
     if(((y-P(i,2))/(x-P(i,1))-(P(i+1,2)-P(i,2))/(P(i+1,1)-P(i,1)))==0) then
     onbound=1 ; out=10
     endif
     
   else										! replacing i+1=1 for i=n
	if (P(i,1)<x .and. x<P(1,1) .or. P(i,1)>x .and. x>P(1,1)) then		! points with x,y coor different from those of the vertex points
	  t= (x-P(1,1))/(P(i,1)-P(1,1))
	  yc= t*P(i,2)+(1-t)*P(1,2)
    
	    if (y==yc) then
	      onbound=1 !return (on boundary)
	      print *, 'onbound'
	    else if (y>yc) then
	      crossings=crossings+1
! 	      print *, '1'
	    else
	      crossings=crossings
	    endif
      
	endif    
! 	
	if (P(i,1)==x .and. P(i,2) <=y) then				! for vertex points or on the vertical line
	    if (P(i,2)==y) then					! vertex point
		 onbound=1!return (on boundary)
		 print *, 'onbound'
! 		 print *,'hello'
	    endif
	    if (P(i+1,1)==x) then
		if ((P(i,2) <= y .and. y <= P(1,2)) .or. (P(i,2) >= y .and. y >= P(1,2))) then 		! on the vertical line joining the vertices
		  onbound=1 !return (on boundary)
		  print *, 'onbound'
		endif
	    else if (P(1,1) >x) then 									! the point just doesnt graze the surface but cuts the curve
	      crossings=crossings+1
! 	      print *, '2'
	    else
	      crossings=crossings
	    endif												! the same case but for the right part of the curve
	    if (i==1) then 	
	      if (P(n,1)>x) then										! assuming the points are in the order of their position on the curve 
		crossings=crossings+1
! 	    print *, '3'
	      endif
	    else
	      if (P(i-1,1)>x) then										
		crossings=crossings+1
! 		print *, 'd4'
	      endif
	    endif
	endif       
!       
         if(((y-P(i,2))/(x-P(i,1))-(P(1,2)-P(i,2))/(P(1,1)-P(i,1)))==0) then
     onbound=1 ; out=10
     endif     
!       
    endif
    
enddo

if (crossings>1) then
 if (mod(crossings,2)==1) then
     !return (inside)
    out=0
 else
     !return (outside)
    out=1
 endif
elseif (crossings==0) then
out=1

else
out=0
endif

 if (onbound==1) then
 out=10
 endif
 
! if (out.eq.0) then 
!  counter=counter+1
! endif
 
!  print *, "out= ,onbound= ", out, onbound, crossings
 
  end subroutine polygon_test

! subroutine polygon_test(x,y,out)
! 
! 
! IMPLICIT NONE
! 
! integer, parameter :: n=140, dims=2		! n- no of points, dims- no of dimensions
! integer :: crossings, COUNT, i, j, onbound, out, ep
! real(4) :: t, x, y , yc
! real(4), dimension(n,dims) :: P 		! POINTS ARRAY	
! 
!  OPEN(10,FILE='NACA0012.dat',FORM='FORMATTED',&
!         RECL=4, STATUS='OLD')
!   
!      DO i=1,n
!            READ(10,*) P(i,1), P(i,2)						! reading points file
!      ENDDO
!     
!   CLOSE(10)
! 
!   ! assuming the points are in a order 
!         
!  crossings=0
!  onbound=0
!  out=0
!  ep=10
!         
! do i=1,n
! 
!     if (i<n) then
! 	
! 	if (P(i,1)<x .and. x<P(i+1,1) .or. P(i,1)>x .and. x>P(i+1,1)) then		! points with x,y coor different from those of the vertex points
! 	  t= (x-P(i+1,1))/(P(i,1)-P(i+1,1))
! 	  yc= t*P(i,2)+(1-t)*P(i+1,2)
!     
! 	    if (y==yc) then
! 	      onbound=1 !return (on boundary)
! 	    else if (y>yc) then
! 	      crossings=crossings+1
! 	     ! print *, '1~'
! 	    else 
! 	      crossings=crossings
! 	    endif
! 	endif
!     
! 	if (P(i,1)==x .and. P(i,2) <=y) then				! for vertex points or on the vertical line
! 	    if (P(i,2)==y) then					! vertex point
! 		 onbound=1!return (on boundary)
! 	    endif
! 	    if (P(i+1,1)==x) then
! 		if ((P(i,2) <= y .and. y <= P(i+1,2)) .or. (P(i,2) >= y .and. y >= P(i+1,2))) then 		! on the vertical line joining the vertices
! 		  onbound=1 !return (on boundary)
! 		endif
! 	    else if (P(i+1,1) >x) then 									! the point just doesnt graze the surface but cuts the curve
! 	      crossings=crossings+1
! 	  
! 	    else
! 	      crossings=crossings
! 	    endif												! the same case but for the right part of the curve
! 	    if (i==1) then 	
! 	      if (P(n,1)>x) then										! assuming the points are in the order of their position on the curve 
! 		crossings=crossings+1
! 	  
! 	      endif
! 	    else
! 	      if (P(i-1,1)>x) then										
! 		crossings=crossings+1
! 		
! 	      endif
! 	    endif
! 	endif       
!      
!    else										! replacing i+1=1 for i=n
! 	if (P(i,1)<x .and. x<P(1,1) .or. P(i,1)>x .and. x>P(1,1)) then		! points with x,y coor different from those of the vertex points
! 	  t= (x-P(1,1))/(P(i,1)-P(1,1))
! 	  yc= t*P(i,2)+(1-t)*P(1,2)
!     
! 	    if (y==yc) then
! 	      onbound=1 !return (on boundary)
! 	    else if (y>yc) then
! 	      crossings=crossings+1
! 	
! 	    else
! 	      crossings=crossings
! 	    endif
!       
! 	endif    
! ! 	
! 	if (P(i,1)==x .and. P(i,2) <=y) then				! for vertex points or on the vertical line
! 	    if (P(i,2)==y) then					! vertex point
! 		 onbound=1!return (on boundary)
! 	
! 	    endif
! 	    if (P(i+1,1)==x) then
! 		if ((P(i,2) <= y .and. y <= P(1,2)) .or. (P(i,2) >= y .and. y >= P(1,2))) then 		! on the vertical line joining the vertices
! 		  onbound=1 !return (on boundary)
! 		endif
! 	    else if (P(1,1) >x) then 									! the point just doesnt graze the surface but cuts the curve
! 	      crossings=crossings+1
! 	
! 	    else
! 	      crossings=crossings
! 	    endif												! the same case but for the right part of the curve
! 	    if (i==1) then 	
! 	      if (P(n,1)>x) then										! assuming the points are in the order of their position on the curve 
! 		crossings=crossings+1
! 	
! 	      endif
! 	    else
! 	      if (P(i-1,1)>x) then										
! 		crossings=crossings+1
! 	
! 	      endif
! 	    endif
! 	endif       
! !       
! !         
! !       
!     Endif
!     
! enddo
! 
! if (crossings>1) then
!  if (mod(crossings,2)==1) then
!     out=0 !return (inside)
!     ep=1
!  else
!     out=1 !return (outside)
!     ep=0
!  endif
! elseif (crossings==0) then
! ep=0
! out=1
! else
! ep=1
! endif
! 
!  if (onbound==1) then
!  ep=1
!  endif
!  
! 
!  
!  end subroutine polygon_test
