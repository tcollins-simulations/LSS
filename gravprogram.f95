program gravity
implicit none 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
type particle 
	real :: m 
	real :: x,y,dx,dy 
	real :: ax,ay,vx,vy,vix,viy,fx,fy 
end type particle
type grav 
	real :: fg, theta, fgx, fgy, r
end type grav
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
integer, parameter :: i=10000 !max is 29584, must be a square
type(particle), dimension(i) :: s
type(grav),dimension(i,i) :: f
real,parameter :: G=6.674e-11, pi=3.1415926
integer :: itr,t,input, n, x
real :: totaltime,tstep, frames
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
open(13,file='gravforce.txt') 
open(14,file='gravxy.txt') 
open(15,file='xyframes.txt')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
call xychoose() 
s%m=5e6
do n=1,i
	print *, n, s(n)%x, s(n)%y
end do
call timedef()
t=0
call xytxt()
call excelinitial()
print *, "Running..."
do t=1,itr
	totaltime=t*tstep
	call gravpotential()
	!call gravtxt()
	do n=1,i
		call objectx()
		call objecty()
	end do
	call excel()
	call xytxt()
end do
print *, "Simulation complete"
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine timedef()
	print *, "How long for each time step?"
	read *, tstep
	print *, "How many iterations to run for?"
	read *, itr
end subroutine timedef
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine xychoose()
	print *, "How would you like to input initial input data?"
	print *, "Enter 1 for direct input " 
	print *, "2 for importing data"
	print *, "3 for a grid defined by the step between particles"
	print *, "4 for a grid defined by area"
	read *, input
	if (input==1) then
		call inputcoordinates()
	else if (input==2) then
		call importcoordinates()
	else if (input==3) then
		call stepgrid()
	else if (input==4) then
		call areagrid()
	else 
		do while (input/=1 .and. input/=2 .and. input/=3 .and. input/=4)
			print *, "Invalid input. Try again:"
			print *, "Enter 1 for direct input " 
			print *, "2 for importing data"
			print *, "3 for a grid defined by the step between particles"
			print *, "4 for a grid defined by area"
			read *, input
			if (input==1) then
				call inputcoordinates()
			else if (input==2) then
				call importcoordinates()
			else if (input==3) then
				call stepgrid()
			else if (input==4) then
				call areagrid()
			else 
				continue
			end if
		end do
	end if
end subroutine xychoose
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine objectx()
	s(n)%ax=s(n)%fx/s(n)%m
	s(n)%vx=s(n)%vix+s(n)%ax*tstep
	call vxcap()
	s(n)%dx=s(n)%vx*tstep
	s(n)%x=s(n)%x+s(n)%dx
	s(n)%vix=s(n)%vx
end subroutine objectx
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine objecty()
	s(n)%ay=s(n)%fy/s(n)%m
	s(n)%vy=s(n)%viy+s(n)%ay*tstep
	call vycap()
	s(n)%dy=s(n)%vy*tstep
	s(n)%y=s(n)%y+s(n)%dy
	s(n)%viy=s(n)%vy
end subroutine objecty
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine vxcap()
	if (s(n)%vx>0.1E-01) then
		s(n)%vx=.1E-01
	else if (s(n)%vx<-0.1E-01) then
		s(n)%vx=-.1E-01
	else
		continue
	end if
end subroutine vxcap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine vycap()
	if (s(n)%vy>0.1E-01) then
		s(n)%vy=.1E-01
	else if (s(n)%vy<-0.1E-01) then
		s(n)%vy=-.1E-01
	else
		continue
	end if
end subroutine vycap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine xytxt()
	integer :: write, finalstep
	finalstep=i-1
	write (14,5,advance='no') totaltime
	5 format(E10.4,2X)
	6 format(2(E10.4,2x))
	do write=1,finalstep
		write(14,6,advance='no') s(write)%x, s(write)%y
	end do
	write (14,6,advance='yes') s(i)%x, s(i)%y
end subroutine xytxt
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine gravpotential()
	do n=1,i
		do x=1,i
			if (n/=x) then
				f(n,x)%r=SQRT(((s(n)%x-s(x)%x)**2)+((s(n)%y-s(x)%y)**2))
				f(n,x)%theta=ATAN(ABS(s(x)%y-s(n)%y)/ABS(s(x)%x-s(n)%x))
				f(n,x)%fg=2*G*s(n)%m*s(x)%m/f(n,x)%r
				f(n,x)%fgx=f(n,x)%fg*COS(f(n,x)%theta)
				f(n,x)%fgy=f(n,x)%fg*SIN(f(n,x)%theta)
				if (s(n)%x>s(x)%x) then
					f(n,x)%fgx=(-1)*f(n,x)%fgx
				else if (s(n)%x==s(x)%x) then
					f(n,x)%fgx=0
				else 
					continue
				end if
				if (s(n)%y>s(x)%y) then
					f(n,x)%fgy=(-1)*f(n,x)%fgy
				else if (s(n)%y==s(x)%y) then
					f(n,x)%fgy=0
				else 
					continue
				end if
			end if
			if (n==x) then
				f(n,x)%r=0
				f(n,x)%theta=0
				f(n,x)%fg=0
				f(n,x)%fgx=0
				f(n,x)%fgy=0
			end if
		end do
	end do
	s%fx=SUM(f%fgx,DIM=2)
	s%fy=SUM(f%fgy,DIM=2)
end subroutine gravpotential
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine gravtxt()
	do n=1,i
		do x=1,i
		write (13,*) n, x, s(n)%x, s(n)%y, s(x)%x, s(x)%y, f(n,x)%r, f(n,x)%theta, f(n,x)%fg, f(n,x)%fgx, f(n,x)%fgy
		end do
	end do
end subroutine gravtxt
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine excelinitial()
	print *, "How many frames to capture?"
	read *, frames
	do n=1,i
		write (15,*) t,s(n)%x,s(n)%y
	end do
end subroutine excelinitial
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine excel()
	integer :: write, finalstep, timetestint
	real :: endtime, timedivide, timetest
	finalstep=i-1
	endtime=tstep*itr
	timedivide=endtime/frames
	timetest=t/timedivide
	timetestint=timetest
	if (timetestint == timetest) then
		do n=1,i
			write (15,*) t, s(n)%x, s(n)%y
		end do
	else
		continue
	end if
end subroutine excel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine inputcoordinates()
	print *, "Enter the coordinate values for object 1:"
	read *, s(1)%x, s(1)%y
	do n=2,i
		print *,"Object",n,":"
		read *, s(n)%x, s(n)%y
	end do
end subroutine inputcoordinates
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine importcoordinates()
	character (len=20) :: fullname, name
	print *, "Enter name of file to read from"
	read *, fullname
	name=TRIM(fullname)
	open (12,file=name)
	do n=1,i
		read (12,*) s(n)%x, s(n)%y
	end do
end subroutine importcoordinates
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine stepgrid()
	integer :: l, p, k, w, u
	real :: d
	l=AINT(SQRT(REAL(i)))
	print *, "Enter spacing between particles"
	read *, d
	do u=1,l
		do p=1,l
			n=l*(u-1)+p
			s(n)%x=(p-1)*d
		end do
	end do
	do w=1,l
		do k=1,l
			n=l*(w-1)+k
			s(n)%y=(w-1)*d
		end do
	end do
end subroutine stepgrid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine areagrid()
	integer :: l, d, p, k, w, u
	real :: A
	l=AINT(SQRT(REAL(i)))
	print *, "Enter area of analysis:"
	read *, A
	d=(SQRT(A))/(SQRT(REAL(i)))
	do u=1,l
		do p=1,l
			n=l*(u-1)+p
			s(n)%x=(p-1)*d
		end do
	end do
	do w=1,l
		do k=1,l
			n=l*(w-1)+k
			s(n)%y=(w-1)*d
		end do
	end do
end subroutine areagrid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++
end program gravity
