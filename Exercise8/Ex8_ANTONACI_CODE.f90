MODULE ex8

IMPLICIT NONE 

integer ::i,j,k,h

contains

function a_1(d,N,coeff) result(p_s_n_i)

integer*4,intent(in) :: d,n
complex,allocatable,dimension(:) :: p_s_n_i
complex,dimension(:),intent(in) :: coeff
double precision :: norm=0.0


allocate(p_s_n_i(1:N*d))

p_s_n_i=coeff


do j=1,n-1

	norm=0.0
	norm=dot_product(p_s_n_i(1+(j-1)*d:j*d), p_s_n_i(1+(j-1)*d:j*d) )
	
	p_s_n_i(1+(j-1)*d:j*d)=p_s_n_i(1+(j-1)*d:j*d)*(1.0/sqrt(norm) )
	
	
end do	
		
WRITE(*,1)"Total wavefunction for ",n," non interacting particles is : "
1 FORMAT (A,I2,A)


do i=1,n
	
	WRITE(*,*) "("
	do j=1,d


	WRITE(*,2)"(",real(p_s_n_i( (i-1)*d +j) ),"+",aimag(p_s_n_i( (i-1)*d +j) ),"*i ) |",j-1,">_particle",i,"+"
	
	2 FORMAT(1X,A1,1X,F10.7,A1,F10.7,A6,I2,A10,I2,A1)
	
	end do
	
	WRITE(*,*) ")"
	
	if (i<=n-1) then 
		WRITE(*,*) "(X)"
	
	end if
	
end do

WRITE (*,*) " "


end function a_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function a_2(d,N,coeff,bool) result(p_s_i)

logical :: bool
integer*4,intent(in) :: d,n
complex,allocatable,dimension(:) :: p_s_i
complex,dimension(:),intent(in) :: coeff
double precision :: norm=0.0


allocate( p_s_i(1:d**N) )

p_s_i=coeff


norm=dot_product(p_s_i,p_s_i)

p_s_i=p_s_i/sqrt(norm)

!print*,1/sqrt(norm)

if (bool .eqv. .TRUE.) then
	WRITE(*,1)"Total wavefunction for ",n," interacting particles is : "
	1 FORMAT (A,I2,A)

	WRITE (*,*) " "

		
		
	do j=1,d**n


		WRITE(*,2)"(",real(p_s_i(j) ),"+",aimag(p_s_i(j) ),"*i )*|",j-1,">"
		
		2 FORMAT(1X,A1,1X,F10.7,A1,F10.7,A6,B2.2,A1)
		
		if (j<=d**n-1) then
		
		WRITE(*,*)"+"
		 
		end if
		
		
	end do

	WRITE (*,*) " "
	
end if	




end function a_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



function b (d,n,psi) result(rho)

complex,dimension(:),intent(in) :: psi
integer*4,intent(in) :: d,n
complex,allocatable,dimension(:,:) :: rho

allocate ( rho(1:d**n,1:d**n) )

rho=(0.0,0.0)

do i=1,d**n
	rho(i,:)=psi(i)*psi

end do

WRITE(*,*)" The density matrix is :"

WRITE (*,*) " "

do i=1,d**n

	WRITE(*,*) "|",rho(i,:),"|"
	

end do

WRITE (*,*) " "


end function b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine c (d,rho)


integer*4,intent(in) :: d
complex,dimension(:,:),intent(in) :: rho
complex,allocatable,dimension(:,:) :: rho_1,rho_2



allocate( rho_1 (1:d,1:d ),rho_2( 1:d,1:d ) )

rho_1=(0.0,0.0)
rho_2=(0.0,0.0)

do k=1,d
	do h=1,d
		do i=1,d
				
				rho_1(k,h)=rho_1(k,h)+rho( i+(k-1)*d, i+(h-1)*d )
				
		end do
	
	end do
	
end do 

do k=1,d
	do h=1,d
		do i=1,d
				
				rho_2(k,h)=rho_2(k,h)+rho( i+(i-1)*(d-1)+(k-1), i+(i-1)*(d-1)+(h-1) )
				
		end do
	
	end do
	
end do


WRITE (*,*) " The reduced density matrix for subsystem 1 is : "
WRITE (*,*) " "
do i=1,d

	WRITE(*,*)"|",rho_1(i,:),"|"
	

end do

WRITE (*,*) " "

WRITE (*,*) " The reduced density matrix for subsystem 2 is : "
WRITE (*,*) " "
do i=1,d

	WRITE(*,*)"|",rho_2(i,:),"|"
	

end do

WRITE (*,*) " "


deallocate(rho_1,rho_2)


end subroutine c


end module ex8




PROGRAM main

Use ex8

IMPLICIT NONE

logical :: bool=.FALSE.

integer*4,parameter :: n=2,d=2

complex,dimension(d*n) :: n_i_states=(/(1.0,0.0),(1.0,0.0),(1.0,0.0),(0.0,0.0)/),p_s_n_i=(0.0,0.0)

complex,dimension(d**n) :: i_states=(/(1.0,0.0),(0.0,0.0),(1.0,0.0),(0.0,0.0)/),p_s_i=(0.0,0.0)

complex,dimension(1:d**n,1:d**n) :: rho=(0.0,0.0)


if (bool .eqv. .FALSE.) then

p_s_n_i= a_1(d,n,n_i_states)

end if


p_s_i=a_2(d,n,i_states,bool)


rho=b(d,n,p_s_i)

call c(d,rho)



END PROGRAM main










