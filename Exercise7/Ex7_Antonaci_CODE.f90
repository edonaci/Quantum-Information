
MODULE harmonic_oscillator_time_indipendent


USE, INTRINSIC :: iso_c_binding

IMPLICIT NONE

INCLUDE 'fftw3.f03'
INCLUDE 'fftw3l.f03'

double precision,parameter :: eps=0.5E-2,L=6,t_max=5,tau=0.005,h_bar=1.0,m=0.5
integer,parameter :: omega=1,n=int(2*L/eps+1)
integer :: i,j,k

CONTAINS

subroutine ground_state(psi0)

!we have used that 2m=1 and h_tagliato =1 so it means x0=sqrt(2/omega) ad V=(1/4)*(x*omega)**2


double precision,allocatable,dimension(:,:) :: H,ID               ! Using double precision type because potential,eigenvalues and eigenvectors are not complex!
double precision,allocatable,dimension(:) :: w,work
double precision,dimension(n) :: psi0


double precision :: xo,start,finish

integer::lwork,lda,info

integer :: a


xo=sqrt( 2.0/real(omega) )




allocate( H(1:n,1:n),ID(1:n,1:n) )


print*,"hermitian matrix's size is: ",size(H,1)


if (L<=xo) then

print*,"System maximum size is less then typical lenght!"

stop

end if



call cpu_time(start)

FORALL (i=1:n ,j=1:n ,i==j)
	
	
	H(i,j)=+2/(eps**2)+ 0.5*m*(omega*(-L+eps*(i-1) ) )**2
		
END FORALL
	
	


FORALL (i=1:n ,j=1:n ,(j==i+1 .OR. j==i-1) ) 

	H(i,j)=-1/(eps**2)
	H(j,i)=-1/(eps**2)
	
END FORALL

	
!print*,H(1:5,1:5)                             to check i'm making a good work

call cpu_time(finish)
	
print*,"In order to define correct hermitian elements ,the code takes ",finish-start,"	seconds"	




lda=n


lwork=3*lda-1
allocate(work(1),w(1:lda) )


lwork=-1				   
call dsyev("V","L",n,H,lda,w,work,lwork,info)

lwork = int(real(work(1)))

deallocate(work)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(	work(max(1,lwork)	)	)

call cpu_time(start )

call dsyev("V","L",n,H,lda,w,work,lwork,info)

call cpu_time(finish)


print*, "Dsyev takes ",finish-start, "	seconds"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (info/=0) then

print*, "Something is gone bad during eigenvalues and eigenvectors extration !"

stop

endif


ID=matmul( transpose(H),H) 

!print*, ID(1:5,1:5)						to check if we obtain an identity matrix

DO i=1,n
	Do j=1,n
		if (i==j .AND. abs(ID(j,i)-1.0)>1E-8) then
		
			print*, "Error in normalization procedure !"
			
			stop
		
		end if
		
		if (i/=j .AND. abs(ID(j,i)-0.0)>1E-8) then
			
			  print*, "Error in normalization procedure !"
			  
			  stop
		
		end if
		
	END DO

END DO


psi0=H(:,1)


END SUBROUTINE ground_state


END MODULE harmonic_oscillator_time_indipendent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main


USE harmonic_oscillator_time_indipendent

IMPLICIT NONE


double precision,dimension(n) :: psi0
double precision,dimension(int(t_max/tau)+1 ):: t_evolution=(/ (j*tau ,j=0,int(t_max/tau) ) /)
double precision,dimension(n) :: space=(/ ( eps*(k-1)-L ,k=1,n ) /)

double precision,dimension(n,int(t_max/tau)+1) :: wf_evolution_norm,wf_evolution_real,wf_evolution_im

double precision :: x_sh,norm,start,finish
complex*16,dimension(n) :: psi,tr_psi,antr_psi

complex*16 :: factor_pot=complex(0.0,-0.5*tau*0.5*m*omega*omega/h_bar)

complex*16 :: factor_mom= complex(0, -tau*h_bar*(2*acos(-1.0)/(2*L) )**2 )
type(C_PTR) :: planf,planb


call ground_state(psi0)                                  !!!!groud state

psi=psi0                                        !!!become complex

call cpu_time(start)

Do i=1,size(t_evolution)

	norm=0.0
	
	do j=1,n
	
		x_sh =space(j)-t_evolution(i)                   
		
		psi(j)=psi(j)*cdexp(factor_pot*(x_sh**2) )     !!exp ( (-i/h_bar)*tau*V(x,t)/2)
		
	end do
		
			call dfftw_plan_dft_1d(planf,n,psi,tr_psi,FFTW_FORWARD,FFTW_ESTIMATE);          !make plan fourier forward
			call dfftw_execute_dft(planf,psi,tr_psi)					   !fourier forward
			call dfftw_destroy_plan(planf)						   !destroy plan
														
	do j=1,int(n/2d0)									!kinetic part
	
		tr_psi(j)=tr_psi(j)*cdexp(factor_mom*(j)**2 )                            !from 0 to +L
		
	end do
	
	  do j=int(N/2d0+1d0),N
	  
                  tr_psi(j)=tr_psi(j)*cdexp(factor_mom*(j-N)**2 )                        !from -L to 0
                  
                  
                  
          end do
		
			call dfftw_plan_dft_1d(planb,n,tr_psi,psi,FFTW_BACKWARD,FFTW_ESTIMATE);          !make plan fourier backward
			call dfftw_execute_dft(planb, tr_psi, psi )					    !fourier backward
	 		call dfftw_destroy_plan(planb)						    !destroy plan
	 do j=1,n
	
		x_sh =space(j)-t_evolution(i)                   
		
		psi(j)=psi(j)*cdexp(factor_pot*(x_sh**2) )     !! exp ( (-i/h_bar)*tau*V(x,t)/2)
		
	end do	
			

			norm=dot_product(psi,psi)*eps
			
			psi=psi/sqrt(norm)
			
			norm=dot_product( psi,psi)*eps                                        !In order to be sure that psi will have norm equal 1
		
			if ( dabs( norm -1.0) >1E-8 ) then
			
			print*,"At",t_evolution(i)," seconds the wavefunction is not normalized!Norm =", norm
			
			
			stop
			
			end if
		
	wf_evolution_real(:,i)=real(psi)
	
	wf_evolution_im(:,i)=aimag(psi)			
			
			FORALL (j=1:n)
			
			psi(j)=psi(j)*conjg( psi(j) )
			
			
			END FORALL
			
			
			
	wf_evolution_norm(:,i)=psi


End Do

call cpu_time(finish)


print*,"Evolution wave function takes",finish-start, "seconds"


OPEN(UNIT=1,FILE='times.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 

DO i=1,size(t_evolution)

	
	WRITE(UNIT=1,FMT=*) t_evolution(i)

	
	

END DO
 
 CLOSE(UNIT=1)
 

OPEN(UNIT=2,FILE='hotd_norm.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 


Do k=1,size(space)
	
	WRITE(UNIT=2,FMT=*) space(k),wf_evolution_norm(k,:)
	

End Do


 CLOSE(UNIT=2) 
 
 
 OPEN(UNIT=3,FILE='hotd_real.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 


Do k=1,size(space)
	
	WRITE(UNIT=3,FMT=*) space(k),wf_evolution_real(k,:)
	

End Do


 CLOSE(UNIT=3)
  
 
OPEN(UNIT=4,FILE='hotd_im.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 


Do k=1,size(space)
	
	WRITE(UNIT=4,FMT=*) space(k),wf_evolution_im(k,:)
	

End Do


 CLOSE(UNIT=4)  



END PROGRAM main







