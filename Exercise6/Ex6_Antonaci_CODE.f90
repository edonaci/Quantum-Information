
MODULE armonic_oscillator

IMPLICIT NONE

CONTAINS

function matr(omega,L,eps) result (H)

!we have used that 2m=1 and h_tagliato =1 so it means x0=sqrt(2/omega) ad V=(1/4)*(x*omega)**2

integer,intent(in) :: omega
double precision,intent(in) :: eps,L

double precision,allocatable,dimension(:,:) :: H                ! In doppia precisione perchè il potenziale,gli autovalori e gli autovettori are not complex!

double precision :: xo,start,finish


integer :: i,j,k,n


xo=sqrt( 2.0/real(omega) )


n=2*L+1



if (eps*L<=xo) then

print*,"System maximum size is less then typical lenght!"

stop

end if



allocate(  H(1:n,1:n) )

call cpu_time(start)

FORALL (i=1:n ,j=1:n ,i==j)
	
	
	H(i,j)=+2/(eps**2)+ 0.25*(omega*eps*(i-L-1 ) )**2
		
END FORALL
	
	


FORALL (i=1:n ,j=1:n ,(j==i+1 .OR. j==i-1) ) 

	H(i,j)=-1/(eps**2)
	H(j,i)=-1/(eps**2)
	
END FORALL

	
!print*,H(1:5,1:5)                             to check i'm making a good work

call cpu_time(finish)
	
print*,"In order to define correct hermitian elements ,the code takes ",finish-start,"	seconds"	


end function matr


END MODULE armonic_oscillator


PROGRAM main


USE armonic_oscillator

implicit none


INTEGER::OMEGA=1E3
DOUBLE PRECISION :: eps=1.0E-4,L=1000,inizio,fine
		
				   
double precision,allocatable,dimension(:,:) ::H,ID
double precision,allocatable,dimension(:) :: w,work



integer::lwork,lda,info

integer :: a,n

integer :: i,j,k


n=2*L+1



allocate( H(1:n,1:n),ID(1:n,1:n) )

H=matr(OMEGA,L,eps)


n=size(H,1)

lda=n

print*,"hermitian matrix's size is: ",size(H,1)

lwork=3*lda-1
allocate(work(1),w(1:lda) )


lwork=-1				   
call dsyev("V","L",n,H,lda,w,work,lwork,info)

lwork = int(real(work(1)))

deallocate(work)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(	work(max(1,lwork)	)	)

call cpu_time(inizio)

call dsyev("V","L",n,H,lda,w,work,lwork,info)

call cpu_time(fine)

w=w/OMEGA

print*, "Dsyev takes ",fine-inizio, "	seconds"

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

print*, "i'm beginnig to write into files :"

call cpu_time(inizio)

OPEN(UNIT=16,FILE='first_three_eigenvector_1E-4.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 

DO i=1,size(H,1)

	
	WRITE(UNIT=16,FMT=*) eps*(i-L-1),H(i,1),H(i,2),H(i,3)
	


END DO

 CLOSE(UNIT=16) 

OPEN(UNIT=1,FILE='eigenvalues_psi_1E-4.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 

DO i=1,size(H,1)

	
	WRITE(UNIT=1,FMT=*) i-1,w(i),( (i-1)+0.5 )

	
	

END DO
 
 CLOSE(UNIT=1)
 
 call cpu_time(fine)
 
 print*,"writing on files takes ", fine-inizio,"	seconds"



















END PROGRAM main

