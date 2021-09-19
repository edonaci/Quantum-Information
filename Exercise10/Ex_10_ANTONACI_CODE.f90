module to_string

implicit none

contains

character(len=20) function strI(k)

!   "Convert an integer to string."
    integer, intent(in) :: k
    write (strI, *) k
    strI = adjustl(strI)
    
end function strI

character(len=20) function strR(k)

!   "Convert an integer to string."
    real, intent(in) :: k
    write (strR, *) k
    strR = adjustl(strR)
    
end function strR

end module to_string


program ex10

USE to_string

implicit none

double precision,parameter :: dlambda=0.1
integer,parameter ::N_max=4,sudd=30,int_max=28

complex,dimension(2,2),parameter :: sigma_z=reshape( (/(1.0,0.0),(0.0,0.0),(0.0,0.0),(-1.0,0.0)/) ,shape=(/2,2/) )
complex,dimension(2,2),parameter :: sigma_x= reshape( (/ (0.0,0.0),(1.0,0.0),(1.0,0.0),(0.0,0.0)/) , shape=(/2,2/) )


integer :: i,j,k,l,n,q,iteration


real,dimension(sudd+1) :: lambda=(/(dlambda*i,i=0,sudd)/)
double precision,dimension(sudd+1,2,N_max) :: lambda_vs_energy=0.0
double precision,dimension(int_max+1,2,sudd+1) :: plateau=0.0




complex,dimension(:,:),allocatable :: H,H_start,H_1,H_2,A,B

complex,allocatable,dimension(:,:) :: I_left,I_right,aux

!! for RSRG algorithm

complex,allocatable,dimension(:,:) :: A_2N,B_2N
complex,allocatable,dimension(:,:) :: H_2N,H_2N_start,I_N,I_N_1,I_2N_1,P,P_adjoint

!!!!!!!!!!!!!!!!!!!!!!for diagonalization procedure !!!!!!!!!!!!!!!!

integer :: lwork,info,lda
real,allocatable,dimension(:) :: w,rwork
complex,allocatable,dimension(:) :: work

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






Do n=2,N_max

print*,"particles at the beginning are ",n
	
	lambda_vs_energy(:,1,n-1)=lambda
	
	Do l=1,size(lambda,1)
	
	
		allocate( H( 2**(n),2**(n) ),H_1 ( 2**(n),2**(n) ),H_2 ( 2**(n),2**(n) ) )
		
		allocate( H_start( size(H,1),size(H,1) ) )
	
		allocate( A( 2**(n),2**(n) ),B( 2**(n),2**(n) ) )
		
		
	
		H=(0.0,0.0)
		H_1=(0.0,0.0)
		H_2=(0.0,0.0)
		A=(0.0,0.0)
		B=(0.0,0.0)

		
		!!FOR SIGMA_z!!!
		
		Do k=0,N-1



			allocate( I_left(2**(k),2**(k) ) , I_right(2**(n-1-k),2**(n-1-k) )  )
			allocate(aux( size(I_right,1)*2,size(I_right,1)*2 ) )
			
			I_left=(1.0,0.0)
			I_right=(1.0,0.0)
			
			FORALL ( i=1:size(I_left,1),j=1:size(I_left,1),i/=j )
			
			
			I_left(i,j)=(0.0,0.0)


			END FORALL
			

			FORALL ( i=1:size(I_right,1),j=1:size(I_right,1),i/=j )
			
			
			I_right(i,j)=(0.0,0.0)


			END FORALL

		

			aux=(0.0,0.0)
			
			FORALL ( i=1:2 ,j=1:2 )
			
			
			aux( 1+(i-1)*size(I_right,1) : i*size(I_right,1) ,1+(j-1)*size(I_right,1) : j*size(I_right,1) )=sigma_z(i,j)*I_right
			
			
			
			END FORALL
			
		
			
			FORALL(i=1:size(I_left,1),j=1:size(I_left,1) )
			
			
			 H_1( 1+(i-1)*size(aux,1) : i*size(aux,1) ,1+(j-1)*size(aux,1) : j*size(aux,1) )= &
			 H_1( 1+(i-1)*size(aux,1) : i*size(aux,1) ,1+(j-1)*size(aux,1) : j*size(aux,1) )+I_left(i,j)*aux
			
			
			
			END FORALL	
			
			
			
			deallocate(I_left,I_right,aux)
			
		END DO
		
		!!!!!!first term adding

		H=H+lambda(l)*H_1

		!!!!!!!!!!!!!!!!!!!!!!
		
		!!! For sigma_x*sigma_x+1

		DO k=0,n-2
			DO q=0,1
			
				allocate( I_left(2**(k+q),2**(k+q) ) , I_right(2**(n-1-k-q),2**(n-1-k-q) )  )
				allocate(aux( size(I_right,1)*2,size(I_right,1)*2 ) )
				
				I_left=(1.0,0.0)
				I_right=(1.0,0.0)
				
				FORALL ( i=1:size(I_left,1),j=1:size(I_left,1),i/=j )
				
				
				I_left(i,j)=(0.0,0.0)


				END FORALL
				

				FORALL ( i=1:size(I_right,1),j=1:size(I_right,1),i/=j )
				
				
				I_right(i,j)=(0.0,0.0)


				END FORALL


				aux=(0.0,0.0)
				
				
				FORALL ( i=1:2 ,j=1:2 )
				
				
				aux( 1+(i-1)*size(I_right,1) : i*size(I_right,1) ,1+(j-1)*size(I_right,1) : j*size(I_right,1) )=sigma_x(i,j)*I_right
				
				
				
				END FORALL
				
				
				
				if (q==0) then 
					FORALL(i=1:size(I_left,1),j=1:size(I_left,1) )
				
				
				 	A( 1+(i-1)*size(aux,1) : i*size(aux,1) ,1+(j-1)*size(aux,1) : j*size(aux,1) )=I_left(i,j)*aux
				
				
				
					END FORALL	
				
				else
				
					FORALL(i=1:size(I_left,1),j=1:size(I_left,1) )
				
				
				 	B( 1+(i-1)*size(aux,1) : i*size(aux,1) ,1+(j-1)*size(aux,1) : j*size(aux,1) )=I_left(i,j)*aux
				
				
				
					END FORALL
			
				end if
				
				
				
				deallocate(I_left,I_right,aux)
				
				
			
			END DO
			
			
			H_2=H_2+matmul(B,A)
			
			
			

		END DO
		
		!second term adding	

		H=H+H_2
		
		!!! Now i'm finished with H matrix definition
		
		H_start=H


  		!!!!!!!!!begin diagonalization procedure!!!!!!!
  		
		lda=size(H,1)


		lwork=2*lda-1
		
		allocate(work(1:lwork),w(1:lda),rwork(1:3*lda-2 ) )


		lwork=-1				   
		
		call cheev("N","U",size(H,1),H,lda,w,work,lwork,rwork,info)

		lwork = int(real(work(1)))

		deallocate(work)

	
		
		allocate(	work(max(1,lwork)	)	)
              


		call cheev("N","U",size(H,1),H,lda,w,work,lwork,rwork,info)

	

		if (info/=0) then

		print*, "Something is gone bad during eigenvalues  extration !"

		stop

		endif
!-----------------------------------------------------------------------------------!		
		if (n==2) then
		
		plateau(1,1,l)=0
		plateau(1,2,l)=w(1)
		
		end if
!-----------------------------------------------------------------------------------!		
		
		
		
		
		deallocate(work,w,rwork)
		
		deallocate(H_1,H_2,H)
		deallocate(A,B)
		
!!!!_______________________________________________________________________________________________________________________________________________________!!!		
		
		! i have finished to compute hamiltonian matrix for given n and lambda but i'm still in lambda loop inner n loop
		
		!Now i'm going to implement RSRG algorithm
		
		allocate(I_N_1( size(H,1)/2, size(H,1)/2 ) )
		
		allocate( I_N(size(H,1),size(H,1) ) )
		
		
		allocate(A(size(H,1),size(H,1) ),B(size(H,1),size(H,1) ) )	
		
		allocate(aux(size(H,1),size(H,1) ) )
		
		allocate( H_2N( size(H,1)**2,size(H,1)**2 ),H_2N_start( size(H,1)**2,size(H,1)**2 ) )
		
		allocate( P(size(H_2N,1),size(H,1) ) ,P_adjoint(size(H,1),size(H_2N,1) )  )
		
		allocate(A_2N(size(H_2N,1),size(H_2N,1) ),B_2N( size(H_2N,1),size(H_2N,1) ) )
		
		allocate( I_2N_1(size(H_2N,1)/2,size(H_2N,1)/2 ) )
		

					
		
		I_N=(1.0,0.0)
		I_N_1=(1.0,0.0)
		I_2N_1=(1.0,0.0)
		
		A=(0.0,0.0)
		B=(0.0,0.0)
		aux=(0.0,0.0)
		
		B_2N=(0.0,0.0)
		A_2N=(0.0,0.0)
		

		
		
		FORALL ( i=1:size(I_N,1),j=1:size(I_N,1),i/=j )
			
			
			I_N(i,j)=(0.0,0.0)


		END FORALL
		
		
		
		
		FORALL ( i=1:size(I_N_1,1),j=1:size(I_N_1,1),i/=j )
			
			
			I_N_1(i,j)=(0.0,0.0)


		END FORALL
		
		
		
		FORALL ( i=1:size(I_2N_1,1),j=1:size(I_2N_1,1),i/=j )
			
			
			I_2N_1(i,j)=(0.0,0.0)


		END FORALL
		
	!!!!!______________________________________________________________________________________________________________________________________________!!!!!!	
				
		FORALL ( i=1:size(I_N_1,1),j=1:size(I_N_1,1) )
			
			
			aux( 1+(i-1)*size(sigma_x,1) : i*size(sigma_x,1) ,1+(j-1)*size(sigma_x,1) : j*size(sigma_x,1) ) = &
			I_N_1(i,j)*sigma_x
			
		END FORALL
		
		
		FORALL( i=1:size(aux,1),j=1:size(aux,1) )
		
			A_2N( 1+(i-1)*size(I_N,1) : i*size(I_N,1) ,1+(j-1)*size(I_N,1) : j*size(I_N,1) ) = &
			aux(i,j)*I_N
			
		
		END FORALL
		
		
		aux=(0.0,0.0) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
		FORALL( i=1:2,j=1:2 )
		
			aux( 1+(i-1)*size(I_N_1,1) : i*size(I_N_1,1) ,1+(j-1)*size(I_N_1,1) : j*size(I_N_1,1) ) = &
			sigma_x(i,j)*I_N_1
		
		END FORALL
		
		
		FORALL ( i=1:size(I_N,1),j=1:size(I_N,1) )
			
			B_2N( 1+(i-1)*size(I_N,1) : i*size(I_N,1) ,1+(j-1)*size(I_N,1) : j*size(I_N,1) ) = &
			I_N(i,j)*aux
			
			
		END FORALL
		
		
		
		
		Do iteration=1,int_max
		
			H_2N=(0.0,0.0)
			H_2N_start=(0.0,0.0)
			
			
			!1° step algorithm : double system
				
				
			FORALL ( i=1:size(I_N,1) ,j=1:size(I_N,1) )
					
					
				H_2N( 1+(i-1)*size(H_start,1) : i*size(H_start,1) ,1+(j-1)*size(H_start,1) : j*size(H_start,1) )= & 
				H_2N( 1+(i-1)*size(H_start,1) : i*size(H_start,1) ,1+(j-1)*size(H_start,1) : j*size(H_start,1) ) + I_N(i,j)*H_start                  !1° tensor product
				
			END FORALL
			
			FORALL ( i=1:size(H_start,1) ,j=1:size(H_start,1) )
				
				H_2N( 1+(i-1)*size(I_N,1) : i*size(I_N,1) ,1+(j-1)*size(I_N,1) : j*size(I_N,1) )= &
				H_2N( 1+(i-1)*size(I_N,1) : i*size(I_N,1) ,1+(j-1)*size(I_N,1) : j*size(I_N,1) ) + H_start(i,j)*I_N            !2° tensor product
			
			END FORALL
			
			
			
			if (iteration==1) then
			
				FORALL ( i=1:size(I_N_1,1),j=1:size(I_N_1,1) )
					
					
					A( 1+(i-1)*size(sigma_x,1) : i*size(sigma_x,1) ,1+(j-1)*size(sigma_x,1) : j*size(sigma_x,1) ) = &
					I_N_1(i,j)*sigma_x
					
				END FORALL
				
				FORALL ( i=1:2,j=1:2 )
				
					B( 1+(i-1)*size(I_N_1,1) : i*size(I_N_1,1) ,1+(j-1)*size(I_N_1,1) : j*size(I_N_1,1) ) = &
					sigma_x(i,j)*I_N_1
					
					
				END FORALL
	
			
			
			else
	
				
				A=matmul(P_adjoint,matmul(A_2N,P) )
				B=matmul(P_adjoint,matmul(B_2N,P) )
			
			end if
			
			
			FORALL (i=1:size(H_start,1),j=1:size(H_start,1) )
			
			
					H_2N( 1+(i-1)*size(B,1) : i*size(B,1) ,1+(j-1)*size(B,1) : j*size(B,1) )= &
					H_2N( 1+(i-1)*size(B,1) : i*size(B,1) ,1+(j-1)*size(B,1) : j*size(B,1) ) + A(i,j)*B              !3° tensor product
				
				
						
			END FORALL
			
			
			H_2N_start=H_2N
			
			! 2° step algorithm : diagonalize matrix and take first eigenvalue
			
			lda=size(H_2N,1)


			lwork=2*lda-1
			
			allocate(work(1:lwork),w(1:lda),rwork(1:3*lda-2 ) )


			lwork=-1				   
			
			call cheev("V","U",size(H_2N,1),H_2N,lda,w,work,lwork,rwork,info)

			lwork = int(real(work(1)))

			deallocate(work)

		
			
			allocate(	work(max(1,lwork)	)	)
		      


			call cheev("V","U",size(H_2N,1),H_2N,lda,w,work,lwork,rwork,info)

		

			if (info/=0) then

			print*, "Something is gone bad during eigenvalues and eigenvectors extration !"

			stop

			endif
			
			!------------------------------------------------------------------------------------------!
				
		       if (n==2) then
			
			plateau(iteration+1,1,l)=iteration 
			plateau(iteration+1,2,l)=w(1)/(n*2**(iteration)-1 ) 
				
			end if
			
			!--------------------------------------------------------------------------------------------!
			
			lambda_vs_energy(l,2,n-1)=w(1)/(n*2**(iteration)-1 ) 
			
			!3° projected the system
			
			
			
			P=H_2N(1:size(H_2N,1),1:size(H,1) )
			
			
			P_adjoint=conjg(transpose(P) )
			
						
			
			H_start=matmul(P_adjoint,matmul(H_2N_start,P) )
			
			
			
			deallocate(work,w,rwork)
			
			
		END DO  !iteration
		
		
			
	deallocate (I_N,I_N_1,I_2N_1)
			
	deallocate(A,A_2N,B,B_2N)
			
	deallocate(H_2N,H_2N_start)
		
	deallocate (aux)
		
	deallocate (P,P_adjoint)
		
	deallocate(H_start)
	
	
	END DO	!lambda
	
END DO !particles

!Writting on file!!

Do l=1,size(lambda,1)

	
	
	OPEN(UNIT=100,FILE='plateau_lambda_'//trim(strR( lambda(l) ))//'.txt',FORM='FORMATTED',&
	     ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 


	Do k=1,size(plateau,1)
		
	WRITE(UNIT=100,FMT=*) plateau(k,1,l),plateau(k,2,l)
		

	End Do


	 CLOSE(UNIT=100) 

end do

Do n=1,N_max

	
	
	OPEN(UNIT=100,FILE='lambda_vs_energy_'//trim(strI(n+1) )//'_particles.txt',FORM='FORMATTED',&
	     ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 


	Do l=1,size(lambda,1)
		
	WRITE(UNIT=100,FMT=*) lambda_vs_energy(l,1,n),lambda_vs_energy(l,2,n)
		

	End Do


	 CLOSE(UNIT=100) 

end do
		

end program ex10









