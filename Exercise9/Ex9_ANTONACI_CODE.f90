module int_to_string

implicit none

contains

character(len=20) function str(k)

!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
    
end function str

end module int_to_string


program ex9

USE int_to_string

implicit none

double precision,parameter :: dlambda=0.1
integer,parameter ::N_max=9,sudd=3.0/dlambda,level=4

complex,dimension(2,2),parameter :: sigma_z=reshape( (/(1.0,0.0),(0.0,0.0),(0.0,0.0),(-1.0,0.0)/) ,shape=(/2,2/) )
complex,dimension(2,2),parameter :: sigma_x= reshape( (/ (0.0,0.0),(1.0,0.0),(1.0,0.0),(0.0,0.0)/) , shape=(/2,2/) )


integer :: i,j,k,l,n,q
double precision :: finish,start


double precision,dimension(sudd+1) :: lambda=(/(dlambda*i,i=0,sudd)/)
double precision,dimension(sudd+1,N_max,level) :: lambda_vs_energy=0.0




complex,dimension(:,:),allocatable :: H,H_1,H_2,A,B

complex,allocatable,dimension(:,:) :: I_left,I_right,aux

!!!!!!!!!!!!!!!!!!!!!!for diagonalization procedure !!!!!!!!!!!!!!!!

integer :: lwork,info,lda
real,allocatable,dimension(:) :: w,rwork
complex,allocatable,dimension(:) :: work

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



Do k=1,level


lambda_vs_energy(:,1,k)=lambda

end do


Do n=2,N_max
	
call cpu_time(start)

	Do l=1,size(lambda,1)
	
	
		allocate( H( 2**(n),2**(n) ),H_1 ( 2**(n),2**(n) ),H_2 ( 2**(n),2**(n) ) )
	
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
		
		!!!!!!!!!!!!!!!!!!!!!!!!!finish diagonalization procedure!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		Do k=1,level

			lambda_vs_energy(l,n,k)=w(k)/(n-1)

		end do
		
		
		
		
		
		
		deallocate(H,H_1,H_2)
		deallocate(A,B)
		deallocate(work,w,rwork)
		

	
	End do
	
call cpu_time(finish)

WRITE(*,1)"For N = ",N," program takes : " ,finish-start," seconds"
1 FORMAT(A,I2,A,F13.7,A)

End Do


!!Writting on file!!

Do k=1,level

	
	
	OPEN(UNIT=100,FILE='state'//trim(str(k-1))//'dlambda_0_1.txt',FORM='FORMATTED',&
	     ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 


	Do l=1,size(lambda,1)
		
		WRITE(UNIT=100,FMT=*) lambda_vs_energy(l,:,k)
		

	End Do


	 CLOSE(UNIT=100) 

end do


end program ex9





