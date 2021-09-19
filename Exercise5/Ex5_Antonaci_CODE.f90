!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EX1
MODULE Ex1

IMPLICIT NONE

CONTAINS

 FUNCTION punto1 (A) result(W)
        
        
        COMPLEX*16,ALLOCATABLE,DIMENSION(:,:),INTENT(IN) :: A
        
        
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION (:) :: W
        
        COMPLEX*16,ALLOCATABLE,DIMENSION (:) :: WORK
        
        INTEGER :: LWORK
        
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: RWORK
        
        INTEGER ::info
        
        INTEGER :: i,j,k,m
        


        m=size(A,1)
        
        LWORK=2*m-1
        
        allocate( W(1:m),WORK(1:LWORK),RWORK(1:3*m-2) )




	CALL zheev('N','L',m,A,m,W,WORK,LWORK,RWORK,info)
	
	if (info/=0) then
	
	WRITE(*,*) "Something is gone bad during diagonalization procedure! "
	
	
	
	end if
	
	

 END FUNCTION punto1
 
 
 FUNCTION punto2(W,suddivisioni) result(S)
 
 IMPLICIT NONE
 

 INTEGER,INTENT(IN) :: suddivisioni
 DOUBLE PRECISION,DIMENSION(:),INTENT(IN) :: W
 
 INTEGER ::j
 
INTEGER,DIMENSION( suddivisioni ) :: levels
 
 DOUBLE PRECISION,DIMENSION(size(W,1)-1,suddivisioni) :: S,R
  
 
 DOUBLE PRECISION::media
 
 INTEGER ::i,k,m,lowest,upest,l_bound,u_bound
 
 levels=(/ ( (size(W)/suddivisioni)*j ,j=1,suddivisioni,1) /)
 
 DO j=1,size(S,2)
 	DO k=1,size(S,1)
 	
 	S(k,j)=W(k+1)-W(k)
 	
 	
 	END DO
 	
 	  m=levels(j)
 		DO i=1,size(S,1)
		 	lowest= i-m
		 	upest=i+m
		 	
		 	if ( (i-m)>0 .AND. i+m<(size(S,1)+1) ) then
		 	   
		 	
		 	 	media=SUM(S(i-m:i+m,j))/size( S(i-m:i+m,1) )
		 	 	
		 	 	S(i-m:i+m,j)=S(i-m:i+m,j)/media
		 	 	
		 	 endif
		 	 
		 	 if ( (lowest <0) .AND. (upest+abs(lowest)<size(S,1)+1) ) then
		 	       
		 	       l_bound=1
		 	       
		 	       u_bound=upest+abs(lowest)
		 	       
		 	       media=SUM( S(l_bound:u_bound,j)  )/size( S(l_bound:u_bound,1 )  )
		 	 	
		 	 	S(l_bound:u_bound,j)=S(l_bound:u_bound,j)/media
		 	
		 	 endif
		 	 
		 	 if ( ( (lowest-abs( size(S,1)-upest ) ) >= 0) .AND. (upest > size(S,1) ) ) then
		 	       
		 	       l_bound=lowest - abs( size(S,1)-upest ) 
		 	       u_bound=size(S,1)
		 	       
		 	       media=SUM( S(l_bound:u_bound,j) )/size( S( l_bound:u_bound,1 ) )
		 	 	
		 	 	S(l_bound:u_bound,j)=S(l_bound:u_bound,j)/media
		 	
		 	 endif
		 	 
		 	 if (lowest <0 .AND. upest >size(S,1) ) then
		 	 
		 	 	l_bound=1 
		 	       u_bound=size(S,1)
		 	       
		 	       media=SUM( S(l_bound:u_bound,j) )/size( S( l_bound:u_bound ,1 ) )
		 	 	
		 	 	S(l_bound:u_bound,j)=S(l_bound:u_bound,j)/media
		 	
		 	 endif
		 	  
		    
		    
 	END DO
 	 
 END DO
 
 
 	
 
 
 
 END FUNCTION punto2
 
 
END MODULE Ex1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Fine EX1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EX2

MODULE Ex2

IMPLICIT NONE

CONTAINS

 FUNCTION punto3(S,input2) result(bool)
  

 LOGICAL :: bool
 INTEGER,PARAMETER :: intervalli=10**(5)
 DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:),INTENT(IN) :: S
 DOUBLE PRECISION,DIMENSION(size(S,1)*size(S,2)) :: s_cons 
 
 DOUBLE PRECISION,DIMENSION(intervalli) :: centroide,counts
 
 DOUBLE PRECISION ::max_min,lower,upper
 
 CHARACTER(1),INTENT(IN) :: input2
 
 INTEGER::i,j,k,N,suddivisione
 
 if (input2=="G" .OR. input2=="g") then
 
 suddivisione =size(s,3)
 
 endif
 
 
  if (input2=="L" .OR. input2=="l") then
  
 suddivisione =size(s,3)/2
 
 endif
 
  if (suddivisione==0) then
 
 write (*,*) "Pay attention suddivision is ",suddivisione ," !" 
 
 endif
	 
 DO i=1,size(S,1)
 		Do k=1,size(S,2)
 			s_cons( k+size(S,2)*(i-1) ) = S(i, k,suddivisione )    

 		end DO
 		
 End Do
 

max_min=maxval( s_cons)-minval(s_cons )

!WRITE(*,*) "il range è ",max_min
DO i=1,size(centroide)
	
	N=0
	
	lower=minval(s_cons ) +(max_min/intervalli)*(i-1)
	
	if (i /= size(centroide) ) then 
	
		upper=minval(s_cons) + (max_min/intervalli)*i
	
	else
	
		upper=maxval( s_cons)
		
	endif
	
	centroide(i)=minval(s_cons) + (max_min/intervalli)/2 + (max_min/intervalli)*(i-1)
		
		DO j=1,size(s_cons )
			
			if ( lower <= s_cons(j) .AND. s_cons(j) <=upper)  then
			
				N=N+1
			endif
			
			
			
		END DO
	
	counts(i)=N	
		
END DO

if ( SUM(counts) /= size(s_cons )  ) then

	WRITE(*,*) SUM(COUNTS),size(s_cons )

	Write(*,*) "I'm losting some date for hist!"

	stop

end if


counts=counts/( SUM(counts)*(max_min/intervalli) )

N=0

Do i=1,size(counts)
	

	if (counts(i) /= 0 ) then
	
		N=N+1
		
	end if
	
end do

WRITE(*,*) "Not zero entries are : ",real(N)*100/size(counts),"%"

if(input2=="G" .OR. input2=="g") then

	OPEN(UNIT=16,FILE='hist_hermitiana_globale.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 
	     
	DO i=1,size(counts)
	     
		WRITE(UNIT=16,FMT=*) centroide(i),counts(i)
		
	END DO
	     
	 CLOSE(UNIT=16) 

endif!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (input2=="L" .OR. input2=="l") then
 
	OPEN(UNIT=16,FILE='hist_hermitiana_locale.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 
	     
	DO i=1,size(counts)
	     
		WRITE(UNIT=16,FMT=*) centroide(i),counts(i)
		
	END DO
	     
	 CLOSE(UNIT=16) 


endif				

bool=.TRUE.

 

 
 
 END FUNCTION punto3
 
 
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 FUNCTION punto4(S,input2) result(bool)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL :: bool
 INTEGER,PARAMETER :: intervalli=10**(2)
 DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:),INTENT(IN) :: S
 DOUBLE PRECISION,DIMENSION(size(S,1)*size(S,2)) :: s_cons 
 
 DOUBLE PRECISION,DIMENSION(intervalli) :: centroide,counts
 
 DOUBLE PRECISION ::max_min,lower,upper
 
 CHARACTER(1),INTENT(IN) :: input2
 
 INTEGER::i,j,k,N,suddivisione
 
 if (input2=="G" .OR. input2=="g") then
 
 suddivisione =size(s,3)
 
 endif
 
 
  if (input2=="L" .OR. input2=="l") then
  
 suddivisione =size(s,3)/2
 
 endif
 
 if (suddivisione==0) then
 
 write (*,*) "Pay attention suddivision is ",suddivisione ," !" 
 
 endif
	 
 DO i=1,size(S,1)
 		Do k=1,size(S,2)
 			s_cons( k+size(S,2)*(i-1) ) = S(i, k,suddivisione )     

 		end DO
 		
 End Do
 

max_min=maxval( s_cons)-minval(s_cons )

if (max_min==0) then

Write(*,*) " you cannot do an histogram with only 2 eigenvalues! "

stop

end if

!WRITE(*,*) "il range è ",max_min
DO i=1,size(centroide)
	
	N=0
	
	lower=minval(s_cons ) +(max_min/intervalli)*(i-1)
	
	if (i /= size(centroide) ) then 
	
		upper=minval(s_cons) + (max_min/intervalli)*i
	
	else
	
		upper=maxval( s_cons)
		
	endif
	
	centroide(i)=minval(s_cons) + (max_min/intervalli)/2 + (max_min/intervalli)*(i-1)
		
		DO j=1,size(s_cons )
			
			if ( lower <= s_cons(j) .AND. s_cons(j) <=upper)  then
			
				N=N+1
			endif
			
			
			
		END DO
	
	counts(i)=N	
		
END DO

if ( SUM(counts) /= size(s_cons )  ) then

	WRITE(*,*) SUM(COUNTS),size(s_cons )

	Write(*,*) "I'm losting some date for hist!"

	stop

end if


counts=counts/( SUM(counts)*(max_min/intervalli) )

N=0

Do i=1,size(counts)
	

	if (counts(i) /= 0 ) then
	
		N=N+1
		
	end if
	
end do

WRITE(*,*) "Not zero entries are : ",real(N)*100/size(counts),"%"

if(input2=="G" .OR. input2=="g") then

	OPEN(UNIT=16,FILE='hist_diagonale_globale.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 
	     
	DO i=1,size(counts)
	     
		WRITE(UNIT=16,FMT=*) centroide(i),counts(i)
		
	END DO
	     
	 CLOSE(UNIT=16) 

endif!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 if (input2=="L" .OR. input2=="l") then
 
	OPEN(UNIT=16,FILE='hist_diagonale_locale.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 
	     
	DO i=1,size(counts)
	     
		WRITE(UNIT=16,FMT=*) centroide(i),counts(i)
		
	END DO
	     
	 CLOSE(UNIT=16) 


endif				

bool=.TRUE.
 
 
 END FUNCTION punto4
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 FUNCTION punto5(S,input2) result(R_average)
 
 IMPLICIT NONE
 
 DOUBLE PRECISION,DIMENSION(:,:,:),INTENT(IN) :: S
 CHARACTER(1),INTENT(IN) ::input2

 
 DOUBLE PRECISION,DIMENSION(size(S,2),size(S,1) ) ::R_matrices
 
 DOUBLE PRECISION,DIMENSION( size(S,1) ) :: R
  
 DOUBLE PRECISION ::R_average  
 INTEGER ::i,j,k,suddivisione
 
 
  if (input2=="G" .OR. input2=="g") then
 
 suddivisione =size(s,3)
 
 endif
 
 
  if (input2=="L" .OR. input2=="l") then
  
 suddivisione =size(s,3)/2
 
 endif
 
 
DO k=1,size(R_matrices,2) 
	DO j=1,size(R_matrices,1)                                      
 	
 		R_matrices(j,k)=min(S(k,j+1,suddivisione),S(k,j,suddivisione))/max(S(k,j+1,suddivisione),S(k,j,suddivisione))
 	END DO
 	
 	R(k)=sum(R_matrices(:,k) )/size( R_matrices (:,k) )
 	
END DO 


R_average=SUM(R)/SIZE(R)

END FUNCTION punto5



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE Ex2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Fine EX2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! inizio program
PROGRAM main

USE Ex1

USE Ex2


IMPLICIT NONE


INTEGER,PARAMETER :: d=1000,suddivisioni=20,matrices=50

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: W_sample
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) ::W,S_sample
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,: ) :: S
COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) :: A_sample
COMPLEX*16,ALLOCATABLE,DIMENSION(:,:,:) :: A

DOUBLE PRECISION ::R,zero

DOUBLE PRECISION:: r1,r2

INTEGER ::i,j,k

CHARACTER(1) ::input1,input2

LOGICAL :: bool

allocate(A(1:matrices,1:d,1:d),S(1:matrices,1:d-1,1:suddivisioni) )
allocate( A_sample(1:d,1:d),S_sample(1:d-1,1:suddivisioni ), W(1:d,1:matrices) )
allocate ( W_sample(1:d) )


Do k=1,matrices

	DO i=1,d
		DO j=1,d
			CALL random_number(r1)
			CALL random_number(r2)	
			
			if (j==i)then
			A(k,i,j)=CMPLX(r1*10.0,zero)
			
			else
			
			A(k,i,j)=CMPLX(r1*10.0,r2*10.0)
			
			end if
			
		END DO
	END DO
	
END DO

WRITE(*,*) "Which matrices you want initialize?[H = hermitian D =diagonal]"
READ(*,*) input1


if (input1=="H" .OR. input1=="h" )  then 
	FORALL(i=1:d,j=1:d,k=1:matrices, i<j)
		

		A(k,i,j)=CONJG(A(k,j,i) )		

	END FORALL


	Write(*,*) "Hermitian matrices initialized !"

end if

if (input1=="D" .OR. input1=="d" )  then 

	FORALL(i=1:d,j=1:d,k=1:matrices, i<j)
		

		A(k,i,j)=CMPLX(zero,zero)
		A(k,j,i)=CMPLX(zero,zero)		

	END FORALL


	Write(*,*) "Diagonal matrices initialized !"

end if

if (input1 /= "H" .AND. input1 /= "D" .AND. input1 /= "h" .AND. input1 /= "d" ) then

	Write(*,*) input1,"value is not permitted!"
	
	stop
endif


WRITE(*,*) "Global or Local?[G = global L=local]"
READ(*,*) input2


if (input2 /= "G" .AND. input2 /= "g" .AND. input2 /= "L" .AND. input2 /= "l" ) then

	Write(*,*) input2,"value is not permitted!"
	
	stop
endif

DO k=1,matrices
	FORALL(i=1:d,j=1:d)
		

		A_sample(i,j)=A(k,i,j)		


	END FORALL

	
	W_sample=punto1(A_sample)
		
	FORALL (i=1:d)
	
		W(i,k)=W_sample(i)
		
	END FORALL
	

	
	
	S_sample=punto2(W_sample,suddivisioni)
	
	DO i=1,d-1
		DO j=1,suddivisioni

			S(k,i,j)=S_sample(i,j)
		
		END DO
	
	END DO
	
	

END DO


			
			R=punto5(S,input2)
			WRITE(*,*) "Average r value is : ", R



if ((input1 == "H" .OR. input1 == "h") .AND.( (input2 == "G" .OR. input2 == "g") .OR. (input2 == "L" .OR. input2 == "l") ) ) then

	bool=punto3(S,input2)
	
end if



if ((input1 == "D" .OR. input1 == "d") .AND.( (input2 == "G" .OR. input2 == "g") .OR. (input2 == "L" .OR. input2 == "l") ) ) then

	bool=punto4(S,input2)
	
end if

if (bool .EQV. .TRUE.) then

	WRITE(*,*) "Program works properly!"
	
endif
 
     

deallocate (A,S)
deallocate(A_sample,S_sample,W)
deallocate(W_sample)





END PROGRAM main

