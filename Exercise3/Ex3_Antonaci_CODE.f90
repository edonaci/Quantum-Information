!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Exercise1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE checkpoint_debug

    implicit none

!Here i define a interface inside a checkpoint_debug MODULE in order to summarize all matrices subroutines checks
    
    INTERFACE checkpoint
    
           MODULE PROCEDURE c_matrix_real
           MODULE PROCEDURE c_matrix_complex
           MODULE PROCEDURE c_matrix_integer
           
    END INTERFACE checkpoint
    
    contains
!################################################################   
!DOCUMENTATION : c_matrix_real takes two real matrices and checks if they have same dimensions and  same elements
!################################################################ 
!ARGUMENTS : A is a double real matrix with no defined dimensions to specify when it's called.
!            B is a double real matrix with no defined dimensions to specify when it's called.
!###############################################################
!CONDITIONS :
! post matrix-matrix calculation: A and B dimensions and elements must be the same

    SUBROUTINE c_matrix_real(A,B)
    
    implicit none
    
    real,parameter :: threshold=10E-5
    real,dimension(:,:),intent(in) :: A,B
    logical:: debug 
    integer :: count=0,i,j
    real :: aux=0.0
     
    c_result :   if(size(A,dim=1)==size(B,dim=1) .AND. size(A,dim=2)==size(B,dim=2)    ) then        ! check if dimensions are equals
    
    WRITE(*,*) "Matrixes have same dimensions" ,shape(source=A)
    
    !WRITE(*,*) "A ",A                                                !i suggest to uncomment these lines if A and B have a not so big rank
    !WRITE(*,*) "B ",B 
    
    else
    
    WRITE(*,*)  "Matrices have different dimensions!",shape (source=A ),"and",shape(source=B)
    
    STOP
    
    end if  c_result
    
    
        
        do i=1,size(A,dim=2)
              do j=1,size(A,dim=1)
                    
                   aux=abs( A(j,i)-B(j,i) )/abs(A(j,i) ) 
                    
                    c_each_value :if ( aux > threshold ) then           !check if elements are almost the same
                            
                                count=count+1
                                !WRITE(*,*) "A(",j,",",i,")" ,A(j,i)                  i suggest to uncomment these lines if A and B have a not so big rank
                                !WRITE(*,*) "B (",j,",",i,")" ,B(j,i)
                                !WRITE(*,*) "il valore è" ,aux
                    end if     c_each_value
                    
               end do
         end do    
                
                    
         if (count==0) then              !debug variable
         
         debug=.true.
         
         WRITE(*,*) "Matrices elements are almost the same?" ,debug
         
         else
         
         debug =.false.
         
          WRITE(*,*) "Matrices elements are almost the same?" ,debug
          
          STOP
          
          end if
        
    
    
    end  SUBROUTINE c_matrix_real
 
!################################################################   
!DOCUMENTATION : c_matrix_complex takes two complex matrices and checks if they have same dimensions and  same elements
!################################################################ 
!ARGUMENTS : A is a double complex matrix with no defined dimensions to specify when it's called.
!            B is a double complex matrix with no defined dimensions to specify when it's called.
!###############################################################
!CONDITIONS :
! post matrix-matrix calculation: A and B dimensions and elements must be the same   
!####################################################
    
    SUBROUTINE c_matrix_complex(A,B)
    
    implicit none
    
    real,parameter :: threshold=10E-5
    complex,dimension(:,:),intent(in) :: A,B
    logical:: debug 
    integer :: count=0,i,j
    real :: aux=0.0
     
    c_result :   if(size(A,dim=1)==size(B,dim=1) .AND. size(A,dim=2)==size(B,dim=2)    ) then               ! check if dimensions are equals
    
    WRITE(*,*) "Matrices have same dimensions" ,shape(source=A)
    
    !WRITE(*,*) "A ",A                                                !i suggest to uncomment these lines if A and B have a not so big rank
    !WRITE(*,*) "B ",B
    else
    
    WRITE(*,*)  "Matrices have different dimensions!",shape (source=A ),"and",shape(source=B)
    
    STOP
    
    end if  c_result
    
    
        
        do i=1,size(A,dim=1)
              do j=1,size(A,dim=2)
                    aux=abs( A(j,i)-B(j,i) )/abs(A(j,i) )
                    
                    c_each_value :if ( aux > threshold ) then                      !check if elements are almost the same
                                !WRITE(*,*) "A(",j,",",i,")" ,A(j,i)                  i suggest to uncomment these lines if A and B have a not so big rank
                                !WRITE(*,*) "B (",j,",",i,")" ,B(j,i)
                                !WRITE(*,*) "il valore è" ,aux
                                count=count+1
                                
                    end if     c_each_value
                    
               end do
         end do    
                
                    
        if (count==0) then              !debug variable
         
         debug=.true.
         
         WRITE(*,*) "Matrices elements are almost the same?" ,debug
         
         else
         
         debug =.false.
         
          WRITE(*,*) "Matrices elements are almost the same?" ,debug
          
          STOP
          
          end if
        
    
    
    end  SUBROUTINE c_matrix_complex
    
 !################################################################   
!DOCUMENTATION : c_matrix_integer takes two integer matrices and checks if they have same dimensions and  same elements
!################################################################ 
!ARGUMENTS : A is a double integer matrix with no defined dimensions to specify when it's called.
!            B is a double integer matrix with no defined dimensions to specify when it's called.
!###############################################################
!CONDITIONS :
! post matrix-matrix calculation: A and B dimensions and elements must be the same  
!####################################################

   
    SUBROUTINE c_matrix_integer(A,B)
    
    implicit none
    
    real,parameter :: threshold=10E-5
    integer,dimension(:,:),intent(in) :: A,B
    logical:: debug 
    integer :: count=0,i,j
    real :: aux=0.0
     
    c_result :   if(size(A,dim=1)==size(B,dim=1) .AND. size(A,dim=2)==size(B,dim=2)    ) then                   ! check if dimensions are equals
    
    WRITE(*,*) "Matrices have same dimensions" ,shape(source=A)
    
     !WRITE(*,*) "A ",A                                                !i suggest to uncomment these lines if A and B have a not so big rank
    !WRITE(*,*) "B ",B
    else
    
    WRITE(*,*)  "Matrices have different dimensions!",shape (source=A ),"and",shape(source=B)
    
    STOP
    
    end if  c_result
    
    
        
        do i=1,size(A,dim=1)
              do j=1,size(A,dim=2)
                    
                    aux=abs( A(j,i)-B(j,i) )/abs(A(j,i) )
                    
                    c_each_value :if ( aux > threshold ) then                !check if elements are almost the same
                                !WRITE(*,*) "A(",j,",",i,")" ,A(j,i)                  i suggest to uncomment these lines if A and B have a not so big rank
                                !WRITE(*,*) "B (",j,",",i,")" ,B(j,i)
                                !WRITE(*,*) "il valore è" ,aux
                                count=count+1
                                
                    end if     c_each_value
                    
               end do
         end do    
                
                    
        if (count==0) then              !debug variable
         
         debug=.true.
         
         WRITE(*,*) "Matrices elements are almost the same?" ,debug
         
         else
         
         debug =.false.
         
          WRITE(*,*) "Matrices elements are almost the same?" ,debug
          
          STOP
          
          end if
        
    
    
    end  SUBROUTINE c_matrix_integer
    
    
    
END MODULE checkpoint_debug

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!2° module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MODULE matrix_mul

    implicit none
    
    INTEGER :: i,j,k
    DOUBLE PRECISION :: summa=0.0
    
    contains
!################################################################   
!DOCUMENTATION : by rows function takes two real matrices and thier dimensions and compute matrix-matrix multiplication by rows
!                farther it checks if every dimension is not zero and if first's 2° dim is different then second's 1° dim
!################################################################ 
!ARGUMENTS : mat1 is a double real matrix with no defined dimensions to specify when it's called.
!            mat2 is a double real matrix with no defined dimensions to specify when it's called.
!            m is an integer stands for A rows 
!            n is an integer stands for A columns and B rows
!            l is an integer stands for B columns 
!###############################################################
!CONDITIONS :
! pre matrix-matrix calculation : m,n,l are not zero and  first's 2° dim is different then second's 1° dim  
!####################################################
    FUNCTION by_rows(mat1,mat2,m,n,l)
    
        implicit none
    
        REAL,DIMENSION(:,:),intent(IN) :: mat1,mat2
        
        INTEGER,intent(in) :: m,n,l
        
        REAL,DIMENSION(l,n) :: mat2_t
        REAL,DIMENSION(m,l) :: mat_mul,by_rows
        
        if(m==0 .OR. n==0 .OR. l==0) then
          
        WRITE(*,*) "Some dimension is 0! "
        
        STOP
          
        else 
          
        continue
          
        endif
          
        
        if (size(mat1,dim=2) /= size(mat2,dim=1) ) then
    
        
        WRITE(*,*)  "Impossible to do matrices multiplication : "
        WRITE(*,*)  "dim2's matrix1 its different than dim1's matrix2"
        
        STOP
        
        end if
        
        mat2_t=transpose(mat2)
    
       	DO i=1,m                                                 
    	  	DO k=1,l
    	  		summa=0.0
    
    			DO j=1,n
    				summa=summa+mat1(i,j)*mat2_t(k,j)            !The first way to do a matrix computation  moving on coloumns
    			
    			END DO
    			
    		mat_mul(i,k)=summa
    			
    		END DO
    	  END DO 
          
          by_rows= mat_mul
    
    	 RETURN 
    
    END FUNCTION by_rows
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!2° function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!################################################################   
!DOCUMENTATION : by columns function takes two real matrices and thier dimensions and compute matrix-matrix multiplication by columns
!                farther it checks if every dimension is not zero and if first's 2° dim is different then second's 1° dim
!################################################################ 
!ARGUMENTS : mat1 is a double real matrix with no defined dimensions to specify when it's called.
!            mat2 is a double real matrix with no defined dimensions to specify when it's called.
!            m is an integer stands for A rows 
!            n is an integer stands for A columns and B rows
!            l is an integer stands for B columns 
!###############################################################
!CONDITIONS :
! pre matrix-matrix calculation : m,n,l are not zero and  first's 2° dim is different then second's 1° dim  
!####################################################
    FUNCTION by_columns(mat1,mat2,m,n,l)
      
          implicit none
      
          REAL,DIMENSION(:,:),intent(IN) :: mat1,mat2
          INTEGER,intent(in) :: m,n,l

          
          REAL,DIMENSION(n,m) :: mat1_t
          REAL,DIMENSION(m,l) :: mat_mul,by_columns
          
          DOUBLE PRECISION :: summa=0.0
          
          if(m==0 .OR. n==0 .OR. l==0) then
          
          WRITE(*,*) "Some dimension is 0! "
          
          STOP
          
          else 
          
          continue
          
          endif
          
          if (size(mat1,dim=2) /= size(mat2,dim=1) ) then
          
          WRITE(*,*)  "Impossible to do matrices multiplication : "
          WRITE(*,*)  "dim2's matrix1 its different than dim1's matrix2"
          
          STOP
          
          end if
          
      
           mat1_t=transpose(mat1)
      
        DO k=1,l                                                 
      	  	DO i=1,m
      	  		summa=0.0
      
      			DO j=1,n
      				summa=summa+mat1_t(j,i)*mat2(j,k)            
      			
      			END DO
      			
      		mat_mul(i,k)=summa
      			
      		END DO
      	END DO 
            
            by_columns=mat_mul
      
      RETURN 
      
      END FUNCTION by_columns	 	 
      
      
END MODULE matrix_mul


PROGRAM test_module

    USE  checkpoint_debug
    USE  matrix_mul
    
    implicit none
    
    INTEGER :: f,m,n,l
    
    INTEGER,PARAMETER :: size=20,r=100,c_1=110,c_2=120 
    
    INTEGER,DIMENSION(size) :: row=(/ (i,i=r,r*size,r) /), col1=(/ (j,j=c_1,c_1*size,c_1) /),col2=(/   (k,k=c_2,c_2*size,c_2)  /)
    
    
    REAL :: start,finish
    
    REAL,DIMENSION(size)    :: T1,T2,T3                            !different array to collect different timing 
    
    REAL,ALLOCATABLE :: A(:,:),B(:,:),C1(:,:),C2(:,:),C3(:,:)
    
    
     DO f=1,size                 ! Do loop on integer arrays 
        
      
        
        
        m=row(f)
        n=col1(f)
        l=col2(f)
        
        allocate(A(1:m,1:n),B(1:n,1:l) )
        allocate( C1(1:m,1:l)  ,C2(1:m,1:l) ,C3(1:m,1:l)   )
        
        
        CALL random_number(A)
  
        CALL random_number(B)
        
        
        print '("   ")'
        WRITE (*,*) "The matrices will have dimensions ",shape(A)," and ",shape(B) 
         
        CALL cpu_time(start)               ! start time consuming first way to do computation
        
        C1=by_rows(A,B,m,n,l)
        
        CALL cpu_time(finish) 
        
        T1(f)=finish-start
        
        WRITE(*,*)"A*B by rows calculation spends ",finish-start," seconds"
        
        
        CALL cpu_time(start)               ! start time consuming second way to do computation
        
        C2=by_columns(A,B,m,n,l)
        
        CALL cpu_time(finish) 
        
        T2(f)=finish-start
        
        WRITE(*,*)"A*B by cloumns calculation spends ",finish-start," seconds"
        
        
        CALL cpu_time(start)               ! start time consuming third way to do computation
        
        C3=matmul(A,B)
        
        CALL cpu_time(finish) 
        
        T3(f)=finish-start
        
        WRITE(*,*)"A*B by matmul function calculation spends ",finish-start," seconds"
        
        

    
    
        WRITE(*,*) "DEBUGGING TEST :" 
       
    
       
         WRITE(*,*) "For by_rows and by_columns :"  
         call checkpoint(C1,C2)
         
         WRITE(*,*) "For by_rows and matmul"  
         call checkpoint(C1,C3)	   
      
         WRITE(*,*) "For by_columns and matmul"  
         call checkpoint(C1,C3)
       
         deallocate(A,B)
         deallocate(C1,C2,C3)
     
     
     END DO
     
     
     !write matrix dimension and time consuming in a txt file	
     OPEN(UNIT=16,FILE='prova.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE') 
     
     DO i=1,size
     
        WRITE(UNIT=16,FMT=*) row(i),col2(i),T1(i),T2(i),T3(i)
        
     END DO
     
     CLOSE(UNIT=16) 
 

END PROGRAM test_module