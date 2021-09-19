  

 !!***Module***!!!
MODULE def_matrix
  IMPLICIT NONE
  SAVE

  INTEGER,PARAMETER :: m=4
  TYPE double_complex_matrix
  
    INTEGER :: matrix_dimension
    COMPLEX,DIMENSION(m,m) :: matrix_elements 
    COMPLEX :: matrix_trace,matrix_determinant
  
  END TYPE double_complex_matrix
  
END MODULE def_matrix           

 !!*****main test program     *****!!
PROGRAM test 
    
      USE def_matrix
      
      IMPLICIT NONE
      
      TYPE (double_complex_matrix) :: inizi_matrix, matrix_adjoint,dcm,dcm1,dcm2
      
     INTEGER ::i,j 
    COMPLEX ::trace
       
      INTERFACE     
      
      
           FUNCTION inizi_matrix(dcm)        !!! ****Inizi_matrix interface***
          
           USE def_matrix
           IMPLICIT NONE
           
           TYPE (double_complex_matrix) :: dcm
           
          END FUNCTION inizi_matrix
          
          FUNCTION  trace(dcm)             !!!***trace interface**!!!
       
          USE def_matrix
          IMPLICIT NONE
          
          TYPE (double_complex_matrix) :: dcm
          
          END FUNCTION trace
          
          
          FUNCTION  matrix_adjoint(dcm)   !!**matrix_adjoint interface**!!
      
          USE def_matrix
          IMPLICIT NONE
          
          TYPE (double_complex_matrix) :: dcm
          
          END FUNCTION matrix_adjoint
          
          
       END INTERFACE
       
 dcm1=inizi_matrix(dcm) 
 
 dcm2=matrix_adjoint(dcm1)
    


 !**** MATRIX  INITIALIZE  *****!
 
 
 
 WRITE (*,*)  "The matrix dimensions are " ,dcm1%matrix_dimension
 
 WRITE(*,*) "The matrix elements are :"
      DO i=1,m
        
            WRITE(*,*)  (dcm1%matrix_elements(i,j), j=1,dcm1%matrix_dimension )       
        
      END DO
 WRITE(*,*) " the matrix trace is ",trace (dcm1)
 
 
 !**** ADJOINT MATRIX ***********!
 
 WRITE (*,*)  "The adjoint matrix dimensions are " , dcm2%matrix_dimension
 
 

 WRITE(*,*) " The adjoint matrix elements are :"
 
      DO i=1,m
        
            WRITE(*,*)  (dcm2%matrix_elements(i,j), j=1,dcm2%matrix_dimension )       
        
      END DO
      
 WRITE(*,*) " the adjoint matrix trace is ",trace (dcm2) 
 
 
 !******FILE******!
 CALL wf (dcm1)
  
END PROGRAM test    


!**********************inizi_matrix***********************!
 
TYPE (double_complex_matrix) FUNCTION  inizi_matrix (dcm)
  
  USE def_matrix
  
  IMPLICIT NONE
  
  SAVE
  
  INTEGER :: i,j
  REAL :: a,b
  COMPLEX :: ini=(0.0,0.0)
  
  TYPE (double_complex_matrix) :: dcm
  
  dcm%matrix_dimension=m
  
  dcm%matrix_elements= ini
  
  
  DO i=1,m
     DO j=1,m
        call random_number(a)
        call random_number(b)
        
        dcm%matrix_elements(i,j)=CMPLX(10.0*a,10.0*b)
    
    END DO
    
  END DO
  
  dcm%matrix_trace = ini
  
  dcm%matrix_determinant = ini

   inizi_matrix=dcm
 
RETURN
   
END FUNCTION inizi_matrix
 
 !***********************trace******************!
 
 COMPLEX FUNCTION  trace(dcm) 
 
 USE def_matrix
 
 IMPLICIT NONE
   
  
  integer :: i
  COMPLEX :: diag_value=(0.0,0.0)
  
  TYPE (double_complex_matrix) :: dcm
   
  
  DO i=1,m
  
     diag_value=diag_value+dcm%matrix_elements(i,i) 
   
  END DO  
   
   
   trace=diag_value

   diag_value=(0.0,0.0)
    
END FUNCTION trace

 !**********************Matrix_adjoint*******************!

TYPE(double_complex_matrix) FUNCTION  matrix_adjoint(dcm)

 USE def_matrix
 
 IMPLICIT NONE
 
 SAVE
  
 TYPE (double_complex_matrix) :: dcm,dcm1
 integer :: i,j
 
 dcm1=dcm
 
 DO i=1,m
    DO j=1,m
    
        dcm1%matrix_elements(i,j)=conjg( dcm%matrix_elements(i,j)   )
        
     END DO   
  END DO
  
  
  dcm1%matrix_elements=transpose(dcm1%matrix_elements)
  
  matrix_adjoint=dcm1
  
  RETURN
  
 END FUNCTION matrix_adjoint 
 
 
 
 !*************Define SUBROUTINE*********!
SUBROUTINE wf(dcm)

USE def_matrix

IMPLICIT NONE

SAVE

INTEGER :: i ,j

TYPE(double_complex_matrix) ::dcm
 COMPLEX :: trace
 

      OPEN(UNIT=16,FILE='matrix_type.txt',FORM='FORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE',ACTION='WRITE',POSITION='REWIND') 
      
      WRITE(UNIT=16,FMT='(A,I20)',ADVANCE='YES') "matrix dimension is : ", dcm%matrix_dimension
      
      WRITE(UNIT=16,FMT='(A)',ADVANCE='YES') "matrix_elements is : "
      
      DO i=1,m
        
            WRITE(UNIT=16,FMT=*)  (dcm%matrix_elements(i,j), j=1,dcm%matrix_dimension )       
        
      END DO 
      
      WRITE(UNIT=16,FMT=*) "matrix_trace is : ", trace (dcm)
      
      
      ENDFILE 16
      
      CLOSE(UNIT=16,STATUS='KEEP') 
      
      WRITE(*,*) 'I have just written on file matrix_type.txt! '
      
END SUBROUTINE wf