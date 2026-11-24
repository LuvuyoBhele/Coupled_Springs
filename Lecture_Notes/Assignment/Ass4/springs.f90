program springs
    implicit none
    integer :: n, row, col, index
    real :: determinanat, trace
    real, allocatable :: D(:,:), transD(:,:), R(:,:), V(:,:)

    interface
        recursive function det(A) result(retval)!
            real :: retval
            real, intent(in) :: A(:,:)
        end function det

        function det2(A) result(retval)
            real :: retval
            real, intent(in) :: A(:,:)
        end function det2

        function det3(A) result(retval)
            real :: retval
            real, intent(in) :: A(:,:)
        end function det3


        function trans(Q,n) result(retval)
            integer, intent(in) :: n
            real, intent(in) :: Q(n,n)
            real :: retval(n,n)
        end function trans

        real function trc(x, n) result(retval)
            implicit none
            integer, intent(in) :: n
            real, intent(in) :: x(n, n) ! Removed 'allocatable'
        end function trc

        subroutine gram_schmidt(Q)
            real, intent(inout) :: Q(:,:)
        end subroutine gram_schmidt

        
        function mtx_multi(mtx1, mtx2, q, n, p) result(retval)
            integer :: row_index, col_index, m, q, p, n 
            real, intent(in) :: mtx1(q,n), mtx2(n,p) !input matrices are adjustable matricies with adjustable dimensions
            real :: mtx3(q,p)  !introduce a matrix with adjustable dimensions
            real :: retval(q,p) 
        end function mtx_multi

        subroutine Q_i(D,V)
            real, intent(inout) :: D(:,:), V(:,:)
        end subroutine Q_i

        subroutine eigenvalue_verification(D,V,n)
            integer, intent(in) :: n
            real, intent(in) :: D(n,n), V(n,n)
        end subroutine eigenvalue_verification

    end interface



    print*, "Write the dimension number of the matrix" !tell user what input you want from them
    read(*,*) n !ask user to input the desired dimensions of matrix
    allocate(D(n,n)) !Decalare teh type of matrix you're working with
    do row = 1, n !for loop to through each row
        do col = 1, n !for loop to go through each column
            if ( row == col ) then !for the entries on the diagonal
                D(row, col) = 2 !declare the valiues in the diagonal
            else if ( row-col == 1 .or. row-col == -1 ) then !if the absolute difference between the column index and row index is 1
                D(row, col) = -1 !those entries must have be defined by -1
            else !all other entries
                D(row,col) = 0 !must be have 0
            end if
        end do
    end do
    print*,'original system'
    do row = 1, n !for every row in the matrix D
        
        print*, D(row,:) !print the row out
    end do

    !allocate(transD(n,n))
    !!Now we use the trans(Q, n) to check if its a transpose
    !transD = trans(D,size(D,1))
    !if (all(D == transD)) then
!
    !!   print*, "matrix is symmetric"
    !end if

    !Now we compute the trace and determinant
    determinanat = det(D)
    !print*, "The determinant of D: Det(D)=", determinanat

    !Now we compute the trace
    
    trace = trc(D,size(D,1))
    !print*, "The trace of trc(D):", trace

    !call gram_schmidt(transD)
   !
    !allocate(R(n,n))
    !R = mtx_multi(transD,D,n,n,n)
    !!print*, 'The triangle matrix is'
    !do row = 1, n !for every row in the matrix D
    !    print*, R(row,:) !print the row out
    !end do 
    allocate(V(n,n))

    do row = 1, n
        do col = 1, n
            if (col == row) then
                V(row,col) = 1.0
            else
                V(row,col) = 0.0
            end if
        end do
    end do
    open(1,file='matrix.dat')
    do index = 1, 100
        write(1,*) 'matrix',index
        call Q_i(D,V)
    end do
    close(1)
    
!    print*,'The trace of the numerical approximation',trace(D)
!    print*,'The determinant of the numerical approximation', det(D)
    print*, "After Q_i subroutine V is"
    do row = 1, n !for every row in the matrix D
        print*, V(row,:) !print the row out
    end do



    !call eigenvalue_verification(D,V,n)
    

end program springs

function trans(Q, n) result(retval)
    integer, intent(in) :: n
    real, intent(in) :: Q(n,n)
    real :: retval(n,n)
    real :: a, b
    integer :: i, j

    retval = Q

    do i = 1, n
        do j = i+1, n
            a = retval(i,j)
            b = retval(j,i)
            retval(i,j) = b
            retval(j,i) = a
        end do
    end do

end function trans


recursive real function det(A) result(retval)
    real, intent(in) :: A(:,:)
    integer :: n, i, row, col, sub_row, sub_col
    !n = rank of original matrix
    !i = row index going through original matrix
    !r = for loop row index going through sub matrix 
    !c = for col index going through sub matrix
    !sub_row = element row index for sub matrix
    !sub_col = element col index for sub matrix

    interface
        function det2(A) result(retval)
            real :: retval
            real, intent(in) :: A(:,:)
        end function det2

        function det3(A) result(retval)
            real :: retval
            real, intent(in) :: A(:,:)
        end function det3
    end interface



    real, allocatable :: sub(:,:) !defining submatrix
    real :: C !C that will determine the sign of multiplication

    n = size(A, 1) !defining rank of original matrix

    if (n == 2) then !if rank is two call det2 function
        retval = det2(A) !get this as a result
    else if (n == 3) then !if rank is three call det3 function
        retval = det3(A) !get this as a result
    else
        retval = 0.0 !initialise the determinant
        !now define the submatrix
        do i = 1, n
            allocate(sub(n-1, n-1))!declare exsistence
            !define each element explicitly
            sub_row = 0!intialise the count of the rows
            do row = 2, n !these are all the rows in A we going to go through
                
                sub_row = sub_row + 1 ! start search
                sub_col = 0 !initialise column 
                do col = 1, n !we start by going through every column of A
                    if (col /= i) then !as long as the column we doing lapace expansion in is exluded is exluded from the 
                        sub_col = sub_col + 1 
                        sub(sub_row, sub_col) = A(row, col) !define the element of submatrix
                    end if
                end do
            end do
            

            if ((i / 2) * 2 == i) then !this line checks if j is even or odd
                C = -1.0             ! Even columns get a negative sign in row expansion
            else
                C = 1.0              ! Odd columns get a positive sign
            end if
            ! definition of determinant
            retval = retval + C * A(1, i) * det(sub)
            
            deallocate(sub)

        end do
    end if

end function det




real function det2(A) result(retval)!this function computes matrix of rank2
    real, intent(in) :: A(:,:)
    !explicit formula to calculate the determinant
    retval = A(1,1)*A(2,2) - A(1,2)*A(2,1)
end function det2


real function det3(A) result(retval)!this function computes matrix of rank2
    real, intent(in) :: A(:,:)
    !explicit formula to calculate the determinant
    retval = A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) - A(1,2) * (A(2,1)*A(3,3) - A(2,3)*A(3,1)) + A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
end function det3

real function trc(x, m) result(retval)
    implicit none
    integer, intent(in) :: m
    real, intent(in) :: x(m, m) ! Removed 'allocatable'
    integer :: i                 ! Declared properly; removed duplicate 'n'
    real :: sum_val             ! Avoid using the reserved keyword 'sum'

    sum_val = 0.0
    do i = 1, m
        sum_val = sum_val + x(i, i) ! Changed x(n,n) to x(i,i)
    end do                         ! Fixed missing 'end do'

    retval = sum_val            ! Assigned the result to your return variable
end function trc


function sol(A,x) result(retval)
    real, intent(in) :: A(:,:), x(:) !declaring A as the system we are solving and x is the vector that satisfies the system
    real :: retval(size(x)) ! decarling the constant variable of system
    real :: sum !this will be the sum for each component of the vector constant
    integer :: i, j, n !i = row index; j=column index;n = dimension of vector 
    n = size(A,1) ! dimension of system of matrix
    do i = 1, n ! row index
        sum = 0.0 !the sum of the element in the row
        do j = 1,n !column index
            sum = sum + A(i,j)*x(j)
        end do!
        retval(i)=sum!
    end do!
end function sol!

real function dot(p,q) result(retval)
    real, intent(in) :: p(:), q(:) !declare vector p and vector q
    real, allocatable :: prod(:) !declare vector prod
    integer :: index, n !declare index for prod vector, declare n the dimension of vectors
    real :: sum !declare the sum for the summation of of the prod elements

    n = size(p) !the dimensions of the vector
    allocate(prod(n)) !define vector of prod
    sum = 0 !the sum for the final product
    prod = p * q !the mulitiplication of each element to each element to form a vector prod
    do index = 1, n! index of prod vector
        sum = sum + prod(index) !take a en elemnt in prod and and it to sum
    end do
    retval = sum !make sum the out 
end function dot


function eigen_vec(A,x) result(lambda) 

    real, intent(in) :: A(:,:), x(:) !declare system of equations A and constant vector x
    real :: lambda !declare the eigenvalue

    interface !here we intro the interface to describe define
        function sol(A,x) result(retval)
           real, intent(in) :: A(:,:), x(:)
            real :: retval(size(x))
        end function sol
    end interface

    integer :: n
    real, allocatable :: sol_vec(:)
    
    n = size(A,1) !define the dimension the of A
    allocate(sol_vec(n)) !declare the dimension of sol_vec
    sol_vec = sol(A,x) !define the vector that is created from multiplying A and x
    print*, "the solution vector", sol_vec
    print*, "element we're dividing by", x(1)
    lambda = sol_vec(1)/x(1) !define the eigen value as the quotient of the constant vector to solution vector
end function eigen_vec


subroutine gram_schmidt(Q) !create a gram_schmidt subroutine 
    interface
        real function dot(p,q) result(retval)
            real, intent(in) :: p(:), q(:)
        end function
    end interface

    real, intent(inout) :: Q(:,:) !declare the matrix Q which we will make orthogonal
    real :: c, norm !declare c 
    integer n, index1, index2 !declare n to define dimension; index1 as the row index; index1 as the row index

    n = size(Q,1) !define the dimensions by the row of Q

    do index1 = 1, n !for the column index we're orthogonalizing
       ! print*, 'Processing column:', index1
        do index2 = 1, index1-1 !all the vectors we're gonna make orthgonal to
            c = dot(Q(:,index1),Q(:,index2))/dot(Q(:,index2),Q(:,index2)) !this is the normalizing constant for the denominator
        !    print*, 'Projection coefficient:', c
            if (abs(c) > 1.0e-12) then !if condition to make sure normaliziation constant isnt zero
                Q(:,index1) = Q(:,index1) - c*Q(:,index2) !this is how we calculate the next vector
        !        print*, 'Vector after projection:'
        !        print*, Q(:,index1)
            end if

        end do

        norm = sqrt(dot_product(Q(:,index1),Q(:,index1)))!normalise the vector the we diagonalised  

        if (norm > 1.0e-12) then !if the normalisation is greater than zero
            Q(:,index1) = Q(:,index1) / norm !normalise the vector
        !    print*, 'Column after normalisation:', index1
        !    print*, Q(:,index1)
        end if

    end do
    !print*, 'the diagonalized matrix is'
    !do index1 = 1, n
    !    print *, Q(index1,:)
    !end do
!print*,'col 1 and col 2 perp verification:', dot(Q(:,index1),Q(:,index1))

end subroutine gram_schmidt


function mtx_multi(mtx1, mtx2, q, n, p) result(retval)
    integer :: row_index, col_index, m, q, p, n 
    real, intent(in) :: mtx1(q,n), mtx2(n,p) !input matrices are adjustable matricies with adjustable dimensions
    real :: mtx3(q,p)  !introduce a matrix with adjustable dimensions
    real :: retval(q,p) 

    !mtx(matrix,1) this give the number of rows
    !mtx(matrix,2) this give the number of col

    ! Matrix multiplication
    do row_index = 1, q       ! for every row
        do col_index = 1, p    ! for every column
            mtx3(row_index, col_index) = 0.0 !initialise the element
            do m = 1, size(mtx1,1)
                mtx3(row_index,col_index) = mtx3(row_index,col_index) + mtx1(row_index,m) * mtx2(m,col_index) !here you're multiply the row
            end do
        end do
    end do

    retval = mtx3

    
end function mtx_multi


subroutine Q_i(D,V)
    implicit none
    

    interface
        subroutine gram_schmidt(D)
            real, intent(inout) :: D(:,:)
        end subroutine gram_schmidt

        function mtx_multi(mtx1, mtx2, q, n, p) result(retval)
            integer :: q, n, p, m
            real, intent(in) :: mtx1(q,n), mtx2(n,p)
            real :: retval(q,p)
        end function mtx_multi

        function trans(Q,n) result(retval)
            integer, intent(in) :: n
            real, intent(in) :: Q(n,n)
            real :: retval(n,n)
        end function trans
    end interface

    real, intent(inout) :: D(:,:), V(:,:)
    real, allocatable :: Q(:,:), trans_Q(:,:), R(:,:), I(:,:), Q_save(:,:)
    integer :: n, row

    n = size(D,1)

    allocate(Q(n,n), trans_Q(n,n), R(n,n), I(n,n), Q_save(n,n))

    Q = D

    call gram_schmidt(Q)

    print*, "Q after Gram Schmidt"
    do row = 1, n
        print*, Q(row,:)
    end do

    Q_save = Q

    trans_Q = trans(Q,n)

    I = mtx_multi(trans_Q,Q,n,n,n)

    print*, "The verification of orthogonality"
    do row = 1, n
        print*, I(row,:)
    end do

    R = mtx_multi(trans_Q,D,n,n,n)

    print*, "R becomes"
    do row = 1, n
        print*, R(row,:)
    end do

    ! Check V before multiplication
    print*, "V before multiplication"
    do row = 1, n
        print*, V(row,:)
    end do

    Q = mtx_multi(R,Q,n,n,n)

    print*, "Q_{i+1} becomes"
    do row = 1, n
        print*, Q(row,:)
    end do

    D = Q

    ! Check the saved orthogonal matrix
    print*, "Q_save used for eigenvectors"
    do row = 1, n
        print*, Q_save(row,:)
    end do

    V = mtx_multi(V,Q_save,n,n,n)

    print*, "V_{i+1} becomes"
    do row = 1, n
        print*, V(row,:)
    end do
end subroutine Q_i


!finish this subroutine
subroutine eigenvalue_verification(D,V,n) 
    !V is the matrix that contains modes
    !D is the original matrix
    !Q is the matrix withe the eigenvalue
    real, intent(in) :: D(n,n), V(n,n)
    integer, intent(in) :: n
    real :: d_eigen_val
    real :: v_i(n)

    interface
        function eigen_vec(A,x) result(lambda) 
            real, intent(in) :: A(:,:), x(:) !declare system of equations A and constant vector x
            real :: lambda !declare the eigenvalue
        end function eigen_vec

        function sol(A,x) result(retval)
            real, intent(in) :: A(:,:), x(:) !declaring A as the system we are solving and x is the vector that satisfies the system
            real :: retval(size(x)) ! decarling the constant variable of system
            real :: sum !this will be the sum for each component of the vector constant
            integer :: i, j, n !i = row index; j=column index;n = dimension of vector 
        end function sol
    end interface

    integer :: row, col
    do col = 1, n
        v_i = V(:,col) !find the mode 
        d_eigen_val = eigen_vec(D,v_i) !find the eigenvalue
        print*, 'the eigenvalue to to col', col
        print*, d_eigen_val
    end do
end subroutine eigenvalue_verification
    