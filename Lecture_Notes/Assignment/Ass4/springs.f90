program springs
    implicit none
    integer :: n, row, col
    real :: determinanat, trace
    real, allocatable :: D(:,:), transD(:,:), R(:,:)

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
            real, intent(inout) :: Q(n,n)
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

    do row = 1, n !for every row in the matrix D
        print*, D(row,:) !print the row out
    end do

    
    !Now we use the trans(Q, n) to check if its a transpose
    transD = trans(D,size(D,1))
    if (all(D == transD)) then
       print*, "matrix is symmetric"
    end if

    !Now we compute the trace and determinant
    determinanat = det(D)
    print*, "The determinant of D: Det(D)=", determinanat

    !Now we compute the trace
    
    trace = trc(D,size(D,1))
    print*, "The trace of trc(D):", trace

    call gram_schmidt(transD)
   
    allocate(R(n,n))
    R = mtx_multi(transD,D,n,n,n)
    print*, 'The triangle matrix is'
    do row = 1, n !for every row in the matrix D
        print*, R(row,:) !print the row out
    end do


end program springs

function trans(Q, n) result(retval) !function is to calculate the transpose of the input matrix
    integer, intent(in) :: n !n is the dimension of matrix 
    real, intent(inout) :: Q(n, n) ! Removed 'allocatable'
    real :: retval(n,n) !define output of function as matrix 
    integer :: i, j   !i and j are the indicies of matrix
    do i = 1, n ! row of original matrix
        do j = 1, n ! column of original matrix
        ! now swap the values that of the elements that have swapped incidees
            a = Q(i,j) 
            b = Q(j,i)
            Q(i,j) = b
            Q(j,i) = a
        end do
    end do                         ! Fixed missing 'end do'

    retval = Q            ! Assigned the result to your return variable
end function trans                 ! Fixed missing 'end do'



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
    real, intent(in) :: A(:,:), x(:)
    real :: retval(size(x))
    real :: sum
    integer :: i, j, n
    n = size(A,1)
    do i = 1, n
        sum = 0.0
        do j = 1,n
            sum = sum + A(i,j)*x(j)
        end do
        retval(i)=sum
    end do

end function sol

real function dot(p,q) result(retval)
    real, intent(in) :: p(:), q(:)
    real, allocatable :: prod(:)
    integer :: index, n
    real :: sum

    n = size(p,1)
    allocate(prod(n))
    sum = 0
    prod = p * q
    do index = 1, n
        sum = sum + prod(index)
    end do
    retval = sum
end function dot


function eigen_vec(A,x) result(lambda) 

    real, intent(in) :: A(:,:), x(:)
    real :: lambda

    interface
        function sol(A,x) result(retval)
           real, intent(in) :: A(:,:), x(:)
            real :: retval(size(x))
        end function sol
    end interface

    integer :: n
    real, allocatable :: sol_vec(:)
    
    n = size(A,1)
    allocate(sol_vec(n))
    sol_vec = sol(A,x)
    lambda = sol_vec(1)/x(1)
end function eigen_vec


subroutine gram_schmidt(Q)
    interface
        real function dot(p,q) result(retval)
            real, intent(in) :: p(:), q(:)
        end function
    end interface

    real, intent(inout) :: Q(:,:)
    real :: c, norm
    integer n, index1, index2

    n = size(Q,1)

    do index1 = 1, n

        do index2 = 1, index1-1

            c = dot(Q(index2,:),Q(index2,:))
            if (c > 1.0e-12) then
                c = dot(Q(index2,:),Q(index1,:)) / c
                Q(index1,:) = Q(index1,:) - c*Q(index2,:)
            end if

        end do

        norm = sqrt(dot_product(Q(index1,:),Q(index1,:)))

        if (norm > 1.0e-12) then
            Q(index1,:) = Q(index1,:) / norm
        end if

    end do

    do index1 = 1, n
        print *, Q(index1,:)
    end do

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
            mtx3(row_index, col_index) = 0.0
            do m = 1, size(mtx1,1)
                mtx3(row_index,col_index) = mtx3(row_index,col_index) + mtx1(row_index,m) * mtx2(m,col_index) !here you're multiply the row
            end do
        end do
    end do

    retval = mtx3

    
end function mtx_multi
