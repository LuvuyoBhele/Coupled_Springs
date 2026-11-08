program springs
    implicit none
    integer :: n, row, col
    real, allocatable :: D(:,:)
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
    transD = trans(D, size(D,1))
    !if all(D == transD) then
    !   print*, "matrix is symmetric"
    
end program springs

real trans(Q, n) result(retval)
    implicit none
    integer, intent(in) :: n
    real, intent(in) :: Q(n, n) ! Removed 'allocatable'
    integer :: i, j              ! Declared properly; removed duplicate 'n'
    do i = 1, n
        do j = 1, n
        a = Q(i,j) ! Changed x(n,n) to x(i,i)
        b = Q(j,i)
        Q(i,j) = b
        Q(j,i) = a
        end do
    end do                         ! Fixed missing 'end do'

    retval = sum_val            ! Assigned the result to your return variable
end function trans


