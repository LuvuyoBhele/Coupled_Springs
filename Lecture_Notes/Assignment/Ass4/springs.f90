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
    
end program springs

