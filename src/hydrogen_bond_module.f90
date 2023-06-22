module hydrogen_bonds

    use distance_module

contains
    subroutine Read_Input(inputfile, criteria, donor_selection, selections, traj_name, idx_name)
    
        implicit none
        ! Input ********************************************************************
        character(len=*), intent(in) :: inputfile
        ! Output *******************************************************************
        real, allocatable, intent(out) :: criteria(:,:)
        character(len=100), dimension(2), intent(out) :: donor_selection
        character(len=100), allocatable, intent(out) :: selections(:)
        character(len=100), intent(out) :: traj_name, idx_name
        ! Local ********************************************************************
        integer :: i, num_acc_selections
        ! **************************************************************************

        open(10, file=trim(inputfile), status='old')
        read(10,*) 
        read(10,*) traj_name, idx_name
        read(10,*)
        read(10,*) donor_selection(1), donor_selection(2)
        read(10,*) 
        read(10,*) num_acc_selections
        read(10,*)
        allocate(criteria(num_acc_selections,3))
        allocate(selections(num_acc_selections))
        criteria = 0.0
        selections = "None"
        do i=1, num_acc_selections
            read(10,*) selections(i), criteria(i,1), criteria(i,2), criteria(i,3)
        end do

    end subroutine Read_Input

    subroutine hbond_criteria(rd_O, rd_H, ra_X, box, criteria, hbonded, rHO, theta)

    implicit none

    ! Input variables *******************************************************
    real, intent(in), dimension(3) :: rd_O(3), rd_H(3), ra_X(3), box(3)
    real, intent(in), dimension(:) :: criteria
    ! Output variables ******************************************************
    integer, intent(out) :: hbonded
    real, intent(out) :: rHO, theta
    ! Local variables *******************************************************
    ! ************************************************************************

    ! Initialize the output variables
    hbonded = 0; rHO = 0.0; theta = 0.0

    ! Calculate the HX distance
    rHO = periodic_distance2(rd_H, ra_X, box)
    if ( rHO < criteria(2) ) then
        theta = angle_between_points(rd_H, rd_O, ra_X, box)
        if ( theta < criteria(3) )then
            hbonded = 1
        endif
    endif

    end subroutine hbond_criteria

    subroutine find_H_bonds(r_don_O, r_don_H, r_acc_X, box, criteria, hbond_values, hbond_donH_idx, hbond_accX_idx, hbond_count)
    ! This is a subroutine to calculate hydrogen bonds from a set of coordinates.
    ! The subroutine takes as input the coordinates of the donor and acceptor atoms,
    ! the box size, the criteria for the hydrogen bond, and the output arrays.
    ! 
    ! This uses the cell_list_distance library, which must be compiled before this library.
    ! 
    ! Input:
    !   r_don_O: 2D array of the coordinates of the donor oxygen atoms
    !   r_don_H: 2D array of the coordinates of the donor hydrogen atoms
    !   r_acc_X: 2D array of the coordinates of the acceptor X=(N,O,etc) atoms
    !   box: 1D array of the box size
    !   criteria: 1D array of the criteria for the hydrogen bond
    !
    ! Output:
    !   hbond_values: sparse 2D array of the hydrogen bond values (OO distance, HO distance, angle)
    !   hbond_donH_idx: 1D array of the donor hydrogen index
    !   hbond_accX_idx: 1D array of the acceptor X index
    !
    
    implicit none

    ! Input variables *******************************************************
    real, intent(in) :: r_don_O(:,:), r_don_H(:,:), r_acc_X(:,:)
    real, intent(in) :: box(:), criteria(:)
    ! Output variables ******************************************************
    real, intent(out) :: hbond_values(:,:)
    integer, intent(out) :: hbond_donH_idx(:), hbond_accX_idx(:)
    integer, intent(out) :: hbond_count
    ! Local variables *******************************************************
    ! Sparse array for the OO distance calculation
    real, allocatable :: dr_OO_values(:)
    integer, allocatable :: dr_OO_don_idx(:), dr_OO_acc_idx(:)
    ! Values needed for the cell_list_distance module
    integer :: pred_max_hbonds
    real :: cell_length
    
    integer :: loop_index, oindex, hindex, aindex, i
    real :: rHO, theta
    integer :: hbonded
    !************************************************************************

    ! Allocate Arrays
    ! Each water can donate to 2 acceptors, so the maximum number of hydrogen bonds is 2*number of waters
    ! Added some buffer to be safe by doubling this.
    pred_max_hbonds = size(r_don_O,1)*2
    allocate(dr_OO_values(pred_max_hbonds*2))
    allocate(dr_OO_don_idx(pred_max_hbonds*2))
    allocate(dr_OO_acc_idx(pred_max_hbonds*2))

    ! Initialize Arrays
    hbond_values = 0.0; hbond_donH_idx = 0; hbond_accX_idx = 0

    ! Calculate the distance between the donor and acceptor atoms
    cell_length = sqrt(criteria(1))/2.0
    call cell_list_distance(r_don_O, r_acc_X, box, cell_length, criteria(1), dr_OO_values, dr_OO_don_idx, dr_OO_acc_idx)

    ! Check the remainder of the hydrogen bonds
    loop_index = 1
    hbond_count = 0
    do while (loop_index < size(dr_OO_values,1) .and. dr_OO_values(loop_index) /= 0)
        ! Check the hydrogen bond criteria, noting that OO already satisfies it.
        oindex = dr_OO_don_idx(loop_index)
        aindex = dr_OO_acc_idx(loop_index)
        do i=1,2
            ! Where rH is in the r_don_H array
            hindex = (oindex-1)*2 + i
            call hbond_criteria(r_don_O(oindex,:), r_don_H(hindex,:), r_acc_X(aindex,:), box, criteria, hbonded, rHO, theta)
            if ( hbonded == 1 ) then
                hbond_count = hbond_count + 1
                hbond_values(hbond_count,1) = dr_OO_values(loop_index)
                hbond_values(hbond_count,2) = rHO
                hbond_values(hbond_count,3) = theta
                hbond_donH_idx(hbond_count) = hindex
                hbond_accX_idx(hbond_count) = aindex
                exit ! skip the rest of the i loop if found a hydrogen bond already
            endif
        enddo
        loop_index = loop_index + 1
    enddo

    end subroutine find_H_bonds

end module hydrogen_bonds