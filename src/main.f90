subroutine Read_Input(inputfile, criteria, donor_selection, selections, traj_name, idx_name)
    
    implicit none
    ! Input ********************************************************************
    character(len=100), intent(in) :: inputfile
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

program hydrogen_bond_analysis

    use gmxfort_trajectory
    use hydrogen_bonds

    implicit none

    type (Trajectory) :: trj
    ! Parameters ***************************************************************
    integer, parameter :: chunk_size = 100

    ! Input ********************************************************************
    character(len=100) :: traj_name, idx_name
    character(len=100), dimension(2) :: donor_selection
    character(len=100), allocatable :: selections(:)
    real, allocatable :: criteria(:,:)


    integer :: number_of_frames, number_of_atoms
    ! Local ********************************************************************
    integer :: chunk, fr_idx, i, hbond, acc
    integer :: number_of_chunks, chunk_stop

    integer :: num_donor_O, num_donor_H
    integer, allocatable :: num_acc_atoms(:)
    real, allocatable :: rdonO(:,:), rdonH(:,:), racc(:,:)

    real, allocatable :: hbond_values(:,:,:,:)
    integer, allocatable :: hbond_donH_idx(:,:,:), hbond_accX_idx(:,:,:)
    integer, allocatable :: hbond_count(:,:)
    integer :: pred_max_hbonds
    ! **************************************************************************
    
    ! Open Output File
    open(10, file="hbonding.out")


    ! Read Input File
    call Read_Input("hbonding.in", criteria, donor_selection, selections, traj_name, idx_name)

    ! Open Trajectory
    call trj%open(trim(traj_name), trim(idx_name))
    number_of_frames = trj%nframes
    number_of_atoms = trj%natoms()

    ! Find number of atoms of each type
    num_donor_O = trj%natoms(trim(donor_selection(1)))
    num_donor_H = trj%natoms(trim(donor_selection(2)))
    do i=1, size(selections,1)
        num_acc_atoms(i) = trj%natoms(trim(selections(i)))
    enddo
    pred_max_hbonds = num_donor_H*2
! ****************************************************************************************************
! ALLOCATE ARRAYS
    allocate(racc(maxval(num_acc_atoms),3))
    allocate(rdonO(num_donor_O,3))
    allocate(rdonH(num_donor_H,3))

    
    allocate(hbond_values( chunk_size,   size(selections,1),    pred_max_hbonds, 3))
    allocate(hbond_donH_idx( chunk_size, size(selections,1),    pred_max_hbonds))
    allocate(hbond_accX_idx( chunk_size, size(selections,1),    pred_max_hbonds))
    allocate(hbond_count( chunk_size,    size(selections,1)))

! ****************************************************************************************************
! Loop over Chunks
! ****************************************************************************************************

    ! Set Calculation Values
    number_of_chunks = ceiling(real(number_of_frames/chunk_size))
    chunks: do chunk=1, number_of_chunks
        ! Zero arrays
        racc = 0.0; rdonO = 0.0; rdonH = 0.0
        hbond_values = 0.0; hbond_donH_idx = 0; hbond_accX_idx = 0; hbond_count = 0

        write(*,*) "Chunk is %d", chunk
        chunk_stop = chunk_size
        ! Grab the frames and bring them into memory
        chunk_stop = min(chunk_size, number_of_frames-(chunk-1)*chunk_size)
        ! Read in chunk of coordinates
        trj%read_next(chunk_stop)

! ****************************************************************************************************
! Calculate H-Bonds
        frames: Do fr_idx=1, chunk_stop
            ! Grab Coordinates Donor
            rdonO = trj%x(fr_idx,  group=trim(donor_selection(1)))
            rdonH = trj%x(fr_idx,  group=trim(donor_selection(2)))

            ! Loop over Acceptor Types
            Do acc=1, size(selections,1)
                ! Grab Coordinates Acceptor
                racc = trj%x(fr_idx,  group=trim(selections(acc)))

                ! Calculate H-Bonds
                call find_H_bonds(rdonO, rdonH, racc, criteria(acc,:) & ! *****
                    , hbond_values(fr_idx, acc, :,:), hbond_donH_idx(fr_idx, acc, :) &                  ! ***** HBOND CALCULATION
                    , hbond_accX_idx(fr_idx, acc, :), hbond_count(fr_idx, acc))                        ! *****
            EndDo
        EndDo frames
! ****************************************************************************************************
        ! Write H-Bonds
        write(*,*) "Writing H-Bonds"
        do fr_idx=1, chunk_stop
            write(10,*) fr_idx, 0
            do acc=1, size(selections,1)
                write(10,*) hbond_count(fr_idx, acc)
                do hbond=1, hbond_count(fr_idx, acc)
                    write(10,*) hbond_donH_idx(fr_idx, acc, hbond), hbond_accX_idx(fr_idx, acc, hbond)
                enddo
            enddo
        enddo
    enddo chunks
! ****************************************************************************************************
! END CHUNKS
! ****************************************************************************************************

    ! Close Output file
    close(10)

end program hydrogen_bond_analysis