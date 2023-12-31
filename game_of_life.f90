program game_of_life

  use BATL_lib
  use BATL_pass_cell, ONLY: message_pass_ng_int1

  implicit none

  integer, allocatable:: iState_IIB(:,:,:), nNei_II(:,:)
  real,    allocatable:: r_II(:,:)
  logical, allocatable:: DoRefine_B(:)
  integer:: i, j, iBlock, iStep
  !----------------------------------------------------------------------------
  call init_mpi

  call init_batl( &
       MaxBlockIn     = 100,            &
       CoordMinIn_D   = [-10.,-10.,-0.5],  &
       CoordMaxIn_D   = [+10.,+10.,+0.5],  &
       IsPeriodicIn_D = [.true., .true., .true.] )

  allocate(iState_IIB(MinI:MaxI,MinJ:MaxJ,MaxBlock), &
       nNei_II(nI,nJ), r_II(nI,nJ), DoRefine_B(MaxBlock))
  iState_IIB = 0

  ! Refine the original block twice, so we get 16 blocks (2D)
  do i = 1, 2
     DoRefine_B = .true.
     call init_grid_batl(DoRefine_B)
  end do

  ! Initialize state
  do iBlock = 1, nBlock
     if(Unused_B(iBlock)) CYCLE
     do i=0, iProc
        call random_number(r_II)
     end do
     where(r_II < 0.2) iState_IIB(1:nI,1:nJ,iBlock) = 1
  end do

  ! Save initial state
  call save_plot_block

  do iStep = 1, 400

     ! Updage ghost cells
     call message_pass_ng_int1(iState_IIB)

     do iBlock = 1, nBlock
        if(Unused_B(iBlock)) CYCLE

        ! Count neighbors
        do j = 1, nJ; do i = 1, nI
           nNei_II(i,j) = sum(iState_IIB(i-1:i+1,j-1:j+1,iBlock)) &
                - iState_IIB(i,j,iBlock)
        end do; end do

        ! Update state
        do j = 1, nJ; do i = 1, nI
           if(nNei_II(i,j) == 3) then
              iState_IIB(i,j,iBlock) = 1
           elseif(nNei_II(i,j) /= 2) then
              iState_IIB(i,j,iBlock) = 0
           end if
        end do; end do

     end do

     ! Save updated state
     call save_plot_block

  end do

  deallocate(iState_IIB, nNei_II, r_II, DoRefine_B)
  call clean_batl
  call clean_mpi

contains
  !============================================================================

  subroutine save_plot_block

    use ModPlotFile, ONLY: save_plot_file
    use BATL_lib, ONLY: iProc, &
         nBlock, Unused_B, CoordMin_DB, CoordMax_DB, CellSize_DB

    integer :: iBlock
    character(len=100):: NameFile
    character (len=10) :: TypePosition = 'rewind'
    !--------------------------------------------------------------------------

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       write(NameFile,'(a,i3.3,a,i5.5,a)') &
            'game_pe',iProc,'_blk',iBlock,'.out'

       call save_plot_file(NameFile,     &
            TypeFileIn='real4',          &
            TypePositionIn=TypePosition, &
            nStepIn = iStep, &
            TimeIn  = 0.0, &
            nDimIn  = nDim, &
            CoordMinIn_D = CoordMin_DB(1:nDim,iBlock)        &
            +              0.5*CellSize_DB(1:nDim,iBlock),   &
            CoordMaxIn_D = CoordMax_DB(1:nDim,iBlock)        &
            -              0.5*CellSize_DB(1:nDim,iBlock),   &
            VarIn_IIV = real(iState_IIB(1:nI,1:nJ,iBlock:iBlock)))
    end do

    TypePosition = 'append'

  end subroutine save_plot_block
  !============================================================================

end program game_of_life
!==============================================================================

