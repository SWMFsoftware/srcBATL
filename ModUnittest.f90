module pass_cell

  use ModMpi
  use BATL_pass_cell, ONLY: message_pass_cell
  use BATL_size, ONLY: MaxDim, nDim, nDimAmr, iRatio, jRatio, kRatio, &
                       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG, nI, nJ, nK, &
                       nBlock, nIJK_D, iRatio_D
  use BATL_grid, ONLY: init_grid, create_grid, clean_grid, Xyz_DGB, &
                       CellSize_DB, CoordMin_DB, CoordMin_D, DomainSize_D
  use BATL_tree, ONLY: init_tree, set_tree_root, find_tree_node, show_tree, &
                       refine_tree_node, distribute_tree, clean_tree, &
                       Unused_B, DiLevelNei_IIIB, iNode_B
  use BATL_geometry, ONLY: IsCartesianGrid, IsRotatedCartesian, IsRoundCube, &
                           IsCylindricalAxis, IsSphericalAxis, IsLatitudeAxis, &
                           Lat_, Theta_, coord_to_xyz, init_geometry, z_, &
                           IsPeriodic_D, rot_to_cart, xyz_to_coord, coord_to_xyz
  use BATL_high_order, ONLY: restriction_high_order_reschange, &
                             prolongation_high_order_amr, &
                             prolongation_high_order_for_face_ghost, &
                             correct_face_ghost_for_fine_block, &
                             limit_interpolation, restriction_high_order_amr
  use BATL_mpi, ONLY: iComm, iProc
  use ModUtilities, ONLY: CON_stop
  use ModMpi, ONLY: MPI_ALLREDUCE
  use omp_lib
  use ModUtilities, ONLY: lower_case
  use ModNumConst, ONLY: cPi, cHalfPi, cTwoPi

  implicit none

contains
  !=============================================================================
  subroutine test_pass_cell

    integer, parameter:: MaxBlockTest = 200
    logical:: IsPeriodicTest_D(MaxDim) = .true.
    integer:: nRootTest_D(MaxDim) = [3, 3, 3]
    real   :: DomainMin_D(MaxDim) = [1.0, 10.0, 100.0]
    real   :: DomainMax_D(MaxDim) = [4.0, 40.0, 400.0]

    real   :: Tolerance = 1e-6

    integer, parameter:: nVar = nDim
    real, allocatable:: State_VGB(:, :, :, :, :)
    real, allocatable:: Scalar_GB(:, :, :, :)
    real, allocatable:: FineGridLocal_III(:, :, :)
    real, allocatable:: FineGridGlobal_III(:, :, :)
    real, allocatable:: XyzCorn_DGB(:, :, :, :, :)
    real :: CourseGridCell_III(iRatio, jRatio, kRatio)

    integer:: nWidth
    integer:: nProlongOrder
    integer:: nCoarseLayer
    integer:: iSendCorner, iRestrictFace
    logical:: DoSendCorner, DoRestrictFace

    real:: Xyz_D(MaxDim)
    integer:: iNode, iBlock, i, j, k, iMin, iMax, jMin, jMax, kMin, kMax, iDim
    integer:: iDir, jDir, kDir, Di, Dj, Dk

    integer:: iOp
    integer, parameter:: nOp = 2
    character(len=4):: NameOperator_I(nOp) = ["min", "max"]
    character(len=4):: NameOperator = "Min"
    real :: FineGridStep_D(MaxDim)
    integer:: iFG, jFG, kFG
    integer:: nFineCell
    integer:: iMpiOperator
    integer:: iError, iTest

    character(len=20):: NameGeometry

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'test_pass_cell'
    !--------------------------------------------------------------------------
    DoTest = iProc == 0

    if (DoTest) write (*, *) 'Starting ', NameSub

    call test_switches
    call test_scalar

    if (nDim == 1) RETURN !------------------------
    call test_non_cartesian

    if (nG < 3) RETURN
    if (nDim > nDimAmr) RETURN
    call test_high_order_cartesian
    call test_high_order_non_cartesian

  contains
    !==========================================================================
    subroutine test_switches

      !------------------------------------------------------------------------
      call init_tree(MaxBlockTest)
      call init_geometry(IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))
      call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
      call set_tree_root(nRootTest_D(1:nDim))

      call find_tree_node([0.5, 0.5, 0.5], iNode)
      if (DoTest) write (*, *) NameSub, ' middle node=', iNode
      call refine_tree_node(iNode)
      call distribute_tree(.true.)
      call create_grid

      if (DoTest) call show_tree(NameSub, .true.)

      allocate (State_VGB(nVar, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlockTest))

      do nProlongOrder = 1, 2; do nCoarseLayer = 1, 2; do nWidth = 1, nG

        ! Second order prolongation does not work with sending multiple
        ! coarse cell layers into the fine cells with their original values.
        if (nProlongOrder == 2 .and. nCoarseLayer == 2) CYCLE

        ! Cannot send more coarse layers than the number of ghost cell layers
        if (nCoarseLayer > nWidth) CYCLE

        if (DoTest) write (*, *) 'testing message_pass_cell with', &
          ' nProlongOrder=', nProlongOrder, &
          ' nCoarseLayer=', nCoarseLayer, &
          ' nWidth=', nWidth

        ! Set the range of ghost cells that should be set
        iMin = 1 - nWidth
        jMin = 1; if (nDim > 1) jMin = 1 - nWidth
        kMin = 1; if (nDim > 2) kMin = 1 - nWidth
        iMax = nI + nWidth
        jMax = nJ; if (nDim > 1) jMax = nJ + nWidth
        kMax = nK; if (nDim > 2) kMax = nK + nWidth

        do iSendCorner = 1, 2; do iRestrictFace = 1, 2

          DoSendCorner = iSendCorner == 2
          DoRestrictFace = iRestrictFace == 2

          ! Second order prolongation does not work with restricting face:
          ! the first order restricted cell cannot be used in the
          ! prolongation.
          if (DoRestrictFace .and. nProlongOrder == 2) CYCLE

          if (DoTest) write (*, *) 'testing message_pass_cell with', &
            ' DoSendCorner=', DoSendCorner, &
            ' DoRestrictFace=', DoRestrictFace

          State_VGB = 0.0

          do iBlock = 1, nBlock
            if (Unused_B(iBlock)) CYCLE
            State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) = &
              Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)
          end do

          call message_pass_cell(nVar, State_VGB, &
                                 nProlongOrderIn=nProlongOrder, &
                                 nCoarseLayerIn=nCoarseLayer, &
                                 nWidthIn=nWidth, &
                                 DoSendCornerIn=DoSendCorner, &
                                 DoRestrictFaceIn=DoRestrictFace)

          do iBlock = 1, nBlock
            if (Unused_B(iBlock)) CYCLE

            ! Loop through all cells including ghost cells
            do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI

              ! The filled in second order accurate ghost cell value
              ! should be the same as the coordinates of the cell center
              Xyz_D = Xyz_DGB(:, i, j, k, iBlock)

              ! Check that no info is sent in the non-used dimensions,
              ! i.e. for all iDim: nDim+1 < iDim < MaxDim
              if (i < iMin .or. i > iMax .or. &
                  j < jMin .or. j > jMax .or. &
                  k < kMin .or. k > kMax) then

                      do iDim = 1, nDim
                        if (abs(State_VGB(iDim, i, j, k, iBlock)) > 1e-6) then
                          write (*, *) 'Face should not be set: ', &
                            'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
                            iProc, iBlock, i, j, k, iDim, &
                            State_VGB(iDim, i, j, k, iBlock), &
                            Xyz_D(iDim)
                        end if
                      end do

                      CYCLE
                    end if

                    ! Get the direction vector
                    iDir = 0; if (i < 1) iDir = -1; if (i > nI) iDir = 1
                    jDir = 0; if (j < 1) jDir = -1; if (j > nJ) jDir = 1
                    kDir = 0; if (k < 1) kDir = -1; if (k > nK) kDir = 1

                    ! If nCoarseLayer==2 and DoSendCorner is true
                    ! the second ghost cells in the corner/edges
                    ! are not well defined (they may contain
                    ! the value coming from the first or second coarse cell).

                    if (nCoarseLayer == 2 .and. DoSendCorner .and. ( &
                        (i < 0 .or. i > nI + 1) .and. (jDir /= 0 .or. kDir /= 0) .or. &
                        (j < 0 .or. j > nJ + 1) .and. (iDir /= 0 .or. kDir /= 0) .or. &
                        (k < 0 .or. k > nK + 1) .and. (iDir /= 0 .or. jDir /= 0) &
                        )) CYCLE

                    ! if we do not send corners and edges, check that the
                    ! State_VGB in these cells is still the unset value
                    if (.not. DoSendCorner .and. ( &
                        iDir /= 0 .and. jDir /= 0 .or. &
                        iDir /= 0 .and. kDir /= 0 .or. &
                        jDir /= 0 .and. kDir /= 0)) then

                      do iDim = 1, nDim
                        if (abs(State_VGB(iDim, i, j, k, iBlock)) > 1e-6) then
                          write (*, *) 'corner/edge should not be set: ', &
                            'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
                            iProc, iBlock, i, j, k, iDim, &
                            State_VGB(iDim, i, j, k, iBlock), &
                            Xyz_D
                        end if
                      end do

                      CYCLE
                    end if

                    ! Shift ghost cell coordinate into periodic domain
                    Xyz_D = CoordMin_D + modulo(Xyz_D - CoordMin_D, DomainSize_D)

                    ! Calculate distance of ghost cell layer
                    Di = 0; Dj = 0; Dk = 0
                    if (i < 1 .and. iRatio == 2) Di = 2*i - 1
                    if (i > nI .and. iRatio == 2) Di = 2*(i - nI) - 1
                    if (j < 1 .and. jRatio == 2) Dj = 2*j - 1
                    if (j > nJ .and. jRatio == 2) Dj = 2*(j - nJ) - 1
                    if (k < 1 .and. kRatio == 2) Dk = 2*k - 1
                    if (k > nK .and. kRatio == 2) Dk = 2*(k - nK) - 1

                    if (DoRestrictFace .and. &
                        DiLevelNei_IIIB(iDir, jDir, kDir, iBlock) == -1) then
                      ! Shift coordinates if only 1 layer of fine cells
                      ! is averaged in the orthogonal direction
                      Xyz_D(1) = Xyz_D(1) - 0.25*Di*CellSize_DB(1, iBlock)
                      Xyz_D(2) = Xyz_D(2) - 0.25*Dj*CellSize_DB(2, iBlock)
                      Xyz_D(3) = Xyz_D(3) - 0.25*Dk*CellSize_DB(3, iBlock)
                    end if

                    if (nProlongOrder == 1 .and. &
                        DiLevelNei_IIIB(iDir, jDir, kDir, iBlock) == 1) then
                      ! Shift depends on the parity of the fine ghost cell
                      ! except when there is no AMR or multiple coarse cell
                      ! layers are sent in that direction
                      if (iRatio == 2 .and. (nCoarseLayer == 1 .or. iDir == 0)) &
                        Di = 2*modulo(i, 2) - 1
                      if (jRatio == 2 .and. (nCoarseLayer == 1 .or. jDir == 0)) &
                        Dj = 2*modulo(j, 2) - 1
                      if (kRatio == 2 .and. (nCoarseLayer == 1 .or. kDir == 0)) &
                        Dk = 2*modulo(k, 2) - 1

                      Xyz_D(1) = Xyz_D(1) + 0.5*Di*CellSize_DB(1, iBlock)
                      Xyz_D(2) = Xyz_D(2) + 0.5*Dj*CellSize_DB(2, iBlock)
                      Xyz_D(3) = Xyz_D(3) + 0.5*Dk*CellSize_DB(3, iBlock)
                    end if

                    do iDim = 1, nDim
                      if (abs(State_VGB(iDim, i, j, k, iBlock) - Xyz_D(iDim)) &
                          > Tolerance) then
                        write (*, *) 'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
                          iProc, iBlock, i, j, k, iDim, &
                          State_VGB(iDim, i, j, k, iBlock), &
                          Xyz_D(iDim)
                      end if
                    end do
                  end do; end do; end do
              end do

            end do; end do; end do
        end do; end do ! test parameters
      deallocate (State_VGB)

      call clean_grid
      call clean_tree
    end subroutine test_switches
    !==========================================================================

    subroutine test_scalar
      !------------------------ Test Scalar -----------------------------

      ! To test the message pass for the cell with min max operators we generate a
      ! fine uniform grid for the whole domain and transfer the cell values from
      ! the block cells to the cells on the fine grid. Then we gather all the data
      ! on the fine grid with the proper operator.
      ! We can then compare the values on the coresponding node after
      ! message_pass_cell_scalar is called with the fine grid values.

      ! rescale the domain to make indexing easier
      !------------------------------------------------------------------------
      DomainSize_D = iRatio_D*nRootTest_D*nIJK_D
      DomainMin_D = [0.0, 0.0, 0.0]
      DomainMax_D = DomainSize_D

      call init_tree(MaxBlockTest)
      call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
      call init_geometry(IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))
      call set_tree_root(nRootTest_D(1:nDim))

      call find_tree_node([0.5, 0.5, 0.5], iNode)
      call refine_tree_node(iNode)
      call distribute_tree(.true.)
      call create_grid

      ! Position of cell corners, for solving problems with round-off
      ! when getting fine grid positions
      allocate (XyzCorn_DGB(MaxDim, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlockTest))
      do iBlock = 1, nBlock
        if (Unused_B(iBlock)) CYCLE
        do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
            XyzCorn_DGB(:, i, j, k, iBlock) = &
              Xyz_DGB(:, i, j, k, iBlock) - &
              0.5*CellSize_DB(:, iBlock)* &
              [min(1, nI - 1), min(1, nJ - 1), min(1, nK - 1)]
          end do; end do; end do
      end do

      allocate (Scalar_GB(MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlockTest))
      Scalar_GB = -7777

      allocate (FineGridLocal_III( &
                nI*iRatio*nRootTest_D(1), &
                nJ*jRatio*nRootTest_D(2), &
                nK*kRatio*nRootTest_D(3)))
      allocate (FineGridGlobal_III( &
                (nI)*iRatio*nRootTest_D(1), &
                (nJ)*jRatio*nRootTest_D(2), &
                (nK)*kRatio*nRootTest_D(3)))

      nFineCell = ((nI)*iRatio*nRootTest_D(1))* &
                  ((nJ)*jRatio*nRootTest_D(2))* &
                  ((nK)*kRatio*nRootTest_D(3))

      FineGridStep_D = DomainSize_D/(DomainMax_D - DomainMin_D)

      do iOp = 1, nOp

        NameOperator = NameOperator_I(iOp)
        select case (NameOperator)
        case ("min")
          FineGridLocal_III(:, :, :) = 1.0e8
          FineGridGlobal_III(:, :, :) = 1.0e8
          iMpiOperator = MPI_MIN
        case ("max")
          FineGridLocal_III(:, :, :) = -1.0e8
          FineGridGlobal_III(:, :, :) = -1.0e8
          iMpiOperator = MPI_MAX
        case default
          call CON_stop(NameSub//': incorrect operator name')
        end select

        if (DoTest) write (*, *) 'testing message_pass_cell_scalar ', &
          'with operator= ', NameOperator

        do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
              Scalar_GB(i, j, k, iBlock) = iNode_B(iBlock) + &
                                           sum(CoordMin_DB(:, iBlock) + &
                                               ([i, j, k])*CellSize_DB(:, iBlock))
            end do; end do; end do
        end do

        do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          do k = 1, nK; do j = 1, nJ; do i = 1, nI

              iFG = nint(XyzCorn_DGB(1, i, j, k, iBlock)*FineGridStep_D(1)) + 1
              jFG = nint(XyzCorn_DGB(2, i, j, k, iBlock)*FineGridStep_D(2)) + 1
              kFG = nint(XyzCorn_DGB(3, i, j, k, iBlock)*FineGridStep_D(3)) + 1

              FineGridLocal_III(iFG, jFG, kFG) = Scalar_GB(i, j, k, iBlock)
            end do; end do; end do
        end do

        call message_pass_cell(Scalar_GB, &
                               nProlongOrderIn=1, nCoarseLayerIn=2, &
                               DoSendCornerIn=.true., DoRestrictFaceIn=.false., &
                               NameOperatorIn=NameOperator_I(iOp))

        call MPI_ALLREDUCE(FineGridLocal_III(1, 1, 1), &
                           FineGridGlobal_III(1, 1, 1), &
                           nFineCell, MPI_REAL, iMpiOperator, iComm, iError)

        ! making sure that we have the center cell along the x=0 side
        ! so the boundary are not tested.
        call find_tree_node([0.0, 0.5, 0.5], iNode)

        do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          if (iNode_B(iBlock) == iNode) then
            do k = MinK, MaxK; do j = MinJ, MaxJ; do i = 1, MaxI

              iFG = nint(XyzCorn_DGB(1, i, j, k, iBlock)*FineGridStep_D(1)) + 1
              jFG = nint(XyzCorn_DGB(2, i, j, k, iBlock)*FineGridStep_D(2)) + 1
              kFG = nint(XyzCorn_DGB(3, i, j, k, iBlock)*FineGridStep_D(3)) + 1
              ! copy cells that are inside the course grid cell
              CourseGridCell_III = FineGridGlobal_III( &
                                   iFG:iFG + min(1, iRatio - 1), &
                                   jFG:jFG + min(1, jRAtio - 1), &
                                   kFG:kFG + min(1, kRatio - 1))
              select case (NameOperator)
              case ("min")
                if (Scalar_GB(i,j,k,iBlock) /= minval(CourseGridCell_III))then
                  write (*, *) "Error for operator, iNode,  iBlock= ", &
                    NameOperator, iNode_B(iBlock), iBlock, ", value=", &
                    minval(CourseGridCell_III), &
                    " should be ", Scalar_GB(i, j, k, iBlock), "index : ", &
                    i, j, k, " ::", iFG, jFG, kFG
                end if
               case ("max")
                if (Scalar_GB(i,j,k,iBlock) /= maxval(CourseGridCell_III)) then
                  write (*, *) "Error for operator, iNode,  iBlock= ", &
                    NameOperator, iNode_B(iBlock), iBlock, ",value=", &
                    maxval(CourseGridCell_III), &
                    " should be ", Scalar_GB(i, j, k, iBlock), "index : ", &
                    i, j, k, " ::", iFG, jFG, kFG
                end if
               end select

            end do; end do; end do
        end if
      end do

    end do
    deallocate (Scalar_GB, FineGridLocal_III, FineGridGlobal_III, XyzCorn_DGB)
    call clean_grid
    call clean_tree

  end subroutine test_scalar
  !==========================================================================
  subroutine test_non_cartesian

    !------------------------------------------------------------------------
    do iTest = 1, 6
      ! The code is quite inaccurate for partial AMR across the pole
      if (nDimAmr < nDim .and. iTest > 3) EXIT
      call init_tree(MaxBlockTest)
      ! Do not test ghost cells in the radial direction
      iMin = 1; iMax = nI
      jMin = MinJ; jMax = MaxJ
      kMin = MinK; kMax = MaxK
      select case (iTest)
      case (1, 4)
        NameGeometry = 'cylindrical'
        ! 0 < r < 10, 0 < phi < 360deg, -5 < z < 5
        DomainMin_D = [0.0, 0.0, -5.0]
        DomainMax_D = [8.0, cTwoPi, +5.0]
        IsPeriodicTest_D = [.false., .true., .true.]
        ! There must be an even number of root blocks in the phi direction
        ! There are 3 root blocks in z so that we can refine the middle
        ! and avoid issues of periodicity in the testing
        nRootTest_D = [2, 4, 3]
        ! Test ghost cells at rMin
        iMin = MinI
      case (2, 5)
        if (nDim < 3) CYCLE
        NameGeometry = 'spherical'
        ! 1 < r < 9, 0 < theta < 180deg, 0 < phi < 360deg
        DomainMin_D = [1.0, 0.0, 0.0]
        DomainMax_D = [9.0, cPi, cTwoPi]
        IsPeriodicTest_D = [.false., .false., .true.]
        ! There must be an even number of root blocks in the phi direction
        ! There are 3 root blocks in r so that we can refine the middle
        ! and avoid issues at inner and outer radial boundaries
          nRootTest_D = [3, 2, 4]

      case (3, 6)
        if (nDim < 3) CYCLE
        NameGeometry = 'rlonlat'
        ! 1 < r < 9, 0 < phi < 360deg, -90 < lat < 90
        DomainMin_D = [1.0, 0.0, -cHalfPi]
        DomainMax_D = [9.0, cTwoPi, cHalfPi]
        IsPeriodicTest_D = [.false., .true., .false.]
        ! There must be an even number of root blocks in the phi direction
        ! There are 3 root blocks in r so that we can refine the middle
        ! and avoid issues at inner and outer radial boundaries
        nRootTest_D = [3, 4, 2]
      end select
      DomainSize_D = DomainMax_D - DomainMin_D
      if (DoTest) then
        if (iTest <= 3) write (*, *) &
          'testing message_pass_cell across '//trim(NameGeometry)// &
          ' pole'
        if (iTest >= 4) write (*, *) &
          'testing message_pass_cell across '//trim(NameGeometry)// &
          ' pole with resolution change'
      end if
      call init_geometry(NameGeometry, &
                         IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))
      call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim), &
                     UseDegreeIn=.false.)
      call set_tree_root(nRootTest_D(1:nDim))
      if (any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
        write (*, *) NameSub, ': IsPeriodic_D=', IsPeriodic_D(1:nDim), &
        ' should agree with ', IsPeriodicTest_D(1:nDim)
      if (iTest > 3) then
        ! Test with refined grid
        if (iTest == 4) then
          ! refine node next to r=0 axis but middle in Z direction
          call refine_tree_node(9)
        else
          ! refine root nodes at min and max (theta/lat/phi) coordinates
          ! but middle in the R direction
          call refine_tree_node(2)
          call refine_tree_node(23)
        end if
        ! Restriction is not linear so there is truncation error
        Tolerance = 0.15
      else
        ! For tests with no AMR the error is round-off only
        Tolerance = 1e-6
      end if

      call distribute_tree(.true.)
      call create_grid
      allocate (State_VGB(nVar, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlockTest))
      State_VGB = 0.0
      do iBlock = 1, nBlock
        if (Unused_B(iBlock)) CYCLE
        State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) = &
          Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)
      end do
      ! Second order
      call message_pass_cell(nVar, State_VGB)
      do iBlock = 1, nBlock
        if (Unused_B(iBlock)) CYCLE
        ! Loop through all cells including ghost cells
        do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
            ! The filled in second order accurate ghost cell value
            ! should be the same as the coordinates of the cell center
            Xyz_D = Xyz_DGB(:, i, j, k, iBlock)
            ! For 3D cylindrical Z coordinate is periodic
            if ((iTest == 1 .or. iTest == 4) .and. nDim == 3) &
              Xyz_D(z_) = DomainMin_D(z_) &
                          + modulo(Xyz_D(z_) - DomainMin_D(z_), DomainSize_D(z_))
            do iDim = 1, nDim
              if (abs(State_VGB(iDim, i, j, k, iBlock) - Xyz_D(iDim)) &
                  /abs(Xyz_D(iDim)) > Tolerance) then
                write (*, *) 'iProc,iBlock,i,j,k,iDim,State,Xyz=', &
                  iProc, iBlock, i, j, k, iDim, &
                  State_VGB(iDim, i, j, k, iBlock), &
                  Xyz_D(iDim)
              end if
            end do
          end do; end do; end do
      end do

      deallocate (State_VGB)

      call clean_grid
      call clean_tree
    end do

    ! In previous do loop, it may cycle some loops without cleaning.
    call clean_grid
    call clean_tree

  end subroutine test_non_cartesian
  !==========================================================================

  subroutine test_high_order_cartesian
    real    :: ExactSolution, Error, ErrorTotal
    integer :: iCount, nCount, nRefineNode, iRefinement, nRefinement
    integer :: iNode_I(8)
    integer :: iNode1_I(8)
    logical :: DoTestMeOnly = .false.
    integer :: nPoly = 3

    !------------------------------------------------------------------------
    if (nDim == 2) then
      nCount = 16; nRefineNode = 4
    else
      nCount = 256; nRefineNode = 8
    endif
    if (nDimAmr < nDim) RETURN
    iMin = MinI; iMax = MaxI
    jMin = MinJ; jMax = MaxJ
    kMin = MinK; kMax = MaxK
    NameGeometry = 'cartesian'
    DomainMin_D = [0.0, 0.0, 0.0]
    DomainMax_D = [8.0, 8.0, 8.0]
    DomainSize_D = DomainMax_D - DomainMin_D

    IsPeriodicTest_D = [.true., .true., .true.]
    nRootTest_D = [4, 4, 4]
    nRefinement = 2
    
    do iRefinement = 1, nRefinement
      ! iRefinement = 1: 1 level refine
      ! iRefinement = 2: 2 level refine
      if (DoTest) then
        write (*, *) &
          'testing message_pass_cell across '//trim(NameGeometry)// &
          ' with high resolution change with refinement level =', &
          iRefinement
      end if
      if (nDim == 2) then
        if (iRefinement == 1) then
          iNode_I = [6, 7, 10, 11, -1, -1, -1, -1]
        else
          iNode1_I = [6, 7, 10, 11, -1, -1, -1, -1]
          iNode_I = [20, 23, 26, 29, -1, -1, -1, -1]
        endif
      else ! 3D
        if (iRefinement == 1) then
          iNode_I = [22, 23, 26, 27, 38, 39, 42, 43]
        else
          iNode1_I = [22, 23, 26, 27, 38, 39, 42, 43]
          iNode_I = [72, 79, 86, 93, 100, 107, 114, 121]
        endif
      endif
      do iCount = 0, nCount - 1
        if (DoTestMeOnly) then
          write (*, *) ''
          write (*, *) 'test_high_order iCount = ', iCount
        endif
        call init_tree(MaxBlockTest)
        call init_geometry(NameGeometry, &
                           IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))

          call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim), &
                         UseDegreeIn=.false.)
          call set_tree_root(nRootTest_D(1:nDim))

          if (any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
            write (*, *) NameSub, ': IsPeriodic_D=', IsPeriodic_D(1:nDim), &
            ' should agree with ', IsPeriodicTest_D(1:nDim)

          if (iRefinement == 2) then
            do iNode = 1, nRefineNode
              call refine_tree_node(iNode1_I(iNode))
            enddo
          endif

          do iNode = 1, nRefineNode
            if (btest(iCount, iNode - 1)) then
              call refine_tree_node(iNode_I(iNode))
              if (DoTestMeOnly) &
                write (*, *) 'iNode IsRefined:', iNode_I(iNode), 'TRUE'
            else
              if (DoTestMeOnly) &
                write (*, *) 'iNode IsRefined:', iNode_I(iNode), 'FALSE'
            endif
          enddo

          Tolerance = 5e-15

          call distribute_tree(.true.)
          call create_grid

          allocate ( &
            State_VGB(nVar, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlockTest))
          State_VGB = 0

          do iBlock = 1, nBlock
            if (Unused_B(iBlock)) CYCLE
            do i = 1, nI; do j = 1, nJ; do k = 1, nK
                State_VGB(1, i, j, k, iBlock) = &
                  exact_solution(Xyz_DGB(:, i, j, k, iBlock), nPolyIn=nPoly)
              enddo; enddo; enddo
          end do

          call message_pass_cell(nVar, nG, State_VGB, nProlongOrderIn=1, &
                                 nCoarseLayerIn=2, DoResChangeOnlyIn=.false., &
                                 UseHighResChangeIn=.true.)
          ErrorTotal = 0
          do iBlock = 1, nBlock
            if (Unused_B(iBlock)) CYCLE

            ! Loop through all cells including ghost cells
            do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
                Xyz_D = Xyz_DGB(:, i, j, k, iBlock)
                if (.not. (all(Xyz_D(1:nDim) < DomainMax_D(1:nDim)) &
                           .and. all(Xyz_D(1:nDim) > DomainMin_D(1:nDim)))) then
                  CYCLE
                endif

                ExactSolution = exact_solution(Xyz_D, nPolyIn=nPoly)
                Error = abs(ExactSolution - State_VGB(1, i, j, k, iBlock))
                ErrorTotal = ErrorTotal + Error
                if (abs(Error)/abs(ExactSolution) > Tolerance) &
                  then
                  write (*, *) &
                    'iProc,iNode,i,j,k,x,y,z,', &
                    'state,exact-solution,error,relative-error='
                  write (*, '(5I5,7e20.12)') &
                    iProc, iNode_B(iBlock), i, j, k, &
                    Xyz_D, State_VGB(1, i, j, k, iBlock), &
                    ExactSolution, Error, abs(Error)/abs(ExactSolution)

                end if
              end do; end do; end do
          end do
          if (DoTestMeOnly) then
            write (*, *) 'Refine level = ', iRefinement
            write (*, *) 'Total error  = ', ErrorTotal
          endif
          deallocate (State_VGB)

          call clean_grid
          call clean_tree
        enddo ! iCount
      enddo ! iRefinement

    end subroutine test_high_order_cartesian
!==========================================================================

    subroutine test_high_order_non_cartesian
      real    :: ErrorTotal, ExactSolution, Error
      real    :: Xyz1_D(3), XyzGeneral_D(3)
      integer :: nPoly
      !------------------------------------------------------------------------
      do iTest = 1, 2

        ! The code is quite inaccurate for partial AMR across the pole
        if (nDimAmr < nDim .and. iTest > 3) EXIT

        call init_tree(MaxBlockTest)

        ! Do not test ghost cells in the radial direction
        iMin = 1; iMax = nI
        jMin = MinJ; jMax = MaxJ
        kMin = MinK; kMax = MaxK

        select case (iTest)
        case (1)

          NameGeometry = 'cylindrical'

          ! 2 < r < 8, 0 < phi < 360deg, -5 < z < 5
          DomainMin_D = [2.0, 0.0, -5.0]
          DomainMax_D = [8.0, cTwoPi, +5.0]
          IsPeriodicTest_D = [.false., .true., .true.]

          ! There must be an even number of root blocks in the phi direction
          ! There are 3 root blocks in z so that we can refine the middle
          ! and avoid issues of periodicity in the testing
          nRootTest_D = [3, 4, 3]

          ! Test ghost cells at rMin
          iMin = MinI

        case (2)
          NameGeometry = 'rotatedcartesian'

          DomainMin_D = [0.0, 0.0, 0.0]
          DomainMax_D = [6.0, 6.0, 6.0]

          IsPeriodicTest_D = [.false., .false., .false.]
          nRootTest_D = [3, 3, 3]
        end select
        DomainSize_D = DomainMax_D - DomainMin_D

        if (DoTest) then
          write (*, *) &
            'testing message_pass_cell across '//trim(NameGeometry)// &
            ' with high resolution change'
        end if

        call init_geometry(NameGeometry, &
                           IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))

        call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim), &
                       UseDegreeIn=.false.)
        call set_tree_root(nRootTest_D(1:nDim))

        if (any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
          write (*, *) NameSub, ': IsPeriodic_D=', IsPeriodic_D(1:nDim), &
          ' should agree with ', IsPeriodicTest_D(1:nDim)

        if (iTest == 1) then
          ! refine node next to r=0 axis but middle in Z direction
          if (nDim == 2) then
            call refine_tree_node(2)
          elseif (nDim == 3) then
            call refine_tree_node(17)
          endif
        elseif (iTest == 2) then
          if (nDim == 2) then
            call refine_tree_node(5)
          elseif (nDim == 3) then
            call refine_tree_node(14)
          endif
        end if

        if (iTest == 2) then
          Tolerance = 1e-14
          nPoly = 3
        else
          Tolerance = 7e-3
          nPoly = 1
        endif

        call distribute_tree(.true.)
        call create_grid

        allocate (State_VGB(nVar, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlockTest))
        State_VGB = 0

        do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          do i = 1, nI; do j = 1, nJ; do k = 1, nK
              State_VGB(1, i, j, k, iBlock) = &
                exact_solution(Xyz_DGB(:, i, j, k, iBlock), nPolyIn=nPoly)
            enddo; enddo; enddo
        end do

        call message_pass_cell(nVar, nG, State_VGB, nProlongOrderIn=1, &
                               nCoarseLayerIn=2, DoResChangeOnlyIn=.false., &
                               UseHighResChangeIn=.true.)

        ! Second order
        ! call message_pass_cell(nVar, State_VGB)

        ErrorTotal = 0
        do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          ! Loop through all cells including ghost cells
          do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax

              Xyz_D = Xyz_DGB(:, i, j, k, iBlock)

              if (iTest == 2) then
                Xyz1_D = rot_to_cart(Xyz_D)
              else
                call Xyz_to_coord(Xyz_D, Xyz1_D)
              endif

              if (.not. (all(Xyz1_D(1:nDim) < DomainMax_D(1:nDim)) &
                         .and. all(Xyz1_D(1:nDim) > DomainMin_D(1:nDim)))) then
                CYCLE
              endif

              ExactSolution = exact_solution(Xyz_D, nPolyIn=nPoly)
              Error = abs(ExactSolution - State_VGB(1, i, j, k, iBlock))
              ErrorTotal = ErrorTotal + Error/abs(ExactSolution)

              if (abs(Error)/abs(ExactSolution) > Tolerance) then
                write (*, *) &
                  'iProc,iNode,i,j,k,x,y,z,', &
                  'state,exact-solution,error,relative-error='
                write (*, '(5I5,7e20.12)') &
                  iProc, iNode_B(iBlock), i, j, k, &
                  Xyz_D, State_VGB(1, i, j, k, iBlock), &
                  ExactSolution, Error, abs(Error)/abs(ExactSolution)
                call Xyz_to_coord(Xyz_D, XyzGeneral_D)
                write (*, *) 'Xyz general = ', XyzGeneral_D
                write (*, *) ''
              end if

            end do; end do; end do
        end do

        deallocate (State_VGB)

        call clean_grid
        call clean_tree
      end do ! iTest

    end subroutine test_high_order_non_cartesian
!==========================================================================

    real function exact_solution(Xyz_D, nPolyIn)
      real, intent(in):: Xyz_D(3)
      integer, optional, intent(in) :: nPolyIn
      integer:: nPoly
      real:: x, y, z

      !------------------------------------------------------------------------
      nPoly = 4
      if (present(nPolyIn)) nPoly = nPolyIn

      x = Xyz_D(1)
      y = Xyz_D(2)
      z = Xyz_D(3)

      ! exact_solution = 4.0+sin(y*2*cpi/8.0)+cos(z*2*cpi/8.0)+sin(x*2*cpi/8.0)
      select case (nPoly)
      case (4)
        exact_solution = x**4 + y**4 + z**4 + x**3*y &
                         + y**3*x + x**3*z + x*z**3 + y**3*z + y*z**3 &
                         + x**2*y**2 + x**2*z**2 + y**2*z**2 &
                         + y*z*x**2 + x*z*y**2 + x*y*z**2
      case (3)
        exact_solution = x**3 + y**3 + z**3 + x*y*z + x**2*y + x*y**2 + &
                         x**2*z + x*z**2 + y**2*z + y*z**2
      case (1)
        exact_solution = x + y + z
      end select
    end function exact_solution
!==========================================================================

  end subroutine test_pass_cell
!============================================================================

end module pass_cell

module pass_face

  use BATL_pass_face
  use BATL_mpi, ONLY: iProc
  use BATL_size, ONLY: MaxDim, nDim, nI, nJ, nK, nBlock
  use BATL_tree, ONLY: init_tree, set_tree_root, find_tree_node, &
                       refine_tree_node, distribute_tree, show_tree, clean_tree, &
                       Unused_B, UseTimeLevel, iTimeLevel_A, nNode, &
                       iTree_IA, Level_, iNode_B, di_level_nei
  use BATL_grid, ONLY: init_grid, create_grid, clean_grid, Xyz_DGB, CellFace_DB
  use BATL_geometry, ONLY: init_geometry
  use ModUtilities, ONLY: CON_stop

  implicit none

contains

  subroutine test_pass_face

    integer, parameter:: MaxBlockTest = 50
    integer, parameter:: nRootTest_D(MaxDim) = [3, 3, 3]
    logical, parameter:: IsPeriodicTest_D(MaxDim) = .true.
    real, parameter:: DomainMin_D(MaxDim) = [1.0, 10.0, 100.0]
    real, parameter:: DomainMax_D(MaxDim) = [4.0, 40.0, 400.0]

    real, parameter:: Tolerance = 1e-6

    integer, parameter:: nVar = nDim
    real, allocatable, dimension(:, :, :, :, :):: &
      Flux_VFD, Flux_VXB, Flux_VYB, Flux_VZB

    integer:: iResChangeOnly
    logical:: DoResChangeOnly

    integer:: iNode, iBlock, i, j, k, iDim, DiLevel
    integer:: iStage, nStage, iTimeLevel
    real:: Flux, FluxGood, FluxUnset = -77.0

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'test_pass_face'
    !--------------------------------------------------------------------------
    DoTest = iProc == 0

    if (DoTest) write (*, *) 'Starting ', NameSub

    call init_tree(MaxBlockTest)
    call init_geometry(IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))
    call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
    call set_tree_root(nRootTest_D(1:nDim))

    call find_tree_node([0.5, 0.5, 0.5], iNode)
    if (DoTest) write (*, *) NameSub, ' middle node=', iNode
    call refine_tree_node(iNode)
    call distribute_tree(.true.)
    call create_grid

    allocate ( &
      Flux_VFD(nVar, nI + 1, nJ + 1, nK + 1, nDim), &
      Flux_VXB(nVar, nJ, nK, 2, MaxBlockTest), &
      Flux_VYB(nVar, nI, nK, 2, MaxBlockTest), &
      Flux_VZB(nVar, nI, nJ, 2, MaxBlockTest), &
      iTimeLevel_A(nNode))

    do iNode = 1, nNode
      iTimeLevel_A(iNode) = modulo(iNode, 4) + 5*iTree_IA(Level_, iNode)
    end do

    if (DoTest) &
      call show_tree(NameSub, DoShowNei=.true.) ! DoShowTimeLevel=.true.)

    do iResChangeOnly = 1, 2

      DoResChangeOnly = iResChangeOnly == 1
      UseTimeLevel = .not. DoResChangeOnly

      if (DoTest) write (*, *) 'testing message_pass_face with', &
        ' DoResChangeOnly, UseTimeLevel=', DoResChangeOnly, UseTimeLevel

      Flux_VFD = 0.0
      Flux_VXB = 0.0
      Flux_VYB = 0.0
      Flux_VZB = 0.0

      do iBlock = 1, nBlock
        if (Unused_B(iBlock)) CYCLE

        ! Set Flux_VFD based on the coordinates
        call set_flux

        if (UseTimeLevel) then
          ! Pretend subsycling nStage times
          nStage = iTimeLevel_A(iNode_B(iBlock))
          do iStage = 1, nStage
            call store_face_flux(iBlock, nVar, Flux_VFD, &
                                 Flux_VXB, Flux_VYB, Flux_VZB, &
                                 DoResChangeOnlyIn=.false., &
                                 DoStoreCoarseFluxIn=.true., &
                                 DtIn=1.0/nStage)
          end do
        else
          call store_face_flux(iBlock, nVar, Flux_VFD, &
                               Flux_VXB, Flux_VYB, Flux_VZB, &
                               DoResChangeOnlyIn=.true., &
                               DoStoreCoarseFluxIn=.true.)
        end if

        iTimeLevel = iTimeLevel_A(iNode_B(iBlock))

        ! Set the fluxes to the unset value except on the fine side
        if (di_level_nei(-1, 0, 0, iBlock, DoResChangeOnly) < 1) &
          Flux_VXB(:, :, :, 1, iBlock) = FluxUnset

        if (di_level_nei(+1, 0, 0, iBlock, DoResChangeOnly) < 1) &
          Flux_VXB(:, :, :, 2, iBlock) = FluxUnset

        if (nDim > 1) then
          if (di_level_nei(0, -1, 0, iBlock, DoResChangeOnly) < 1) &
            Flux_VYB(:, :, :, 1, iBlock) = FluxUnset

          if (di_level_nei(0, +1, 0, iBlock, DoResChangeOnly) < 1) &
            Flux_VYB(:, :, :, 2, iBlock) = FluxUnset

        end if

        if (nDim > 2) then
          if (di_level_nei(0, 0, -1, iBlock, DoResChangeOnly) < 1) &
            Flux_VZB(:, :, :, 1, iBlock) = FluxUnset

          if (di_level_nei(0, 0, +1, iBlock, DoResChangeOnly) < 1) &
            Flux_VZB(:, :, :, 2, iBlock) = FluxUnset

        end if

      end do

      ! Message pass from fine to coarse grid/time levels
      ! Zero out the fine side.
      call message_pass_face(nVar, Flux_VXB, Flux_VYB, Flux_VZB, &
                             DoSubtractIn=.false., DoResChangeOnlyIn=DoResChangeOnly)

      do iBlock = 1, nBlock
        if (Unused_B(iBlock)) CYCLE

        ! Calculate the coordinate based Flux_VFD again
        call set_flux

        ! Check min X face
        DiLevel = di_level_nei(-1, 0, 0, iBlock, DoResChangeOnly)
        do k = 1, nK; do j = 1, nJ; do iDim = 1, nDim
            Flux = Flux_VXB(iDim, j, k, 1, iBlock)
            FluxGood = flux_good(DiLevel, 1, Flux_VFD(iDim, 1, j, k, 1))

            if (abs(Flux - FluxGood) > Tolerance) then
              write (*, *) 'Error at min X face: ', &
                'iNode,DiLevel,j,k,iDim,Flux,Good=', &
                iNode_B(iBlock), DiLevel, j, k, iDim, Flux, FluxGood, &
                CellFace_DB(1, iBlock)
            end if
          end do; end do; end do

        ! Check max X face
        DiLevel = di_level_nei(+1, 0, 0, iBlock, DoResChangeOnly)
        do k = 1, nK; do j = 1, nJ; do iDim = 1, nDim
            Flux = Flux_VXB(iDim, j, k, 2, iBlock)
            FluxGood = flux_good(DiLevel, 1, Flux_VFD(iDim, nI + 1, j, k, 1))
            if (abs(Flux - FluxGood) > Tolerance) &
              write (*, *) 'Error at max X face: ', &
              'iNode,DiLevel,j,k,iDim,Flux,Good=', &
              iNode_B(iBlock), DiLevel, j, k, iDim, Flux, FluxGood
          end do; end do; end do

        if (nDim > 1) then
          ! Check min Y face
          DiLevel = di_level_nei(0, -1, 0, iBlock, DoResChangeOnly)
          do k = 1, nK; do i = 1, nI; do iDim = 1, nDim
              Flux = Flux_VYB(iDim, i, k, 1, iBlock)
              FluxGood = flux_good(DiLevel, 2, Flux_VFD(iDim, i, 1, k, 2))
              if (abs(Flux - FluxGood) > Tolerance) &
                write (*, *) 'Error at min Y face: ', &
                'iNode,DiLevel,i,k,iDim,Flux,Good=', &
                iNode_B(iBlock), DiLevel, i, k, iDim, Flux, FluxGood
            end do; end do; end do

          ! Check max Y face
          DiLevel = di_level_nei(0, +1, 0, iBlock, DoResChangeOnly)
          do k = 1, nK; do i = 1, nI; do iDim = 1, nDim
              Flux = Flux_VYB(iDim, i, k, 2, iBlock)
              FluxGood = flux_good(DiLevel, 2, Flux_VFD(iDim, i, nJ + 1, k, 2))
              if (abs(Flux - FluxGood) > Tolerance) &
                write (*, *) 'Error at max Y face: ', &
                'iNode,DiLevel,i,k,iDim,Flux,Good=', &
                iNode_B(iBlock), DiLevel, i, k, iDim, Flux, FluxGood
            end do; end do; end do
        end if

        if (nDim > 2) then
          ! Check min Z face
          DiLevel = di_level_nei(0, 0, -1, iBlock, DoResChangeOnly)
          do j = 1, nJ; do i = 1, nI; do iDim = 1, nDim
              Flux = Flux_VZB(iDim, i, j, 1, iBlock)
              FluxGood = flux_good(DiLevel, 3, Flux_VFD(iDim, i, j, 1, 3))
              if (abs(Flux - FluxGood) > Tolerance) &
                write (*, *) 'Error at min Z face: ', &
                'iNode,DiLevel,i,j,iDim,Flux,Good=', &
                iNode_B(iBlock), DiLevel, i, j, iDim, Flux, FluxGood
            end do; end do; end do

          ! Check max Z face
          DiLevel = di_level_nei(0, 0, +1, iBlock, DoResChangeOnly)
          do j = 1, nJ; do i = 1, nI; do iDim = 1, nDim
              Flux = Flux_VZB(iDim, i, j, 2, iBlock)
              FluxGood = flux_good(DiLevel, 3, Flux_VFD(iDim, i, j, nK + 1, 3))
              if (abs(Flux - FluxGood) > Tolerance) &
                write (*, *) 'Error at max Z face: ', &
                'iNode,DiLevel,i,j,iDim,Flux,Good=', &
                iNode_B(iBlock), DiLevel, i, j, iDim, Flux, FluxGood
            end do; end do; end do
        end if

      end do ! iBlock

    end do ! test parameters
    deallocate (Flux_VFD, Flux_VXB, Flux_VYB, Flux_VZB)

    call clean_grid
    call clean_tree

  contains
    !==========================================================================
    subroutine set_flux

      ! Fill in Flux_VFD with coordinates of the face center

      !------------------------------------------------------------------------
      Flux_VFD(:, :, 1:nJ, 1:nK, 1) = 0.5*CellFace_DB(1, iBlock)* &
                                      (Xyz_DGB(1:nDim, 0:nI, 1:nJ, 1:nK, iBlock) &
                                       + Xyz_DGB(1:nDim, 1:nI + 1, 1:nJ, 1:nK, iBlock))

      if (nDim > 1) &
        Flux_VFD(:, 1:nI, :, 1:nK, 2) = 0.5*CellFace_DB(2, iBlock)* &
                                        (Xyz_DGB(1:nDim, 1:nI, 0:nJ, 1:nK, iBlock) &
                                         + Xyz_DGB(1:nDim, 1:nI, 1:nJ + 1, 1:nK, iBlock))

      if (nDim > 2) &
        Flux_VFD(:, 1:nI, 1:nJ, :, 3) = 0.5*CellFace_DB(3, iBlock)* &
                                        (Xyz_DGB(1:nDim, 1:nI, 1:nJ, 0:nK, iBlock) &
                                         + Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK + 1, iBlock))

    end subroutine set_flux
    !==========================================================================

    real function flux_good(DiLevel, iDir, XyzCellFace)

      integer, intent(in):: DiLevel, iDir
      real, intent(in):: XyzCellFace

      ! Calculate the correct flux solution in direction iDir based on the
      ! grid/time level change DiLevel at the given face.
      ! The coordinate XyzCellFace is the correct answer on the coarse side.
      ! Zero is the correct value on the fine side.
      ! FluxUnset is the correct value everywhere else.

      real:: CellFace
      !------------------------------------------------------------------------
      select case (DiLevel)
      case (0)
        flux_good = FluxUnset
      case (1)
        ! The flux is zeroed out on the fine side
        flux_good = 0.0
      case (-1)
        ! The flux is set to coordinate value on the coarse side
        flux_good = XyzCellFace

        ! Fix periodic coordinates. Flux = Coordinate*CellArea
        if (iDim == iDir) then
          CellFace = CellFace_DB(iDim, iBlock)
          if (abs(XyzCellFace/CellFace - DomainMin_D(iDim)) < Tolerance) &
            flux_good = DomainMax_D(iDim)*CellFace
          if (abs(XyzCellFace/CellFace - DomainMax_D(iDim)) < Tolerance) &
            flux_good = DomainMin_D(iDim)*CellFace
        endif
      case default
        call CON_stop(NameSub//' invalid value for DiLevel')
      end select

    end function flux_good
    !==========================================================================

  end subroutine test_pass_face
  !============================================================================

end module pass_face

module amr

  use BATL_amr, ONLY: init_amr, do_amr
  use BATL_tree, ONLY: iAmrChange_B, &
                       AmrRemoved_, AmrMoved_, AmrRefined_, AmrCoarsened_
  use BATL_mpi, ONLY: iProc, nProc, barrier_mpi
  use BATL_size, ONLY: MaxDim, nDim, nIJK_D, iDimAmr_D, &
                       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, nBlock, &
                       iRatio, jRatio, kRatio
  use BATL_tree, ONLY: init_tree, set_tree_root, refine_tree_node, &
                       coarsen_tree_node, distribute_tree, move_tree, show_tree, clean_tree, &
                       iProcNew_A, Unused_B, nNode, iNode_B, iTree_IA, Proc_, Block_, &
                       Child1_, unset_, nNode
  use BATL_grid, ONLY: init_grid, create_grid, clean_grid, Xyz_DGB, CellSize_DB
  use BATL_geometry, ONLY: init_geometry
  use ModUtilities, ONLY: CON_stop

  implicit none

  public :: test_amr

  ! Private data for test_amr
  integer, parameter:: nExtraData = 2
  real, allocatable:: ExtraData_IB(:, :)

contains
  !=============================================================================
  subroutine test_amr

    integer, parameter:: MaxBlockTest = 100
    integer, parameter:: nRootTest_D(MaxDim) = [3, 3, 3]
    logical, parameter:: IsPeriodicTest_D(MaxDim) = .true.
    real, parameter:: DomainMin_D(MaxDim) = [1.0, 10.0, 100.0]
    real, parameter:: DomainMax_D(MaxDim) = [10.0, 100.0, 1000.0]
    real, parameter:: DomainSize_D(MaxDim) = DomainMax_D - DomainMin_D

    integer, parameter:: nVar = nDim
    real, allocatable:: State_VGB(:, :, :, :, :), Dt_B(:)
    real, allocatable:: TestState_VC(:, :, :, :)
    logical, allocatable:: Used_GB(:, :, :, :)
    integer, allocatable:: iEffectedNode_A(:)
    integer:: iBlock, iDim, iNode, iVar, i, j, k
    integer:: iChild

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'test_amr'
    !--------------------------------------------------------------------------
    DoTest = iProc == 0

    if (DoTest) write (*, *) 'Starting ', NameSub

    call init_tree(MaxBlockTest)
    call init_geometry(IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))
    call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
    call set_tree_root(nRootTest_D(1:nDim))
    call distribute_tree(.true.)
    call create_grid
    call init_amr

    if (DoTest) call show_tree('after create_grid')

    allocate (State_VGB(nVar, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlockTest), &
              Dt_B(MaxBlockTest))

    do iBlock = 1, nBlock
      if (Unused_B(iBlock)) CYCLE
      State_VGB(:, :, :, :, iBlock) = Xyz_DGB(1:nDim, :, :, :, iBlock)
      ! set the time step to the cell size in the first AMR direction
      iDim = iDimAmr_D(1)
      Dt_B(iBlock) = DomainSize_D(iDim)/(nIjk_D(iDim)*nRootTest_D(iDim))
    end do

    if (DoTest) write (*, *) 'test prolong and balance'
    call refine_tree_node(1)
    if (DoTest) call show_tree('after refine_tree_node')
    call distribute_tree(.false.)
    if (DoTest) then
      call show_tree('after distribute_tree(.false.)')
      write (*, *) 'iProc, iProcNew_A=', iProc, iProcNew_A(1:nNode)
    end if

    call do_amr(nVar, State_VGB, Dt_B)
    if (DoTest) call show_tree('after do_amr')
    call move_tree
    if (DoTest) call show_tree('after move_tree')
    if (DoTest) write (*, *) 'iAmrChange_B=', iAmrChange_B(1:nBlock)

    call check_state

    if (DoTest) write (*, *) 'test restrict and balance'
    call coarsen_tree_node(1)
    if (DoTest) call show_tree('after coarsen_tree_node')
    call distribute_tree(.false.)
    if (DoTest) then
      call show_tree('after distribute_tree(.false.)')
      write (*, *) 'iProc, iProcNew_A=', iProc, iProcNew_A(1:nNode)
    end if
    call do_amr(nVar, State_VGB, Dt_B)
    if (DoTest) call show_tree('after do_amr')
    call move_tree
    if (DoTest) call show_tree('after move_tree')
    if (DoTest) write (*, *) 'iAmrChange_B=', iAmrChange_B(1:nBlock)

    call check_state

    ! tests with mask
    if (DoTest) write (*, *) 'test masked cells and extra data'

    allocate (Used_GB(MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlockTest), &
              TestState_VC(nVar, nI, nJ, nK), ExtraData_IB(nExtraData, MaxBlockTest))

    Used_GB = .true.

    do iBlock = 1, nBlock
      if (Unused_B(iBlock)) CYCLE
      State_VGB(:, :, :, :, iBlock) = Xyz_DGB(1:nDim, :, :, :, iBlock)
      ! set the time step to the cells size in the first AMR direction
      iDim = iDimAmr_D(1)
      Dt_B(iBlock) = DomainSize_D(iDim)/(nIjk_D(iDim)*nRootTest_D(iDim))
      ! Set extra data
      ExtraData_IB(:, iBlock) = &
        [real(iNode_B(iBlock)), Xyz_DGB(1, 1, 1, 1, iBlock)]
    end do

    call refine_tree_node(1)
    call distribute_tree(.false.)

    ! Mask out the middle cell of node 1
    if (iTree_IA(Proc_, 1) == iProc) then
      iBlock = iTree_IA(Block_, 1)
      State_VGB(:, (MinI + MaxI)/2, (MinJ + MaxJ)/2, (MinK + MaxK)/2, iBlock) = -7777
      Used_GB((MinI + MaxI)/2, (MinJ + MaxJ)/2, (MinK + MaxK)/2, iBlock) = .false.
    end if

    call do_amr(nVar, State_VGB, Dt_B, Used_GB=Used_GB, &
                nExtraData=nExtraData, &
                pack_extra_data=test_pack, &
                unpack_extra_data=test_unpack)

    call move_tree

    call check_extra_data

    ! Nodes that needs special care in the testing
    allocate (iEffectedNode_A(nNode))
    iEffectedNode_A = unset_

    ! 0: the masked cell, 1: x+ neighbor, 2: y+ neighbor
    ! 3: z+ neighbor
    ! The first child of Node of the refined node
    iChild = Child1_
    iEffectedNode_A(iTree_IA(iChild, 1)) = 0
    iChild = iChild + 1
    ! As masked cell is on the corner of the child nodes
    ! it will also affect the neighboring cell in i-direction
    ! The running "iChild" counter is using the Morton indexing
    ! to give the right neighbor independent of which directions
    ! are refined.
    if (iRatio == 2) then
      iEffectedNode_A(iTree_IA(iChild, 1)) = 1
      iChild = iChild + 1
    end if
    ! Neighbor in j direction
    if (jRatio == 2) then
      iEffectedNode_A(iTree_IA(iChild, 1)) = 2
      iChild = iChild + iRatio
    end if
    ! and k direction
    if (kRatio == 2) iEffectedNode_A(iTree_IA(iChild, 1)) = 3

    call check_state_mask

    !--------------- END refine --------------------
    iEffectedNode_A = unset_
    Used_GB = .true.

    if (iTree_IA(Proc_, iTree_IA(Child1_, 1)) == iProc) then
      iBlock = iTree_IA(Block_, iTree_IA(Child1_, 1))
      if (iRatio == 2) then
        i = nI
      else
        i = (MinI + MaxI)/2
      end if
      if (jRatio == 2) then
        j = nJ
      else
        j = (MinJ + MaxJ)/2
      end if
      if (kRatio == 2) then
        k = nK
      else
        k = (MinK + MaxK)/2
      end if

      State_VGB(:, i, j, k, iBlock) = -7777
      Used_GB(i, j, k, iBlock) = .false.
    end if

    call coarsen_tree_node(1)
    call distribute_tree(.false.)

    call do_amr(nVar, State_VGB, Dt_B, Used_GB=Used_GB, &
                nExtraData=nExtraData, &
                pack_extra_data=test_pack, &
                unpack_extra_data=test_unpack)

    call move_tree

    call check_extra_data

    call check_state_mask

    deallocate (State_VGB, Dt_B, Used_GB, iEffectedNode_A, TestState_VC, &
                ExtraData_IB)

    call clean_grid
    call clean_tree

  contains
    !==========================================================================
    subroutine check_extra_data

      !------------------------------------------------------------------------
      do iBlock = 1, nBlock
        if (iAmrChange_B(iBlock) /= AmrMoved_) CYCLE
        if (abs(ExtraData_IB(1, iBlock) - iNode_B(iBlock)) > 1e-6 .or. &
            abs(ExtraData_IB(2, iBlock) - Xyz_DGB(1, 1, 1, 1, iBlock)) > 1e-6) &
          write (*, *) NameSub, ' error for iProc,iBlock,ExtraData,iNode,x=', &
          iProc, iBlock, ExtraData_IB(:, iBlock), &
          iNode_B(iBlock), Xyz_DGB(1, 1, 1, 1, iBlock)
      end do

    end subroutine check_extra_data
    !==========================================================================
    subroutine check_state

      !------------------------------------------------------------------------
      do iBlock = 1, nBlock
        if (Unused_B(iBlock)) CYCLE

        if (any(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) &
                    - Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)) > 1e-6)) then
          write (*, *) NameSub, ' error for iProc,iBlock,maxloc=', iProc, iBlock, &
            maxloc(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) &
                       - Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)))
        end if

        iDim = iDimAmr_D(1)
        if (abs(Dt_B(iBlock) - CellSize_DB(iDim, iBlock)) > 1e-6) &
          write (*, *) NameSub, ' error for iProc,iBlock,dt,dx=', &
          iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim, iBlock)

      end do

    end subroutine check_state
    !==========================================================================
    subroutine check_state_mask

      integer :: iDn, iUp, jDn, jUp, kDn, kUp
      integer :: Di, Dj, Dk

      !------------------------------------------------------------------------
      do iBlock = 1, nBlock
        if (Unused_B(iBlock)) CYCLE
        iNode = iNode_B(iBlock)
        select case (iEffectedNode_A(iNode))
        case (unset_)
          if (any(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) &
                      - Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)) > 1e-6)) then
            write (*, *) NameSub, ' error for iProc,iBlock,maxloc=', &
              iProc, iBlock, &
              maxloc(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) &
                         - Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)))
          end if

          iDim = iDimAmr_D(1)
          if (abs(Dt_B(iBlock) - CellSize_DB(iDim, iBlock)) > 1e-6) &
            write (*, *) NameSub, ' error for iProc,iBlock,dt,dx=', &
            iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim, iBlock)
        case (0) ! the masked block

          if (iRatio == 2) then
            iDn = nI - iRatio + 1
            iUp = nI
            Di = iRatio
          else
            iDn = (MinI + MaxI)/2
            iUp = iDn
            Di = 0
          end if

          if (jRatio == 2) then
            jDn = nJ - jRatio + 1
            jUp = nJ
            Dj = jRatio
          else
            jDn = (MinJ + MaxJ)/2
            jUp = jDn
            Dj = 0
          end if

          if (kRatio == 2) then
            kDn = nK - kRatio + 1
            kUp = nK
            Dk = kRatio
          else
            kDn = (MinK + MaxK)/2
            kUp = kDn
            Dk = 0
          end if

          TestState_VC = Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)

          ! the average of the masked Cell
          do iVar = 1, nVar
            TestState_VC(iVar, iDn:iUp, jDn:jUp, kDn:kUp) = &
              sum(Xyz_DGB(iVar, iDn:iUp, jDn:jUp, kDn:kUp, iBlock))/ &
              (iRatio*jRatio*kRatio)
          end do

          ! The neighbor cell in i direction will only have
          ! 1st order prolongation
          TestState_VC(1, iDn - Di:iUp - Di, jDn:jUp, kDn:kUp) = &
            sum(Xyz_DGB(1, iDn - Di:iUp - Di, jDn:jUp, kDn:kUp, iBlock))/ &
            (iRatio*jRatio*kRatio)

          ! The neighbor cell in j direction will only have
          ! 1st order prolongation
          if (nDim > 1) then
            TestState_VC(2, iDn:iUp, jDn - Dj:jUp - Dj, kDn:kUp) = &
              sum(Xyz_DGB(2, iDn:iUp, jDn - Dj:jUp - Dj, kDn:kUp, iBlock))/ &
              (iRatio*jRatio*kRatio)
          end if

          ! The neighbor cell in k direction will only have
          ! 1st order prolongation
          if (nDim > 2) then
            TestState_VC(3, iDn:iUp, jDn:jUp, kDn - Dk:kUp - Dk) = &
              sum(Xyz_DGB(3, iDn:iUp, jDn:jUp, kDn - Dk:kUp - Dk, iBlock))/ &
              (iRatio*jRatio*kRatio)
          end if

          if (any(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) &
                      - TestState_VC) > 1e-6)) then
            write (*, *) NameSub, ' case 0 error for iProc,iBlock,maxloc=', &
              iProc, iBlock, &
              maxloc(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) - TestState_VC))
          end if

          iDim = iDimAmr_D(1)
          if (abs(Dt_B(iBlock) - CellSize_DB(iDim, iBlock)) > 1e-6) &
            write (*, *) NameSub, ' error for iProc,iBlock,dt,dx=', &
            iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim, iBlock)

        case (1) ! i neigbor of mask block

          if (iRatio == 2) then
            iDn = 1
            iUp = iRatio
            Di = iRatio
          else
            iDn = (MinI + MaxI)/2
            iUp = iDn
            Di = 0
          end if

          if (jRatio == 2) then
            jDn = nJ - jRatio + 1
            jUp = nJ
            Dj = jRatio
          else
            jDn = (MinJ + MaxJ)/2
            jUp = jDn
            Dj = 0
          end if

          if (kRatio == 2) then
            kDn = nK - kRatio + 1
            kUp = nK
            Dk = kRatio
          else
            kDn = (MinK + MaxK)/2
            kUp = kDn
            Dk = 0
          end if

          TestState_VC = Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)

          ! The neighbor cell in i direction will only have
          ! 1st order prolongation
          TestState_VC(1, iDn:iUp, jDn:jUp, kDn:kUp) = &
            sum(Xyz_DGB(1, iDn:iUp, jDn:jUp, kDn:kUp, iBlock))/ &
            (iRatio*jRatio*kRatio)

          if (any(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) &
                      - TestState_VC) > 1e-6)) then
            write (*, *) NameSub, ' case 1 error for iProc,iBlock,maxloc=', &
              iProc, iBlock, &
              maxloc(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) - TestState_VC))
          end if

          iDim = iDimAmr_D(1)
          if (abs(Dt_B(iBlock) - CellSize_DB(iDim, iBlock)) > 1e-6) &
            write (*, *) NameSub, ' error for iProc,iBlock,dt,dx=', &
            iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim, iBlock)

        case (2) ! j neighbore of masked block

          if (jRatio == 2) then
            jDn = 1
            jUp = jRatio
            Dj = jRatio
          else
            jDn = (MinJ + MaxJ)/2
            jUp = jDn
            Dj = 0
          end if

          if (iRatio == 2) then
            iDn = nI - iRatio + 1
            iUp = nI
            Di = iRatio
          else
            iDn = (MinI + MaxI)/2
            iUp = iDn
            Di = 0
          end if

          if (kRatio == 2) then
            kDn = nK - kRatio + 1
            kUp = nK
            Dk = kRatio
          else
            kDn = (MinK + MaxK)/2
            kUp = kDn
            Dk = 0
          end if

          TestState_VC = Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)

          ! The neighbor cell in j direction will only have
          ! 1st order prolongation
          TestState_VC(2, iDn:iUp, jDn:jUp, kDn:kUp) = &
            sum(Xyz_DGB(2, iDn:iUp, jDn:jUp, kDn:kUp, iBlock))/ &
            (iRatio*jRatio*kRatio)

          if (any(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) &
                      - TestState_VC) > 1e-6)) then
            write (*, *) NameSub, ' case 2 error for iProc,iBlock,maxloc=', &
              iProc, iBlock, &
              maxloc(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) - TestState_VC))
          end if

          iDim = iDimAmr_D(1)
          if (abs(Dt_B(iBlock) - CellSize_DB(iDim, iBlock)) > 1e-6) &
            write (*, *) NameSub, ' error for iProc,iBlock,dt,dx=', &
            iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim, iBlock)

        case (3) ! k neighbore of masked block

          if (kRatio == 2) then
            kDn = 1
            kUp = kRatio
            Dk = kRatio
          else
            kDn = (MinK + MaxK)/2
            kUp = kDn
            Dk = 0
          end if

          if (iRatio == 2) then
            iDn = nI - iRatio + 1
            iUp = nI
            Di = iRatio
          else
            iDn = (MinI + MaxI)/2
            iUp = iDn
            Di = 0
          end if

          if (jRatio == 2) then
            jDn = nJ - jRatio + 1
            jUp = nJ
            Dj = jRatio
          else
            jDn = (MinJ + MaxJ)/2
            jUp = jDn
            Dj = 0
          end if

          TestState_VC = Xyz_DGB(1:nDim, 1:nI, 1:nJ, 1:nK, iBlock)

          ! The neighbor cell in k direction will only have
          ! 1st order prolongation
          TestState_VC(3, iDn:iUp, jDn:jUp, kDn:kUp) = &
            sum(Xyz_DGB(3, iDn:iUp, jDn:jUp, kDn:kUp, iBlock))/ &
            (iRatio*jRatio*kRatio)

          if (any(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) &
                      - TestState_VC) > 1e-6)) then
            write (*, *) NameSub, ' case 3 error for iProc,iBlock,maxloc=', &
              iProc, iBlock, &
              maxloc(abs(State_VGB(:, 1:nI, 1:nJ, 1:nK, iBlock) - TestState_VC))
          end if

          iDim = iDimAmr_D(1)
          if (abs(Dt_B(iBlock) - CellSize_DB(iDim, iBlock)) > 1e-6) &
            write (*, *) NameSub, ' error for iProc,iBlock,dt,dx=', &
            iProc, iBlock, Dt_B(iBlock), CellSize_DB(iDim, iBlock)

        end select
      end do

    end subroutine check_state_mask
!==========================================================================
    subroutine show_state

      integer :: iProcShow
      !------------------------------------------------------------------------
      call barrier_mpi

      do iProcShow = 0, nProc - 1
        if (iProc == iProcShow) then
          do iBlock = 1, nBlock
            if (Unused_B(iBlock)) CYCLE
            write (*, '(a,2i4,100f8.4)') &
              'iProc, iBlock, State(1,1:nI,1,1,iBlock)=', &
              iProc, iBlock, State_VGB(1, 1:nI, 1, 1, iBlock)

            if (nDim > 1) write (*, '(a,2i4,100f8.3)') &
              'iProc, iBlock, State(2,1,1:nJ,1,iBlock)=', &
              iProc, iBlock, State_VGB(2, 1, 1:nJ, 1, iBlock)

            if (nDim > 2) write (*, '(a,2i4,100f8.2)') &
              'iProc, iBlock, State(3,1,1,1:nK,iBlock)=', &
              iProc, iBlock, State_VGB(3, 1, 1, 1:nK, iBlock)
          end do
        end if
        call barrier_mpi
      end do

    end subroutine show_state
!==========================================================================
    subroutine show_dt

      integer:: iProcShow
      !------------------------------------------------------------------------
      call barrier_mpi

      do iProcShow = 0, nProc - 1
        if (iProc == iProcShow) then
          do iBlock = 1, nBlock
            if (Unused_B(iBlock)) CYCLE
            write (*, *) 'iProc, iBlock, Dx, Dt =', &
              iProc, iBlock, CellSize_DB(1, iBlock), Dt_B(iBlock)
          end do
        end if
        call barrier_mpi
      end do

    end subroutine show_dt
!==========================================================================

  end subroutine test_amr
!============================================================================

  subroutine test_pack(iBlock, nBuffer, Buffer_I)

    integer, intent(in) :: iBlock
    integer, intent(in) :: nBuffer
    real, intent(out):: Buffer_I(nBuffer)
    !--------------------------------------------------------------------------

    if (nBuffer /= nExtraData) write (*, *) 'ERROR in test_pack: ', &
      'nBuffer, nExtraData=', nBuffer, nExtraData

    Buffer_I = ExtraData_IB(:, iBlock)

  end subroutine test_pack
!============================================================================
  subroutine test_unpack(iBlock, nBuffer, Buffer_I)

    integer, intent(in):: iBlock
    integer, intent(in):: nBuffer
    real, intent(in):: Buffer_I(nBuffer)
    !--------------------------------------------------------------------------

    if (nBuffer /= nExtraData) write (*, *) 'ERROR in test_unpack: ', &
      'nBuffer, nExtraData=', nBuffer, nExtraData

    ExtraData_IB(:, iBlock) = Buffer_I

  end subroutine test_unpack
!============================================================================

end module amr