program unit_test

  use ModMpi
  use BATL_lib
  use ModUtilities, ONLY: CON_stop
  use ModRandomNumber, ONLY: random_real

  use pass_cell, ONLY: test_pass_cell
  use pass_face, ONLY: test_pass_face
  use amr, ONLY: test_amr

  implicit none

  integer:: iError

  !----------------------------------------------------------------------------
  call MPI_init(iError)
  call init_mpi(MPI_COMM_WORLD)

  call test_tree
  call test_geometry
  call test_grid
  call test_pass_cell
  call test_pass_face
  call test_pass_node
  call test_amr
  call test_amr_criteria

  call MPI_finalize(iError)

contains
  subroutine test_tree ! unit test

    integer, parameter:: MaxBlockTest = 50
    integer, parameter:: nRootTest_D(MaxDim) = [3, 2, 1]
    logical, parameter:: IsPeriodicTest_D(MaxDim) = [.true., .true., .false.]
    real, parameter:: CoordTest_D(MaxDim) = 0.99

    integer :: iNode, iCoord_D(MaxDim), iCell_D(MaxDim)
    real    :: Distance_D(MaxDim), DistanceGood_D(MaxDim), CellSize_D(MaxDim)

    integer :: iNodeCell_II(0:nDim, 2**nDim), iNodeCellGood_II(0:3, 8)
    real    :: Weight_I(2**nDim), WeightGood_I(8)

    integer, allocatable:: iTypeNode_I(:)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'test_tree'
    !--------------------------------------------------------------------------
    DoTest = iProc == 0

    if (DoTest) write (*, *) 'Starting ', NameSub
    if (DoTest) write (*, *) 'Testing init_tree'
    call init_tree(MaxBlockTest)
    if (MaxBlock /= MaxBlockTest) &
         write (*, *) 'init_tree failed, MaxBlock=', &
         MaxBlock, ' should be ', MaxBlockTest
    if (MaxNode /= 2*ceiling(MaxBlockTest*nProc*(1 + 1.0/(2**nDimAmr - 1)))) &
         write (*, *) 'init_tree failed, MaxNode=', MaxNode, &
         ' should be', 2*ceiling(50*nProc*(1 + 1.0/(2**nDimAmr - 1)))

    if (DoTest) write (*, *) 'Testing init_geometry'
    call init_geometry('cartesian', IsPeriodicTest_D(1:nDim))
    if (any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
         write (*, *) 'init_geometry failed, IsPeriodic_D=', &
         IsPeriodic_D(1:nDim), ' should be ', IsPeriodicTest_D(1:nDim)

    if (DoTest) write (*, *) 'Testing i_node_new()'
    iNode = i_node_new()
    if (iNode /= 1) &
         write (*, *) 'i_node_new() failed, iNode=', iNode, ' should be 1'

    if (DoTest) write (*, *) 'Testing set_tree_root'
    call set_tree_root(nRootTest_D(1:nDim))

    if (DoTest) call show_tree('after set_tree_root')

    if (any(nRoot_D(1:nDim) /= nRootTest_D(1:nDim))) &
         write (*, *) 'set_tree_root failed, nRoot_D=', nRoot_D(1:nDim), &
         ' should be ', nRootTest_D(1:nDim)

    iCoord_D = [3, 1, 1]

    if (any(iTree_IA(Coord1_:Coord0_ + nDim, 3) /= iCoord_D(1:nDim))) &
         write (*, *) 'set_tree_root failed, coordinates of node four=', &
         iTree_IA(Coord1_:Coord0_ + nDim, 3), ' should be ', iCoord_D(1:nDim)

    if (DoTest) write (*, *) 'Testing find_tree_cell'
    call find_tree_cell(CoordTest_D, iNode, iCell_D, Distance_D)
    if (iNode /= nRoot) write (*, *) 'ERROR: Test find point failed, iNode=', &
         iNode, ' instead of', nRoot

    if (.not. is_point_inside_node(CoordTest_D, iNode)) &
         write (*, *) 'ERROR: Test find point failed'

    if (any(iCell_D(1:nDim) /= nIJK_D(1:nDim))) &
         write (*, *) 'ERROR: Test find point failed, iCell_D=', &
         iCell_D(1:nDim), ' instead of', nIjk_D(1:nDim)

    ! Cell size in units where the whole domain is 1.0
    CellSize_D = 1.0/(nRoot_D*nIJK_D)
    ! Distance to the last grid cell, normalized to the cell size
    DistanceGood_D = (CoordTest_D - (1.0 - CellSize_D/2))/CellSize_D
    if (any(abs(Distance_D(1:nDim) - DistanceGood_D(1:nDim)) > 1e-6)) &
         write (*, *) 'ERROR: Test find point failed, Distance_D=', &
         Distance_D(1:nDim), ' instead of ', DistanceGood_D(1:nDim)

    if (DoTest) write (*, *) 'Testing interpolate_tree'
    call interpolate_tree(CoordTest_D, iNodeCell_II, Weight_I)
    select case (nDim)
    case (1)
       iNodeCellGood_II(0:1, 1:2) = reshape([nRoot, nI, nRoot, nI + 1], &
            [2, 2])
       WeightGood_I(1:2) = [1 - Distance_D(1), Distance_D(1)]
    case (2)
       iNodeCellGood_II(0:2, 1:4) = reshape( &
            [nRoot, nI, nJ, nRoot, nI + 1, nJ, nRoot, nI, nJ + 1, nRoot, nI + 1, nJ + 1], &
            [3, 4])
       WeightGood_I(1:4) = [ &
            (1 - Distance_D(1))*(1 - Distance_D(2)), &
            Distance_D(1)*(1 - Distance_D(2)), &
            (1 - Distance_D(1))*Distance_D(2), &
            Distance_D(1)*Distance_D(2) &
            ]
    case (3)
       iNodeCellGood_II(0:3, 1:8) = reshape([ &
            nRoot, nI, nJ, nK, &
            nRoot, nI + 1, nJ, nK, nRoot, nI, nJ + 1, nK, nRoot, nI + 1, nJ + 1, nK, &
            nRoot, nI, nJ, nK + 1, &
            nRoot, nI + 1, nJ, nK + 1, nRoot, nI, nJ + 1, nK + 1, nRoot, nI + 1, nJ + 1, nK + 1], &
            [4, 8])
       WeightGood_I(1:8) = [ &
            (1 - Distance_D(1))*(1 - Distance_D(2))*(1 - Distance_D(3)), &
            Distance_D(1)*(1 - Distance_D(2))*(1 - Distance_D(3)), &
            (1 - Distance_D(1))*Distance_D(2)*(1 - Distance_D(3)), &
            Distance_D(1)*Distance_D(2)*(1 - Distance_D(3)), &
            (1 - Distance_D(1))*(1 - Distance_D(2))*Distance_D(3), &
            Distance_D(1)*(1 - Distance_D(2))*Distance_D(3), &
            (1 - Distance_D(1))*Distance_D(2)*Distance_D(3), &
            Distance_D(1)*Distance_D(2)*Distance_D(3) &
            ]
    end select

    if (any(iNodeCell_II /= iNodeCellGood_II(0:nDim, 1:2**nDim))) &
         write (*, *) 'ERROR: Test interpolate_tree failed, iNodeCell_II=', &
         iNodeCell_II, ' instead of ', iNodeCellGood_II

    if (any(abs(Weight_I - WeightGood_I(1:2**nDim)) > 1e-6)) &
         write (*, *) 'ERROR: Test interpolate_tree failed, Weight_I=', &
         Weight_I, ' instead of ', WeightGood_I

    if (DoTest) write (*, *) 'Testing distribute_tree 1st'
    call distribute_tree(.true.)
    if (DoTest) call show_tree('after distribute_tree 1st', .true.)

    if (DoTest) write (*, *) 'Testing refine_tree_node'
    ! Refine the node where the point was found and find it again
    call refine_tree_node(iNode)

    if (DoTest) write (*, *) 'Testing distribute_tree 2nd'

    ! Set node type to level+1
    allocate (iTypeNode_I(MaxNode))
    iTypeNode_I = 1 + iTree_IA(Level_, :)
    call distribute_tree(.true., iTypeNode_I)
    if (DoTest) call show_tree('after distribute_tree with type=level', .true.)

    ! Set node type to the second coordinate index
    iTypeNode_I = iTree_IA(Coord2_, :)
    call distribute_tree(.true., iTypeNode_I)
    if (DoTest) call show_tree('after distribute_tree with type=Coord2', .true.)

    ! Use default (single type)
    call distribute_tree(.true.)
    if (DoTest) call show_tree('after distribute_tree 2nd', .true.)

    call find_tree_node(CoordTest_D, iNode)
    if (.not. is_point_inside_node(CoordTest_D, iNode)) &
         write (*, *) 'ERROR: Test find point failed for iNode=', iNode

    ! Refine another node
    if (DoTest) write (*, *) 'nRoot=', nRoot
    call refine_tree_node(2)

    if (DoTest) call show_tree('after another refine_tree_node')

    if (DoTest) write (*, *) 'Testing coarsen_tree_node'

    ! Coarsen back the last root node and find point again
    call coarsen_tree_node(nRoot)
    if (DoTest) call show_tree('after coarsen_tree_node')

    ! Distribute the new tree
    if (DoTest) write (*, *) 'Testing distribute_tree 3rd'
    call distribute_tree(.true.)
    if (DoTest) call show_tree('after distribute_tree 3rd', .true.)

    call find_tree_node(CoordTest_D, iNode)
    if (iNode /= nRoot) write (*, *) &
         'ERROR: coarsen_tree_node+compact failed, iNode=', &
         iNode, ' instead of', nRoot
    if (.not. is_point_inside_node(CoordTest_D, iNode)) &
         write (*, *) 'ERROR: is_point_inside_node failed'

    if (iTree_IA(Status_, nNode + 1) /= Unset_) &
         write (*, *) 'ERROR: compact_tree failed, nNode=', nNode, &
         ' but status of next node is', iTree_IA(Status_, nNode + 1), &
         ' instead of ', Unset_
    if (any(iTree_IA(Status_, 1:nNode) == Unset_)) &
         write (*, *) 'ERROR: compact_tree failed, nNode=', nNode, &
         ' but iTree_IA(Status_, 1:nNode)=', &
         iTree_IA(Status_, 1:nNode), ' contains unset=', Unset_
    call find_tree_node(CoordTest_D, iNode)
    if (iNode /= nRoot) write (*, *) 'ERROR: compact_tree faild, iNode=', &
         iNode, ' instead of', nRoot
    if (.not. is_point_inside_node(CoordTest_D, iNode)) &
         write (*, *) 'ERROR: is_point_inside_node failed'

    if (DoTest) write (*, *) 'Testing write_tree_file'
    call write_tree_file('tree.rst')

    if (DoTest) write (*, *) 'Testing read_tree_file'
    iTree_IA = Unset_
    nRoot_D = 0
    call read_tree_file('tree.rst')
    if (DoTest) call show_tree('after read_tree')

    call find_tree_node(CoordTest_D, iNode)
    if (iNode /= nRoot) write (*, *) 'ERROR: compact_tree failed, iNode=', &
         iNode, ' instead of', nRoot

    if (DoTest) write (*, *) 'Testing distribute_tree 4th'
    call distribute_tree(.true.)
    if (DoTest) call show_tree('after distribute_tree 4th', .true.)

    if (DoTest) write (*, *) 'Testing clean_tree'
    call clean_tree
    if (DoTest) write (*, *) 'MaxNode=', MaxNode

  end subroutine test_tree
  !============================================================================

  subroutine test_geometry ! unit test

    logical:: IsPeriodicTest_D(MaxDim)
    real:: Xyz_D(MaxDim), Coord_D(MaxDim), Good_D(MaxDim)

    real:: r, GenR
    real:: Rgen_I(5) = [1.0, 1.2, 5.0, 25.0, 100.0]

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'test_geometry'
    !--------------------------------------------------------------------------
    DoTest = iProc == 0

    if (DoTest) write (*, *) 'Starting ', NameSub
    if (DoTest) write (*, *) 'Testing init_geometry for Cartesian'
    call init_geometry

    if (TypeGeometry /= 'cartesian') &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be Cartesian by default'

    if (.not. IsCartesian .or. IsRotatedCartesian .or. &
         IsRzGeometry .or. IsSpherical .or. IsCylindrical) &
         write (*, *) 'ERROR: init_geometry failed for Cartesian grid, ', &
         'IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical=', &
         IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

    if (IsLogRadius .or. IsGenRadius) &
         write (*, *) 'ERROR: init_geometry failed for Cartesian grid, ', &
         'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

    if (any(IsPeriodic_D)) &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'IsPeriodic_D =', IsPeriodic_D, ' should be all false'

    if (DoTest) write (*, *) 'Testing xyz_to_coord for Cartesian'
    Xyz_D = [1., 2., 3.]
    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(Coord_D /= Xyz_D)) &
         write (*, *) 'ERROR: xyz_to_coord failed for Cartesian grid, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for Cartesian'
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(Coord_D /= Xyz_D)) &
         write (*, *) 'ERROR: coord_to_xyz failed for Cartesian grid, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D

    if (nDim == 2) then
       if (DoTest) write (*, *) 'Testing init_geometry for RZ geometry'
       IsPeriodicTest_D = [.true., .false., .false.]

       call init_geometry('rz', IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))

       if (TypeGeometry /= 'rz') &
            write (*, *) 'ERROR: init_geometry failed, ', &
            'TypeGeometry=', TypeGeometry, ' should be rz'

       if (.not. IsRzGeometry .or. IsCartesian .or. IsSpherical .or. IsCylindrical) &
            write (*, *) 'ERROR: init_geometry failed for RZ grid, ', &
            'IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical=', &
            IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

       if (IsLogRadius .or. IsGenRadius) &
            write (*, *) 'ERROR: init_geometry failed for RZ grid, ', &
            'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

       if (any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
            write (*, *) 'ERROR: init_geometry failed, ', &
            'for TypeGeometry=', TypeGeometry, &
            'IsPeriodic_D =', IsPeriodic_D(1:nDim), &
            ' should be ', IsPeriodicTest_D(1:nDim)

       if (DoTest) write (*, *) 'Testing xyz_to_coord for RZ geometry'
       Xyz_D = [1., 2., 3.]
       call xyz_to_coord(Xyz_D, Coord_D)
       if (any(Coord_D /= Xyz_D)) &
            write (*, *) 'ERROR: xyz_to_coord failed for RZ grid, ', &
            'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D

       if (DoTest) write (*, *) 'Testing coord_to_xyz for RZ geometry'
       call coord_to_xyz(Coord_D, Xyz_D)
       if (any(Coord_D /= Xyz_D)) &
            write (*, *) 'ERROR: coord_to_xyz failed for RZ grid, ', &
            'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D

    end if

    if (nDim == 1) RETURN

    if (DoTest) write (*, *) 'Testing init_geometry for cylindrical_lnr'
    IsPeriodicTest_D = [.false., .true., .true.]

    call init_geometry('cylindrical_lnr', &
         IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))

    if (TypeGeometry /= 'cylindrical_lnr') &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be cylindrical_lnr'

    if (.not. IsCylindrical .or. IsCartesian .or. IsRzGeometry .or. IsSpherical) &
         write (*, *) 'ERROR: init_geometry failed for cylindrical_lnr, ', &
         'IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical=', &
         IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

    if (.not. IsLogRadius .or. IsGenRadius) &
         write (*, *) 'ERROR: init_geometry failed for cylindrical_lnr, ', &
         'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

    if (any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry, &
         'IsPeriodic_D =', IsPeriodic_D(1:nDim), &
         ' should be ', IsPeriodicTest_D(1:nDim)

    if (DoTest) write (*, *) 'Testing xyz_to_coord for cylindrical_lnr'
    Xyz_D = [1., 2., 3.]
    Good_D = [log(sqrt(5.)), atan2(2., 1.), 3.]
    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for cylindrical_lnr, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing init_geometry for roundcube'
    IsPeriodicTest_D = [.false., .false., .false.]

    call init_geometry('roundcube', &
         IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))

    if (TypeGeometry /= 'roundcube') &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be roundcube'

    if (.not. IsRoundCube .or. IsCartesian .or. IsRzGeometry &
         .or. IsCylindrical .or. IsSpherical) &
         write (*, *) 'ERROR: init_geometry failed for roundcube, ', &
         'IsRoundCube,IsCartesian,IsRzGeometry,IsSpherical,IsCylindrical=', &
         IsRoundCube, IsCartesian, IsRzGeometry, IsSpherical, IsCylindrical

    if (IsLogRadius .or. IsGenRadius) &
         write (*, *) 'ERROR: init_geometry failed for roundcube, ', &
         'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

    if (any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry, &
         'IsPeriodic_D =', IsPeriodic_D(1:nDim), &
         ' should be ', IsPeriodicTest_D(1:nDim)

    if (DoTest) write (*, *) 'Testing roundcube with rRound0=200 rRound1=320'
    rRound0 = 200.0
    rRound1 = 320.0

    if (DoTest) write (*, *) 'Testing xyz_to_coord for roundcube along X axis'

    ! points along main axes with L1 = rRound1 are most distorted
    Good_D = [320., 0., 0.]
    Xyz_D = sqrt(real(nDim))*Good_D

    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for roundcube, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for roundcube along X axis'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(abs(Xyz_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: coord_to_xyz failed for roundcube, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing xyz_to_coord for roundcube along diagonal'
    if (nDim == 3) then
       Xyz_D = [300., 300., 300.]
    elseif (nDim == 2) then
       Xyz_D = [300., 300., 0.]
    end if
    Good_D = Xyz_D

    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for roundcube, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for roundcube along diagonal'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(abs(Xyz_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: coord_to_xyz failed for roundcube, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing xyz_to_coord for arbitrary point'
    if (nDim == 3) then
       Xyz_D = [397.1825374147, 264.7883582764, 132.394179138]
       Good_D = [300., 200., 100.]
    elseif (nDim == 2) then
       Xyz_D = [344.1742027, 229.4494684, 0.]
       Good_D = [300., 200., 0.]
    end if

    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for roundcube, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for arbitrary point'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(abs(Xyz_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: coord_to_xyz failed for roundcube, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing xyz_to_coord for roundcube inside rRound0'
    Xyz_D = [100., 90., 0.]    ! Inside rRound0, points are not distorted
    Good_D = Xyz_D
    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for roundcube, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for roundcube inside rRound0'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(abs(Xyz_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: coord_to_xyz failed for roundcube, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing roundcube with rRound0=1, rRound1=0.6'
    rRound0 = 1.0
    rRound1 = 0.6

    if (DoTest) write (*, *) 'Testing xyz_to_coord for roundcube'
    if (nDim == 2) then
       Xyz_D = [0.0964809, 0.1929618, 0.]
       Good_D = [0.1, 0.2, 0.]
    else if (nDim == 3) then
       Xyz_D = [0.09008918, 0.18017837, 0.2702675]
       Good_D = [0.1, 0.2, 0.3]
    end if

    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for roundcube, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for roundcube'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(abs(Xyz_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: coord_to_xyz failed for roundcube, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing xyz_to_coord for roundcube'
    if (nDim == 2) then
       Xyz_D = [0.5736097, 0.4916654, 0.]
       Good_D = [0.7, 0.6, 0.]
    else if (nDim == 3) then
       Xyz_D = [0.52539750154, 0.450340715612, 0.37528392967]
       Good_D = [0.7, 0.6, 0.5]
    endif

    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for roundcube, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for roundcube'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(abs(Xyz_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: coord_to_xyz failed for roundcube, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing xyz_to_coord for roundcube along X axis'
    if (nDim == 2) then
       Xyz_D = [0.7, 0., 0.]
    else if (nDim == 3) then
       Xyz_D = [0.3, 0., 0.]
    endif
    Good_D = Xyz_D

    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for roundcube, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for roundcube along X axis'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(abs(Xyz_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: coord_to_xyz failed for roundcube, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D, ' should be ', Good_D

    if (nDim < 3) RETURN

    if (DoTest) write (*, *) 'Testing init_geometry for spherical_genr'
    IsPeriodicTest_D = [.false., .true., .false.]

    call init_geometry('spherical_genr', IsPeriodicIn_D=IsPeriodicTest_D, &
         RgenIn_I=Rgen_I)

    if (TypeGeometry /= 'spherical_genr') &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'TypeGeometry=', TypeGeometry, ' should be spherical_genr'

    if (.not. IsSpherical .or. IsCartesian .or. IsRzGeometry .or. IsCylindrical) &
         write (*, *) 'ERROR: init_geometry failed for spherical_genr, ', &
         'IsCartesian, IsRzGeometry, IsCylindrical, IsSpherical=', &
         IsCartesian, IsRzGeometry, IsCylindrical, IsSpherical

    if (.not. IsGenRadius .or. IsLogRadius) &
         write (*, *) 'ERROR: init_geometry failed for spherical_genr, ', &
         'IsLogRadius, IsGenRadius =', IsLogRadius, IsGenRadius

    if (any(IsPeriodic_D(1:nDim) .neqv. IsPeriodicTest_D(1:nDim))) &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry, &
         'IsPeriodic_D =', IsPeriodic_D(1:nDim), &
         ' should be ', IsPeriodicTest_D(1:nDim)

    if (nRgen /= size(Rgen_I)) &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry, &
         'nRgen=', nRgen, ' should be ', size(Rgen_I)

    if (.not. allocated(LogRgen_I)) &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry, &
         'LogRgen_I is not allocated'

    if (any(abs(exp(LogRgen_I) - Rgen_I) > 1e-6)) &
         write (*, *) 'ERROR: init_geometry failed, ', &
         'for TypeGeometry=', TypeGeometry, &
         'exp(LogRgen_I) =', exp(LogRgen_I), ' should be ', Rgen_I

    if (DoTest) write (*, *) 'Testing radius_to_gen and gen_to_radius'
    r = sqrt(Rgen_I(2)*Rgen_I(3))
    GenR = r
    call radius_to_gen(GenR)
    if (abs(GenR - 1.5/4) > 1e-6) &
         write (*, *) 'ERROR: radius_to_gen failed for spherical_genr, ', &
         'r=', r, ' GenR =', GenR, ' should be ', 1.5/4

    ! Test conversion back
    call gen_to_radius(GenR)
    if (abs(GenR - r) > 1e-6) &
         write (*, *) 'ERROR: gen_to_radius failed for spherical_genr, ', &
         'Orig r=', r, ' new r =', GenR

    r = 1.0/Rgen_I(2)**2
    GenR = r
    call radius_to_gen(GenR)
    if (abs(GenR + 2.0/4) > 1e-6) &
         write (*, *) 'ERROR: radius_to_gen failed for spherical_genr, ', &
         'r=', r, ' GenR =', GenR, ' should be ', -2.0/4

    ! Test conversion back
    call gen_to_radius(GenR)
    if (abs(GenR - r) > 1e-6) &
         write (*, *) 'ERROR: gen_to_radius failed for spherical_genr, ', &
         'Orig r=', r, ' new r =', GenR

    r = 1600.0
    GenR = r
    call radius_to_gen(GenR)
    if (abs(GenR - (1 + 2./4)) > 1e-6) &
         write (*, *) 'ERROR: radius_to_gen failed for spherical_genr, ', &
         'r=', r, ' GenR =', GenR, ' should be ', 1 + 2./4.

    ! Test conversion back
    call gen_to_radius(GenR)
    if (abs(GenR - r) > 1e-6) &
         write (*, *) 'ERROR: gen_to_radius failed for spherical_genr, ', &
         'Orig r=', r, ' new r =', GenR

    if (DoTest) write (*, *) 'Testing xyz_to_coord for spherical_genr'
    Xyz_D = [9., 12., 20.]
    Good_D = [0.75, atan2(15., 20.), atan2(12., 9.)]
    call xyz_to_coord(Xyz_D, Coord_D)
    if (any(abs(Coord_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: xyz_to_coord failed for spherical_genr, ', &
         'Xyz_D =', Xyz_D, ' Coord_D =', Coord_D, ' should be ', Good_D

    if (DoTest) write (*, *) 'Testing coord_to_xyz for spherical_genr'
    Good_D = Xyz_D
    call coord_to_xyz(Coord_D, Xyz_D)
    if (any(abs(Xyz_D - Good_D) > 1e-6)) &
         write (*, *) 'ERROR: coord_to_xyz failed for spherical, ', &
         'Coord_D =', Coord_D, ' Xyz_D =', Xyz_D, ' should be ', Good_D

  end subroutine test_geometry
  !============================================================================

  subroutine test_grid

    use ModNumConst, ONLY: i_DD, cPi

    integer :: iBlock

    integer, parameter:: MaxBlockTest = 50
    integer, parameter:: nRootTest_D(MaxDim) = [3, 3, 3]
    logical, parameter:: IsPeriodicTest_D(MaxDim) = [.true., .true., .false.]
    real:: DomainMin_D(MaxDim) = [3.0, 2.0, 1.0]
    real:: DomainMax_D(MaxDim) = [9.0, 6.0, 4.0]

    ! number of points in each dimension to test interpolation
    integer, parameter:: &
         nPointI = 91, &
         nPointJ = 1 + 90*min(1, nDim - 1), &
         nPointK = 1 + 90*max(0, nDim - 2)

    integer, parameter:: nPoint_D(MaxDim) = [nPointI, nPointJ, nPointK]
    real, allocatable:: Point_VIII(:, :, :, :), PointAll_VIII(:, :, :, :)
    integer, parameter:: nVarPoint = nDim
    real:: XyzPoint_D(MaxDim), Point_V(nVarPoint), Weight
    integer:: iPoint, jPoint, kPoint, iPoint_D(MaxDim), iCell, nCell, iError
    integer:: iCell_II(0:nDim, 2**nDim)
    logical:: IsSecondOrder
    integer:: iDiscr_D(MaxDim)
    real   :: Weight_I(2**nDim)

    real:: Tolerance

    integer:: i, j, k, Di, Dj, Dk, iDim, iBlockOut, iProcOut, iCell_D(MaxDim)
    integer:: iNodeCenter
    real:: Radius, Phi, Xyz_D(MaxDim), Coord_D(MaxDim), Distance_D(MaxDim)
    real:: Good, Good_D(MaxDim)
    real, allocatable:: CellVolumeCart_B(:), CellFaceCart_DB(:, :)
    real:: Volume, VolumeAll

    ! testing check_interpolate_amr_gc
    integer, parameter:: nPointCheck = 10000
    integer, parameter:: nRootCheck = 2**nDim
    integer, parameter:: nBlockCheckMax = nRootCheck**2
    integer, parameter:: nCaseCheck = 2**nRootCheck
    integer:: iCase, iNodeCheck, iBlockCheck, iProcCheck, iSeed = 1

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'test_grid'
    !--------------------------------------------------------------------------
    DoTest = iProc == 0

    if (DoTest) then
       write (*, *) 'Starting ', NameSub
       write (*, *) 'Testing init_grid'
       write (*, *) 'nDimAmr, nIJK_D=', nDimAmr, nIJK_D
    end if
    ! Set Cartesian grid geometry before initializing tree and grid
    call init_geometry(IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))
    call init_tree(MaxBlockTest)
    call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
    call set_tree_root(nRootTest_D(1:nDim))

    call find_tree_node([0.5, 0.5, 0.5], iNodeCenter)
    call refine_tree_node(iNodeCenter)
    call distribute_tree(.true.)
    if (DoTest) call show_tree('After distribute_tree')

    if (DoTest) write (*, *) 'Testing create_grid'
    call create_grid

    if (iProc == 0) call show_grid_proc

    if (DoTest) write (*, *) 'Testing find_grid_block'
    Xyz_D = 0.0
    Xyz_D(1:nDim) = DomainMin_D(1:nDim)
    call find_grid_block(Xyz_D, iProcOut, iBlockOut, &
         iCell_D, Distance_D, UseGhostCell=.true.)
    if (iProc == iProcOut) then
       Xyz_D = Xyz_DGB(:, iCell_D(1), iCell_D(2), iCell_D(3), iBlockOut) &
            + 0.5*CellSize_DB(:, iBlockOut)
       if (any(abs(DomainMin_D(1:nDim) - Xyz_D(1:nDim)) > 1e-6)) then
          write (*, *) 'Error: DomainMin_D, Xyz_D=', &
               DomainMin_D, Xyz_D
          write (*, *) 'iProcOut, iBlockOut, iCell_D = ', &
               iProcOut, iBlockOut, iCell_D
       end if
    end if

    if (any(iCell_D(1:nDim) /= 0)) then
       write (*, *) 'Error: iCell_D=', iCell_D(1:nDim), ' should be 0'
       write (*, *) 'iProcOut, iBlockOut, Distance_D = ', &
            iProcOut, iBlockOut, Distance_D
    end if

    if (any(abs(Distance_D(1:nDim) - 0.5) > 1e-6)) then
       write (*, *) 'Error: Distance_D=', Distance_D(1:nDim), ' should be -0.5'
       write (*, *) 'iProcOut, iBlockOut, iCell_D = ', &
            iProcOut, iBlockOut, iCell_D
    end if

    Xyz_D = 0.0
    Xyz_D(1:nDim) = DomainMax_D(1:nDim)
    call find_grid_block(Xyz_D, iProcOut, iBlockOut, &
         iCell_D, Distance_D, UseGhostCell=.true.)
    if (iProc == iProcOut) then
       Xyz_D = Xyz_DGB(:, iCell_D(1), iCell_D(2), iCell_D(3), iBlockOut) &
            + 0.5*CellSize_DB(:, iBlockOut)
       if (any(abs(DomainMax_D(1:nDim) - Xyz_D(1:nDim)) > 1e-6)) then
          write (*, *) 'Error: DomainMax_D, Xyz_D=', &
               DomainMax_D, Xyz_D
          write (*, *) 'iProcOut, iBlockOut, iCell_D = ', &
               iProcOut, iBlockOut, iCell_D
       end if
    end if

    if (any(iCell_D(1:nDim) /= nIJK_D(1:nDim))) then
       write (*, *) 'Error: iCell_D=', iCell_D(1:nDim), &
            ' should be ', nIJK_D(1:nDim)
       write (*, *) 'iProcOut, iBlockOut, Distance_D = ', &
            iProcOut, iBlockOut, Distance_D
    end if

    if (any(abs(Distance_D(1:nDim) - 0.5) > 1e-6)) then
       write (*, *) 'Error: Distance_D=', Distance_D(1:nDim), ' should be +0.5'
       write (*, *) 'iProcOut, iBlockOut, iCell_D = ', &
            iProcOut, iBlockOut, iCell_D
    end if

    if (DoTest) write (*, *) 'Testing interpolate_grid'

    Xyz_D = 0.0
    if (.not. allocated(Point_VIII)) &
         allocate (Point_VIII(0:nVarPoint, nPointI, nPointJ, nPointK))
    Point_VIII = 0.0
    do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
       iPoint_D = [iPoint, jPoint, kPoint]
       XyzPoint_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim) - 0.5) &
            *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

       call interpolate_grid(XyzPoint_D, nCell, iCell_II, Weight_I)

       do iCell = 1, nCell
          Point_VIII(0, iPoint, jPoint, kPoint) = &
               Point_VIII(0, iPoint, jPoint, kPoint) + Weight_I(iCell)
          iBlock = iCell_II(0, iCell)
          iCell_D = 1
          iCell_D(1:nDim) = iCell_II(1:nDim, iCell)

          ! Interpolate the coordinates to check order of accuracy
          ! Note: Using array syntax in combination with the indirect
          ! iCell_D index fails with optimization for NAG v5.1
          do iDim = 1, nDim
             Xyz_D(iDim) = &
                  Xyz_DGB(iDim, iCell_D(1), iCell_D(2), iCell_D(3), iBlock)
          end do

          ! Take care of periodic dimensions: shift coordinates as necessary
          do iDim = 1, nDim
             if (.not. IsPeriodicTest_D(iDim)) CYCLE
             if (XyzPoint_D(iDim) < Xyz_D(iDim) - 2*CellSize_DB(iDim, iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) - DomainSize_D(iDim)

             if (XyzPoint_D(iDim) > Xyz_D(iDim) + 2*CellSize_DB(iDim, iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) + DomainSize_D(iDim)
          end do

          Point_VIII(1:nDim, iPoint, jPoint, kPoint) = &
               Point_VIII(1:nDim, iPoint, jPoint, kPoint) &
               + Weight_I(iCell)*Xyz_D(1:nDim)

       end do
    end do; end do; end do

    ! Collect contributions from all processors to proc 0
    if (nProc > 1) then
       allocate (PointAll_VIII(0:nVarPoint, nPointI, nPointJ, nPointK))
       call MPI_reduce(Point_VIII, PointAll_VIII, size(PointAll_VIII), &
            MPI_REAL, MPI_SUM, 0, iComm, iError)
       Point_VIII = PointAll_VIII
       deallocate (PointAll_VIII)
    end if

    if (iProc == 0) then
       ! Check interpolated coordinate values against point coordinates
       do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
          iPoint_D = [iPoint, jPoint, kPoint]
          Xyz_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim) - 0.5) &
               *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

          Weight = Point_VIII(0, iPoint, jPoint, kPoint)
          Point_V = Point_VIII(1:nDim, iPoint, jPoint, kPoint)/Weight

          if (abs(Weight - 1.0) < 1e-6) then
             Tolerance = 1e-6
          else
             Tolerance = 3e-2
          end if

          if (any(abs(Xyz_D(1:nDim) - Point_V) > Tolerance)) then
             write (*, *) 'ERROR: Point_V=', Point_V(1:nDim), &
                  ' should be ', Xyz_D(1:nDim)
             write (*, *) 'Total weight=', Weight
             write (*, *) 'i,j,kPoint=', iPoint_D(1:nDim)
             write (*, *) 'CoordMin,Max=', CoordMin_D(1:nDim), CoordMax_D(1:nDim)
          end if
       end do; end do; end do
    end if

    if (DoTest) write (*, *) 'Testing interpolate_grid_amr'
    Xyz_D = 0.0
    if (.not. allocated(Point_VIII)) &
         allocate (Point_VIII(0:nVarPoint, nPointI, nPointJ, nPointK))
    Point_VIII = 0.0
    do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
       iPoint_D = [iPoint, jPoint, kPoint]
       XyzPoint_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim) - 0.5) &
            *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

       call interpolate_grid_amr(XyzPoint_D, nCell, iCell_II, Weight_I, &
            IsSecondOrder)

       do iCell = 1, nCell
          Point_VIII(0, iPoint, jPoint, kPoint) = &
               Point_VIII(0, iPoint, jPoint, kPoint) + Weight_I(iCell)
          iBlock = iCell_II(0, iCell)
          iCell_D = 1
          iCell_D(1:nDim) = iCell_II(1:nDim, iCell)

          ! Interpolate the coordinates to check order of accuracy
          ! Note: Using array syntax in combination with the indirect
          ! iCell_D index fails with optimization for NAG v5.1
          do iDim = 1, nDim
             Xyz_D(iDim) = &
                  Xyz_DGB(iDim, iCell_D(1), iCell_D(2), iCell_D(3), iBlock)
          end do

          ! Take care of periodic dimensions: shift coordinates as necessary
          do iDim = 1, nDim
             if (.not. IsPeriodicTest_D(iDim)) CYCLE
             if (XyzPoint_D(iDim) < Xyz_D(iDim) - 2*CellSize_DB(iDim, iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) - DomainSize_D(iDim)

             if (XyzPoint_D(iDim) > Xyz_D(iDim) + 2*CellSize_DB(iDim, iBlock)) &
                  Xyz_D(iDim) = Xyz_D(iDim) + DomainSize_D(iDim)
          end do

          ! if point is close to boundary then interpolation isn't 2nd order
          where (.not. (IsSecondOrder .or. IsPeriodicTest_D(1:nDim))) &
               Xyz_D(1:nDim) = XyzPoint_D(1:nDim)

          Point_VIII(1:nDim, iPoint, jPoint, kPoint) = &
               Point_VIII(1:nDim, iPoint, jPoint, kPoint) &
               + Weight_I(iCell)*Xyz_D(1:nDim)

       end do

    end do; end do; end do

    ! Collect contributions from all processors to proc 0
    if (nProc > 1) then
       allocate (PointAll_VIII(0:nVarPoint, nPointI, nPointJ, nPointK))
       call MPI_reduce(Point_VIII, PointAll_VIII, size(PointAll_VIII), &
            MPI_REAL, MPI_SUM, 0, iComm, iError)
       Point_VIII = PointAll_VIII
       deallocate (PointAll_VIII)
    end if

    if (iProc == 0) then
       ! Check interpolated coordinate values against point coordinates
       do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
          iPoint_D = [iPoint, jPoint, kPoint]
          Xyz_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim) - 0.5) &
               *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

          Weight = Point_VIII(0, iPoint, jPoint, kPoint)
          Point_V = Point_VIII(1:nDim, iPoint, jPoint, kPoint)/Weight

          if (abs(Weight - 1.0) < 1e-6) then
             Tolerance = 1e-6
          else
             Tolerance = 3e-2
          end if

          if (any(abs(Xyz_D(1:nDim) - Point_V) > Tolerance)) then
             write (*, *) 'ERROR: Point_V=', Point_V(1:nDim), &
                  ' should be ', Xyz_D(1:nDim)
             write (*, *) 'Total weight=', Weight
             write (*, *) 'i,j,kPoint=', iPoint_D(1:nDim)
             write (*, *) 'CoordMin,Max=', CoordMin_D(1:nDim), CoordMax_D(1:nDim)
          end if
       end do; end do; end do
    end if

    if (nDim == nDimAmr) then
       if (DoTest) write (*, *) 'Testing interpolate_grid_amr_gc'
       Xyz_D = 0.0
       if (.not. allocated(Point_VIII)) &
            allocate (Point_VIII(0:nVarPoint, nPointI, nPointJ, nPointK))
       Point_VIII = 0.0
       do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
          iPoint_D = [iPoint, jPoint, kPoint]
          XyzPoint_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim) - 0.5) &
               *DomainSize_D(1:nDim)/nPoint_D(1:nDim)
          call interpolate_grid_amr_gc(XyzPoint_D, nCell, iCell_II, Weight_I, &
               IsSecondOrder)

          do iCell = 1, nCell
             Point_VIII(0, iPoint, jPoint, kPoint) = &
                  Point_VIII(0, iPoint, jPoint, kPoint) + Weight_I(iCell)
             iBlock = iCell_II(0, iCell)
             iCell_D = 1
             iCell_D(1:nDim) = iCell_II(1:nDim, iCell)

             ! Interpolate the coordinates to check order of accuracy
             ! Note: Using array syntax in combination with the indirect
             ! iCell_D index fails with optimization for NAG v5.1
             do iDim = 1, nDim
                Xyz_D(iDim) = &
                     Xyz_DGB(iDim, iCell_D(1), iCell_D(2), iCell_D(3), iBlock)
             end do

             ! Take care of periodic dimensions: shift coordinates as necessary
             do iDim = 1, nDim
                if (.not. IsPeriodicTest_D(iDim)) CYCLE
                if (XyzPoint_D(iDim) < Xyz_D(iDim) - 2*CellSize_DB(iDim, iBlock)) &
                     Xyz_D(iDim) = Xyz_D(iDim) - DomainSize_D(iDim)

                if (XyzPoint_D(iDim) > Xyz_D(iDim) + 2*CellSize_DB(iDim, iBlock)) &
                     Xyz_D(iDim) = Xyz_D(iDim) + DomainSize_D(iDim)
             end do

             ! Take care of ghost cells: shift coordinates for coarser neighbor
             iDiscr_D = 0
             where (iCell_D(1:nDim) < 1)
                iDiscr_D(1:nDim) = -1
             elsewhere(iCell_D(1:nDim) > nIJK_D(1:nDim))
                iDiscr_D(1:nDim) = 1
             end where
             ! check that neighbor is coarser
             if (DiLevelNei_IIIB(iDiscr_D(1), iDiscr_D(2), iDiscr_D(3), iBlock) == 1) &
                  Xyz_D = Xyz_D + &
                  0.5*(2*modulo(iCell_D, 2) - 1)*CellSize_DB(:, iBlock)

             ! if point is close to boundary then interpolation isn't 2nd order
             where (.not. (IsSecondOrder .or. IsPeriodicTest_D(1:nDim))) &
                  Xyz_D(1:nDim) = XyzPoint_D(1:nDim)

             Point_VIII(1:nDim, iPoint, jPoint, kPoint) = &
                  Point_VIII(1:nDim, iPoint, jPoint, kPoint) &
                  + Weight_I(iCell)*Xyz_D(1:nDim)
          end do
       end do; end do; end do
       ! Collect contributions from all processors to proc 0
       if (nProc > 1) then
          allocate (PointAll_VIII(0:nVarPoint, nPointI, nPointJ, nPointK))
          call MPI_reduce(Point_VIII, PointAll_VIII, size(PointAll_VIII), &
               MPI_REAL, MPI_SUM, 0, iComm, iError)
          Point_VIII = PointAll_VIII
          deallocate (PointAll_VIII)
       end if

       if (iProc == 0) then
          ! Check interpolated coordinate values against point coordinates
          do kPoint = 1, nPointK; do jPoint = 1, nPointJ; do iPoint = 1, nPointI
             iPoint_D = [iPoint, jPoint, kPoint]
             Xyz_D(1:nDim) = CoordMin_D(1:nDim) + (iPoint_D(1:nDim) - 0.5) &
                  *DomainSize_D(1:nDim)/nPoint_D(1:nDim)

             Weight = Point_VIII(0, iPoint, jPoint, kPoint)
             if (Weight < 1e-6) CYCLE
             Point_V = Point_VIII(1:nDim, iPoint, jPoint, kPoint)/Weight

             if (abs(Weight - 1.0) < 1e-6) then
                Tolerance = 1e-6
             else
                Tolerance = 3e-2
             end if

             if (any(abs(Xyz_D(1:nDim) - Point_V) > Tolerance)) then
                write (*, *) 'ERROR: Point_V=', Point_V(1:nDim), &
                     ' should be ', Xyz_D(1:nDim)
                write (*, *) 'Total weight=', Weight
                write (*, *) 'i,j,kPoint=', iPoint_D(1:nDim)
                write (*, *) 'CoordMin,Max=', CoordMin_D(1:nDim), CoordMax_D(1:nDim)
             end if
          end do; end do; end do
       end if
    end if

    if (nDim == 2) then
       if (DoTest) write (*, *) 'Testing create_grid in RZ geometry'

       ! Store Cartesian values for checking
       allocate (CellVolumeCart_B(MaxBlock), CellFaceCart_DB(MaxDim, MaxBlock))
       CellFaceCart_DB = CellFace_DB
       CellVolumeCart_B = CellVolume_B

       ! Clean Cartesian grid
       call clean_grid

       ! Initialize RZ grid
       call init_geometry(TypeGeometryIn='rz')
       call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
       call create_grid
       if (iProc == 0) call show_grid_proc

       ! Check relative to Cartesian
       Tolerance = 1e-6
       do iBlock = 1, nBlock
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             if (abs(CellVolume_GB(i, j, k, iBlock) &
                  - abs(Xyz_DGB(2, i, j, k, iBlock))*CellVolumeCart_B(iBlock)) &
                  < Tolerance) CYCLE
             write (*, *) NameSub, ' ERROR: incorrect cell volume=', &
                  CellVolume_GB(i, j, k, iBlock), ' should be', &
                  abs(Xyz_DGB(2, i, j, k, iBlock))*CellVolumeCart_B(iBlock), &
                  ' at i,j,k,iBlock,iProc=', i, j, k, iBlock, iProc
          end do; end do; end do
          do iDim = 1, nDim
             Di = i_DD(1, iDim); Dj = i_DD(2, iDim)
             do k = 1, nK; do j = 1, nJ + Dj; do i = 1, nI + Di
                Radius = 0.5*sum(abs(Xyz_DGB(2, i - Di:i, j - Dj:j, k, iBlock)))
                if (abs(CellFace_DFB(iDim, i, j, k, iBlock) - &
                     Radius*CellFaceCart_DB(iDim, iBlock)) &
                     < Tolerance) CYCLE
                write (*, *) NameSub, ' ERROR: incorrect face area=', &
                     CellFace_DFB(iDim, i, j, k, iBlock), ' should be', &
                     Radius*CellFaceCart_DB(iDim, iBlock), &
                     ' at iDim,i,j,k,iBlock,iProc=', &
                     iDim, i, j, k, iBlock, iProc

             end do; end do; end do
          end do
       end do
    end if

    if (nDim >= 2) then
       if (DoTest) write (*, *) 'Testing create_grid in cylindrical geometry'

       ! Clean  grid
       call clean_grid

       ! Initialize cylindrical grid
       call init_geometry(TypeGeometryIn='cylindrical')

       DomainMin_D = [1., 0., -0.5]
       DomainMax_D = [3., 90., 0.5]

       ! This is temporary solution to keep the test working
       IsNodeBasedGrid = .false.
       call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
       call create_grid

       if (iProc == 0) call show_grid_proc

       ! Check relative to generalized coordinate volumes and areas
       Tolerance = 1e-6
       do iBlock = 1, nBlock
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             Good = sqrt(sum(Xyz_DGB(1:2, i, j, k, iBlock)**2)) &
                  *CellVolume_B(iBlock)
             if (abs(CellVolume_GB(i, j, k, iBlock) - Good) < Tolerance) CYCLE
             write (*, *) NameSub, ' ERROR: incorrect cell volume=', &
                  CellVolume_GB(i, j, k, iBlock), ' should be', Good, &
                  ' at i,j,k,iBlock,iProc=', i, j, k, iBlock, iProc
          end do; end do; end do
          do iDim = 1, nDim
             Di = i_DD(1, iDim); Dj = i_DD(2, iDim); Dk = i_DD(3, iDim)
             do k = 1, nK + Dk; do j = 1, nJ + Dj; do i = 1, nI + Di
                ! Calculate face center in generalized coordinates
                Coord_D = CoordMin_DB(:, iBlock) + CellSize_DB(:, iBlock) &
                     *[i - 0.5*(1 + Di), j - 0.5*(1 + Dj), k - 0.5*(1 + Dk)]

                Good = CellFace_DB(iDim, iBlock)
                if (iDim /= 2) Good = Good*Coord_D(1)
                if (abs(CellFace_DFB(iDim, i, j, k, iBlock) - Good) > Tolerance) &
                     write (*, *) NameSub, ' ERROR: incorrect face area=', &
                     CellFace_DFB(iDim, i, j, k, iBlock), ' should be', Good, &
                     ' at iDim,i,j,k,iBlock,iProc=', &
                     iDim, i, j, k, iBlock, iProc

                Phi = Coord_D(2)
                if (iDim == 1) then
                   Good_D = [cos(Phi), sin(Phi), 0.0]
                elseif (iDim == 2) then
                   Good_D = [-sin(Phi), cos(Phi), 0.0]
                else
                   Good_D = [0.0, 0.0, 1.0]
                end if
                ! Multiply by area (for now)
                Good_D = Good_D*CellFace_DFB(iDim, i, j, k, iBlock)

                if (any(Tolerance < abs(FaceNormal_DDFB(:, iDim, i, j, k, iBlock) &
                     - Good_D(1:nDim)))) &
                     write (*, *) NameSub, ' ERROR: incorrect face area=', &
                     FaceNormal_DDFB(:, iDim, i, j, k, iBlock), ' should be', Good_D, &
                     ' at iDim,i,j,k,iBlock,iProc=', &
                     iDim, i, j, k, iBlock, iProc

             end do; end do; end do
          end do
       end do
    end if

    if (nDim == 3) then
       if (DoTest) write (*, *) 'Testing create_grid in spherical geometry'

       ! Clean  grid
       call clean_grid

       ! Initialize spherical grid
       call init_geometry(TypeGeometryIn='spherical')

       DomainMin_D = [1., 45.0, 0.0]
       DomainMax_D = [3., 135.0, 90.0]

       ! This is temporary solution to keep the test working
       IsNodeBasedGrid = .false.
       call init_grid(DomainMin_D, DomainMax_D)
       call create_grid

       if (iProc == 0) call show_grid_proc

       if (DoTest) write (*, *) 'Testing create_grid in rlonlat geometry'

       ! Clean grid
       call clean_grid

       ! Initialize r-lon-lat grid
       call init_geometry(TypeGeometryIn='rlonlat')

       DomainMin_D = [1., 0.0, -45.0]
       DomainMax_D = [3., 90.0, 45.0]

       ! This is temporary solution to keep the test working
       IsNodeBasedGrid = .false.
       call init_grid(DomainMin_D, DomainMax_D)
       call create_grid

       if (iProc == 0) call show_grid_proc

       if (DoTest) write (*, *) 'Testing create_grid in roundcube geometry'

       ! Clean  grid
       call clean_grid

       ! Initialize roundcube grid
       call init_geometry(TypeGeometryIn='roundcube')

       rRound0 = 2.0
       rRound1 = 3.2

       DomainMin_D = [-3.2, -3.2, -3.2]
       DomainMax_D = [3.2, 3.2, 3.2]

       IsNodeBasedGrid = .true.
       call init_grid(DomainMin_D, DomainMax_D)
       call create_grid

       if (iProc == 0) call show_grid_proc

       ! Check total volume.
       ! It should be approximately the volume of a sphere of radius 3.2
       Volume = sum(CellVolume_GB(1:nI, 1:nJ, 1:nK, 1:nBlock))
       if (nProc > 1) then
          call MPI_reduce(Volume, VolumeAll, 1, MPI_REAL, MPI_SUM, 0, iComm, iError)
          Volume = VolumeAll
       end if

       if (iProc == 0) then
          ! Analytic volume of the sphere
          VolumeAll = 4./3.*cPi*(sqrt(3.)*rRound1)**3

          if (abs(VolumeAll - Volume)/VolumeAll > 0.02) &
               write (*, *) 'ERROR: total volume numerical vs analytic:', &
               Volume, VolumeAll
       end if
    end if

    if (DoTest) write (*, *) 'Testing clean_grid'
    call clean_grid
    call clean_tree

    if (nDim == nDimAmr .and. nDim > 1) then
       if (DoTest) write (*, *) 'Testing check_interpolate_amr_gc'

       DomainMin_D = [0.0, 0.0, 0.0]
       DomainMax_D = [1.0, 1.0, 1.0]
       Coord_D = [0.0, 0.0, 0.0]

       ! go over all geometry cases
       do iCase = 1, nCaseCheck
          if (DoTest) &
               write (*, '(a,i3)') '  Testing case ', iCase
          ! Set Cartesian grid geometry before initializing tree and grid
          call init_geometry(IsPeriodicIn_D=spread([.false.], 1, nDim))
          call init_tree(nBlockCheckMax)
          call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
          call set_tree_root(spread([2], 1, nDim))

          do iNodeCheck = 1, nRootCheck
             if (BTEST(iCase - 1, iNodeCheck - 1)) &
                  call refine_tree_node(iNodeCheck)
          end do

          call distribute_tree(.true.)

          call create_grid

          ! generate random points and test check_interpolate_amr_gc
          do iPoint = 1, nPointCheck
             do iDim = 1, nDim
                Coord_D(iDim) = &
                     0.75*DomainMin_D(iDim) + 0.25*DomainMin_D(iDim) + &
                     0.5*(DomainMax_D(iDim) - DomainMin_D(iDim))* &
                     random_real(iSeed)
             end do
             call check_interpolate_amr_gc(Coord_D, 1, iProcCheck, iBlockCheck)
             do iBlock = 2, nBlock
                call check_interpolate_amr_gc( &
                     Coord_D, iBlock, iProcOut, iBlockOut)
                if (iProcOut /= iProcCheck .or. iBlockOut /= iBlockCheck) then
                   write (*, '(a,3f11.8,a,i3,a,2i3,a,i3,a,2i3)') &
                        'ERROR: for point ', Coord_D, &
                        ' result from iBlock =', 1, ' is iProcOut, iBlockOut =', &
                        iProcCheck, iBlockCheck, ' but from iBlock =', iBlock, &
                        ' is iProcOut, iBlockOut =', iProcOut, iBlockOut
                end if
             end do
          end do

          call clean_grid
          call clean_tree
       end do
    end if

  end subroutine test_grid
  !============================================================================

  subroutine test_pass_node ! unit test

    ! To test the message pass for the node we generate a fine uniform grid
    ! for the whole domain and transfer the node values from the
    ! block nodes to the nodes on the fine grid. Then we gather all the
    ! data on the fine grid with the proper operator.
    ! We can then compare the values on the coresponding node after
    ! message_pass_node is called with the fine grid values.

    integer, parameter:: MaxBlockTest = 50
    integer, parameter:: nRootTest_D(MaxDim) = [3, 3, 3]
    logical, parameter:: IsPeriodicTest_D(MaxDim) = .false.
    real, parameter:: DomainMin_D(MaxDim) = [0.0, 0.0, 0.0]
    real, parameter:: DomainMax_D(MaxDim) = [48.0, 48.0, 48.0]

    integer, parameter:: nVar = 3! nDim
    integer :: iVar

    character(len=4) :: NameOperator = "Min"

    integer ::iOp
    integer, parameter :: nOp = 3
    character(len=4) :: NameOperator_I(nOp) = ["mean", "min ", "max "]

    ! Variable on the nodes
    real, allocatable:: State_VNB(:, :, :, :, :)

    integer, allocatable:: i_NB(:, :, :, :)

    real, allocatable, dimension(:, :, :, :, :) :: Xyz_DNB

    real, allocatable, dimension(:, :, :, :) :: FineGridLocal_IIIV
    real, allocatable, dimension(:, :, :, :) :: FineGridGlobal_IIIV
    real :: FineGridStep_D(MaxDim)
    integer :: iFG, jFG, kFG
    integer :: nFineCell
    integer :: iMpiOperator
    integer :: iError

    integer:: iNode, iBlock, i, j, k

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'test_pass_node'
    !--------------------------------------------------------------------------
    DoTest = iProc == 0

    if (DoTest) write (*, *) 'Starting ', NameSub

    call init_tree(MaxBlockTest)
    call init_geometry(IsPeriodicIn_D=IsPeriodicTest_D(1:nDim))
    call init_grid(DomainMin_D(1:nDim), DomainMax_D(1:nDim))
    call set_tree_root(nRootTest_D(1:nDim))

    call find_tree_node([0.5, 0.5, 0.5], iNode)
    call refine_tree_node(iNode)
    call distribute_tree(.true.)
    call create_grid

    allocate (Xyz_DNB(MaxDim, nINode, nJNode, nKNode, MaxBlockTest))
    allocate (State_VNB(nVar, nI + 1, nJ + 1, nK + 1, MaxBlockTest))
    allocate (i_NB(nI + 1, nJ + 1, nK + 1, MaxBlockTest))

    allocate (FineGridLocal_IIIV( &
         nI*iRatio*nRootTest_D(1) + 1, &
         nJ*jRatio*nRootTest_D(2) + 1, &
         nK*kRatio*nRootTest_D(3) + 1, &
         nVar + 1))
    allocate (FineGridGlobal_IIIV( &
         nI*iRatio*nRootTest_D(1) + 1, &
         nJ*jRatio*nRootTest_D(2) + 1, &
         nK*kRatio*nRootTest_D(3) + 1, &
         nVar + 1))

    nFineCell = product(nIJK_D*iRatio_D*nRootTest_D + 1)

    FineGridStep_D = (nIJK_D*iRatio_D*nRootTest_D) &
         /(DomainMax_D - DomainMin_D)

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE
       do k = 1, nKNode; do j = 1, nJNode; do i = 1, nINode
          Xyz_DNB(:, i, j, k, iBlock) = CoordMin_DB(:, iBlock) + &
               ([i, j, k] - 1.0)*CellSize_DB(:, iBlock)
       end do; end do; end do
    end do

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE
       iNode = (iNode_B(iBlock) - 1)*(nI + 1)*(nK + 1)*(nJ + 1)
       do k = 1, nK + 1; do j = 1, nJ + 1; do i = 1, nI + 1
          iNode = iNode + 1
          i_NB(i, j, k, iBlock) = iNode
       end do; end do; end do
    end do

    do iOp = 1, nOp

       NameOperator = NameOperator_I(iOp)

       if (DoTest) write (*, *) 'testing message_pass_node with operator=', &
            NameOperator

       select case (NameOperator)
       case ("mean")
          FineGridLocal_IIIV(:, :, :, :) = 0.0
          FineGridGlobal_IIIV(:, :, :, :) = 0.0
          iMpiOperator = MPI_SUM
       case ("min")
          FineGridLocal_IIIV(:, :, :, :) = 1.0e8
          FineGridGlobal_IIIV(:, :, :, :) = 1.0e8
          iMpiOperator = MPI_MIN
       case ("max")
          FineGridLocal_IIIV(:, :, :, :) = -1.0e8
          FineGridGlobal_IIIV(:, :, :, :) = -1.0e8
          iMpiOperator = MPI_MAX
       case default
          call CON_stop(NameSub//' incorrect operator name='//NameOperator)
       end select

       State_VNB(1, :, :, :, :) = i_NB(:, :, :, :)
       do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          do k = 1, nKNode; do j = 1, nJNode; do i = 1, nINode
             iNode = iNode_B(iBlock) - 1
             State_VNB(2:3, i, j, k, iBlock) = [1.0, real(iNode)]
          end do; end do; end do
       end do

       do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          do k = 1, nKNode; do j = 1, nJNode; do i = 1, nINode

             iFG = int(Xyz_DNB(1, i, j, k, iBlock)*FineGridStep_D(1)) + 1
             jFG = int(Xyz_DNB(2, i, j, k, iBlock)*FineGridStep_D(2)) + 1
             kFG = int(Xyz_DNB(3, i, j, k, iBlock)*FineGridStep_D(3)) + 1

             select case (NameOperator)
             case ("mean")
                FineGridLocal_IIIV(iFG, jFG, kFG, 1:nVar) = &
                     FineGridLocal_IIIV(iFG, jFG, kFG, 1:nVar) + &
                     State_VNB(:, i, j, k, iBlock)
                FineGridLocal_IIIV(iFG, jFG, kFG, nVar + 1) = &
                     FineGridLocal_IIIV(iFG, jFG, kFG, nVar + 1) + 1
             case ("min")
                do iVar = 1, nVar
                   FineGridLocal_IIIV(iFG, jFG, kFG, iVar) = min( &
                        FineGridLocal_IIIV(iFG, jFG, kFG, iVar), &
                        State_VNB(iVar, i, j, k, iBlock))
                end do
             case ("max")
                do iVar = 1, nVar
                   FineGridLocal_IIIV(iFG, jFG, kFG, iVar) = max( &
                        FineGridLocal_IIIV(iFG, jFG, kFG, iVar), &
                        State_VNB(iVar, i, j, k, iBlock))
                end do
             end select
          end do; end do; end do
       end do

       call message_pass_node(nVar, State_VNB, NameOperator_I(iOp))

       call MPI_ALLREDUCE(FineGridLocal_IIIV(1, 1, 1, 1), &
            FineGridGlobal_IIIV(1, 1, 1, 1), &
            nFineCell*(nVar + 1), &
            MPI_REAL, iMpiOperator, iComm, iError)

       do iBlock = 1, nBlock
          if (Unused_B(iBlock)) CYCLE
          do k = 1, nKNode; do j = 1, nJNode; do i = 1, nINode
             iFG = nint(Xyz_DNB(1, i, j, k, iBlock)*FineGridStep_D(1)) + 1
             jFG = nint(Xyz_DNB(2, i, j, k, iBlock)*FineGridStep_D(2)) + 1
             kFG = nint(Xyz_DNB(3, i, j, k, iBlock)*FineGridStep_D(3)) + 1
             do iVar = 1, nVar
                select case (NameOperator)
                case ("mean")
                   if (FineGridGlobal_IIIV(iFG, jFG, kFG, iVar)/ &
                        FineGridGlobal_IIIV(iFG, jFG, kFG, nVar + 1) /= &
                        State_VNB(iVar, i, j, k, iBlock)) then
                      write (*, *) "Error for operator, variable, iBlock= ", &
                           NameOperator, iVar, iBlock, ", value=", &
                           FineGridGlobal_IIIV(iFG, jFG, kFG, iVar)/ &
                           FineGridGlobal_IIIV(iFG, jFG, kFG, nVar + 1), &
                           " should be ", State_VNB(iVar, i, j, k, iBlock)
                   end if
                case ("min", "max")
                   if (FineGridGlobal_IIIV(iFG, jFG, kFG, iVar) /= &
                        State_VNB(iVar, i, j, k, iBlock)) then
                      write (*, *) "Error for operator, variable, iBlock= ", &
                           NameOperator, iVar, iBlock, ", value=", &
                           FineGridGlobal_IIIV(iFG, jFG, kFG, iVar), &
                           " should be ", State_VNB(iVar, i, j, k, iBlock)
                   end if
                end select
             end do
          end do; end do; end do
       end do

    end do

    deallocate (FineGridLocal_IIIV)
    deallocate (FineGridGlobal_IIIV)
    deallocate (Xyz_DNB)
    call clean_grid
    call clean_tree

  end subroutine test_pass_node
  !============================================================================
  subroutine test_amr_criteria ! unit test

!!$    use BATL_size, ONLY : MaxBlock,nBlock, iRatio, jRatio, kRatio
!!$    use BATL_mpi, ONLY: iProc, nProc
!!$    use BATL_tree, ONLY: init_tree, set_tree_root, find_tree_node, &
!!$         refine_tree_node, distribute_tree, clean_tree, Unused_B, iNode_B, &
!!$         nNode, show_tree, iStatusNew_A, &
!!$         Coarsen_, Unset_, adapt_tree, move_tree, distribute_tree
!!$    use BATL_grid, ONLY: init_grid, create_grid, clean_grid, CellSize_DB
!!$    use BATL_geometry, ONLY: init_geometry
!!$    use BATL_amr, ONLY: do_amr, init_amr
!!$    use BATL_size, ONLY: MaxDim, nDim, MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
!!$         nI, nJ, nK
!!$    ! For Random generation
!!$    integer :: jSeed
!!$    logical :: IsFirst
!!$    integer, parameter :: iMPLIER=16807, &
!!$         iMODLUS=2147483647, &
!!$         iMOBYMP=127773, &
!!$         iMOMDMP=2836
!!$
!!$    integer, parameter:: MaxBlockTest            = 50
!!$    integer, parameter:: nRootTest_D(MaxDim)     = (/3,3,3/)
!!$    logical, parameter:: IsPeriodicTest_D(MaxDim)= .false.
!!$    real, parameter:: DomainMin_D(MaxDim) = (/ -24.0, -24.0, -24.0 /)
!!$    real, parameter:: DomainMax_D(MaxDim) = (/ 24.0, 24.0, 24.0 /)
!!$    integer :: iNode, iBlock
!!$
!!$    real, allocatable :: Criterias_IB(:,:),AllCriterias_IBP(:,:,:)
!!$    real, allocatable :: PreCriterias_IB(:,:)
!!$    real, allocatable :: RefineLevel_I(:), CoursenLevel_I(:)
!!$    logical, allocatable :: Used_GB(:,:,:,:)
!!$
!!$    integer :: nCritExt = 4
!!$    integer, allocatable :: iA_I(:)
!!$    real, allocatable :: TestState_VGB(:,:,:,:,:)
!!$    integer :: nVar =3
!!$    integer :: i,j,k,iVar
!!$    logical:: DoTestMe
!!$    character(len=*), parameter :: NameSub = 'test_amr_criteria'
!!$    !-----------------------------------------------------------------------
!!$    DoTestMe = iProc == 0
!!$
!!$    write(*,*) " Temural not testing :: test_amr_criteria"
!!$    RETURN
!!$
!!$    if(DoTestMe) write(*,*) 'Starting ',NameSub
!!$
!!$    call init_tree(MaxBlockTest)
!!$    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
!!$    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
!!$    call set_tree_root( nRootTest_D(1:nDim))
!!$
!!$    call find_tree_node( (/0.5,0.5,0.5/), iNode)
!!$    call refine_tree_node(iNode)
!!$    call distribute_tree(.true.)
!!$    call create_grid
!!$    call init_amr
!!$
!!$    call srand(123456789+iProc)
!!$
!!$    allocate(Criterias_IB(nCritExt,nBlock), &
!!$         AllCriterias_IBP(nCritExt,nBlock,nProc),&
!!$         PreCriterias_IB(nCritExt,nBlock))
!!$    allocate(RefineLevel_I(nCritExt),CoursenLevel_I(nCritExt))
!!$
!!$
!!$    !===================== Begin Test  =======================
!!$    DoSortAmrCrit = .true.
!!$    !--------------------- internal --------------------------
!!$    nIntCrit = 1
!!$    nExtCritUsed = 0
!!$    allocate(CoarsenCritAll_I(nVar), &
!!$         RefineCritAll_I(nVar),iVarCritAll_I(nVar),iMapToStateVar_I(nVar))
!!$
!!$    allocate(TestState_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
!!$    CoarsenCritAll_I = -1.0
!!$    RefineCritAll_I  =  1.0
!!$    TestState_VGB = 1.0
!!$    iVarCritAll_I(nIntCrit)=1
!!$    iMapToStateVar_I(nIntCrit)=1
!!$    nAmrCritUsed = 1
!!$    nAmrCrit = nIntCrit + nExtCritUsed
!!$    ! allocate(AmrCrit_IB(nIntCrit+nExtCritUsed,nBlock))
!!$    ! AmrCrit_IB = 0.0
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB)
!!$
!!$    ! do iBlock = 1, nBlock
!!$    !   if(Unused_B(iBlock)) CYCLE
!!$    !   write(*,'(i6,f16.12)') iNode_B(iBlock),AmrCrit_IB(1,iBlock)
!!$    ! end do
!!$    !
!!$    ! do iNode = 1, nNode-1
!!$    !   print *,iNode, iRank_A(iNode)
!!$    ! end do
!!$
!!$    do iNode = 2, nNode-1
!!$       if(iRank_A(iNode-1) > iRank_A(iNode)) then
!!$          write(*,*) " ERROR in ",NameSub, " Internal test"
!!$       end if
!!$    end do
!!$
!!$    ! deallocate(AmrCrit_IB)
!!$
!!$
!!$
!!$    !-------------------- external --------------------------
!!$    Criterias_IB = 0.0
!!$    RefineLevel_I =  1.0
!!$    CoursenLevel_I = -1.0
!!$    ! external criteria copyed into the global criteia
!!$    CoarsenCritAll_I = CoursenLevel_I
!!$    RefineCritAll_I  = RefineLevel_I
!!$
!!$    nExtCritUsed = 1
!!$    nIntCrit = 0
!!$    nCritExt = nExtCritUsed
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) then
!!$          Criterias_IB(1,iBlock) = 10.0
!!$       else
!!$          Criterias_IB(1,iBlock) = AmrCrit_IB(1,iBlock)
!!$       end if
!!$    end do
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB, &
!!$         nCritExt, Criterias_IB)
!!$
!!$    do iNode = 2, nNode-1
!!$       if(iRank_A(iNode-1)>iRank_A(iNode)) then
!!$          write(*,*) " ERROR in ",NameSub, " External=Intenal test"
!!$       end if
!!$    end do
!!$
!!$    !--------------------- internal with masked cells --------
!!$    nIntCrit = 1
!!$    nExtCritUsed = 0
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed
!!$
!!$    CoarsenCritAll_I = -1.0
!!$    RefineCritAll_I  =  1.0
!!$    TestState_VGB = 1.0
!!$    iVarCritAll_I(nIntCrit)=1
!!$    iMapToStateVar_I(nIntCrit)=1
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    UseAmrMask = .true.
!!$    nAmrBox = 1
!!$    allocate(AmrBox_DII(3,2,nAmrBox))
!!$    AmrBox_DII(1,1,1) = DomainMin_D(1)
!!$    AmrBox_DII(1,2,1) = DomainMin_D(1) + CellSize_DB(1,1)*nI
!!$    AmrBox_DII(2,1,1) = DomainMin_D(2)
!!$    AmrBox_DII(2,2,1) = DomainMin_D(2) + CellSize_DB(2,1)*nJ
!!$    AmrBox_DII(3,1,1) = DomainMin_D(3)
!!$    AmrBox_DII(3,2,1) = DomainMin_D(3) + CellSize_DB(3,1)*nK
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB)
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       if( iNode_B(iBlock) == 1) then
!!$          if( AmrCrit_IB(1,iBlock)  == 0.0) &
!!$               write (*,*) " ERROR in ",NameSub, &
!!$               " in  Internal test masked cells", &
!!$               " AmrCrit_IB of Node == 1 shoud be none zero"
!!$       else
!!$          if( AmrCrit_IB(1,iBlock)  /= 0.0) &
!!$               write (*,*) " ERROR in ",NameSub, &
!!$               " in  Internal test masked cells", &
!!$               " AmrCrit_IB of Node /= 1 shoud be zero"
!!$       end if
!!$    end do
!!$
!!$   ! Using any becouse of the ghost cells
!!$   ! do iBlock = 1, nBlock
!!$   !    if(Unused_B(iBlock)) CYCLE
!!$   !    write(*,*) iNode_B(iBlock), &
!!$   !         any(DoAmr_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,iBlock))
!!$   ! end do
!!$
!!$
!!$    UseAmrMask = .false.
!!$    deallocate(AmrBox_DII)
!!$    deallocate(DoAmr_GB)
!!$    !--------------------- internal with masked cells --------
!!$
!!$    !--------------------- internal with masked body -------
!!$    nIntCrit = 1
!!$    nExtCritUsed = 0
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed
!!$
!!$    allocate(Used_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
!!$    Used_GB = .true.
!!$
!!$    CoarsenCritAll_I = -1.0
!!$    RefineCritAll_I  =  1.0
!!$    TestState_VGB = 1.0
!!$    iVarCritAll_I(nIntCrit)=1
!!$    iMapToStateVar_I(nIntCrit)=1
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    do iBlock = 1, nBlock
!!$       if(iNode_B(iBlock) /= nNode)  CYCLE
!!$
!!$       Used_GB(1,1,1,iBlock) = .false.
!!$
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = dexp(0.1*(i*0.5))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$
!!$       TestState_VGB(:,1,1,1,iBlock) = 1.0e18
!!$    end do
!!$
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB,Used_GB=Used_GB)
!!$
!!$    if(iRank_A(1) /= nNode) &
!!$         write(*,*) " ERROR in ",NameSub, " Internal test masked body"
!!$
!!$    do iNode = 3, nNode-1
!!$       if(iRank_A(iNode-1) > iRank_A(iNode)) then
!!$          write(*,*) " ERROR in ",NameSub, " Internal test masked body"
!!$       end if
!!$    end do
!!$
!!$
!!$
!!$    deallocate(Used_GB)
!!$    !-------------------- internal x 2 -------------------
!!$
!!$    if(iRatio > 1 .and. jRatio > 1 .and. kRatio >1 ) then
!!$
!!$       nCritExt = 0
!!$       nIntCrit = 2
!!$       nExtCritUsed = nCritExt
!!$       nAmrCritUsed = nIntCrit + nExtCritUsed
!!$       iVarCritAll_I(1:nIntCrit)=(/ 1,2 /)
!!$       iMapToStateVar_I =  iVarCritAll_I
!!$       ! external criteria copyed into the global criteia
!!$       CoarsenCritAll_I = CoursenLevel_I
!!$       RefineCritAll_I  = RefineLevel_I
!!$
!!$       do iBlock=1,nBlock
!!$          if(Unused_B(iBlock)) CYCLE
!!$          do k = MinK, MaxK
!!$             do j = MinJ, MaxJ
!!$                do i = MinI,MaxI
!!$                   TestState_VGB(1,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(35-iNode_B(iBlock)+1)))
!!$
!!$                   TestState_VGB(2,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$
!!$       do iBlock = 1, nBlock
!!$          if(Unused_B(iBlock)) then
!!$             Criterias_IB(1,iBlock) = 10.0
!!$          else
!!$             Criterias_IB(1,iBlock) = AmrCrit_IB(1,nBlock-iBlock+1)
!!$          end if
!!$
!!$       end do
!!$
!!$
!!$       call set_amr_criteria(nVar,TestState_VGB, &
!!$            nCritExt, Criterias_IB)
!!$
!!$       allocate(iA_I(nNode-1))
!!$
!!$       iA_I =(/ 1,35, 2, 34, 33,  3,  4, 32, 31,  5,  6, 30, 29,  7, &
!!$            8, 28, 27,  9, 26, 10, 25, 11, 12, 24, 13, 23, 22, 15, 21, &
!!$            16, 17, 20, 19, 18 /)
!!$
!!$       ! do iBlock = 1, nBlock
!!$       !   if(Unused_B(iBlock)) CYCLE
!!$       !   print *,"Criterias_IB(1,iBlock) = ", &
!!$       !        AmrCrit_IB(1:2,iBlock), iNode_B(iBlock)
!!$       ! end do
!!$
!!$
!!$       ! do iNode = 1, nNode-1
!!$       !   print *,iRank_A(iNode)," :: ", iA_I(iNode)
!!$       ! end do
!!$
!!$       do iNode = 1, nNode-1, 2
!!$          if(iRank_A(iNode) /= iA_I(iNode)) &
!!$               write(*,*) " ERROR in ",NameSub, "2 x Intenal test"
!!$       end do
!!$
!!$       deallocate(iA_I)
!!$    end if
!!$    !-------------------- testing levels ---------------------
!!$    ! intenals
!!$    CoarsenCritAll_I(1) = -1.0
!!$    RefineCritAll_I(1)  =  1.0
!!$    ! externals
!!$    RefineLevel_I(1) =  0.030
!!$    CoursenLevel_I(1) = 0.020
!!$
!!$
!!$    ! external criteria copyed into the global criteia
!!$    CoarsenCritAll_I(2) = CoursenLevel_I(1)
!!$    RefineCritAll_I(2)  = RefineLevel_I(1)
!!$
!!$    PercentRefine  = 30.0
!!$    PercentCoarsen = 10.0
!!$    nCritExt = 1
!!$    nIntCrit = 1
!!$
!!$    nExtCritUsed = nCritExt
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed
!!$
!!$    ! init Internals
!!$    iVarCritAll_I(nIntCrit) = 1
!!$    iMapToStateVar_I(nIntCrit) = 1
!!$    do iBlock=1,nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                ! do iVar=1,nVar
!!$                ! print *,nrand()
!!$                TestState_VGB(2,i,j,k,iBlock) = &
!!$                     dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                ! end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    ! init Externals
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       Criterias_IB(1,iBlock) = 0.001*(nNode - iNode_B(iBlock)+1)
!!$    end do
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB, &
!!$         nCritExt, Criterias_IB)
!!$
!!$    do k = nNodeCoarsen+1, nNode-1-nNodeRefine
!!$       if( iRank_A(k) < nNodeRefine .and. &
!!$            iRank_A(k) > nNode-1 - nNodeCoarsen ) &
!!$            write(*,*) "Error in seting refinment by criteria"
!!$    end do
!!$
!!$    ! test unused block
!!$    Criterias_IB = 0.0
!!$    RefineLevel_I =  1.0
!!$    CoursenLevel_I = -1.0
!!$    nCritExt = 0
!!$    nIntCrit = 1
!!$
!!$    iStatusNew_A = Unset_
!!$    nExtCritUsed = nCritExt
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed
!!$
!!$    do iBlock=1,nBlock
!!$       if(iNode_B(iBlock) > 3**nDim) then
!!$          iStatusNew_A(iNode_B(iBlock)) =  Coarsen_
!!$       end if
!!$    end do
!!$    ! call show_tree(NameSub,.true.)
!!$    call adapt_tree
!!$    call distribute_tree(DoMove=.false.)
!!$    call do_amr(nVar,TestState_VGB)
!!$    call move_tree
!!$    ! call show_tree(NameSub,.true.)
!!$
!!$    do iBlock=1,nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB)
!!$
!!$    do iNode = 2, nNode
!!$       if(iRank_A(iNode-1)>iRank_A(iNode)) then
!!$          write(*,*) " ERROR in ",NameSub, " unused  test"
!!$       end if
!!$    end do
!!$
!!$
!!$    !===================== End Test  =======================
!!$    deallocate(CoarsenCritAll_I,RefineCritAll_I,iVarCritAll_I,iMapToStateVar_I)
!!$    deallocate(TestState_VGB)
!!$
!!$    deallocate(Criterias_IB,PreCriterias_IB,AllCriterias_IBP)
!!$    deallocate(RefineLevel_I,CoursenLevel_I)
!!$    call clean_grid
!!$    call clean_tree
!!$
!!$  contains
!!$
!!$    ! The saudo random number generator is for testing performense in
!!$    ! parallel sorting
!!$    subroutine srand(iSeed)
!!$      integer, intent(in) :: iSeed
!!$      jSeed = iSeed
!!$      IsFirst = .true.
!!$    end subroutine srand
!!$
!!$    real function rand()
!!$      !  A pseudo-random number generator implemented to make sure that
!!$      !  all platform reproduce the same sequence for testing and compering.
!!$      !  The algorithm is based on "Integer Version 2" given in :
!!$      !
!!$      !       Park, Steven K. and Miller, Keith W., "Random Number Generators:
!!$      !       Good Ones are Hard to Find", Communications of the ACM,
!!$      !       October, 1988.
!!$
!!$      integer :: nHvalue,nLvalue,nTestv
!!$      integer, save :: nExtn
!!$
!!$      if(IsFirst) then
!!$         nExtn=jSeed
!!$         IsFirst = .false.
!!$      end if
!!$
!!$      nHvalue = nExtn/iMOBYMP
!!$      nLvalue = mod(nExtn,iMOBYMP)
!!$      nTestv = iMPLIER*nLvalue - iMOMDMP*nHvalue
!!$      if(nTestv > 0) then
!!$         nExtn = nTestv
!!$      else
!!$         nExtn = nTestv + iMODLUS
!!$      end if
!!$
!!$      rand = real(nExtn)/real(iMODLUS)
!!$
!!$    end function rand

    !--------------------------------------------------------------------------
  end subroutine test_amr_criteria
  !============================================================================

end program unit_test

include 'external_routines.f90'
