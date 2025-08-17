!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_lib

  ! Collection of all public methods and data that an application can access

  use BATL_size
  use BATL_mpi
  use BATL_tree
  use BATL_geometry
  use BATL_grid
  use BATL_region
  use BATL_amr
  use BATL_amr_criteria
  use BATL_pass_cell
  use BATL_pass_face
  use BATL_pass_node
  use BATL_high_order
  use BATL_particles
  use BATL_test

  implicit none

  private ! except

  ! Public methods and variables of this module
  public:: init_batl
  public:: clean_batl
  public:: init_grid_batl
  public:: regrid_batl

  integer, public, allocatable:: iVectorVar_I(:)

  logical, public:: IsBatlInitialized = .false.

  !--------------------
  ! Variables inherited
  !--------------------

  ! Inherited from BATL_size
  public:: MaxDim, nDim, Dim1_, Dim2_, Dim3_, iDim_, jDim_, kDim_
  public:: iRatio, jRatio, kRatio
  public:: nDimAmr, iDimAmr_D
  public:: MaxBlock, nBlock
  public:: nI, nJ, nK, nIJK, nIJK_D, MinIJK_D, MaxIJK_D
  public:: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG
  public:: j0_, j2_, nJp1_, nJm1_, k0_, k2_, nKp1_, nKm1_, i0_, nIp1_
  public:: nINode, nJNode, nKNode
  public:: iRatio_D

  ! Inherited from BATL_amr_criteria
  public:: AmrCrit_IB, nAmrCrit, DoCritAmr, DoAutoAmr, DoStrictAmr

  ! Inherited from BATL_test
  public:: lVerbose, StringTest
  public:: UseTestCell, UseTestXyz
  public:: iProcTest, iBlockTest, iTest, jTest, kTest
  public:: XyzTestCell_D, xTest, yTest, zTest
  public:: UseTest2Cell, UseTest2Xyz
  public:: iProcTest2, iBlockTest2, iTest2, jTest2, kTest2
  public:: XyzTestCell2_D, xTest2, yTest2, zTest2
  public:: iVarTest, NameVarTest, iDimTest, iSideTest

  ! Inherited from BATL_region
  public:: nInitialAmrLevel

  ! Inherited from BATL_particles
  public:: Particle_I

  ! Inherited from BATL_amr
  public:: BetaProlong

  ! Inherited from BATL_mpi
  public:: init_mpi, clean_mpi, barrier_mpi
  public:: iComm, nProc, iProc, nThread, iThread

  ! Inherited from BATL_tree
  public:: MaxNode, nNode, nNodeUsed, nRoot_D, nRoot
  public:: MaxLevel, nLevelMin, nLevelMax, MaxCoord_I
  public:: Unused_B, Unused_BP
  public:: iNode_B, iMortonNode_A, iNodeMorton_I
  public:: DiLevelNei_IIIB, iNodeNei_IIIB, IsNeighbor_P
  public:: iStatusNew_A, Refine_, Coarsen_, Unset_
  public:: iTree_IA, Status_,  Level_, MinLevel_, MaxLevel_, Block_, Proc_
  public:: Coord0_, Coord1_, Coord2_, Coord3_, Used_
  public:: UseTimeLevel, nTimeLevel, iTimeLevel_A
  public:: IsNewDecomposition, IsNewTree
  public:: iAmrChange_B
  public:: AmrRemoved_, AmrUnchanged_, AmrMoved_, AmrRefined_, AmrCoarsened_

  ! Inherited from BATL_geometry
  public:: TypeGeometry, IsCartesianGrid, IsCartesian, IsRzGeometry
  public:: IsRotatedCartesian, GridRot_DD
  public:: IsSpherical, IsRLonLat, IsCylindrical
  public:: IsCylindricalAxis, IsSphericalAxis, IsLatitudeAxis, IsAnyAxis
  public:: x_, y_, z_, r_, Phi_, Theta_, Lon_, Lat_
  public:: IsLogRadius, IsGenRadius, nRgen, LogRgen_I
  public:: IsPeriodic_D, IsPeriodicCoord_D
  public:: Xi_, Eta_, Zeta_
  public:: UseHighFDGeometry
  public:: rRound0, rRound1, IsRoundCube, SqrtNDim
  public:: IsCubedSphere, iCubedSphere

  ! Inherited from BATL_grid
  public:: IsNodeBasedGrid
  public:: CoordMin_D, CoordMax_D, DomainSize_D
  public:: CoordMin_DB, CoordMax_DB, CellSize_DB, CellSizeRoot
  public:: Xyz_DGB, Xyz_DNB, Used_GB
  public:: CellFace_DB, CellFace_DFB, FaceNormal_DDFB
  public:: CellVolume_B, CellVolume_GB
  public:: find_grid_block, interpolate_grid, average_grid_node
  public:: CellMetrice_DDG, CellCoef_DDGB
  public:: check_interpolate_amr_gc
  public:: create_grid, show_grid_cell
  public:: integrate_grid, minval_grid, maxval_grid
  public:: interpolate_state_vector ! Interpolate state vactor to a given loc,
                                    ! Need all ghost cells with ProlongOrder=1

  ! Inherited from BATL_pass_cell
  public:: UseSimpleProlongation

  !--------------------
  ! Functions inherited
  !--------------------

  ! Inherited from BATL_amr
  public:: sync_cpu_gpu_amr

  ! Inherited from BATL_amr_criteria
  public:: set_amr_criteria, clean_amr_criteria, read_amr_criteria
  public:: calc_error_amr_criteria, set_amr_geometry
  public:: is_masked_amr_criteria, init_amr_criteria

  ! Inherited from BATL_region
  public:: read_region_param, get_region_indexes, block_inside_regions

  ! Inherited from BATL_pass_cell
  public:: message_pass_cell

  ! Inherited from BATL_pass_face
  public:: message_pass_face
  public:: store_face_flux, correct_face_flux
  public:: apply_flux_correction, apply_flux_correction_block

  ! Inherited from BATL_pass_node
  public:: message_pass_node

  ! Inherited from BATL_high_order
  public:: correct_face_value, calc_center_first_derivate, calc_face_value

  ! Inherited from BATL_particles
  public:: check_particle_location
  public:: put_particles, trace_particles
  public:: message_pass_particles, remove_undefined_particles, mark_undefined

  ! Inherited from BATL_tree
  public:: init_tree        ! initialize tree
  public:: clean_tree       ! clean tree data
  public:: set_tree_param   ! set parameters like UseUniformAxis
  public:: set_tree_root    ! set root nodes
  public:: refine_tree_node ! refined one tree node
  public:: coarsen_tree_node! coarsen one tree node
  public:: adapt_tree       ! refine and coarsen the whole tree
  public:: distribute_tree  ! distribute tree nodes over processors
  public:: move_tree        ! finish load balance and compact tree
  public:: get_tree_position! get node position in the domain
  public:: find_tree_node   ! find tree node containing a given point
  public:: find_tree_cell   ! find tree cell containing a given point
  public:: interpolate_tree ! not yet complete
  public:: write_tree_file  ! save tree info for restart
  public:: read_tree_file   ! read tree info for restart
  public:: min_tree_level   ! return min level for a subcycling stage
  public:: set_tree_periodic! switch periodicity info on/off as needed
  public:: show_tree        ! show info for debugging
  public:: find_neighbor_for_anynode
  public:: i_node_new
  public:: is_point_inside_node

  ! Inherited from BATL_grid
  public :: init_grid           ! initializa module
  public :: clean_grid          ! clean module
  public :: interpolate_grid_amr
  public :: interpolate_grid_amr_gc
  public :: show_grid_proc      ! only used in test?

  ! Inherited from BATL_geometry
  public:: init_geometry  ! initialize the module
  public:: clean_geometry ! clean up storage
  public:: xyz_to_coord, coord_to_xyz, radius_to_gen, gen_to_radius
  public:: set_high_geometry, rot_to_cart, cart_to_rot

  ! Inherited from BATL_test
  public:: read_test_param, find_test_cell
  public:: test_start, test_stop, test_cell

contains
  !============================================================================
  subroutine init_batl(CoordMinIn_D, CoordMaxIn_D, MaxBlockIn, &
       TypeGeometryIn, IsPeriodicIn_D, nRootIn_D, UseRadiusIn, UseDegreeIn, &
       rGenIn_I, UseUniformAxisIn, UseFDFaceFluxIn, iVectorVarIn_I)

    integer, intent(in):: MaxBlockIn         ! max number of blocks/processor
    real,    intent(in):: CoordMinIn_D(nDim) ! min (gen) coordinates of domain
    real,    intent(in):: CoordMaxIn_D(nDim) ! max (gen) coordinates of domain

    ! Number of blocks per dimension on the tree root level
    integer,          optional, intent(in):: nRootIn_D(nDim)

    ! Grid geometry type (cartesian, spherical, etc.)
    character(len=*), optional, intent(in):: TypeGeometryIn

    ! Periodicity of grid boundaries per dimension
    logical,          optional, intent(in):: IsPeriodicIn_D(nDim)

    ! Use true radius or generalized radial coordinate
    logical,          optional, intent(in):: UseRadiusIn

    ! Use radian or degrees for angle coordinate
    logical,          optional, intent(in):: UseDegreeIn

    ! Array describing the stretched radial coordinate
    real,             optional, intent(in):: rGenIn_I(:)

    ! Logical for uniform grid in the Phi direction around the axis
    logical,          optional, intent(in):: UseUniformAxisIn
    logical,          optional, intent(in):: UseFDFaceFluxIn

    ! Indexes of first components of vector variables in State_VGB
    integer,          optional, intent(in):: iVectorVarIn_I(:)

    ! Initialize the block-adaptive tree and the domain.
    !
    ! CoordMin_D and CoordMax_D are the domain boundaries. For non-Cartesian
    ! grids this is meant in generalized coordinates.
    !
    ! MaxBlockIn gives the maximum number of blocks per processor
    ! during the whole run.
    ! At the base level the domain is decomposed into the root blocks.
    !
    ! The optional nRootIn_D argument provides the number of root blocks
    ! in each dimension. The default is a single root block at the base level.
    !
    ! The optional TypeGeometry describes the geometry of the coordinate
    ! system. Default is "Cartesian". Currently the only other option is
    ! RZ geometry.
    !
    ! The optional IsPeriodicIn_D argument tells if a certain direction
    ! is periodic or not.
    !
    ! At the completion of this subroutine, the domain and the block-tree
    ! are initialized, the root blocks are distributed and their coordinates
    ! cell volumes, face areas, etc. are all set.
    !
    !--------------------------------------------------------------------------
    if(IsBatlInitialized) RETURN

    call init_tree(MaxBlockIn)
    call init_geometry(TypeGeometryIn, IsPeriodicIn_D, rGenIn_I, &
         UseFDFaceFluxIn)
    call init_grid(CoordMinIn_D, CoordMaxIn_D, UseRadiusIn, UseDegreeIn)
    if(present(UseUniformAxisIn))then
       ! IsAnyAxis is set by init_grid.
       if(IsAnyAxis .and. UseUniformAxisIn) &
            call set_tree_param(UseUniformAxisIn=.true.)
    end if
    call set_tree_root(nRootIn_D)
    call distribute_tree(DoMove=.true.)
    call create_grid
    call init_amr
    call init_amr_criteria

    if(present(iVectorVarIn_I))then
       allocate(iVectorVar_I(size(iVectorVarIn_I)))
       iVectorVar_I = iVectorVarIn_I
    end if

    IsBatlInitialized = .true.

  end subroutine init_batl
  !============================================================================
  subroutine clean_batl

    !--------------------------------------------------------------------------
    ! Free up memory
    call clean_amr_criteria
    call clean_grid
    call clean_geometry
    call clean_tree

    IsBatlInitialized = .false.

  end subroutine clean_batl
  !============================================================================

  subroutine init_grid_batl(DoRefine_B)

    logical, optional, intent(in):: DoRefine_B(MaxBlock)

    ! Initialize the grid. The grid does not contain any data at this point.
    ! The optional DoRefine_B argument allows initial refinement of the grid.
    ! After completion, the grid nodes are distributed over the processors,
    ! using Morton ordering, and the coordinates and related variables
    ! (cell volume, face areas, etc.) are all set.

    integer:: iBlock

    !--------------------------------------------------------------------------
    if(present(DoRefine_B))then
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(DoRefine_B(iBlock)) iStatusNew_A(iNode_B(iBlock)) = Refine_
       end do
    end if
    call adapt_tree
    call distribute_tree(DoMove=.true.)
    call create_grid
    call init_amr

  end subroutine init_grid_batl
  !============================================================================

  subroutine regrid_batl(nVar, State_VGB, Dt_B, DoRefine_B, DoCoarsen_B, &
       DoBalanceEachLevelIn, iTypeBalance_A, iTypeNode_A, &
       Used_GB, DoBalanceOnlyIn, DoTestIn, &
       nExtraData, pack_extra_data, unpack_extra_data, &
       UseHighOrderAMRIn, DefaultStateIn_V)

    integer, intent(in)   :: nVar                         ! number of variables
    real,    intent(inout):: &                            ! state variables
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)

    real, intent(inout), optional:: Dt_B(MaxBlock)        ! time step limit

    logical, intent(in), optional:: DoRefine_B(MaxBlock)  ! request to refine
    logical, intent(in), optional:: DoCoarsen_B(MaxBlock) ! request to coarsen
    logical, intent(in), optional:: DoBalanceEachLevelIn  ! balance per level?
    integer,intent(in),  optional:: iTypeBalance_A(MaxNode)! balance by types
    integer,intent(inout),optional::iTypeNode_A(MaxNode)  ! adapt node types
    logical, intent(in), optional:: &
         Used_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)  ! used cells
    logical, intent(in), optional:: DoBalanceOnlyIn       ! balance with ghosts
    logical, intent(in), optional:: DoTestIn              ! print test info
    integer, intent(in), optional:: nExtraData            ! size of extra data
    ! Optional methods to send extra information
    interface
       subroutine pack_extra_data(iBlock, nBuffer, Buffer_I)

         ! Pack extra data into Buffer_I

         integer, intent(in) :: iBlock            ! block index
         integer, intent(in) :: nBuffer           ! size of buffer
         real,    intent(out):: Buffer_I(nBuffer) ! buffer

       end subroutine pack_extra_data

       subroutine unpack_extra_data(iBlock, nBuffer, Buffer_I)

         ! Unpack extra data from Buffer_I

         integer, intent(in) :: iBlock            ! block index
         integer, intent(in) :: nBuffer           ! size of buffer
         real,    intent(in) :: Buffer_I(nBuffer) ! buffer

       end subroutine unpack_extra_data
    end interface
    optional:: pack_extra_data, unpack_extra_data

    ! If UseHighOrderAMRIn is true, 5th (6th) order accuracy will be achieved
    ! for refined (coarsened) blocks.
    logical, intent(in), optional:: UseHighOrderAMRIn

    ! Default values of each states. Used for high order AMR. If the default
    ! vaule is positive (like density, pressure), the value after AMR should
    ! be positive too.
    real,intent(in), optional:: DefaultStateIn_V(nVar)

    ! Refine, coarsen and load balance the blocks containing the nVar
    ! state variables in State_VGB. Use second order accurate conservative
    ! restriction and prolongation operators. Load balance by Morton ordering.
    !
    ! If iTypeBalance_A is present, it contains positive integers corresponding
    ! to different groups of blocks. Each group is load balanced separately.
    !
    ! If iTypeNode_A is present, it contains integers corresponding
    ! to different block types. The block types are inherited during AMR.
    !
    ! If DoBalanceEachLevelIn is true, balance each AMR level independently.
    !
    ! If DoBalnceOnlyIn is true, no AMR is performed, the blocks are simply
    ! moved among the processors together with ghost cells.
    !
    ! The Used_GB array can describe unused cells (e.g. inside an internal
    ! boundary) which cannot be used for prolongation or restriction.
    !
    ! The nExtraData, pack_extra_data and unpack_extra_data arguments
    ! allow sending nExtraData real numbers together with the blocks.
    !
    ! Refinement and coarsening is primarily based on the iStatusNew_A array
    ! (available via the BATL_lib module) that is indexed by nodes
    ! (so a processor can request refinement for a non-local block).
    ! It should be set to the values Refine_ and Coarsen_. The AMR algorithm
    ! checks whether the requests can be done while keeping resolution
    ! changes at most a factor of 2 (proper nesting). Refinement requests
    ! take precedence over coarsening requests when multiple processors
    ! set the iStatusNew_A for the same node.
    !
    ! The optional arguments DoRefine_B and DoCoarsen_B are provided
    ! for convenience. They can be used to request refinement and
    ! coarsening for local blocks only.
    !
    ! The optional Dt_B argument contains the time step limit per block.
    ! This information is moved together with the block during load balance,
    ! The time step limit is divided by 2 for prolonged blocks, and
    ! multiplied by 2 (and minimum is taken over the children) for coarsened
    ! blocks.
    !
    ! Note that this estimate of the time step limit for prolonged/restricted
    ! blocks is only an approximate solution, and for partial AMR
    ! it is not correct at all. One can always recalculate the time step
    ! limit for the coarsened and prolonged blocks.
    !
    ! The iAmrChange_B array provides information about the changes due to AMR
    ! in the local blocks. Its value is set to AmrRemoved_ for blocks that
    ! became unused, to AmrUnchanged_ for no change, to AmrMoved_ for blocks
    ! that moved from one processor to another, AmrRefined_ for blocks that
    ! were just refined, and AmrCoarsened_ for blocks that were just coarsened.

    logical:: DoBalanceEachLevel
    integer:: iBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'regrid_batl'
    !--------------------------------------------------------------------------
    DoTest = .false.
    if(present(DoTestIn)) DoTest = DoTestIn

    if(DoTest)write(*,*) NameSub,' starting with nVar=', nVar

    DoBalanceEachLevel = .false.
    if(present(DoBalanceEachLevelIn)) DoBalanceEachLevel = DoBalanceEachLevelIn

    if(present(DoCoarsen_B))then
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(DoCoarsen_B(iBlock)) iStatusNew_A(iNode_B(iBlock)) = Coarsen_
       end do
    end if

    if(present(DoRefine_B))then
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(DoRefine_B(iBlock)) iStatusNew_A(iNode_B(iBlock)) = Refine_
       end do
    end if

    ! Coarsen and refine the tree nodes
    if(DoTest)write(*,*) NameSub,' call adapt_tree'
    call adapt_tree(iTypeNode_A)

    ! Load balance the tree
    if(DoTest)write(*,*) NameSub, &
         ' call distribute_tree with DoBalanceEachLevel=', DoBalanceEachLevel

    if(DoBalanceEachLevel)then
       call distribute_tree(DoMove=.false., &
            iTypeBalance_A=iTree_IA(Level_,:)+1)
    else
       if(present(iTypeBalance_A).or.IsNewTree)&
            call distribute_tree(DoMove=.false., &
            iTypeBalance_A=iTypeBalance_A)
    end if

    ! Initialize iAmrChange
    iAmrChange_B(1:nBlock) = AmrUnchanged_

    ! No grid changes, no need for do_amr
    ! IsNewTree  == .true. also implies IsNewDecomposition == .true.
    if(.not.IsNewDecomposition) RETURN

    ! Coarsen, refine and load balance the flow variables, and set Dt_B.
    if(DoTest)write(*,*) NameSub,' call do_amr'
    call do_amr(nVar, State_VGB, Dt_B, &
         Used_GB=Used_GB, &
         DoBalanceOnlyIn=DoBalanceOnlyIn, &
         DoTestIn=DoTestIn, &
         nExtraData=nExtraData, &
         pack_extra_data=pack_extra_data, &
         unpack_extra_data=unpack_extra_data, &
         UseHighOrderAMRIn=UseHighOrderAMRIn, &
         DefaultStateIn_V=DefaultStateIn_V)

    ! This logical tells find_neighbor (called by move_tree) to check
    ! if the neighbor levels of a block (otherwise not affected by AMR) changed
    DoCheckResChange = nDim == 3 .and. IsNodeBasedGrid &
         .and. .not.IsCartesianGrid

    ! Finalize the tree information and compact the tree
    if(DoTest)write(*,*) NameSub,' call move_tree'
    call move_tree(iTypeNode_A)

    ! Fix 3D curvilinear grid at resolution changes so that faces match
    if(DoCheckResChange) call fix_grid_res_change
    DoCheckResChange = .false.

    if(DoTest)write(*,*) NameSub,' finished'

  end subroutine regrid_batl
  !============================================================================
end module BATL_lib
!==============================================================================
