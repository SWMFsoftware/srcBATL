!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_tree

  use BATL_size, ONLY: &
       MaxBlock, nBlock, MaxDim, nDim, iRatio_D, nDimAmr, iDimAmr_D, nIJK_D
  use BATL_geometry, ONLY: &
       IsPeriodic_D, IsPeriodicCoord_D, Phi_, Theta_, &
       IsCylindricalAxis, IsSphericalAxis, IsLatitudeAxis, IsAnyAxis

  use ModUtilities, ONLY: CON_stop

  implicit none
  save

  private ! except

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
  public:: i_node_new       ! only used in test?
  public:: is_point_inside_node ! only used in test?

  ! DoStrictAmr :  true if we want the program to stop because we can not
  ! refine/coarsen all the blocks we want
  logical, public :: DoStrictAmr = .true.

  ! nDesiredRefine and nDesiredCoarsen set the number of blocks that we
  ! wish to refine/coarsen but may not happen (e.g. DoStrictAmr = .false.)
  ! nNodeRefine and nNodeCoarsen set the number of blocks that has to be
  ! refined/coarsened or the program will stop.
  ! nNodeSort is the number of nodes in the ranking list iRank_A
  ! sorted by one or more AMR criteria

  integer, public :: nDesiredRefine,nNodeRefine, &
       nDesiredCoarsen, nNodeCoarsen, nNodeSort

  ! If the difference in the criteria is less than DiffRange for two blocks,
  ! then the blocks are refined/coarsened together (preserves symmetry).
  real, public, parameter :: DiffRange = 1.0e-6

  ! We can also specify the percentage of blocks we want to refine. For doing
  ! this we need to sort them into a priority list iRank_A. This priority list
  ! can also be used for refining/coarsening blocks to the point where we
  ! have no more blocks available.
  ! Rank_A stores the criteria values for iRank_A
  integer, public, allocatable :: iRank_A(:)

  ! Large Rank_A value means high priority for refinement, while
  ! low Rank_A value means high priority for coarsening.
  real, public, allocatable :: Rank_A(:)

  ! Maximun number of try to refine/coarsen the grid based on iRank_A
  integer, public :: iMaxTryAmr = 100

  ! Input parameter
  integer, public :: MaxTotalBlock = 0

  ! Number of children per node
  integer, public, parameter :: nChild = 2**nDimAmr

  ! Global tree information
  integer, public, allocatable :: iTree_IA(:,:)

  ! Named indexes of iTree_IA
  integer, public, parameter :: &
       Status_   =  1, &
       Level_    =  2, & ! grid level
       Proc_     =  3, & ! processor index
       Block_    =  4, & ! block index
       MinLevel_ =  5, & ! minimum level allowed
       MaxLevel_ =  6, & ! maximum level allowed
       Coord0_   =  6, & ! equal to Coord1_-1
       Coord1_   =  7, & ! coordinate of node in 1st dimension
       Coord2_   =  8, & ! coordinate of node in 2nd dimension
       Coord3_   =  9, & ! coordinate of node in 3rd dimension
       CoordLast_=  9, & ! Coord0_ + MaxDim (?)
       Parent_   = 10, & ! Parent_ must be
       Child0_   = 10, & ! equal to Child0_
       Child1_   = Child0_ + 1,      &
       ChildLast_= Child0_ + nChild

  ! Tell if the grid has changed (refined/coarsened) or blocks were moved.
  ! If IsNewGrid == .true. IsNewDecomposition should also be .true.
  logical, public :: IsNewDecomposition, IsNewTree

  ! Number of items stored in iTree_IA
  integer, parameter :: nInfo = ChildLast_

  character(len=10), parameter:: NameTreeInfo_I(Child0_+8) = [ &
       'Status   ', &
       'Level    ', &
       'Proc     ', &
       'Block    ', &
       'MinLevel ', &
       'MaxLevel ', &
       'Coord1   ', &
       'Coord2   ', &
       'Coord3   ', &
       'Parent   ', &
       'Child1   ', &
       'Child2   ', &
       'Child3   ', &
       'Child4   ', &
       'Child5   ', &
       'Child6   ', &
       'Child7   ', &
       'Child8   ' ]

  ! New status (refine, coarsen etc) requested for nodes
  integer, public, allocatable :: iStatusNew_A(:)
  integer, public, allocatable :: iStatusAll_A(:) ! needed for MPI_allreduce

  ! New processor index of a given node after next load balance
  integer, public, allocatable :: iProcNew_A(:)

  ! Mapping from local block index to global node index
  integer, public, allocatable :: iNode_B(:)

  logical, public, allocatable :: &
       Unused_B(:), Unused_BP(:,:) ! Unused blocks on local/all processors
       

  ! Target is useful to allow pointer from C wrapper
  integer, public, allocatable, target :: &
       DiLevelNei_IIIB(:,:,:,:),  &  ! Level difference relative to neighbors
       iNodeNei_IIIB(:,:,:,:)        ! Node index of neighboring blocks

  ! True for procs that contain blocks adjacent with blocks of this processor
  logical, public, allocatable :: IsNeighbor_P(:)

  ! Index for unset values (that are otherwise larger)
  integer, public, parameter :: Unset_ = -100

  ! Possible values for the status variable
  integer, public, parameter :: &
       Unused_      = -1, & ! unused block (not a leaf)
       Refine_      = -2, & ! parent block to be refined
       DontCoarsen_ = -3, & ! block not to be coarsened
       Coarsen_     = -4, & ! child block to be coarsened
       Used_        =  1, & ! currently used block (leaf)
       RefineNew_   =  2, & ! child block to be refined
       Refined_     =  3, & ! refined child block
       CoarsenNew_  =  4, & ! parent block to be coarsened
       Coarsened_   =  5    ! coarsened parent block

  ! Number of maximum, total and used nodes (leaves of the node tree)
  integer, public :: MaxNode, nNode = 0, nNodeUsed = 0

  ! Ordering along the Morton-Hilbert space filling curve
  integer, public, allocatable :: iNodeMorton_I(:), iMortonNode_A(:)

  ! Levels for the finest and coarsest nodes in the current grid
  integer, public:: nLevelMin = 0, nLevelMax = 0

  ! Deepest AMR level relative to root nodes (limited by 32 bit integers)
  integer, parameter, public :: MaxLevel = 30

  ! The maximum integer coordinate for a given level below root nodes
  ! Implied do loop was not understooed by the pgf90 compiler, so list them
  integer, parameter, public :: MaxCoord_I(0:MaxLevel) = &
       [ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, &
       16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, &
       4194304, 8388608, 16777216, 33554432, 67108864, 134217728, &
       268435456, 536870912, 1073741824 ]

  ! The number of root nodes in all dimensions, and altogether
  integer, public :: nRoot_D(MaxDim) = 0, nRoot = 0

  ! Status change due to AMR is registered in this array
  integer, public, allocatable:: iAmrChange_B(:)
  integer, public, parameter  :: AmrRemoved_ = -1, AmrUnchanged_ = 0, &
       AmrNeiChanged_ = 1, AmrMoved_ = 2, AmrRefined_ = 3, AmrCoarsened_ = 4

  ! Check for changes in resolution change
  logical, public:: DoCheckResChange = .false.

  ! Time level information for subcycling algorithm
  logical, public:: UseTimeLevel = .false.
  integer, public:: nTimeLevel = 0
  integer, public, allocatable:: iTimeLevel_A(:)

  ! Function generalizing the DiLevelNei_IIIB array for time levels
  public:: di_level_nei

  !$acc declare create(MaxNode, nNode, nNodeUsed, nRoot_D, nRoot)
  !$acc declare create(MaxLevel, nLevelMin, nLevelMax, MaxCoord_I)
  !$acc declare create(Unused_B, Unused_BP, Used_GB)
  !$acc declare create(iNode_B, iMortonNode_A, iNodeMorton_I)
  !$acc declare create(DiLevelNei_IIIB, iNodeNei_IIIB, IsNeighbor_P)
  !$acc declare create(iStatusNew_A, Refine_, Coarsen_, Unset_)
  !$acc declare create(iTree_IA, Status_, Level_, MinLevel_, MaxLevel_, Used_)
  !$acc declare create(Block_, Proc_, Coord0_, Coord1_, Coord2_, Coord3_)
  !$acc declare create(UseTimeLevel, nTimeLevel, iTimeLevel_A)
  !$acc declare create(IsNewDecomposition, IsNewTree)
  !$acc declare create(iAmrChange_B)
  !$acc declare create(AmrRemoved_, AmrUnchanged_, AmrMoved_, AmrRefined_)
  !$acc declare create(AmrCoarsened_)

  ! Local variables -----------------------------------------------

  integer :: iNodeNew = 0

  ! The index along the Morton curve is global so that it can be used by the
  ! recursive subroutine order_children
  integer :: iMorton

  ! Needed for compact_tree
  integer, allocatable:: iNodeNew_A(:)

  ! Neighbor in the +phi direction for nodes around the poles
  ! All neighbors can be found by going from node to node
  integer, allocatable:: iNodeAxisNei_A(:)

  ! Use uniform resolution around axis
  logical:: UseUniformAxis = .false.

  character(len=*), parameter:: NameMod = "BATL_tree"

contains
  !============================================================================
  subroutine init_tree(MaxBlockIn)

    use BATL_mpi, ONLY: nProc

    ! Initialize the tree assuming MaxBlockIn blocks per processor

    integer, intent(in) :: MaxBlockIn ! Max number of blocks per processor
    !--------------------------------------------------------------------------
    if(allocated(iTree_IA)) RETURN

    ! Store tree size and maximum number of blocks/processor
    MaxBlock = MaxBlockIn

    ! Initialize max number of blocks for all processors.
    if(MaxTotalBlock <= 0)then
       MaxTotalBlock = nProc*MaxBlock
    else
       ! MaxTotalBlock was set by the #AMRLIMIT command in
       ! BATL_amr_criteria::read_amr_criteria
       MaxTotalBlock = min(MaxTotalBlock, nProc*MaxBlock)
    end if

    ! During AMR we may need extra nodes, because coarsening does not
    ! free up nodes immediately (Coarsen_ status instead of Unset_).
    ! So we have a factor 2 in front, which is the worst case scenario.
    MaxNode = 2*ceiling(nProc*MaxBlock*(1 + 1.0/(nChild - 1)))

    IsNewDecomposition = .true.
    IsNewTree = .true.

    ! Allocate and initialize all elements of tree as unset
    allocate(iTree_IA(nInfo, MaxNode));                 iTree_IA       = Unset_
    allocate(iNodeMorton_I(MaxNode));                   iNodeMorton_I  = Unset_
    allocate(iMortonNode_A(MaxNode));                   iMortonNode_A  = Unset_
    allocate(iStatusNew_A(MaxNode));                    iStatusNew_A   = Unset_
    allocate(iStatusAll_A(MaxNode));                    iStatusAll_A   = Unset_
    allocate(iProcNew_A(MaxNode));                      iProcNew_A     = Unset_
    allocate(iNodeNew_A(MaxNode));                      iNodeNew_A     = Unset_
    allocate(iNode_B(MaxBlock));                        iNode_B        = Unset_
    allocate(Unused_B(MaxBlock));                       Unused_B       = .true.
    allocate(Unused_BP(MaxBlock,0:nProc-1));            Unused_BP      = .true.
    allocate(iNodeNei_IIIB(0:3,0:3,0:3,MaxBlock));      iNodeNei_IIIB  = Unset_
    allocate(DiLevelNei_IIIB(-1:1,-1:1,-1:1,MaxBlock)); DiLevelNei_IIIB= Unset_
    allocate(IsNeighbor_P(0:nProc-1));                  IsNeighbor_P   = .true.
    allocate(iAmrChange_B(MaxBlock));                   iAmrChange_B   = Unset_

    ! Initialize minimum and maximum levels of refinement
    iTree_IA(MinLevel_,:) = 0
    iTree_IA(MaxLevel_,:) = MaxLevel

    ! Initialize min and max refinement levels
    nLevelMin = 0
    nLevelMax = 0

  end subroutine init_tree
  !============================================================================
  subroutine set_tree_param(UseUniformAxisIn)

    logical, optional:: UseUniformAxisIn
    !--------------------------------------------------------------------------
    if(.not.present(UseUniformAxisIn)) RETURN
    UseUniformAxis = UseUniformAxisIn
    if(UseUniformAxis)then
       if(MaxNode == 0)call CON_stop('Set UseUniformAxis after init_tree!')
       if(.not.allocated(iNodeAxisNei_A))then
          allocate(iNodeAxisNei_A(MaxNode)); iNodeAxisNei_A = Unset_
       end if
    else
       if(allocated(iNodeAxisNei_A)) deallocate(iNodeAxisNei_A)
    end if

  end subroutine set_tree_param
  !============================================================================
  subroutine clean_tree
    !--------------------------------------------------------------------------
    if(.not.allocated(iTree_IA)) RETURN
    deallocate(iTree_IA, iNodeMorton_I, iMortonNode_A, &
         iStatusNew_A, iStatusAll_A, &
         iProcNew_A, iNodeNew_A, &
         iNode_B, Unused_B, Unused_BP, &
         iNodeNei_IIIB, DiLevelNei_IIIB, IsNeighbor_P, iAmrChange_B)

    if(allocated(iRank_A)) deallocate(iRank_A)

    call set_tree_param(UseUniformAxisIn=.false.)

    MaxNode = 0
    iNodeNew = 0

  end subroutine clean_tree
  !============================================================================
  integer function i_node_new()

    ! Find a skipped element in the iTree_IA array
    !--------------------------------------------------------------------------
    ! Try next node first
    if(iNodeNew < MaxNode)then
       iNodeNew = iNodeNew + 1
       if(iTree_IA(Status_, iNodeNew) == Unset_)then
          i_node_new = iNodeNew
          RETURN
       end if
    end if

    ! Search from beginning
    do iNodeNew = 1, MaxNode
       if(iTree_IA(Status_, iNodeNew) == Unset_)then
          i_node_new = iNodeNew
          RETURN
       end if
    end do

    ! Could not find any usable node
    call CON_stop('i_node_new: ran out of nodes')

  end function i_node_new
  !============================================================================
  subroutine set_tree_root(nRootIn_D)

    integer, optional, intent(in) :: nRootIn_D(nDim)

    integer :: iRoot, jRoot, kRoot, iNode, iRoot_D(MaxDim)

    character(len=*), parameter:: NameSub = 'set_tree_root'
    !--------------------------------------------------------------------------

    ! Set number of root blocks: default or input arguments
    nRoot_D = 1
    if(present(nRootIn_D)) nRoot_D(1:nDim) = nRootIn_D
    nRoot   = product(nRoot_D)

    ! Check for even number of root blocks in phi direction around the axis
    if(IsAnyAxis)then
       if(modulo(nRoot_D(Phi_),2) /= 0) call CON_stop(NameSub // &
            ': there must be an even number of root blocks around the axis')
    end if

    ! Use the first product(nRoot_D) nodes as root nodes in the tree
    iNode = 0
    do kRoot = 1, nRoot_D(3)
       do jRoot = 1, nRoot_D(2)
          do iRoot = 1, nRoot_D(1)

             iRoot_D = [ iRoot, jRoot, kRoot ]

             iNode = iNode + 1
             iTree_IA(Status_, iNode)            = Used_
             iTree_IA(Parent_, iNode)            = Unset_
             iTree_IA(Child1_:ChildLast_, iNode) = Unset_
             iTree_IA(Level_ , iNode)            = 0
             iTree_IA(Coord1_:CoordLast_, iNode) = iRoot_D

          end do
       end do
    end do

    nNodeUsed = nRoot
    nNode     = nRoot
    !$acc update device(nRoot_D)
  end subroutine set_tree_root
  !============================================================================
  subroutine refine_tree_node(iNode, iTypeNode_A)

    integer, intent(in) :: iNode
    integer, intent(inout), optional:: iTypeNode_A(MaxNode)

    integer :: iChild, DiChild, iLevelChild, iProc, iBlock
    integer :: iCoord_D(nDim)
    integer :: iDim, iDimAmr, iNodeChild

    character(len=*), parameter:: NameSub = 'refine_tree_node'
    !--------------------------------------------------------------------------

    if(iTree_IA(Status_, iNode) == Unused_) &
         call CON_stop(NameSub//' trying to refine and unused block')

    iTree_IA(Status_, iNode) = Refine_

    iLevelChild = iTree_IA(Level_, iNode) + 1
    iProc       = iTree_IA(Proc_,  iNode)
    iBlock      = iTree_IA(Block_, iNode)

    ! Check levels
    if(iLevelChild > MaxLevel) &
         call CON_stop('Error in refine_tree_node: too many levels')

    iCoord_D = 2*iTree_IA(Coord1_:Coord0_+nDim, iNode) - 1

    do iChild = Child1_, ChildLast_

       iNodeChild = i_node_new()

       ! Increase nNode if necessary
       nNode = max(nNode, iNodeChild)

       iTree_IA(iChild, iNode) = iNodeChild

       iTree_IA(Status_,   iNodeChild) = RefineNew_
       iTree_IA(Level_,    iNodeChild) = iLevelChild
       iTree_IA(MinLevel_, iNodeChild) = iTree_IA(MinLevel_,iNode)
       iTree_IA(MaxLevel_, iNodeChild) = iTree_IA(MaxLevel_,iNode)
       iTree_IA(Parent_,   iNodeChild) = iNode
       iTree_IA(Child1_:ChildLast_, iNodeChild) = Unset_

       ! Data will come from the parent's proc/block
       iTree_IA(Proc_,     iNodeChild) = iProc
       iTree_IA(Block_,    iNodeChild) = iBlock

       ! Calculate the coordinates of the child node
       DiChild = iChild - Child1_

       iDimAmr = 0
       do iDim = 1, MaxDim
          if(iRatio_D(iDim) == 2)then
             iDimAmr = iDimAmr + 1
             iTree_IA(Coord0_+iDim,iNodeChild) = &
                  iCoord_D(iDim) + ibits(DiChild, iDimAmr-1, 1)
          else
             ! The non-AMR coordinates remain the same as for the parent node
             iTree_IA(Coord0_+iDim,iNodeChild) = iTree_IA(Coord0_+iDim,iNode)
          endif
       end do

       ! Set node type to be the same as the parent's
       if(present(iTypeNode_A)) iTypeNode_A(iNodeChild) = iTypeNode_A(iNode)
    end do

    ! Keep track of used nodes in the future tree
    nNodeUsed = nNodeUsed + nChild - 1

  end subroutine refine_tree_node
  !============================================================================
  subroutine coarsen_tree_node(iNode, iTypeNode_A)

    integer, intent(in) :: iNode
    integer, intent(inout), optional:: iTypeNode_A(MaxNode)

    integer :: iChild, iNodeChild
    !--------------------------------------------------------------------------

    do iChild = Child1_, ChildLast_
       iNodeChild = iTree_IA(iChild, iNode)

       ! Set the status of the child node
       iTree_IA(Status_, iNodeChild) = Coarsen_
    end do

    ! Make this node used with no children
    iTree_IA(Status_, iNode) = CoarsenNew_

    ! Keep track of used nodes in the future tree
    nNodeUsed = nNodeUsed - nChild + 1

    ! Increase nNode if necessary
    nNode = max(nNode, iNode)

    ! Set node type to be the largest of all children
    if(present(iTypeNode_A)) iTypeNode_A(iNode) = &
         maxval(iTypeNode_A(iTree_IA(Child1_:ChildLast_,iNode)))

  end subroutine coarsen_tree_node
  !============================================================================
  subroutine adapt_tree(iTypeNode_A)

    use BATL_size, ONLY: iRatio, jRatio, kRatio
    use BATL_mpi, ONLY: iProc, nProc, iComm
    use ModMpi, ONLY: MPI_allreduce, MPI_INTEGER, MPI_MAX

    ! Optional node type:
    ! refined nodes inherit the type of the parent node
    ! coarsened node receive the largest integer type of the children
    integer, intent(inout), optional:: iTypeNode_A(MaxNode)

    ! All processors can request some status changes in iStatusNew_A.
    ! Here we collect requests, check for proper nesting,
    ! limitations on level, number of blocks, etc,
    ! modify iStatusNew_A and set iTree_IA.

    integer:: nNodeUsedNow, iMorton, iBlock, iChild, iStatus, iError
    integer:: iNode, iNodeParent, iNodeChild, iNodeChild_I(nChild)
    integer:: jNode, jNodeParent
    integer:: iLevel, iLevelNew, iLevelMax, iLevelMin
    integer:: jLevel, jLevelNew
    integer:: iSide, iSideMin, iSideMax
    integer:: jSide, jSideMin, jSideMax
    integer:: kSide, kSideMin, kSideMax
    integer:: iTryAmr, iRank, iRankLast
    real:: RankLimit

    ! Expected number of nodes minus availabla number of nodes
    integer:: DnNodeUsed = 0

    logical, parameter :: DoTest = .false., DoTestNei = .false.
    integer, allocatable :: iStatusNew0_A(:)

    character(len=*), parameter:: NameSub = 'adapt_tree'
    !--------------------------------------------------------------------------

    ! Collect the local status requests into a global request
    if(nProc > 1)then
       call MPI_allreduce(iStatusNew_A, iStatusAll_A, nNode, MPI_INTEGER, &
            MPI_MAX, iComm, iError)
       iStatusNew_A(1:nNode) = iStatusAll_A(1:nNode)
    end if

    ! if(iProc == 0) then
    !   do iNode=1,nNode
    !      if(iStatusNew_A(iNode) == Refine_) then
    !         print *," Want to Refine node =", iNode
    !      end if
    !   end do
    !   do iNode=1,nNode
    !      if(iStatusNew_A(iNode) == Coarsen_) then
    !         print *," Want to Coarsen node =", iNode
    !      end if
    !   end do
    ! end if

    ! store the initall list that will not be changed by the proper nesting
    if(.not.DoStrictAmr) then
       allocate(iStatusNew0_A(MaxNode))
       iStatusNew0_A(1:nNode) = iStatusNew_A(1:nNode)
    end if

    LOOPTRY: do iTryAmr = 1, iMaxTryAmr

       ! Check max and min levels and coarsening of all siblings
       iLevelMin = nLevelMax
       iLevelMax = 1
       do iMorton = 1, nNodeUsed
          iNode   = iNodeMorton_I(iMorton)

          iLevel    = iTree_IA(Level_,iNode)
          iLevelMin = min(iLevelMin, iLevel)

          ! Check MaxLevel_ of node to see if it can be refined
          if(iStatusNew_A(iNode) == Refine_)then
             if(iLevel >= iTree_IA(MaxLevel_,iNode))then
                iStatusNew_A(iNode) = Unset_
                CYCLE
             end if
             iLevelMax = max(iLevelMax, iLevel)

             if(UseUniformAxis)then
                ! All axis neighbors have to be refined
                if(iNodeAxisNei_A(iNode) /= Unset_) call refine_axis(iNode)
             end if
          end if

          ! Only nodes to be coarsened need further checking
          if(iStatusNew_A(iNode) /= Coarsen_) CYCLE

          ! Check MinLevel_ of node to see if it can be coarsened
          if(iLevel <= iTree_IA(MinLevel_,iNode))then
             iStatusNew_A(iNode) = Unset_
             CYCLE
          end if
          iLevelMax = max(iLevelMax, iLevel)

          ! Cancel coarsening if any axis neighbor is not to be coarsened
          if(UseUniformAxis)then
             if(iNodeAxisNei_A(iNode) /= Unset_) &
                  call check_axis_coarsening(iNode)
          end if

          ! Check if all siblings want to be coarsened
          iNodeParent = iTree_IA(Parent_,iNode)
          iNodeChild_I = iTree_IA(Child1_:ChildLast_,iNodeParent)
          if(.not.all(iStatusNew_A(iNodeChild_I) == Coarsen_))then
             ! Cancel coarsening requests for all siblings
             where(iStatusNew_A(iNodeChild_I) == Coarsen_) &
                  iStatusNew_A(iNodeChild_I) = Unset_

             if(UseUniformAxis)then
                ! Cancel coarsening around axis if any child is at the axis
                do iChild = 1, nChild
                   iNodeChild = iNodeChild_I(iChild)
                   if(iNodeAxisNei_A(iNodeChild) /= Unset_) &
                        call cancel_axis_coarsening(iNodeChild)
                end do
             end if

          end if

       end do

       ! Check proper nesting. Go down level by level. No need to
       ! check base level. Changes in the requests will be applied
       ! to all siblings immediately.
       LOOPLEVEL: do iLevel = iLevelMax, max(iLevelMin, 1), -1

          ! Parallel processing of nodes (blocks)
          LOOPBLOCK: do iBlock = 1, nBlock

             if(Unused_B(iBlock)) CYCLE LOOPBLOCK
             iNode = iNode_B(iBlock)
             if(iTree_IA(Level_,iNode) /= iLevel) CYCLE LOOPBLOCK

             ! Calculate requested level
             if(iStatusNew_A(iNode) == Refine_)then
                iLevelNew = iLevel + 1
             elseif(iStatusNew_A(iNode) == Coarsen_)then
                iLevelNew = iLevel - 1
             else
                CYCLE LOOPBLOCK
             end if

             ! Check neighbors around this corner of the parent block

             ! Loop from 0 to 1 or from 2 to 3 in the side index
             ! depending on which side this node
             ! is relative to its parent. If there is no refinement in some
             ! direction, then loop from 0 to 3 (and skip 2).

             if(iRatio==1)then
                iSideMin = 0; iSideMax = 3
             else
                iSideMin = 2*modulo(iTree_IA(Coord1_,iNode)-1, 2)
                iSideMax = iSideMin + 1
             end if
             if(nDim < 2)then
                ! 1D, no neighbors in j direction
                jSideMin = 1; jSideMax = 1
             elseif(jRatio == 1)then
                ! 2D or 3D but no AMR, check neighbors in both directions
                jSideMin = 0; jSideMax = 3
             else
                ! 2D or 3D and AMR: check only the directions
                ! corresponding to the corner occupied by this child
                jSideMin = 2*modulo(iTree_IA(Coord2_,iNode)-1, 2)
                jSideMax = jSideMin + 1
             end if
             if(nDim < 3)then
                ! 1 or 2D
                kSideMin = 1; kSideMax = 1
             elseif(kRatio == 1)then
                ! 3D but but no AMR, check neighbors in both directions
                kSideMin = 0; kSideMax = 3
             else
                ! 3D and AMR, check only the directions
                ! corresponding to the corner occupied by this child
                kSideMin = 2*modulo(iTree_IA(Coord3_,iNode)-1, 2)
                kSideMax = kSideMin + 1
             end if

             ! Loop through the at most seven neighbors
             do kSide = kSideMin, kSideMax
                if(kRatio == 1 .and. kSide == 2) CYCLE
                do jSide = jSideMin, jSideMax
                   if(jRatio == 1 .and. jSide == 2) CYCLE
                   do iSide = iSideMin, iSideMax
                      if(iRatio == 1 .and. iSide == 2) CYCLE

                      jNode = iNodeNei_IIIB(iSide,jSide,kSide,iBlock)

                      ! Don't check the node itself
                      if(iNode == jNode) CYCLE

                      ! Don't check if neighbor is outside the domain
                      if(jNode == Unset_) CYCLE

                      ! Get the current and requested levels for the neighbor
                      jLevel = iTree_IA(Level_,jNode)

                      jLevelNew = jLevel
                      if(iStatusNew_A(jNode) == Refine_) jLevelNew = jLevel + 1
                      if(iStatusNew_A(jNode) == Coarsen_)jLevelNew = jLevel - 1

                      ! Nothing to worry about
                      if(iLevelNew == jLevelNew) CYCLE

                      ! Fix levels if difference is too much
                      if(iLevelNew >= jLevelNew + 2)then

                         if(jLevelNew < jLevel)then
                            ! Neighbor and its siblings cannot be coarsened
                            jNodeParent = iTree_IA(Parent_,jNode)
                            do iChild = Child1_, ChildLast_
                               iNodeChild = iTree_IA(iChild,jNodeParent)
                               if(iStatusNew_A(iNodeChild) /= Coarsen_) CYCLE
                               iStatusNew_A(iNodeChild) = DontCoarsen_
                               if(.not.UseUniformAxis) CYCLE
                               if(iNodeAxisNei_A(iNodeChild) == Unset_) CYCLE
                               call cancel_axis_coarsening(iNodeChild)
                            end do
                         endif

                         ! If neighbor was coarser it has to be refined
                         if(jLevel < iLevel)then
                            iStatusNew_A(jNode) = Refine_
                            if(UseUniformAxis)then
                               if(iNodeAxisNei_A(jNode) /= Unset_) &
                                    call refine_axis(jNode)
                            end if
                         end if

                      elseif(iLevelNew <= jLevelNew - 2)then
                         ! Cannot coarsen this node
                         iNodeParent = iTree_IA(Parent_,iNode)
                         do iChild = Child1_, ChildLast_
                            iNodeChild = iTree_IA(iChild,iNodeParent)
                            iStatusNew_A(iNodeChild) = DontCoarsen_
                            if(.not.UseUniformAxis) CYCLE
                            if(iNodeAxisNei_A(iNodeChild) == Unset_) CYCLE
                            call cancel_axis_coarsening(iNodeChild)
                         end do
                         CYCLE LOOPBLOCK
                      end if

                   end do ! iSide
                end do ! jSide
             end do ! kSide
          end do LOOPBLOCK

          ! Collect the local status requests into a global request
          if(nProc > 1)then
             call MPI_allreduce(iStatusNew_A, iStatusAll_A, nNode, &
                  MPI_INTEGER, MPI_MAX, iComm, iError)
             iStatusNew_A(1:nNode) = iStatusAll_A(1:nNode)
          end if

       end do LOOPLEVEL! levels

       ! For strict AMR (or if there is no refinement) don't check anything
       if(DoStrictAmr .or. count(iStatusNew_A(1:nNode) == Refine_) == 0) &
            EXIT LOOPTRY

       ! Estimate the difference between the expected number of blocks
       ! after AMR and the number of available/allowed blocks.
       DnNodeUsed = (nNodeUsed - count(iStatusNew_A(1:nNode) == Coarsen_) &
            + count(iStatusNew_A(1:nNode) == Coarsen_)/nChild &
            + count(iStatusNew_A(1:nNode) == Refine_)*(nChild-1)) &
            - min(nProc*(MaxBlock-1), MaxTotalBlock)

       ! If we have geometry based AMR only and we ran out of blocks,
       ! then we abort refining/coarsening
       if(.not.allocated(iRank_A) .and. DnNodeUsed > 0) then
          if(iProc==0) write(*,*)'!!! WARNING in ',NameSub, &
               ': Need ', DnNodeUsed, &
               ' more blocks in total! Skipping refinement!'
          where(iStatusNew_A(1:nNode) == Refine_ )
            iStatusNew_A(1:nNode) = Unset_
          end where
          EXIT LOOPTRY
       end if

       ! exit if we have enough space for all the new blocks
       if(DnNodeUsed <= 0) then
          ! if(iProc==0)then
          !   write(*,*)'!!! DnNodeUsed, nNodeUsed=', DnNodeUsed, nNodeUsed
          !   write(*,*)'!!! nCoarsen=',count(iStatusNew_A(1:nNode)==Coarsen_)
          !   write(*,*)'!!! nRefine =',count(iStatusNew_A(1:nNode)==Refine_)
          !   write(*,*)'!!! nProc*MaxBlock, MaxTotalBlock=', &
          !        nProc*MaxBlock, MaxTotalBlock
          ! end if
          EXIT LOOPTRY
       end if

       ! Number of blocks we want to remove from the refinement list.
       ! Increase with the number of iterations
       DnNodeUsed = (iTryAmr-1)*(nChild-1) + DnNodeUsed

       ! Reset iStatusNew_A to its inital value
       iStatusNew_A = iStatusNew0_A

       ! nodes to be refined are indexed
       ! from nNodeSort-nDesiredRefine to nNodeSort in the sorted list.
       ! The lowest priority corresponds to the first element.
       ! Remove blocks starting from the lowest priority until the
       ! number of blocks in the new grid will not exceed the maximum.
       LOOPREMOVE: do iRank = nNodeSort - nDesiredRefine + 1, nNodeSort
          iNode = iRank_A(iRank)

          ! if a block has reached the max level, it has no impact
          iLevel = iTree_IA(Level_,iNode)
          if(iLevel >= iTree_IA(MaxLevel_,iNode)) CYCLE LOOPREMOVE

          ! Remove the block from the to-be-refined list and modify the
          ! estimate of the number of blocks above the available max
          if(iStatusNew_A(iNode) == Refine_ ) then
             iStatusNew_A(iNode) = Unset_
             DnNodeUsed = DnNodeUsed - (nChild-1)

             ! Check if we have reached the goal
             if(DnNodeUsed < 0) EXIT LOOPREMOVE
          end if
       end do LOOPREMOVE

       ! make sure that the refinment list starts at a point that has
       ! a jump in the criteria larger than DiffRange. This makes
       ! sure that blocks with identical criteria are all refined together.
       ! This helps preserving symmetry.

       ! iRank is the last indirect index that was removed from refinement list
       ! in LOOPREMOVE. Blocks with Rank below RankLimit should be removed too
       iRankLast = min(iRank, nNodeSort)
       RankLimit = Rank_A(iRankLast) + DiffRange
       do iRank = iRankLast + 1, nNodeSort
          iNode = iRank_A(iRank)
          if( Rank_A(iRank) >= RankLimit ) CYCLE LOOPTRY
          iStatusNew_A(iNode) = Unset_
       end do

    end do LOOPTRY

    ! Geometry criteria is not coverd by sorting but can also
    ! demand to many blocks
    if(iTryAmr > iMaxTryAmr) then
       where(iStatusNew_A(1:nNode) == Refine_ )
          iStatusNew_A(1:nNode) = Unset_
       end where
       write(*,*)"   BATL_tree::adapt_tree: No refinment done"
       if(.not. any(iStatusNew_A(1:nNode) == Coarsen_))RETURN
    end if

    if(.not.DoStrictAmr) deallocate(iStatusNew0_A)

    nNodeUsedNow = nNodeUsed

    IsNewTree          = .false.
    IsNewDecomposition = .false.

    ! Coarsen first to reduce number of nodes and used blocks
    do iMorton = 1, nNodeUsedNow
       iNode   = iNodeMorton_I(iMorton)
       iStatus = iStatusNew_A(iNode)

       if(iStatus /= Coarsen_) CYCLE
       IsNewTree          = .true.
       IsNewDecomposition = .true.

       iNodeParent = iTree_IA(Parent_,iNode)

       ! Coarsen the parent node based on the request stored in the first child
       if(iTree_IA(Child1_,iNodeParent) /= iNode) CYCLE

       call coarsen_tree_node(iNodeParent, iTypeNode_A)
    end do

    ! Refine next
    do iMorton = 1, nNodeUsedNow
       iNode   = iNodeMorton_I(iMorton)
       iStatus = iStatusNew_A(iNode)

       if(iStatus /= Refine_) CYCLE
       IsNewTree          = .true.
       IsNewDecomposition = .true.

       ! Refine tree node
       call refine_tree_node(iNode, iTypeNode_A)

       if(nNodeUsed > MaxBlock*nProc) EXIT

    end do

    iStatusNew_A(1:nNode) = Unset_

  contains
    !==========================================================================
    subroutine refine_axis(iNode)

      ! Refine all axis neighbors of iNode

      integer, intent(in):: iNode
      integer:: jNode
      !------------------------------------------------------------------------
      jNode = iNode
      do
         jNode = iNodeAxisNei_A(jNode)
         if(jNode == iNode) EXIT
         iStatusNew_A(jNode) = Refine_
      end do
    end subroutine refine_axis
    !==========================================================================
    subroutine check_axis_coarsening(iNode)

      ! Check if any of the axis neighbors are not to be coarsened.
      ! If not then cancel coarsening for all axis neighbors.

      integer, intent(in):: iNode
      integer:: jNode
      !------------------------------------------------------------------------
      jNode = iNode
      do
         jNode = iNodeAxisNei_A(jNode)
         if(jNode == iNode) EXIT
         if(iStatusNew_A(jNode) /= Coarsen_) EXIT
      end do
      ! No problem if all nodes around the axis are to be coarsened
      if(jNode == iNode) RETURN
      ! Cancel coarsening of the axis
      call cancel_axis_coarsening(iNode)

    end subroutine check_axis_coarsening
    !==========================================================================
    subroutine cancel_axis_coarsening(iNode)

      ! Cancel coarsening for all axis neighbors

      integer, intent(in):: iNode
      integer:: jNode, jNodeParent, jNodeChild_I(nChild)
      !------------------------------------------------------------------------
      jNode = iNode
      do
         jNode = iNodeAxisNei_A(jNode)
         jNodeParent = iTree_IA(Parent_,jNode)
         jNodeChild_I = iTree_IA(Child1_:ChildLast_,jNodeParent)
         iStatusNew_A(jNodeChild_I) = &
              max(iStatusNew_A(jNodeChild_I), DontCoarsen_)
         if(jNode == iNode) EXIT
      end do
    end subroutine cancel_axis_coarsening
    !==========================================================================
  end subroutine adapt_tree
  !============================================================================
  subroutine get_tree_position(iNode, PositionMin_D, PositionMax_D)
    !$acc routine seq

    integer, intent(in) :: iNode
    real,    intent(out):: PositionMin_D(MaxDim), PositionMax_D(MaxDim)

    ! Calculate normalized position of the edges of node inode.
    ! Zero is at the minimum boundary of the grid, one is at the max boundary

    integer :: iLevel
    integer :: MaxIndex_D(MaxDim)
    !--------------------------------------------------------------------------
    iLevel = iTree_IA(Level_, iNode)

    ! For non-AMR directions MaxIndex_D = nRoot_D
    ! For AMR     directions MaxIndex_D = nRoot_D*MaxCoord_I(iLevel)
    MaxIndex_D = ((MaxCoord_I(iLevel)-1)*(iRatio_D-1) + 1)*nRoot_D

    ! Convert to real by adding -1.0 or 0.0 for the two edges, respectively
    PositionMin_D = (iTree_IA(Coord1_:CoordLast_,iNode) - 1.0)/MaxIndex_D
    PositionMax_D = (iTree_IA(Coord1_:CoordLast_,iNode) + 0.0)/MaxIndex_D

  end subroutine get_tree_position
  !============================================================================
  subroutine find_tree_node(CoordIn_D, iNode)
    !$acc routine seq

    ! Find the node that contains a point. The point coordinates should
    ! be given in generalized coordinates normalized to the domain size:
    ! CoordIn_D = (CoordOrig_D - CoordMin_D)/DomainSize_D

    real, intent(in):: CoordIn_D(MaxDim)
    integer, intent(out):: iNode

    real :: Coord_D(MaxDim)
    integer :: iLevel, iChild
    integer :: iRoot_D(MaxDim), iCoord_D(nDimAmr), iBit_D(nDimAmr)

    ! Scale coordinates so that 1 <= Coord_D <= nRoot_D+1
    !--------------------------------------------------------------------------
    Coord_D = 1.0 + nRoot_D*max(0.0, min(1.0, CoordIn_D))

    ! Get root node index
    iRoot_D = min(int(Coord_D), nRoot_D)

    ! Root node indexes are ordered
    iNode = &
         iRoot_D(1) + nRoot_D(1)*((iRoot_D(2)-1) + nRoot_D(2)*(iRoot_D(3)-1))

    if(iTree_IA(Status_,iNode) == Used_) RETURN

    ! Get normalized coordinates within root node and scale it up
    ! to the largest resolution: 0 <= iCoord_D <= MaxCoord_I(nLevelMax)-1
    iCoord_D = min(MaxCoord_I(nLevelMax) - 1, &
         int((Coord_D(iDimAmr_D) - iRoot_D(iDimAmr_D))*MaxCoord_I(nLevelMax)))

    ! Go down the tree using bit information
    do iLevel = nLevelMax-1,0,-1
       ! Get the binary bits based on the coordinates
       iBit_D = ibits(iCoord_D, iLevel, 1)
       ! Construct child index as iChild = Sum Bit_i*2**i
       ! The powers of 2 are stored in MaxCoord_I
       iChild = sum(iBit_D*MaxCoord_I(0:nDimAmr-1)) + Child1_
       iNode  = iTree_IA(iChild,iNode)

       if(iTree_IA(Status_,iNode) == Used_) RETURN
    end do

    ! Did not find the point so set iNode as unset
    iNode = Unset_

  end subroutine find_tree_node
  !============================================================================
  subroutine find_tree_cell(Coord_D, iNode, iCell_D, CellDistance_D, &
       UseGhostCell)
    !$acc routine seq

    ! Find the node that contains a point. The point coordinates should
    ! be given in generalized coordinates normalized to the domain size:
    ! CoordIn_D = (CoordOrig_D - CoordMin_D)/DomainSize_D
    ! If UseGhostCell is not present or false
    !    then iCell_D returns the cell that contains the point.
    ! If UseGhostCell is present and true
    !    then iCell_D will contain the cell indexes to the left of the point.
    ! If CellDistance_D is present, return the signed distances per dimension
    ! normalized to the cell size. This can be used as interpolation weight.
    !

    real,           intent(in) :: Coord_D(MaxDim)
    integer,        intent(out):: iNode
    integer,        intent(out):: iCell_D(MaxDim)
    real,    optional, intent(out):: CellDistance_D(MaxDim)
    logical, optional, intent(in) :: UseGhostCell

    real:: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    real:: CellCoord_D(MaxDim)
    !--------------------------------------------------------------------------
    call find_tree_node(Coord_D, iNode)

    if(iNode == Unset_)then
       iCell_D = Unset_
       if(present(CellDistance_D)) CellDistance_D = Unset_
       RETURN
    end if

    call get_tree_position(iNode, PositionMin_D, PositionMax_D)
    CellCoord_D = 0.5 + &
         nIJK_D*(Coord_D - PositionMin_D)/(PositionMax_D - PositionMin_D)

    if(present(UseGhostCell))then
       if(UseGhostCell) then
          iCell_D = floor(CellCoord_D)
       else
          iCell_D = max(1, min(nIJK_D, nint(CellCoord_D)))
       end if
    else
       iCell_D = max(1, min(nIJK_D, nint(CellCoord_D)))
    end if
    if(present(CellDistance_D)) CellDistance_D = CellCoord_D - iCell_D

  end subroutine find_tree_cell
  !============================================================================
  subroutine interpolate_tree(Coord_D, iNodeCell_II, Weight_I)

    integer, parameter:: nPoint = 2**nDim

    real, intent(in)    :: Coord_D(MaxDim)
    integer, intent(out):: iNodeCell_II(0:nDim,nPoint)
    real,    intent(out):: Weight_I(nPoint)

    ! Find the nPoint=2**nDim cell centers that surround point Coord_D
    ! given in normalized coordinates (0<Coord_D<1).
    ! The cells are described by the node index and nDim cell indexes.
    ! Also provide the proper weights for a second order interpolation.

    integer:: iCell_D(MaxDim), jCell_D(MaxDim)
    real:: CellDistance_D(MaxDim), Weight_D(MaxDim)
    ! real:: CellSize_D(MaxDim), CoordShifted_D(MaxDim)
    integer:: iNode, i, j, k, iPoint, iDim
    !--------------------------------------------------------------------------
    call find_tree_cell(Coord_D, iNode, iCell_D, CellDistance_D)
    if(iNode == Unset_)then
       iNodeCell_II = Unset_
       Weight_I     = Unset_
       RETURN
    end if

    ! Initialize the cell indexes
    jCell_D = iCell_D
    Weight_D = 1.0

    ! In the non-ignored directions the point is between iCell_D and jCell_D
    ! Calculate interpolation weights for iCell
    do iDim = 1, nDim
       if(CellDistance_D(iDim) > 0.0 .or. &
            (CellDistance_D(iDim) == 0 .and. iCell_D(iDim) == 1) )then
          jCell_D(iDim)  = iCell_D(iDim)+1
          Weight_D(iDim) = 1.0 - CellDistance_D(iDim)
       else
          iCell_D(iDim)  = iCell_D(iDim) - 1
          Weight_D(iDim) = abs(CellDistance_D(iDim))
       end if
    end do

    iPoint = 0
    do k = iCell_D(3), jCell_D(3)
       do j = iCell_D(2), jCell_D(2)
          do i = iCell_D(1), jCell_D(1)
             iPoint = iPoint + 1
             iNodeCell_II(0,iPoint)                        = iNode
             iNodeCell_II(1,iPoint)                        = i
             if(nDim > 1) iNodeCell_II(min(2,nDim),iPoint) = j
             if(nDim > 2) iNodeCell_II(min(3,nDim),iPoint) = k
             Weight_I(iPoint) = product(Weight_D(1:nDim))

             ! Flip weight for the other cell
             Weight_D(1) = 1.0 - Weight_D(1)
          end  do
          Weight_D(2) = 1.0 - Weight_D(2)
       end  do
       Weight_D(3) = 1.0 - Weight_D(3)
    end  do

  end subroutine interpolate_tree
  !============================================================================
  logical function is_point_inside_node(Position_D, iNode)

    ! Check if position is inside node or not

    real,    intent(in):: Position_D(MaxDim)
    integer, intent(in):: iNode

    real    :: PositionMin_D(MaxDim), PositionMax_D(MaxDim)
    !--------------------------------------------------------------------------
    call get_tree_position(iNode, PositionMin_D, PositionMax_D)

    ! Include min edge but exclude max edge for sake of uniqueness
    is_point_inside_node = &
         all(Position_D >= PositionMin_D) .and. &
         all(Position_D <  PositionMax_D)

  end function is_point_inside_node
  !============================================================================
  subroutine find_neighbor_for_anynode(iNode, DiLevelNei_III)
    ! Find neighbours for any node in this processor or not.

    use BATL_size, ONLY: iRatio_D

    integer, intent(in):: iNode
    integer, intent(inout):: DiLevelNei_III(-1:1,-1:1,-1:1)

    integer :: iLevel, i, j, k, Di, Dj, Dk, jNode
    real :: Scale_D(MaxDim), x, y, z, y0, z0
    integer:: iNodeNei_III(0:3,0:3,0:3)

    logical:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'find_neighbor_for_anynode'
    !--------------------------------------------------------------------------
    if(DoTest)write(*,*)'Starting find neighbors for node ',iNode

    ! Get AMR level of the node
    iLevel = iTree_IA(Level_,iNode)

    ! Calculate scaling factor from integer index to 0<x,y,z<1 real coordinates
    Scale_D = 1.0/nRoot_D
    where(iRatio_D == 2) &
         Scale_D = Scale_D/MaxCoord_I(iLevel)

    if(DoTest)then
       write(*,*)'iNode, iLevel, Scale_D=', iNode, iLevel, Scale_D
       write(*,*)'scaled coordinates=', &
            iTree_IA(Coord1_:CoordLast_, iNode)*Scale_D
    end if

    ! Fill in self-referring info
    iNodeNei_III(1:2,1:2,1:2) = iNode
    DiLevelNei_III(0,0,0)     = 0

    ! Loop through neighbors
    do k=0,3
       Dk = nint((k - 1.5)/1.5)
       if(nDim < 3)then
          if(k/=1) CYCLE
          z = 0.3
       else
          z = (iTree_IA(Coord3_, iNode) + 0.4*k - 1.1)*Scale_D(3)
          if(z > 1.0 .or. z < 0.0)then
             if(IsPeriodic_D(3))then
                z = modulo(z, 1.0)
             elseif(.not.IsLatitudeAxis)then
                iNodeNei_III(:,:,k) = Unset_
                DiLevelNei_III(:,:,Dk) = Unset_
                CYCLE
             end if
          end if
       end if
       ! store z for spherical axis
       z0 = z
       do j=0,3
          z = z0
          Dj = nint((j - 1.5)/1.5)
          if(nDim < 2)then
             if(j/=1) CYCLE
             y = 0.3
          else
             y = (iTree_IA(Coord2_, iNode) + 0.4*j - 1.1)*Scale_D(2)
             if(y > 1.0 .or. y < 0.0)then
                if(IsPeriodic_D(2))then
                   y = modulo(y, 1.0)
                elseif(IsSphericalAxis)then
                   ! Push back theta and go around half way in phi
                   y = max(0.0, min(1.0, y))
                   z = modulo(z0+0.5, 1.0)
                else
                   iNodeNei_III(:,j,k) = Unset_
                   DiLevelNei_III(:,Dj,Dk) = Unset_
                   CYCLE
                end if
             end if
             if(z0 > 1.0 .or. z0 < 0.0)then
                ! Push back latitude and go around half way in longitude
                z = max(0.0, min(1.0, z0))
                y = modulo(y+0.5, 1.0)
             end if
          end if
          ! store y for cylindrical axis case
          y0 = y
          do i=0,3
             ! Exclude inner points
             if(0<i.and.i<3.and.0<j.and.j<3.and.0<k.and.k<3) CYCLE

             Di = nint((i - 1.5)/1.5)

             ! If neighbor is not finer, fill in the i=2 or j=2 or k=2 elements
             if(DiLevelNei_III(Di,Dj,Dk) >= 0)then
                if(i==2)then
                   iNodeNei_III(i,j,k) = iNodeNei_III(1,j,k)
                   CYCLE
                end if
                if(j==2)then
                   iNodeNei_III(i,j,k) = iNodeNei_III(i,1,k)
                   CYCLE
                end if
                if(k==2)then
                   iNodeNei_III(i,j,k) = iNodeNei_III(i,j,1)
                   CYCLE
                end if
             end if

             x = (iTree_IA(Coord1_, iNode) + 0.4*i - 1.1)*Scale_D(1)
             y = y0
             if(x > 1.0 .or. x < 0.0)then
                if(IsPeriodic_D(1))then
                   x = modulo(x, 1.0)
                elseif(IsCylindricalAxis .and. x < 0.0)then
                   ! Push back radius and go around half way in phi direction
                   x = 0.0
                   y = modulo(y0+0.5, 1.0)
                else
                   iNodeNei_III(i,j,k) = Unset_
                   DiLevelNei_III(Di,Dj,Dk) = Unset_
                   CYCLE
                end if
             end if

             call find_tree_node( [x, y, z], jNode)

             iNodeNei_III(i,j,k) = jNode
             DiLevelNei_III(Di,Dj,Dk) = &
                  iLevel - iTree_IA(Level_, jNode)

             if(DoTest) write(*,'(a,3i2,3f6.3,i4)') &
                  'i,j,k,x,y,z,jNode=',i,j,k,x,y,z,jNode

          end do
       end do
    end do

  end subroutine find_neighbor_for_anynode
  !============================================================================
  subroutine find_neighbor(iBlock)

    use BATL_size, ONLY: iRatio_D

    integer, intent(in):: iBlock

    integer :: iNode, iLevel, i, j, k, Di, Dj, Dk, jNode
    real :: Scale_D(MaxDim), x, y, z, y0, z0, x_D(3)

    integer:: DiLevelNeiOld_III(-1:1,-1:1,-1:1)

    logical, parameter :: DoTest = .false.
    !--------------------------------------------------------------------------
    iNode = iNode_B(iBlock)
    if(DoTest)write(*,*)'Starting find neighbors for node ',iNode

    ! Get AMR level of the node
    iLevel = iTree_IA(Level_,iNode)

    ! Calculate scaling factor from integer index to 0<x,y,z<1 real coordinates
    Scale_D = 1.0/nRoot_D
    where(iRatio_D == 2) &
         Scale_D = Scale_D/MaxCoord_I(iLevel)

    if(DoTest)then
       write(*,*)'iNode, iLevel, Scale_D=', iNode, iLevel, Scale_D
       write(*,*)'scaled coordinates=', &
            iTree_IA(Coord1_:CoordLast_, iNode)*Scale_D
    end if

    if(DoCheckResChange) DiLevelNeiOld_III = DiLevelNei_IIIB(:,:,:,iBlock)

    ! Fill in self-referring info
    iNodeNei_IIIB(1:2,1:2,1:2,iBlock) = iNode
    DiLevelNei_IIIB(0,0,0,iBlock)     = 0

    ! Loop through neighbors
    do k=0,3
       Dk = nint((k - 1.5)/1.5)
       if(nDim < 3)then
          if(k/=1) CYCLE
          z = 0.3
       else
          z = (iTree_IA(Coord3_, iNode) + 0.4*k - 1.1)*Scale_D(3)
          if(z > 1.0 .or. z < 0.0)then
             if(IsPeriodic_D(3))then
                z = modulo(z, 1.0)
             elseif(.not.IsLatitudeAxis)then
                iNodeNei_IIIB(:,:,k,iBlock) = Unset_
                DiLevelNei_IIIB(:,:,Dk,iBlock) = Unset_
                CYCLE
             end if
          end if
       end if
       ! store z for spherical axis
       z0 = z
       do j=0,3
          z = z0
          Dj = nint((j - 1.5)/1.5)
          if(nDim < 2)then
             if(j/=1) CYCLE
             y = 0.3
          else
             y = (iTree_IA(Coord2_, iNode) + 0.4*j - 1.1)*Scale_D(2)
             if(y > 1.0 .or. y < 0.0)then
                if(IsPeriodic_D(2))then
                   y = modulo(y, 1.0)
                elseif(IsSphericalAxis)then
                   ! Push back theta and go around half way in phi
                   y = max(0.0, min(1.0, y))
                   z = modulo(z0+0.5, 1.0)
                else
                   iNodeNei_IIIB(:,j,k,iBlock) = Unset_
                   DiLevelNei_IIIB(:,Dj,Dk,iBlock) = Unset_
                   CYCLE
                end if
             end if
             if(z0 > 1.0 .or. z0 < 0.0)then
                ! Push back latitude and go around half way in longitude
                z = max(0.0, min(1.0, z0))
                y = modulo(y+0.5, 1.0)
             end if
          end if
          ! store y for cylindrical axis case
          y0 = y
          do i=0,3
             ! Exclude inner points
             if(0<i.and.i<3.and.0<j.and.j<3.and.0<k.and.k<3) CYCLE

             Di = nint((i - 1.5)/1.5)

             ! If neighbor is not finer, fill in the i=2 or j=2 or k=2 elements
             if(DiLevelNei_IIIB(Di,Dj,Dk,iBlock) >= 0)then
                if(i==2)then
                   iNodeNei_IIIB(i,j,k,iBlock) = iNodeNei_IIIB(1,j,k,iBlock)
                   CYCLE
                end if
                if(j==2)then
                   iNodeNei_IIIB(i,j,k,iBlock) = iNodeNei_IIIB(i,1,k,iBlock)
                   CYCLE
                end if
                if(k==2)then
                   iNodeNei_IIIB(i,j,k,iBlock) = iNodeNei_IIIB(i,j,1,iBlock)
                   CYCLE
                end if
             end if

             x = (iTree_IA(Coord1_, iNode) + 0.4*i - 1.1)*Scale_D(1)
             y = y0
             if(x > 1.0 .or. x < 0.0)then
                if(IsPeriodic_D(1))then
                   x = modulo(x, 1.0)
                elseif(IsCylindricalAxis .and. x < 0.0)then
                   ! Push back radius and go around half way in phi direction
                   x = 0.0
                   y = modulo(y0+0.5, 1.0)
                else
                   iNodeNei_IIIB(i,j,k,iBlock) = Unset_
                   DiLevelNei_IIIB(Di,Dj,Dk,iBlock) = Unset_
                   CYCLE
                end if
             end if

             x_D = [x,y,z]

             call find_tree_node( x_D, jNode)

             iNodeNei_IIIB(i,j,k,iBlock) = jNode
             DiLevelNei_IIIB(Di,Dj,Dk,iBlock) = &
                  iLevel - iTree_IA(Level_, jNode)
             IsNeighbor_P(iTree_IA(Proc_, jNode)) = .true.
             if(DoTest) write(*,'(a,3i2,3f6.3,i4)') &
                  'i,j,k,x,y,z,jNode=',i,j,k,x,y,z,jNode

          end do
       end do
    end do

    ! Check here if there is a need to redo faces of this block
    if(DoCheckResChange)then
       ! Blocks that were created or moved are fine
       if(iAmrChange_B(iBlock) >= AmrMoved_) RETURN
       do k=-1,1; do j=-1,1; do i=-1,1
          ! Ignore corners
          if(abs(i) + abs(j) + abs(k) == 3) CYCLE

          ! Check if a coarser neighbor has been created or eliminated
          if(DiLevelNei_IIIB(i,j,k,iBlock) /= DiLevelNeiOld_III(i,j,k) &
              .and. (DiLevelNei_IIIB(i,j,k,iBlock) == 1 &
              .or.   DiLevelNeiOld_III(i,j,k) == 1)) then
             iAmrChange_B(iBlock) = AmrNeiChanged_
             RETURN
          end if
       end do; end do; end do
    end if

  end subroutine find_neighbor
  !============================================================================
  subroutine find_axis_neighbor

    ! Set iNodeAxisNei_A array for axis neighbor in +phi direction.
    ! There is no MPI communication. There are few blocks next to the axis.

    integer:: iNode, jNode, iLevel, iCoord, MaxCoord
    real:: PositionMin_D(MaxDim), PositionMax_D(MaxDim), Coord_D(MaxDim)
    !--------------------------------------------------------------------------
    if(.not.IsAnyAxis) RETURN

    do iNode = 1, nNode

       iNodeAxisNei_A(iNode) = Unset_
       if(iTree_IA(Status_,iNode) /= Used_) CYCLE

       ! Check if node is next to the pole
       if(IsCylindricalAxis)then
          if(iTree_IA(Coord1_,iNode) > 1) CYCLE
       else
          iLevel = iTree_IA(Level_,iNode)
          MaxCoord = nRoot_D(Theta_)
          if(iRatio_D(Theta_) == 2) MaxCoord = MaxCoord*MaxCoord_I(iLevel)
          iCoord = iTree_IA(Coord0_+Theta_,iNode)
          if(iCoord > 1 .and. iCoord < MaxCoord) CYCLE
       end if

       ! Get node position
       call get_tree_position(iNode, PositionMin_D, PositionMax_D)

       ! Calculate normalized coordinates for the node center
       Coord_D = 0.5*(PositionMax_D + PositionMin_D)

       ! Shift the Phi coordinate in the positive direction to the next node
       Coord_D(Phi_) = &
            modulo(Coord_D(Phi_)&
            +      0.6*(PositionMax_D(Phi_) - PositionMin_D(Phi_)), 1.0)

       call find_tree_node(Coord_D, jNode)
       iNodeAxisNei_A(iNode) = jNode
    end do

  end subroutine find_axis_neighbor
  !============================================================================
  subroutine compact_tree(iTypeNode_A)

    integer, intent(inout), optional:: iTypeNode_A(MaxNode)

    ! Eliminate holes from the tree.
    ! If iTypeNode_A is present, move the integer type with the node.

    ! Amount of shift for each node
    integer :: iNode, iNodeSkipped, iNodeOld, iNodeNew, i, iBlock
    character(len=*), parameter:: NameSub = 'compact_tree'
    !--------------------------------------------------------------------------
    ! Set impossible initial value
    iNodeSkipped = MaxNode + 1

    do iNode = 1, nNode

       if(iTree_IA(Status_, iNode) == Unset_)then
          ! Store the first skipped position
          iNodeSkipped = min(iNodeSkipped, iNode)
       elseif(iNodeSkipped < iNode)then
          ! Move node to the first skipped position
          iTree_IA(:,iNodeSkipped) = iTree_IA(:,iNode)
          if(present(iTypeNode_A)) &
               iTypeNode_A(iNodeSkipped) = iTypeNode_A(iNode)

          iTree_IA(Status_, iNode) = Unset_
          ! Store new node index
          iNodeNew_A(iNode) = iNodeSkipped
          ! Advance iNodeSkipped
          iNodeSkipped = iNodeSkipped + 1
       else
          ! The node did not move
          iNodeNew_A(iNode) = iNode
       endif
    end do

    ! Apply shifts
    do iNode = 1, MaxNode
       if(iTree_IA(Status_, iNode) == Unset_) EXIT
       do i = Parent_, ChildLast_
          iNodeOld = iTree_IA(i, iNode)
          if(iNodeOld /= Unset_) &
               iTree_IA(i, iNode) = iNodeNew_A(iNodeOld)
       end do
    end do

    ! Set number of nodes and starting point for new nodes (note EXIT above)
    nNode = iNode - 1
    iNodeNew = nNode

    ! Fix the node indexes along the Morton curve
    do iMorton = 1, nNodeUsed
       iNodeOld = iNodeMorton_I(iMorton)
       iNodeNew = iNodeNew_A(iNodeOld)
       iNodeMorton_I(iMorton) = iNodeNew
       iMortonNode_A(iNodeNew)= iMorton
    end do

    ! Fix iNode_B indexes
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       iNodeOld = iNode_B(iBlock)
       iNode_B(iBlock) = iNodeNew_A(iNodeOld)
    end do

    ! Reset iNodeNew_A
    iNodeNew_A(1:nNode) = Unset_

  end subroutine compact_tree
  !============================================================================
  subroutine write_tree_file(NameFile, IsFormattedIn)

    ! Write tree information into a file

    use ModUtilities, ONLY: UnitTmp_, open_file, close_file
    use BATL_mpi, ONLY: iProc, barrier_mpi

    character(len=*), intent(in):: NameFile
    logical, optional:: IsFormattedIn

    logical:: IsFormatted
    integer:: iNode, iInfo
    !--------------------------------------------------------------------------
    call compact_tree

    if(iProc == 0)then
       IsFormatted = .false.
       if(present(IsFormattedIn)) IsFormatted = IsFormattedIn
       if(IsFormatted)then
          call open_file(file=NameFile, status='replace')
          write(UnitTmp_,'(a)') 'BATL tree information after #START'
          write(UnitTmp_,'(a)') 'Line 1: nDim, nInfo, nNode'
          write(UnitTmp_,'(a)') 'Line 2: iRatio_D(nDim)  AMR ratio'
          write(UnitTmp_,'(a)') 'Line 3: nRoot_D(nDim)   Number of root blocks'
          write(UnitTmp_,'(a)') 'Rest  : iTree_IA(nInfo,iNode)   in nNode rows'
          do iInfo = 1, nInfo
             write(UnitTmp_,'(a,i2.2,a)') &
                  'Column', iInfo, ': '//NameTreeInfo_I(iInfo)
          end do
          write(UnitTmp_,'(a)') '#START'
          write(UnitTmp_,'(3i10)') nDim, nInfo, nNode
          write(UnitTmp_,'(3i10)') iRatio_D(1:nDim)
          write(UnitTmp_,'(3i10)') nRoot_D(1:nDim)
          do iNode = 1, nNode
             write(UnitTmp_,'(100i10)') iTree_IA(:,iNode)
          end do
       else
          call open_file(file=NameFile, status='replace', form='unformatted')
          write(UnitTmp_) nDim, nInfo, nNode
          write(UnitTmp_) iRatio_D(1:nDim)
          write(UnitTmp_) nRoot_D(1:nDim)
          write(UnitTmp_) iTree_IA(1:nInfo,1:nNode)
       end if
       call close_file
    end if
    call barrier_mpi

  end subroutine write_tree_file
  !============================================================================
  subroutine read_tree_file(NameFile)

    ! Read tree information from a file

    use ModUtilities, ONLY: UnitTmp_, open_file, close_file

    character(len=*),  intent(in):: NameFile

    logical:: IsFormatted
    integer:: iNode, iLine
    character(len=100):: StringLine
    integer:: nDimIn, nInfoIn, nNodeIn, iRatioIn_D(nDim), nRootIn_D(nDim)
    character(len=*), parameter:: NameSub = 'read_tree_file'
    !--------------------------------------------------------------------------
    ! Try reading the file as formatted first
    call open_file(file=NameFile, status='old')
    read(UnitTmp_,*) StringLine
    call close_file
    ! Check the line: if it reads BATL, it is a text file.
    IsFormatted = StringLine(1:4) == 'BATL'

    if(IsFormatted)then
       call open_file(file=NameFile, status='old')
       ! Skip header
       do iLine = 1, 6 + nInfo
          read(UnitTmp_, *) StringLine
          if(StringLine == '#START') EXIT
       end do
       read(UnitTmp_, *) nDimIn, nInfoIn, nNodeIn
    else
       call open_file(file=NameFile, status='old', form='unformatted')
       read(UnitTmp_) nDimIn, nInfoIn, nNodeIn
    end if
    if(nDimIn /= nDim)then
       write(*,*) NameSub,' nDimIn, nDim=',nDimIn, nDim
       call CON_stop(NameSub//' nDim is different in tree file!')
    end if
    if(nInfoIn /= nInfo)then
       write(*,*) NameSub,' nInfoIn, nInfo=',nInfoIn, nInfo
       call CON_stop(NameSub//' nInfo is different in tree file!')
    end if
    if(nNodeIn > MaxNode)then
       write(*,*) NameSub,' nNodeIn, MaxNode=',nNodeIn, MaxNode
       call CON_stop(NameSub//' too many nodes in tree file!')
    end if
    if(IsFormatted)then
       read(UnitTmp_, *) iRatioIn_D
       read(UnitTmp_, *) nRootIn_D
    else
       read(UnitTmp_) iRatioIn_D
       read(UnitTmp_) nRootIn_D
    end if
    if( any(iRatioIn_D /= iRatio_D(1:nDim)) )then
       write(*,*) NameSub, &
            ' iRatioIn_D=', iRatioIn_D,' iRatio_D=', iRatio_D(1:nDim)
       call CON_stop(NameSub//' iRatio_D is different in tree file!')
    end if

    call set_tree_root(nRootIn_D)

    if(IsFormatted)then
       do iNode = 1, nNodeIn
          read(UnitTmp_,*) iTree_IA(:,iNode)
       end do
    else
       read(UnitTmp_) iTree_IA(:,1:nNodeIn)
    end if
    call close_file

    ! Set nLevelMin and nLevelMax
    call order_tree

  end subroutine read_tree_file
  !============================================================================
  subroutine distribute_tree(DoMove, iTypeBalance_A, iTypeNode_A)

    ! Order tree with the space filling curve then
    ! - if DoMove=T, assign tree nodes to processors and blocks immediately
    ! - if DoMove=F, set iProcNew_A only with future processor index
    ! - if iTypeBalance_A is present, it contains block types 1, 2, .., nType
    !   each type is balanced separately. The total is also balanced.
    ! - if iTypeNode_A is present, it contains node types that should be
    !   moved together with the nodes (only used if DoMove=T).

    use BATL_mpi, ONLY: nProc

    ! Are nodes moved immediately or just assigned new processor/node
    logical, intent(in):: DoMove

    ! Optional block type. Each type is balanced separately
    integer, intent(in),    optional:: iTypeBalance_A(MaxNode)

    ! Optional node type that should be moved together with the node.
    integer, intent(inout), optional:: iTypeNode_A(MaxNode)

    integer :: iMorton, iNode, iBlockTo, iProcTo

    integer :: iType, nType, iProcStart, iProcStop, iProcExtraBlock
    integer, allocatable :: iNodeType_I(:), nNodeType_I(:), iProcType_I(:), &
         iBlock_P(:), nBlockType_PI(:,:)

    character(len=*), parameter:: NameSub = 'distribute_tree'
    !--------------------------------------------------------------------------
    ! DoMove is only true when we initialize the grid, so this is done only
    ! a few times.
    if(DoMove)Unused_BP = .true.

    ! Initialize processor and block indexes
    iProcNew_A(1:nNode) = Unset_

    ! Set iNodeMorton_I and iMortonNode_A
    call order_tree

    ! Check if there are multiple node types that need separate balancing
    if(present(iTypeBalance_A))then
       ! Find number of types and allocate arrays
       nType = maxval(iTypeBalance_A(1:nNode), &
            MASK=iTree_IA(Status_,1:nNode)>=Used_)
    else
       nType = 1
    end if

    ! Allocate load balance tables:
    ! number of blocks per type and per processor and type
    allocate(&
         iNodeType_I(nType), nNodeType_I(nType), &
         iProcType_I(nType), iBlock_P(0:nProc-1), &
         nBlockType_PI(0:nProc-1,nType))

    ! Initialize number of type, counter, processor index and block index
    ! for various node types
    nNodeType_I = 0
    iNodeType_I = 0
    iProcType_I = 0
    iBlock_P    = 0

    if(present(iTypeBalance_A))then
       ! Count number of nodes for each type.
       do iNode = 1, nNode
          if(iTree_IA(Status_,iNode)<=0) CYCLE
          iType = iTypeBalance_A(iNode)
          if(iType > 0) nNodeType_I(iType) = nNodeType_I(iType) + 1
       end do
    else
       nNodeType_I(1) = nNodeUsed
    end if

    ! Construct load balance table for various types
    do iType = 1, nType
       ! minimum number of blocks of type iType for each processor
       nBlockType_PI(:,iType) = nNodeType_I(iType)/nProc

       ! The processors with extra blocks are filled in from nProc-1 backwards
       iProcStart = nProc - modulo(sum(nNodeType_I(1:iType)),nProc)
       iProcStop  = iProcStart + modulo(nNodeType_I(iType),nProc) - 1
       do iProcExtraBlock = iProcStart, iProcStop
          iProcTo = modulo(iProcExtraBlock,nProc)
          nBlockType_PI(iProcTo,iType) = nBlockType_PI(iProcTo,iType) + 1
       end do

       ! convert nBlockType_PI to cummulative load table for easier use
       do iProcTo = 1, nProc-1
          nBlockType_PI(iProcTo,iType) = nBlockType_PI(iProcTo,iType) &
               + nBlockType_PI(iProcTo-1,iType)
       end do

    end do

    ! Distribute the nodes over the processors
    do iMorton = 1, nNodeUsed

       ! Get the node index and type
       iNode = iNodeMorton_I(iMorton)
       if(present(iTypeBalance_A))then
          iType = iTypeBalance_A(iNode)
       else
          iType = 1
       end if

       ! Increase the index for this node type
       iNodeType_I(iType) = iNodeType_I(iType) + 1

       ! Select target processor.
       ! Use iProcType_I to remember last proc. used for the given type
       do iProcTo = iProcType_I(iType), nProc-1
          if(iNodeType_I(iType) <= nBlockType_PI(iProcTo,iType))then
             iProcType_I(iType) = iProcTo
             EXIT
          end if
       end do

       ! Assign future processor index for node
       iProcNew_A(iNode) = iProcTo

       if(iProcTo /= iTree_IA(Proc_,iNode)) IsNewDecomposition = .true.

       if(DoMove)then
          ! Assign block index right away
          iBlock_P(iProcTo) = iBlock_P(iProcTo) + 1
          iBlockTo = iBlock_P(iProcTo)
          if(iBlockTo > MaxBlock) &
               call CON_stop(NameSub//' too many blocks per processor')
          iTree_IA(Block_,iNode) = iBlockTo
          Unused_BP(iBlockTo,iProcTo) = .false.
       end if

    end do

    deallocate(iNodeType_I, nNodeType_I, iProcType_I, iBlock_P, nBlockType_PI)

    if(DoMove) call move_tree(iTypeNode_A)

    !$acc update device(Unused_BP)
  end subroutine distribute_tree
  !============================================================================
  subroutine move_tree(iTypeNode_A)

    ! Finish the load balancing (with or without data movement)
    ! Set status for newly used and unused/unset nodes.
    ! Then compact the tree and find all the neighbors
    ! If iTypeNode_A is present, compact_tree moves the type with the node.

    use BATL_mpi, ONLY: iProc

    integer, intent(inout), optional:: iTypeNode_A(MaxNode)

    integer:: iMorton, iNode, iNodeChild, iNodeParent, iChild, iBlock

    character(len=*), parameter:: NameSub = 'move_tree'
    !--------------------------------------------------------------------------
    ! Update local Unused_B array
    Unused_B(:) = Unused_BP(:,iProc)

    ! Update nBlock too as we loop through the used blocks
    nBlock = 0
    do iMorton = 1, nNodeUsed
       iNode = iNodeMorton_I(iMorton)

       ! Move the node to new processor/node
       iTree_IA(Proc_,iNode) = iProcNew_A(iNode)

       if(       iTree_IA(Status_,iNode) == CoarsenNew_ &
            .or. iTree_IA(Status_,iNode) == Coarsened_) then

          ! Remove the children of newly coarsened blocks from the tree
          do iChild = Child1_, ChildLast_
             iNodeChild = iTree_IA(iChild, iNode)
             iTree_IA(:,iNodeChild) = Unset_
          end do
          iTree_IA(Child1_:ChildLast_, iNode) = Unset_

       elseif(   iTree_IA(Status_,iNode) == RefineNew_ &
            .or. iTree_IA(Status_,iNode) == Refined_)then

          ! Make the parent of newly refined blocks unused
          iNodeParent = iTree_IA(Parent_, iNode)
          iTree_IA(Proc_:Block_,iNodeParent) = Unset_
          iTree_IA(Status_,iNodeParent)      = Unused_
       end if

       ! Now newly formed blocks are simply used
       iTree_IA(Status_,iNode) = Used_

       ! Set local information for this processor
       if(iProc == iTree_IA(Proc_,iNode))then
          iBlock = iTree_IA(Block_,iNode)
          iNode_B(iBlock)  = iNode
          if(.not.Unused_B(iBlock)) nBlock = max(nBlock, iBlock)
       end if

       !$acc update device(nBlock)

    end do
    ! Now that we removed children of coarsened blocks, compact the tree
    call compact_tree(iTypeNode_A)

    ! Inintialize list of neighboring processors:
    IsNeighbor_P = .false.
    ! Set neighbor info
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call find_neighbor(iBlock)
    end do

    if(UseUniformAxis) call find_axis_neighbor
    ! iProc is not a neighbor to itself
    IsNeighbor_P(iProc) = .false.

    !$acc update device(iTree_IA, iNode_B)
    !$acc update device(iNodeNei_IIIB, DiLevelNei_IIIB)
  end subroutine move_tree
  !============================================================================
  subroutine order_tree

    ! Set iNodeMorton_I and iMortonNode_A indirect index arrays according to
    ! 1. root node order
    ! 2. Morton ordering for each root node

    integer :: iNode, iRoot, jRoot, kRoot, iLevel
    !--------------------------------------------------------------------------
    nNode = nRoot
    iNode = 0
    iMorton = 0
    iNodeMorton_I(1:nNodeUsed) = Unset_
    iMortonNode_A(1:nNodeUsed) = Unset_
    do kRoot = 1, nRoot_D(3)
       do jRoot = 1, nRoot_D(2)
          do iRoot = 1, nRoot_D(1)
             ! Root nodes are the first ones
             iNode = iNode + 1

             ! All root nodes are handled as if they were first child
             call order_children(iNode)
          end do
       end do
    end do

    nNodeUsed = iMorton

    ! Set min and max refinement levels
    nLevelMin = MaxLevel
    nLevelMax = 0
    do iMorton = 1, nNodeUsed
       iNode = iNodeMorton_I(iMorton)
       iLevel = iTree_IA(Level_,iNode)
       nLevelMin = min(iLevel, nLevelMin)
       nLevelMax = max(iLevel, nLevelMax)
    end do

   !$acc update device(nLevelMin, nLevelMax)
  end subroutine order_tree
  !============================================================================
  recursive subroutine order_children(iNode)

    ! Recursively apply Morton ordering for nodes below a root block.
    ! Store result into iNodeMorton_I and iMortonNode_A using the global
    ! iMorton index.

    integer, intent(in) :: iNode
    integer :: iChild
    !--------------------------------------------------------------------------
    nNode = max(nNode, iNode)

    if(iTree_IA(Status_, iNode) >= Used_)then
       iMorton = iMorton + 1
       iNodeMorton_I(iMorton) = iNode
       iMortonNode_A(iNode)   = iMorton
    else
       do iChild = Child1_, ChildLast_
          call order_children(iTree_IA(iChild, iNode))
       end do
    end if

  end subroutine order_children
  !============================================================================

  integer function min_tree_level(iStage)

    integer, intent(in):: iStage

    ! Usage: in the iStage stage of the subcycling scheme the
    ! grid blocks at or above min_tree_level(iStage) are advanced.

    ! Theory: if iStage-1 contains 2^n in its prime factorization
    ! then grid blocks with grid/time levels equal or above Maxlevel-n
    ! are advanced.

    integer:: i, n
    !--------------------------------------------------------------------------
    if(iStage == 1)then
       ! All blocks are advanced in the first stage
       if(UseTimeLevel)then
          min_tree_level = 0
       else
          min_tree_level = nLevelMin
       end if
       RETURN
    end if

    i = iStage - 1
    n = 0
    do
       if(mod(i, 2) == 0)then
          i = i/2
          n = n + 1
       else
          EXIT
       end if
    end do
    if(UseTimeLevel)then
       min_tree_level = max(0, nTimeLevel - n)
    else
       min_tree_level = max(nLevelMin, nLevelMax - n)
    end if

  end function min_tree_level
  !============================================================================
  function di_level_nei(iDir, jDir, kDir, iBlock, DoResChangeOnly) &
       RESULT(DiLevel)

    integer, intent(in):: iDir, jDir, kDir, iBlock
    logical, intent(in), optional:: DoResChangeOnly

    integer:: DiLevel

    ! Generalize the DiLevelNei_IIIB array for time levels
    ! If DoReschangeOnly is present and false and the grid levels
    ! are the same (DiLevelNei_IIIB == 0) then set
    ! DiLevel according to the time level difference.

    integer:: iTimeLevel, iNodeNei, iTimeLevelNei

    character(len=*), parameter:: NameSub = 'di_level_nei'
    !--------------------------------------------------------------------------
    DiLevel = DiLevelNei_IIIB(iDir,jDir,kDir,iBlock)

    if(DiLevel /= 0) RETURN
    if(.not.present(DoResChangeOnly)) RETURN
    if(DoResChangeOnly) RETURN

    if(.not.allocated(iTimeLevel_A)) call CON_stop(NameSub// &
         ' called with DoResChangeOnly=F while iTimeLevel_A is not allocated')

    iTimeLevel    = iTimeLevel_A(iNode_B(iBlock))

    ! Get the neighbor node index and time level
    iNodeNei = iNodeNei_IIIB((3*iDir+3)/2,(3*jDir+3)/2,(3*kDir+3)/2,iBlock)
    iTimeLevelNei = iTimeLevel_A(iNodeNei)

    if(iTimeLevel == iTimeLevelNei) RETURN

    ! Set DiLevel = 1 if time level of the block is larger than its neighbor's
    ! otherwise set -1.
    DiLevel = sign(1, iTimeLevel - iTimeLevelNei)

  end function di_level_nei
  !============================================================================
  subroutine set_tree_periodic(IsOn)

    logical, intent(in):: IsOn

    ! Switch on or off the periodic connectivity as needed

    integer:: iSign, iBlock, iNode, iLevel, iCoord, MaxIndex_D(MaxDim)

    ! The only reason to do this if there is at least one
    ! periodic direction that is NOT a true periodic coordinate.
    !--------------------------------------------------------------------------
    if(.not.any(IsPeriodic_D(1:nDim) .and. .not. IsPeriodicCoord_D(1:nDim)))&
         RETURN

    if(IsOn)then
       ! Subtract Unset_ value from DiLevelNei
       iSign = -1
    else
       ! Add Unset_ value to DiLevelNei
       iSign = +1
    end if

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       iNode  = iNode_B(iBlock)
       iLevel = iTree_IA(Level_, iNode)

       ! Largest possible value for the node coordinate
       ! See get_tree_position for an explanation
       MaxIndex_D = ((MaxCoord_I(iLevel)-1)*(iRatio_D-1) + 1)*nRoot_D

       if(IsPeriodic_D(1) .and. .not.IsPeriodicCoord_D(1))then
          iCoord = iTree_IA(Coord1_, iNode)
          if(iCoord == 1) &
               DiLevelNei_IIIB(-1,:,:,iBlock) = &
               DiLevelNei_IIIB(-1,:,:,iBlock) + iSign*Unset_
          if(iCoord == MaxIndex_D(1)) &
               DiLevelNei_IIIB(+1,:,:,iBlock) = &
               DiLevelNei_IIIB(+1,:,:,iBlock) + iSign*Unset_
       end if

       if(nDim == 1) CYCLE

       if(IsPeriodic_D(2) .and. .not.IsPeriodicCoord_D(2))then
          iCoord = iTree_IA(Coord2_, iNode)
          if(iCoord == 1) &
               DiLevelNei_IIIB(:,-1,:,iBlock) = &
               DiLevelNei_IIIB(:,-1,:,iBlock) + iSign*Unset_
          if(iCoord == MaxIndex_D(2)) &
               DiLevelNei_IIIB(:,+1,:,iBlock) = &
               DiLevelNei_IIIB(:,+1,:,iBlock) + iSign*Unset_
       end if

       if(nDim == 2) CYCLE

       if(IsPeriodic_D(3) .and. .not.IsPeriodicCoord_D(3))then
          iCoord = iTree_IA(Coord3_, iNode)
          if(iCoord == 1) &
               DiLevelNei_IIIB(:,:,-1,iBlock) = &
               DiLevelNei_IIIB(:,:,-1,iBlock) + iSign*Unset_
          if(iCoord == MaxIndex_D(3)) &
               DiLevelNei_IIIB(:,:,+1,iBlock) = &
               DiLevelNei_IIIB(:,:,+1,iBlock) + iSign*Unset_
       end if
    end do

    !$acc update device(DiLevelNei_IIIB)
  end subroutine set_tree_periodic
  !============================================================================
  subroutine show_tree(String, DoShowNei, DoShowTimeLevel)

    character(len=*), intent(in):: String
    logical, optional,intent(in):: DoShowNei, DoShowTimeLevel

    ! Show complete tree information. Also write out string as an identifier.

    character(len=10) :: Name
    character(len=200):: StringHeader
    logical:: DoShowTime
    integer:: iInfo, iNode, iBlock
    !--------------------------------------------------------------------------
    DoShowTime = .false.
    if(present(DoShowTimeLevel)) DoShowTime = DoShowTimeLevel

    StringHeader = 'iNode'
    do iInfo = 1, nInfo
       Name = NameTreeInfo_I(iInfo)
       StringHeader(7*iInfo+1:7*(iInfo+1)-1) = Name(1:6)
    end do

    if(DoShowTime) StringHeader = trim(StringHeader)//' iTimeLevel'

    write(*,*) String
    write(*,*) trim(StringHeader)
    do iNode = 1, MaxNode
       if(iTree_IA(Status_, iNode) == Unset_) CYCLE
       if(DoShowTime)then
          write(*,'(100i7)') iNode, iTree_IA(:, iNode), iTimeLevel_A(iNode)
       else
          write(*,'(100i7)') iNode, iTree_IA(:, iNode)
       end if
    end do

    if(.not.present(DoShowNei)) RETURN
    if(.not.DoShowNei) RETURN

    write(*,*)'nNode, nNodeUsed, nBlock=',nNode, nNodeUsed, nBlock
    write(*,*)'iNodeMorton_I =', iNodeMorton_I(1:nNodeUsed)
    write(*,*)'iMortonNode_A =', iMortonNode_A(1:nNode)
    write(*,*)'IsPeriodic_D =', IsPeriodic_D

    iNode = iNodeMorton_I(1)
    iBlock = iTree_IA(Block_,iNode)
    write(*,*)'DiLevelNei_IIIB(:,0,0,First)=', DiLevelNei_IIIB(:,0,0,iBlock)
    write(*,*)'DiLevelNei_IIIB(0,:,0,First)=', DiLevelNei_IIIB(0,:,0,iBlock)
    write(*,*)'DiLevelNei_IIIB(0,0,:,First)=', DiLevelNei_IIIB(0,0,:,iBlock)
    write(*,*)'iNodeNei_IIIB(:,1,1,  First)=',   iNodeNei_IIIB(:,1,1,iBlock)
    write(*,*)'iNodeNei_IIIB(1,:,1,  First)=',   iNodeNei_IIIB(1,:,1,iBlock)
    write(*,*)'iNodeNei_IIIB(1,1,:,  First)=',   iNodeNei_IIIB(1,1,:,iBlock)

    if(UseUniformAxis)write(*,*)'iNodeAxisNei_A=',iNodeAxisNei_A(1:nNode)

  end subroutine show_tree
  !============================================================================
end module BATL_tree
!==============================================================================
