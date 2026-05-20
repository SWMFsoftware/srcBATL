!  Copyright (C) 2024 Regents of the University of Michigan
!  C-compatible wrapper for BATL for Julia integration

module BATL_c
  use iso_c_binding
  use BATL_lib
  ! Additional modules for internal functions not in BATL_lib
  use BATL_tree, ONLY: init_tree, set_tree_root, refine_tree_node, &
                       coarsen_tree_node, distribute_tree, &
                       find_tree_node, find_tree_cell, interpolate_tree, &
                       i_node_new, is_point_inside_node, show_tree
  use BATL_geometry, ONLY: init_geometry, xyz_to_coord, coord_to_xyz

  implicit none

contains
  !============================================================================

  ! --- BATL_size parameters ---

  integer(c_int) function c_batl_get_ni() bind(c)
    c_batl_get_ni = int(nI, c_int)
  end function c_batl_get_ni
!==============================================================================

  integer(c_int) function c_batl_get_nj() bind(c)
    c_batl_get_nj = int(nJ, c_int)
  end function c_batl_get_nj
!==============================================================================

  integer(c_int) function c_batl_get_nk() bind(c)
    c_batl_get_nk = int(nK, c_int)
  end function c_batl_get_nk
!==============================================================================

  integer(c_int) function c_batl_get_ng() bind(c)
    c_batl_get_ng = int(nG, c_int)
  end function c_batl_get_ng
!==============================================================================

  integer(c_int) function c_batl_get_ndim() bind(c)
    c_batl_get_ndim = int(nDim, c_int)
  end function c_batl_get_ndim
!==============================================================================

  integer(c_int) function c_batl_get_maxblock() bind(c)
    c_batl_get_maxblock = int(MaxBlock, c_int)
  end function c_batl_get_maxblock
!==============================================================================

  ! --- BATL_geometry variables & accessors ---

  integer(c_int) function c_batl_is_cartesian() bind(c)
    if (IsCartesian) then
       c_batl_is_cartesian = 1
    else
       c_batl_is_cartesian = 0
    end if
  end function c_batl_is_cartesian
!==============================================================================

  integer(c_int) function c_batl_is_rz_geometry() bind(c)
    if (IsRzGeometry) then
       c_batl_is_rz_geometry = 1
    else
       c_batl_is_rz_geometry = 0
    end if
  end function c_batl_is_rz_geometry
!==============================================================================

  ! --- BATL_tree variables ---

  integer(c_int) function c_batl_get_nroot() bind(c)
    c_batl_get_nroot = int(nRoot, c_int)
  end function c_batl_get_nroot
!==============================================================================

  integer(c_int) function c_batl_get_nnode() bind(c)
    c_batl_get_nnode = int(nNode, c_int)
  end function c_batl_get_nnode
!==============================================================================

  ! --- Lifecycle & Core API ---

  subroutine c_init_mpi() bind(c)
  !----------------------------------------------------------------------------------------------
    call init_mpi()
  end subroutine c_init_mpi
!==============================================================================

  subroutine c_init_batl(CoordMinIn_D, CoordMaxIn_D, MaxBlockIn) bind(c)
    real(c_double), intent(in) :: CoordMinIn_D(nDim)
    real(c_double), intent(in) :: CoordMaxIn_D(nDim)
    integer(c_int), intent(in), value :: MaxBlockIn
  !----------------------------------------------------------------------------------------------
    call init_batl(real(CoordMinIn_D), real(CoordMaxIn_D), int(MaxBlockIn))
  end subroutine c_init_batl
!==============================================================================

  subroutine c_clean_batl() bind(c)
  !----------------------------------------------------------------------------------------------
    call clean_batl()
  end subroutine c_clean_batl
!==============================================================================

  subroutine c_init_grid_batl() bind(c)
  !----------------------------------------------------------------------------------------------
    call init_grid_batl()
  end subroutine c_init_grid_batl
!==============================================================================

  subroutine c_regrid_batl(nVar, State_VGB) bind(c)
    integer(c_int), intent(in), value :: nVar
    real(c_double), intent(inout) :: State_VGB(nVar, MinI:MaxI, MinJ:MaxJ, MinK:MaxK, MaxBlock)
  !----------------------------------------------------------------------------------------------
    call regrid_batl(int(nVar), State_VGB)
  end subroutine c_regrid_batl
!==============================================================================

  ! --- BATL_tree subroutines ---

  subroutine c_init_tree(MaxBlockIn) bind(c)
    integer(c_int), intent(in), value :: MaxBlockIn
  !----------------------------------------------------------------------------------------------
    call init_tree(int(MaxBlockIn))
  end subroutine c_init_tree
!==============================================================================

  subroutine c_set_tree_root(nRootIn_D) bind(c)
    integer(c_int), intent(in) :: nRootIn_D(nDim)
    integer :: nRoot_F(MaxDim)
  !----------------------------------------------------------------------------------------------
    nRoot_F = 1
    nRoot_F(1:nDim) = int(nRootIn_D)
    call set_tree_root(nRoot_F)
  end subroutine c_set_tree_root
!==============================================================================

  integer(c_int) function c_i_node_new() bind(c)
    c_i_node_new = int(i_node_new(), c_int)
  end function c_i_node_new
!==============================================================================

  subroutine c_refine_tree_node(iNode) bind(c)
    integer(c_int), intent(in), value :: iNode
  !------------------------------------------------------------------------------------------------
    call refine_tree_node(int(iNode))
  end subroutine c_refine_tree_node
!==============================================================================

  subroutine c_coarsen_tree_node(iNode) bind(c)
    integer(c_int), intent(in), value :: iNode
  !------------------------------------------------------------------------------------------------
    call coarsen_tree_node(int(iNode))
  end subroutine c_coarsen_tree_node
!==============================================================================

  subroutine c_distribute_tree(DoBalanceOnlyIn) bind(c)
    integer(c_int), intent(in), value :: DoBalanceOnlyIn
  !------------------------------------------------------------------------------------------------
    call distribute_tree(DoBalanceOnlyIn /= 0)
  end subroutine c_distribute_tree
!==============================================================================

  subroutine c_find_tree_node(Coord_D, iNode) bind(c)
    real(c_double), intent(in) :: Coord_D(nDim)
    integer(c_int), intent(out) :: iNode
    integer :: iNode_F
    real :: Coord_F(MaxDim)
  !------------------------------------------------------------------------------------------------
    Coord_F = 0.0
    Coord_F(1:nDim) = real(Coord_D)
    call find_tree_node(Coord_F, iNode_F)
    iNode = int(iNode_F, c_int)
  end subroutine c_find_tree_node
!==============================================================================

  integer(c_int) function c_is_point_inside_node(Coord_D, iNode) bind(c)
    real(c_double), intent(in) :: Coord_D(nDim)
    integer(c_int), intent(in), value :: iNode
    real :: Coord_F(MaxDim)
    Coord_F = 0.0
    Coord_F(1:nDim) = real(Coord_D)
    if (is_point_inside_node(Coord_F, int(iNode))) then
       c_is_point_inside_node = 1
    else
       c_is_point_inside_node = 0
    end if
  end function c_is_point_inside_node
!==============================================================================

  subroutine c_show_tree(StringIn) bind(c)
    character(kind=c_char), intent(in) :: StringIn(*)
    character(len=100) :: String_F
    integer :: i
  !--------------------------------------------------------------------------------------------------
    i = 1
    String_F = ''
    do while (StringIn(i) /= c_null_char)
       String_F(i:i) = StringIn(i)
       i = i + 1
    end do
    if (i > 1) then
       call show_tree(trim(String_F))
    else
       call show_tree('')
    end if
  end subroutine c_show_tree
!==============================================================================

  ! --- BATL_geometry subroutines ---

  subroutine c_init_geometry(TypeGeometryIn) bind(c)
    character(kind=c_char), intent(in) :: TypeGeometryIn(*)
    character(len=100) :: Type_F
    integer :: i
  !--------------------------------------------------------------------------------------------------
    i = 1
    Type_F = ''
    do while (TypeGeometryIn(i) /= c_null_char)
       Type_F(i:i) = TypeGeometryIn(i)
       i = i + 1
    end do
    call init_geometry(TypeGeometryIn=trim(Type_F))
  end subroutine c_init_geometry
!==============================================================================

  subroutine c_xyz_to_coord(Xyz_D, Coord_D) bind(c)
    real(c_double), intent(in)  :: Xyz_D(nDim)
    real(c_double), intent(out) :: Coord_D(nDim)
    real :: Xyz_F(MaxDim), Coord_F(MaxDim)
  !--------------------------------------------------------------------------------------------------
    Xyz_F = 0.0
    Xyz_F(1:nDim) = real(Xyz_D)
    call xyz_to_coord(Xyz_F, Coord_F)
    Coord_D = real(Coord_F(1:nDim), c_double)
  end subroutine c_xyz_to_coord
!==============================================================================

  subroutine c_coord_to_xyz(Coord_D, Xyz_D) bind(c)
    real(c_double), intent(in)  :: Coord_D(nDim)
    real(c_double), intent(out) :: Xyz_D(nDim)
    real :: Coord_F(MaxDim), Xyz_F(MaxDim)
  !--------------------------------------------------------------------------------------------------
    Coord_F = 0.0
    Coord_F(1:nDim) = real(Coord_D)
    call coord_to_xyz(Coord_F, Xyz_F)
    Xyz_D = real(Xyz_F(1:nDim), c_double)
  end subroutine c_coord_to_xyz
!==============================================================================

  ! --- High-level Unit Tests ---

  subroutine c_test_tree() bind(c)
    use BATL_unit_test, ONLY: test_tree
  !--------------------------------------------------------------------------------------------------
    call test_tree()
  end subroutine c_test_tree
!==============================================================================

  subroutine c_test_geometry() bind(c)
    use BATL_unit_test, ONLY: test_geometry
  !--------------------------------------------------------------------------------------------------
    call test_geometry()
  end subroutine c_test_geometry
!==============================================================================

  ! --- Pointer Accessors ---

  function c_batl_get_xyz_ptr() bind(c) result(res)
    type(c_ptr) :: res
  !--------------------------------------------------------------------------------------------------
    res = c_loc(Xyz_DGB)
  end function c_batl_get_xyz_ptr
!==============================================================================

end module BATL_c
!==============================================================================
