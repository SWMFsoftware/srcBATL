!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_high_order
  implicit none

  save
  
  private ! except

  ! Make (f(j+1/2) - f(j-1/2))/dx high order accurate. 
  public:: correct_face_value

  ! Input:  f(i-3)...f(i+3), dx
  ! Output: (f(i+1/2) - f(i-1/2))/dx, which is 6th order accurate.
  public:: calc_center_first_derivate

  ! Input:  f(i-2)...f(i+3)
  ! Output: f(i+1/2), which is 6th order accurate. 
  public:: calc_face_value

  ! 5th order prolongation for resolution change.
  public:: prolongation_high_order_for_face_ghost

  ! 5th order restriction for resolution change. 
  public:: restriction_high_order_reschange

  ! Make sure all the ghost cells are high order accurate.
  public:: correct_face_ghost_for_fine_block

  ! 5th order AMR.
  public:: prolongation_high_order_amr
  public:: restriction_high_order_amr

  public:: limit_interpolation
contains

  !======================================================================
  real function correct_face_value(FaceValue, CellValue_I)

    ! FaceValue is at cell face. CellValue_I are cell centered.
    ! Return 6th order approximation

    real, intent(in):: CellValue_I(4), FaceValue
    real:: Der2, Der4
    real, parameter:: c1over6 = 1./6, c1over180 = 1./180
    !--------------------------------------------------------------------
    Der2 = c1over6*(CellValue_I(2) - 2*FaceValue + CellValue_I(3))
    Der4 = c1over180*(16*FaceValue - &
         9*(CellValue_I(2) + CellValue_I(3)) + &
         CellValue_I(1) + CellValue_I(4))
    correct_face_value = FaceValue - Der2 + Der4

  end function correct_face_value
  !======================================================================

  real function calc_center_first_derivate(CellValue_I, DxIn)

    ! Calculate df/dx at x=x_i with f(k), where k = i-3,i-2 ... i+3.
    ! Directly combine CellValue_I to get df/dx can save some 
    ! computation compare to the approach in this subroutine.

    real, intent(in):: CellValue_I(7)
    real, optional, intent(in):: DxIn
    real:: Dx
    real:: FaceL, FaceR ! FaceL = F_{i-1/2}, FaceR = F_{i+1/2}
    real:: CorrectedFaceL, CorrectedFaceR
    real, parameter:: c3over256 = 3./256, c25over256 = 25./256, &
         c150over256 = 150./256
    integer, parameter:: i = 4
    ! ----------------------------------------------------------------------

    Dx = 1.0
    if(present(DxIn)) Dx = DxIn

    FaceL = calc_face_value(CellValue_I(1:6))
    FaceR = calc_face_value(CellValue_I(2:7))

    CorrectedFaceL = correct_face_value(FaceL, CellValue_I(i-2:i+1))
    CorrectedFaceR = correct_face_value(FaceR, CellValue_I(i-1:i+2)) 

    calc_center_first_derivate = &
         (CorrectedFaceR - CorrectedFaceL)/Dx
  end function calc_center_first_derivate
  !======================================================================

  real function calc_face_value(CellValue_I, DoLimitIn)
    ! Calculate f_{i+1/2} with f(k), where k = i-2,i-1 ... i+3
    real, intent(in):: CellValue_I(6)
    logical, optional, intent(in):: DoLimitIn
    logical:: DoLimit
    real, parameter:: c3over256 = 3./256, c25over256 = 25./256, &
         c150over256 = 150./256
    real:: FaceValue
    real:: Distance_I(4) = (/-1.5, -0.5, 0.5, 1.5/)
    !----------------------------------------------------------------------

    DoLimit = .false. 
    if(present(DoLimitIn)) DoLimit = DoLimitIn

    FaceValue = c3over256*(CellValue_I(1) + CellValue_I(6)) - &
         c25over256*(CellValue_I(2) + CellValue_I(5)) + &
         c150over256*(CellValue_I(3) + CellValue_I(4))

    calc_face_value = FaceValue
    if(DoLimit) calc_face_value = &
         limit_interpolation(FaceValue, CellValue_I(2:5), Distance_I)
  end function calc_face_value

  !===========================================================================

  real function limit_interpolation(FaceOrig, CellValue_I, Distance_I, &
       MaxValueIn, MinValueIn)
    ! This Limiter works for ununiform grid interpolation. 
    ! See (2.18) in 'Accurate Monotonicity-Preserving Schemes with 
    ! Runge-Kutta Time Stepping' by A. Suresh & H. T. Huynh (1997)

    real, intent(in):: FaceOrig, CellValue_I(4), Distance_I(4)

    ! The return value should less/larger than MaxValueIn/MinValueIn if
    ! they present. Both of them should present at the same time.
    ! This option maybe useful, but have been not used at anywhere so far.
    real, optional, intent(in):: MaxValueIn, MinValueIn

    real:: FaceL, FaceR, FaceAV, FaceMD, FaceMin, FaceMax, MaxValue, MinValue
    real:: MP5Result
    !----------------------------------------------------------------------

    FaceAV = two_points_interpolation(CellValue_I(2:3), Distance_I(2:3)) 
    FaceL = two_points_interpolation(CellValue_I(1:2), Distance_I(1:2))
    FaceR = two_points_interpolation(CellValue_I(3:4), Distance_I(3:4))
    
    FaceMD = median(FaceAV, FaceL, FaceR)
    FaceMin = min(FaceMD, CellValue_I(2), CellValue_I(3))  
    FaceMax = max(FaceMD, CellValue_I(2), CellValue_I(3))

    MP5Result = median(FaceOrig, FaceMin, FaceMax)

    MaxValue = MP5Result
    if(present(MaxValueIn)) MaxValue = MaxValueIn
    MinValue = MP5Result
    if(present(MinValueIn)) MinValue = MinValueIn

    limit_interpolation = median(MP5Result, MaxValue, MinValue)
  end function limit_interpolation

  !======================================================================

  real function median(a,b,c)
    real, intent(in):: a, b, c
    median = max(min(a, max(b,c)), min(b,c))
  end function median
  !======================================================================

  real function two_points_interpolation(Cell_I, Distance_I)
    ! Cell_I(i) is at xi, calculate the value at x=0. 
    ! Distance_I(i) = xi - x0
    real, intent(in):: Cell_I(2), Distance_I(2)
    real:: c1, c2
    !----------------------------------------------------------------------
    c1 = Distance_I(2)/(Distance_I(2) - Distance_I(1))
    c2 = -Distance_I(1)/(Distance_I(2) - Distance_I(1))
    two_points_interpolation = c1*Cell_I(1) + c2*Cell_I(2)
  end function two_points_interpolation

  !======================================================================
  subroutine restriction_high_order_reschange(CoarseCell, FineCell_III, &
       Ghost_I, DoSymInterpIn)
    ! For 2D: 
    !         _________________________________
    !         | u1|   |   |   |   |   |   |   |
    !         |___|___|___|___|___|___|___|___|
    !         | u2|   |   |   |   |   |   |   |
    !  _______|___|___|___|___|___|___|___|___|
    ! |       | u3| u7|   |   |   |   |   |   |
    ! |Coarse |__G1___|__G2___|__G3___|___|___|
    ! | u0    | u4| u8|   |   |   |   |   |   |
    ! |_______|___|___|___|___|___|___|___|___|
    !         | u5|   |   |   |   |   |   |   |
    !         |___|___|___|___|___|___|___|___|
    !         | u6|   |   |   |   |   |   |   |
    !         |___|___|___|___|___|___|___|___|


    ! First calculate the face value (f1) between u3 and u4 with u1, u2...u6.
    ! Face value between u7 and u8 (f2) and more face value can be got in the
    ! same way. 8 face values are needed. 

    ! Use u0, f1, f2, f3, f4 to calculate the ghost cell G1. 
    ! Use f1, f2, f3, f4, f5, f6 to interpolate G2
    ! Use f3, f4...f7,f8 to interpolate the third ghost cell (G3). But if 
    ! there are only 6 cells, then use f1..f6 to interpolate G3 (not 
    ! DoSymInterp).

    use BATL_size, ONLY: nJ, nk

    ! In 2D there is only 1 cell in the K direction
    integer, parameter:: k6_ = min(nK, 6), j6_ = min(nJ,6)

    real, intent(in) :: CoarseCell              ! value of coarse neighbor cell
    real, intent(in) :: FineCell_III(8,j6_,k6_) ! value of local fine cells
    real, intent(out):: Ghost_I(3)              ! coarse ghost cells for neighbor
    logical, optional, intent(in) :: DoSymInterpIn

    ! Local variables
    real:: FineCell_I(8), FaceValue_I(6)
    real:: Ghost, Cell_I(4)
    real:: CellInput_I(6), CellInput_II(6,6)
    logical:: DoSymInterp

    ! Distances of coarse neighbor and fine cells from 1st ghost cell
    real:: Distance_I(4)=(/-4,-1,1,3/)
    real:: Distance1_I(4) = (/-3,-1,1,3/)

    ! Interpolation coefficients for G1
    real, parameter:: c1=-1./63, c2=5./12, c3 = 3./4, c4=-5./28, c5=1./36
    
    real:: Coef_I(6) = ([7./256, -45./256, 63./128, -105./128, &
         315./256, 63./256])

    logical, parameter:: DoLimit = .true. ! change only for debugging
    integer:: i, j, k
    character(len=*), parameter:: NameSub = 'restriction_high_order_reschange'
    !----------------------------------------------------------------------
    DoSymInterp = .true. ! Use f3...f8 to interpolate G3.
    if(present(DoSymInterpIn)) DoSymInterp = DoSymInterpIn

    ! Integerpolate fine cell centers to line connecting coarse cell centers
    if(nK == 1) then  ! 2D resolution change       
       do i = 1, 8          
          ! Use a temporary variable CellInput_I to avoid compile error. 
          CellInput_I(1:j6_) = FineCell_III(i,:,1)
          FineCell_I(i) = calc_face_value(CellInput_I, DoLimit)
       enddo
    else   
       ! 3D resolution change. Need more tests to make sure it works!!
       do i = 1, 8
          CellInput_II(1:j6_,1:k6_) = FineCell_III(I,:,:)
          FineCell_I(i) = &
               calc_edge_value(CellInput_II,Dolimit)
       enddo
    endif

    ! High order interpolation for first ghost cell using coarse neighbor
    Ghost = c1*CoarseCell + c2*FineCell_I(1) + c3*FineCell_I(2) + &
         c4*FineCell_I(3) + c5*FineCell_I(4)

    ! Limit value if requested
    if(DoLimit) then
       Cell_I(1)   = CoarseCell
       Cell_I(2:4) = FineCell_I(1:3)
       Ghost_I(1)  = limit_interpolation(Ghost, Cell_I, Distance_I)
    else
       Ghost_I(1) = Ghost
    endif

    ! Interpolate G2 and G3 from fine edge values
    Ghost_I(2) = calc_face_value(FineCell_I(1:6), DoLimit)
    if(DoSymInterp) then
       Ghost_I(3) = calc_face_value(FineCell_I(3:8), DoLimit)
    else
       Ghost = sum(FineCell_I(1:6)*Coef_I)
       if(DoLimit) then
          Cell_I(1:3) = FineCell_I(4:6)
          Cell_I(4)  = FineCell_I(6)
          Ghost_I(3) = limit_interpolation(Ghost, Cell_I, Distance1_I)
       else
          Ghost_I(3) = Ghost
       endif
    endif

  end subroutine restriction_high_order_reschange
  !======================================================================
  real function calc_edge_value(CellValue_II,DoLimitIn)
    ! For 3D, need more tests. 
    real, intent(in) :: CellValue_II(6,6)
    logical, optional, intent(in):: DoLimitIn
    logical:: DoLimit
    real:: CellValue_I(6)
    integer:: i, j
    integer:: iBegin=1, iEnd=6
    !----------------------------------------------------------------------
    DoLimit = .true. 
    if(present(DoLimitIn)) DoLimit = DoLimitIn

    do i = iBegin, iEnd       
       CellValue_I(i) = calc_face_value(CellValue_II(i,:),DoLimit)
    enddo
    calc_edge_value = calc_face_value(CellValue_I,DoLimit)
  end function calc_edge_value
  !======================================================================

  subroutine get_ghost_for_fine_blk(CoarseCell_III, FineCell_I, Ghost_I, &
       Use4thOrderIn)
    ! 2D: 
    ! __________________________
    ! |        |       |       |
    ! |        |       |  jmm  |
    ! |        |       |       |
    ! |________|_______|_______|
    ! |        |       |       |
    ! |        |       |  jm   |
    ! |        |       |       |
    ! |________|_______|_______|_____________
    ! |        |       |       |   |    |   |
    ! |        |       |  j2   |___|____|___|
    ! |   *    |  *  ? |?  * ? |u0 | u1 |u2 |
    ! |________|_______|_______|___|____|___|
    ! |        |       |       |
    ! |        |       |  jp   |
    ! |        |       |       |
    ! |________|_______|_______|
    ! |        |       |       | 
    ! |        |       | jpp   |
    ! |        |       |       |
    ! | _______|_______|_______|

    ! The 3*5 cells represent CoarseCell_II. u0, u1, and u2 are FineCell_I.
    ! '?' are the ghost values needed. 
    ! First calculate the values represented by '*', and then use these star 
    ! values and u0, u1, u2 to calculate the ghost cells. 
    use BATL_size, ONLY: nK

    real, intent(in):: CoarseCell_III(5,5,3), FineCell_I(3)
    real, intent(out):: Ghost_I(3)
    logical, optional, intent(in):: Use4thOrderIn

    real:: CoarseCell_I(3), CoarseCell_II(5,3)
    integer:: i
    logical:: DoLimit = .true., Use4thOrder = .true.     
    !----------------------------------------------------------------------

    DoLimit = .true. 
    Use4thOrder = .false. 
    if(present(Use4thOrderIn)) Use4thOrder = Use4thOrderIn

    if(nK == 1) then
       CoarseCell_II = CoarseCell_III(:,1,:)
       do i = 1, 3       
          CoarseCell_I(i) = &
               interpolate_in_coarse_blk_1d(CoarseCell_II(:,i), &
               DoLimitIn=DoLimit)
       enddo
    else
       do i = 1, 3
          CoarseCell_I(i) = interpolate_in_coarse_blk_2d(&
               CoarseCell_III(:,:,i), DoLimitIn=DoLimit,&
               Use4thOrderIn=Use4thOrder)
       enddo
    endif

    call interpolate_ghost_for_fine_blk

  contains
    !----------------------------------------------------------------------
    subroutine interpolate_ghost_for_fine_blk
      real, parameter:: c11=-4./231, c12=4./7,c13=5./7, c14=-1./3, c15=5./77
      real, parameter:: c21=-9./572, c22=1./6,c23=1.05, c24=-3./11, c25=14./195
      real, parameter:: c31=-9./286, c32=5./7,c33=0.5, c34=-20./77, c35=1./13
      real:: Ghost, Cell_I(4), Distance_I(4)

      !-----------------
      Ghost = c11*CoarseCell_I(2) + c12*CoarseCell_I(1) + &
           c13*FineCell_I(1) + c14*FineCell_I(2) + c15*FineCell_I(3)

      if(DoLimit) then
         Distance_I(1) = -5; Distance_I(2) = -1
         Distance_I(3) = 2; Distance_I(4) = 4
         Cell_I(1) = CoarseCell_I(2); Cell_I(2) = CoarseCell_I(1)
         Cell_I(3) = FineCell_I(1); Cell_I(4) = FineCell_I(2)

         Ghost_I(1) = limit_interpolation(Ghost, Cell_I, Distance_I)
      else
         Ghost_I(1) = Ghost
      endif
      !---------------------

      Ghost = c21*CoarseCell_I(3) + c22*CoarseCell_I(2) + &
           c23*CoarseCell_I(1) + c24*FineCell_I(1) + c25*FineCell_I(2)

      if(DoLimit) then
         Distance_I(1) = -7; Distance_I(2) = -3
         Distance_I(3) = 1; Distance_I(4) = 4
         Cell_I(1) = CoarseCell_I(3); Cell_I(2) = CoarseCell_I(2)
         Cell_I(3) = CoarseCell_I(1); Cell_I(4) = FineCell_I(1)

         Ghost_I(2) = limit_interpolation(Ghost, Cell_I, Distance_I)
      else
         Ghost_I(2) = Ghost
      endif
      !----------------------

      Ghost = c31*CoarseCell_I(3) + c32*CoarseCell_I(2) + &
           c33*CoarseCell_I(1) + c34*FineCell_I(1) + c35*FineCell_I(2)

      if(DoLimit) then
         Distance_I(1) = -5; Distance_I(2) = -1
         Distance_I(3) = 3; Distance_I(4) = 6
         Cell_I(1) = CoarseCell_I(3); Cell_I(2) = CoarseCell_I(2)
         Cell_I(3) = CoarseCell_I(1); Cell_I(4) = FineCell_I(1)

         Ghost_I(3) = limit_interpolation(Ghost, Cell_I, Distance_I)
      else
         Ghost_I(3) = Ghost
      endif

    end subroutine interpolate_ghost_for_fine_blk

  end subroutine get_ghost_for_fine_blk
  !======================================================================
  real function interpolate_in_coarse_blk_1d(Cell_I, DoLimitIn, Use4thOrderIn)
    real, intent(in):: Cell_I(5)
    logical, optional, intent(in):: DoLimitIn, Use4thOrderIn
    logical:: DoLimit, Use4thOrder
    real   :: cpp, cp, c0, cm, cmm

    integer, parameter:: Ipp_=1, Ip_=2, I_=3, Im_=4, Imm_=5
    real:: Temp, Distance_I(4)=(/-7,-3,1,5/)
    !----------------------------------------------------------------------

    DoLimit = .true.    
    if(present(DoLimitIn)) DoLimit = DoLimitIn

    Use4thOrder = .false. 
    if(present(Use4thOrderIn)) Use4thOrder = Use4thOrderIn
    if(Use4thOrder) then 
       ! Cell_I(5) is not used. 
       cpp = -5./128;  cp = 35./128
       c0  = 105./128; cm = -7./128
       cmm = 0
    else
       cpp=-45./2048; cp=105./512
       c0=945./1024 ; cm=-63./512 
       cmm=35./2048
    endif

    Temp = cpp*Cell_I(Ipp_) + cp*Cell_I(Ip_) + &
         c0*Cell_I(I_) + cm*Cell_I(Im_) + cmm*Cell_I(Imm_)        
    if(DoLimit) then
       interpolate_in_coarse_blk_1d = limit_interpolation(Temp, Cell_I(Ipp_:Im_), Distance_I)
    else
       interpolate_in_coarse_blk_1d = Temp
    endif
  end function interpolate_in_coarse_blk_1d
  !======================================================================
  real function interpolate_in_coarse_blk_2d(Cell_II, DoLimitIn, Use4thOrderIn)
    real, intent(in):: Cell_II(5,5)
    logical, optional,intent(in):: DoLimitIn, Use4thOrderIn
    logical:: DoLimit, Use4thOrder
    real:: Cell_I(5)
    integer:: i
    
    DoLimit = .true.
    if(present(DoLimitIn)) DoLimit = DoLimitIn

    Use4thOrder = .false. 
    if(present(Use4thOrderIn)) Use4thOrder = Use4thOrderIn

    do i = 1, 5 ! Eliminate j dimension
       Cell_I(i) = &
            interpolate_in_coarse_blk_1d(Cell_II(i,:), DoLimit, Use4thOrder)
    enddo
    interpolate_in_coarse_blk_2d = &
         interpolate_in_coarse_blk_1d(Cell_I, DoLimit, Use4thOrder)    
  end function interpolate_in_coarse_blk_2d
  !======================================================================

  subroutine prolongation_high_order_for_face_ghost(&
       iBlock, nVar, Field1_VG, Field_VG, Do5thFace_G)
    ! High order prolongation for simple resolution change (resolution 
    ! change in only one direction).

    use BATL_tree, ONLY: DiLevelNei_IIIB, iNodeNei_IIIB
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         nI, nJ, nK, i0_,i2_,j0_, j2_, nJp1_, nJm1_, k0_, k2_, nKp1_, nKm1_, &
         jm2_,jm1_,nJm2_,nJm1_,nJp2_,nJp3_,km2_,km1_,nKm2_,nKm1_,nKp2_,&
         nKp3_,im2_,im1_,nIm2_,nIm1_,nIp1_,nIp2_,nIp3_, i3_,j3_,k3_

    integer, intent(in) :: iBlock, nVar
    real, intent(inout) :: Field1_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(inout) :: Field_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    logical, intent(in) :: Do5thFace_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    integer :: i1, j1, k1, i2, j2, k2
    integer :: ip, im, jp, jm, kp, km, iVar, jpp, jmm, ipp, imm, kpp, kmm

    real :: FieldCoarse_VII(nVar,5,3), FieldFine_VI(nVar,3),Ghost_I(3)
    real :: FieldCoarse_VIII(nVar,5,5,3)
    integer, parameter:: Ipp_=1, Ip_=2, I_=3, Im_=4, Imm_=5

    integer:: iNode1, iNode2, iNode3, iNode0

    integer:: Index1_I(5), Index2_I(5)
    integer:: iPara1,iPara2, i, j, k

    ! If it is non-simple resolution change, some face ghost cells does not
    ! have enough information to reach 5th order accuracy, then use 4th order
    ! interpolation. These 4th order face ghost cells will be sent to neighbour 
    ! block as edge/corner ghost cells. As face ghost cell value, it will 
    ! be overwritten by remote prolongation (iSendStage == 3). 
    logical:: Use4thOrder
    !-------------------------------------------------------------------------

    Field1_VG = Field_VG

    ! Do six faces
    if(DiLevelNei_IIIB(-1,0,0,iBlock) == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2;
          do k2 = k1,k1+min(1,nK-1); do j2 = j1,j1+1
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             jpp = 7*j2 - 6*j1 -3; jmm = 8*j1 - 7*j2 + 4

             iNode0 = iNodeNei_IIIB(0,0,1,iBlock)
             iNode1 = iNodeNei_IIIB(0,1,1,iBlock)
             iNode2 = iNodeNei_IIIB(0,2,1,iBlock)
             iNode3 = iNodeNei_IIIB(0,3,1,iBlock)

             if(DiLevelNei_IIIB(-1,-1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,-1,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                ! Be careful what are stored in the corner ghost cells!!
                jpp = max(-1, jpp)
                jmm = max(-1, jmm)
             endif

             if(DiLevelNei_IIIB(-1,1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,1,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                jpp = min(jpp, nJ+2)
                jmm = min(jmm, nJ+2)
             endif

             if(nK == 1)then
                kp  = 1; km  = 1
                kpp = 1; kmm = 1
             else
                kp  = 3*k2 - 2*k1 -1; km  = 4*k1 -3*k2 +2
                kpp = 7*k2 - 6*k1 -3; kmm = 8*k1 -7*k2 +4
             end if

             iNode0 = iNodeNei_IIIB(0,1,0,iBlock)
             iNode1 = iNodeNei_IIIB(0,1,1,iBlock)
             iNode2 = iNodeNei_IIIB(0,1,2,iBlock)
             iNode3 = iNodeNei_IIIB(0,1,3,iBlock)

             if(DiLevelNei_IIIB(-1,0,-1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,0,-1,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                kpp = max(-1, kpp)
                kmm = max(-1, kmm)
             endif

             if(DiLevelNei_IIIB(-1,0,1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,0,1,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                kpp = min(kpp, nK+2)
                kmm = min(kmm, nK+2)
             endif

             ! Works for both 2D and 3D.
             Index1_I(1) = jpp
             Index1_I(2) = jp
             Index1_I(3) = j2
             Index1_I(4) = jm
             Index1_I(5) = jmm

             Index2_I(1) = kpp
             Index2_I(2) = kp
             Index2_I(3) = k2
             Index2_I(4) = km 
             Index2_I(5) = kmm

             do k = 1, min(5,nK)
                iPara2 = Index2_I(k)
                do j = 1, 5
                   iPara1 = Index1_I(j)
                   FieldCoarse_VIII(:,j,k,1:abs(im2_-i0_)+1) = &
                        Field1_VG(:,&
                        i0_:im2_:sign(1,im2_-i0_),&
                        iPara1,iPara2)
                enddo
             enddo

             FieldFine_VI(:,1) = Field1_VG(:,1,  j2,k2)
             FieldFine_VI(:,2) = Field1_VG(:,i2_,j2,k2)
             FieldFine_VI(:,3) = Field1_VG(:,i3_,j2,k2)

             Use4thOrder = .not.Do5thFace_G(i0_,j2,k2)
             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VIII(iVar,:,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I, Use4thOrder)
                Field_VG(iVar,i0_ ,j2,k2) = Ghost_I(1)
                Field_VG(iVar,im1_,j2,k2) = Ghost_I(2)
                Field_VG(iVar,im2_,j2,k2) = Ghost_I(3)
             enddo

          end do; end do
       end do; end do
    end if

    if(DiLevelNei_IIIB(1,0,0,iBlock) == 1)then
       do k1=1, nK, 2; do j1=1, nJ, 2
          do k2 = k1,k1+min(1,nK-1); do j2 = j1,j1+1
             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             jpp = 7*j2 - 6*j1 -3; jmm = 8*j1 - 7*j2 + 4

             iNode0 = iNodeNei_IIIB(3,0,1,iBlock)
             iNode1 = iNodeNei_IIIB(3,1,1,iBlock)
             iNode2 = iNodeNei_IIIB(3,2,1,iBlock)
             iNode3 = iNodeNei_IIIB(3,3,1,iBlock)

             if(DiLevelNei_IIIB(1,-1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,-1,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                jpp = max(-1, jpp)
                jmm = max(-1, jmm)
             endif

             if(DiLevelNei_IIIB(1,1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,1,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                jpp = min(jpp, nJ+2)
                jmm = min(jmm, nJ+2)
             endif

             if(nK == 1)then
                kp = 1; km = 1
                kpp= 1; kmm= 1
             else
                kp  = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
                kpp = 7*k2 - 6*k1 -3; kmm = 8*k1 -7*k2 +4
             end if

             iNode0 = iNodeNei_IIIB(3,1,0,iBlock)
             iNode1 = iNodeNei_IIIB(3,1,1,iBlock)
             iNode2 = iNodeNei_IIIB(3,1,2,iBlock)
             iNode3 = iNodeNei_IIIB(3,1,3,iBlock)

             if(DiLevelNei_IIIB(1,0,-1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,0,-1,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                kpp = max(-1, kpp)
                kmm = max(-1, kmm)
             endif

             if(DiLevelNei_IIIB(1,0,1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,0,1,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                kpp = min(kpp, nK+2)
                kmm = min(kmm, nK+2)
             endif

             Index1_I(1) = jpp
             Index1_I(2) = jp
             Index1_I(3) = j2
             Index1_I(4) = jm
             Index1_I(5) = jmm

             Index2_I(1) = kpp
             Index2_I(2) = kp
             Index2_I(3) = k2
             Index2_I(4) = km 
             Index2_I(5) = kmm

             do k = 1, min(5,nK)
                iPara2 = Index2_I(k)
                do j = 1, 5
                   iPara1 = Index1_I(j)
                   FieldCoarse_VIII(:,j,k,1:abs(nIp3_-nIp1_)+1) = &
                        Field1_VG(:,nIp1_:nIp3_,iPara1,iPara2)
                enddo
             enddo

             FieldFine_VI(:,1) = Field1_VG(:,nI,j2,k2)
             FieldFine_VI(:,2) = Field1_VG(:,nIm1_,j2,k2)
             FieldFine_VI(:,3) = Field1_VG(:,nIm2_,j2,k2)

             Use4thOrder = .not.Do5thFace_G(nIp1_,j2,k2)
             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VIII(iVar,:,:,:),&
                     FieldFine_VI(iVar,:), Ghost_I, Use4thOrder)

                Field_VG(iVar,nIp1_,j2,k2) = Ghost_I(1)
                Field_VG(iVar,nIp2_,j2,k2) = Ghost_I(2)
                Field_VG(iVar,nIp3_,j2,k2) = Ghost_I(3)
             enddo

          end do; end do
       end do; end do
    end if

    if(DiLevelNei_IIIB(0,-1,0,iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2
          do k2 = k1,k1+min(1,nK-1); do i2 = i1,i1+1
             ip  = 3*i2 - 2*i1 -1; im  = 4*i1 - 3*i2 + 2
             ipp = 7*i2 - 6*i1 -3; imm = 8*i1 - 7*i2 + 4

             iNode0 = iNodeNei_IIIB(0,0,1,iBlock)
             iNode1 = iNodeNei_IIIB(1,0,1,iBlock)
             iNode2 = iNodeNei_IIIB(2,0,1,iBlock)
             iNode3 = iNodeNei_IIIB(3,0,1,iBlock)

             if(DiLevelNei_IIIB(-1,-1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(-1,0,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                ipp = max(-1, ipp)
                imm = max(-1, imm)
             endif

             if(DiLevelNei_IIIB(1,-1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(1,0,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                ipp = min(ipp, nI+2)
                imm = min(imm, nI+2)
             endif

             if(nK == 1)then
                kp  = 1; km  = 1
                kpp = 1; kmm = 1
             else
                kp  = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
                kpp = 7*k2 - 6*k1 -3; kmm = 8*k1 -7*k2 +4
             end if


             iNode0 = iNodeNei_IIIB(1,0,0,iBlock)
             iNode1 = iNodeNei_IIIB(1,0,1,iBlock)
             iNode2 = iNodeNei_IIIB(1,0,2,iBlock)
             iNode3 = iNodeNei_IIIB(1,0,3,iBlock)

             if(DiLevelNei_IIIB(0,-1,-1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,0,-1,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                kpp = max(-1, kpp)
                kmm = max(-1, kmm)
             endif

             if(DiLevelNei_IIIB(0,-1,1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,0,1,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                kpp = min(kpp, nK+2)
                kmm = min(kmm, nK+2)
             endif

             Index1_I(1) = ipp
             Index1_I(2) = ip
             Index1_I(3) = i2
             Index1_I(4) = im
             Index1_I(5) = imm

             Index2_I(1) = kpp
             Index2_I(2) = kp
             Index2_I(3) = k2
             Index2_I(4) = km 
             Index2_I(5) = kmm

             do k = 1, min(5,nK)
                iPara2 = Index2_I(k)
                do i = 1, 5
                   iPara1 = Index1_I(i)
                   FieldCoarse_VIII(:,i,k,1:abs(jm2_-j0_)+1) = &
                        Field1_VG(:,iPara1,&
                        j0_:jm2_:sign(1,jm2_-j0_),&
                        iPara2)
                enddo
             enddo

             FieldFine_VI(:,1) = Field1_VG(:,i2,1,k2)
             FieldFine_VI(:,2) = Field1_VG(:,i2,j2_,k2)
             FieldFine_VI(:,3) = Field1_VG(:,i2,j3_,k2)

             Use4thOrder = .not.Do5thFace_G(i2,j0_,k2)
             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VIII(iVar,:,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I, Use4thOrder)
                Field_VG(iVar,i2,j0_ ,k2) = Ghost_I(1)
                Field_VG(iVar,i2,jm1_,k2) = Ghost_I(2)
                Field_VG(iVar,i2,jm2_,k2) = Ghost_I(3)
             enddo

          end do; end do
       end do; end do
    end if

    if(DiLevelNei_IIIB(0,1,0,iBlock) == 1)then
       do k1=1, nK, 2; do i1=1, nI, 2
          do k2 = k1,k1+min(1,nK-1); do i2 = i1,i1+1
             ip = 3*i2 - 2*i1 -1 ; im = 4*i1 -3*i2 +2
             ipp = 7*i2 - 6*i1 -3; imm = 8*i1 - 7*i2 + 4

             iNode0 = iNodeNei_IIIB(0,3,1,iBlock)
             iNode1 = iNodeNei_IIIB(1,3,1,iBlock)
             iNode2 = iNodeNei_IIIB(2,3,1,iBlock)
             iNode3 = iNodeNei_IIIB(3,3,1,iBlock)

             if(DiLevelNei_IIIB(-1,1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(-1,0,0,iBlock) == 1 .or. &
                  iNode0 /=iNode1)) then
                ipp = max(-1, ipp)
                imm = max(-1, imm)
             endif

             if(DiLevelNei_IIIB(1,1,0,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(1,0,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                ipp = min(ipp, nI+2)
                imm = min(imm, nI+2)
             endif

             if(nK == 1)then
                kp = 1; km = 1
                kpp = 1; kmm = 1
             else
                kp = 3*k2 - 2*k1 -1 ; km = 4*k1 -3*k2 +2
                kpp = 7*k2 - 6*k1 -3; kmm = 8*k1 -7*k2 +4
             end if

             iNode0 = iNodeNei_IIIB(1,3,0,iBlock)
             iNode1 = iNodeNei_IIIB(1,3,1,iBlock)
             iNode2 = iNodeNei_IIIB(1,3,2,iBlock)
             iNode3 = iNodeNei_IIIB(1,3,3,iBlock)

             if(DiLevelNei_IIIB(0,1,-1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,0,-1,iBlock) == 1 .or. &
                  iNode0 /=iNode1)) then
                kpp = max(-1, kpp)
                kmm = max(-1, kmm)
             endif

             if(DiLevelNei_IIIB(0,1,1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,0,1,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                kpp = min(kpp, nK+2)
                kmm = min(kmm, nK+2)
             endif

             Index1_I(1) = ipp
             Index1_I(2) = ip
             Index1_I(3) = i2
             Index1_I(4) = im
             Index1_I(5) = imm

             Index2_I(1) = kpp
             Index2_I(2) = kp
             Index2_I(3) = k2
             Index2_I(4) = km 
             Index2_I(5) = kmm

             do k = 1, min(5,nK)
                iPara2 = Index2_I(k)
                do i = 1, 5
                   iPara1 = Index1_I(i)
                   FieldCoarse_VIII(:,i,k,1:abs(nJp3_-nJp1_)+1) = &
                        Field1_VG(:,iPara1,&
                        nJp1_:nJp3_:sign(1,nJp3_-nJp1_),&
                        iPara2)
                enddo
             enddo

             FieldFine_VI(:,1) = Field1_VG(:,i2,nJ,   k2)
             FieldFine_VI(:,2) = Field1_VG(:,i2,nJm1_,k2)
             FieldFine_VI(:,3) = Field1_VG(:,i2,nJm2_,k2)

             Use4thOrder = .not.Do5thFace_G(i2,nJp1_,k2)


             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VIII(iVar,:,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I, Use4thOrder)
                Field_VG(iVar,i2,nJp1_, k2) = Ghost_I(1)
                Field_VG(iVar,i2,nJp2_, k2) = Ghost_I(2)
                Field_VG(iVar,i2,nJp3_, k2) = Ghost_I(3)
             enddo
          end do; end do
       end do; end do
    end if


    if(DiLevelNei_IIIB(0,0,-1,iBlock) == 1) then
       do j1 = 1, nJ, 2; do i1 = 1, nI, 2
          do j2 = j1, j1+1; do i2 = i1, i1+1
             ip  = 3*i2 - 2*i1 -1; im  = 4*i1 - 3*i2 + 2
             ipp = 7*i2 - 6*i1 -3; imm = 8*i1 - 7*i2 + 4

             iNode0 = iNodeNei_IIIB(0,1,0,iBlock)
             iNode1 = iNodeNei_IIIB(1,1,0,iBlock)
             iNode2 = iNodeNei_IIIB(2,1,0,iBlock)
             iNode3 = iNodeNei_IIIB(3,1,0,iBlock)

             if(DiLevelNei_IIIB(-1,0,-1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(-1,0,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                ipp = max(-1, ipp)
                imm = max(-1, imm)
             endif

             if(DiLevelNei_IIIB(1,0,-1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(1,0,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                ipp = min(ipp, nI+2)
                imm = min(imm, nI+2)
             endif


             jp = 3*j2 - 2*j1 -1 ; jm = 4*j1 -3*j2 +2
             jpp = 7*j2 - 6*j1 -3; jmm = 8*j1 - 7*j2 + 4

             iNode0 = iNodeNei_IIIB(1,0,0,iBlock)
             iNode1 = iNodeNei_IIIB(1,1,0,iBlock)
             iNode2 = iNodeNei_IIIB(1,2,0,iBlock)
             iNode3 = iNodeNei_IIIB(1,3,0,iBlock)

             if(DiLevelNei_IIIB(0,-1,-1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,-1,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                jpp = max(-1, jpp)
                jmm = max(-1, jmm)
             endif

             if(DiLevelNei_IIIB(0,1,-1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,1,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                jpp = min(jpp, nJ+2)
                jmm = min(jmm, nJ+2)
             endif

             Index1_I(1) = ipp
             Index1_I(2) = ip
             Index1_I(3) = i2
             Index1_I(4) = im
             Index1_I(5) = imm

             Index2_I(1) = jpp
             Index2_I(2) = jp             
             Index2_I(3) = j2
             Index2_I(4) = jm 
             Index2_I(5) = jmm

             do j = 1, 5
                iPara2 = Index2_I(j)
                do i = 1, 5
                   iPara1 = Index1_I(i)
                   FieldCoarse_VIII(:,i,j,1:min(nk,3)) = &
                        Field1_VG(:,iPara1,iPara2,&
                        k0_:km2_:sign(1,km2_-k0_))
                enddo
             enddo

             FieldFine_VI(:,1) = Field1_VG(:,i2,j2,1)
             FieldFine_VI(:,2) = Field1_VG(:,i2,j2,k2_)
             FieldFine_VI(:,3) = Field1_VG(:,i2,j2,k3_)

             Use4thOrder = .not.Do5thFace_G(i2,j2,k0_) 
             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VIII(iVar,:,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I, Use4thOrder)
                Field_VG(iVar,i2,j2,k0_ ) = Ghost_I(1)
                Field_VG(iVar,i2,j2,km1_) = Ghost_I(2)
                Field_VG(iVar,i2,j2,km2_) = Ghost_I(3)
             enddo
          enddo; enddo
       enddo; enddo

    endif

    if(DiLevelNei_IIIB(0,0,1,iBlock) == 1) then
       do j1 = 1, nJ, 2; do i1 = 1, nI, 2
          do j2 = j1, j1+1; do i2 = i1, i1+1
             ip  = 3*i2 - 2*i1 -1; im  = 4*i1 - 3*i2 + 2
             ipp = 7*i2 - 6*i1 -3; imm = 8*i1 - 7*i2 + 4

             iNode0 = iNodeNei_IIIB(0,1,3,iBlock)
             iNode1 = iNodeNei_IIIB(1,1,3,iBlock)
             iNode2 = iNodeNei_IIIB(2,1,3,iBlock)
             iNode3 = iNodeNei_IIIB(3,1,3,iBlock)

             if(DiLevelNei_IIIB(-1,0,1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(-1,0,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                ipp = max(-1, ipp)
                imm = max(-1, imm)
             endif

             if(DiLevelNei_IIIB(1,0,1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(1,0,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                ipp = min(ipp, nI+2)
                imm = min(imm, nI+2)
             endif

             jp  = 3*j2 - 2*j1 -1; jm  = 4*j1 - 3*j2 +2
             jpp = 7*j2 - 6*j1 -3; jmm = 8*j1 - 7*j2 + 4

             iNode0 = iNodeNei_IIIB(1,0,3,iBlock)
             iNode1 = iNodeNei_IIIB(1,1,3,iBlock)
             iNode2 = iNodeNei_IIIB(1,2,3,iBlock)
             iNode3 = iNodeNei_IIIB(1,3,3,iBlock)

             if(DiLevelNei_IIIB(0,-1,1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,-1,0,iBlock) == 1 .or. &
                  iNode0 /= iNode1)) then
                jpp = max(-1, jpp)
                jmm = max(-1, jmm)
             endif

             if(DiLevelNei_IIIB(0,1,1,iBlock) == 1 .and. &
                  (DiLevelNei_IIIB(0,1,0,iBlock) == 1 .or. &
                  iNode2 /= iNode3)) then
                jpp = min(jpp, nJ+2)
                jmm = min(jmm, nJ+2)
             endif

             Index1_I(1) = ipp
             Index1_I(2) = ip
             Index1_I(3) = i2
             Index1_I(4) = im
             Index1_I(5) = imm

             Index2_I(1) = jpp
             Index2_I(2) = jp             
             Index2_I(3) = j2
             Index2_I(4) = jm 
             Index2_I(5) = jmm

             do j = 1, 5
                iPara2 = Index2_I(j)
                do i = 1, 5
                   iPara1 = Index1_I(i)
                   FieldCoarse_VIII(:,i,j,1:min(nk,3)) = &
                        Field1_VG(:,iPara1,iPara2,&
                        nKp1_:nKp3_:sign(1,nKp3_-nKp1_))
                enddo
             enddo

             FieldFine_VI(:,1) = Field1_VG(:,i2,j2,nk   )
             FieldFine_VI(:,2) = Field1_VG(:,i2,j2,nKm1_)
             FieldFine_VI(:,3) = Field1_VG(:,i2,j2,nKm2_)

             Use4thOrder = .not.Do5thFace_G(i2,j2,nKp1_)
             do iVar = 1, nVar
                call get_ghost_for_fine_blk(FieldCoarse_VIII(iVar,:,:,:), &
                     FieldFine_VI(iVar,:), Ghost_I, Use4thOrder)
                Field_VG(iVar,i2,j2,nKp1_) = Ghost_I(1)
                Field_VG(iVar,i2,j2,nKp2_) = Ghost_I(2)
                Field_VG(iVar,i2,j2,nKp3_) = Ghost_I(3)
             enddo
          enddo; enddo
       enddo; enddo


    endif
  end subroutine prolongation_high_order_for_face_ghost
  !======================================================================

  subroutine correct_face_ghost_for_fine_block(iBlock, nVar, Field_VG)

    use BATL_tree, ONLY: DiLevelNei_IIIB, iNode_B
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         nI, nJ, nK, nDim
    integer, intent(in):: iBlock, nVar
    real, intent(inout):: Field_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    logical:: IsCorrected_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real,parameter::  Distance_II(4,4) = reshape((/&
         -2,-1,4,5, &
         -3,-2,3,4, &
         -4,-3,2,3, &
         -5,-4,1,2 &
         /), (/4,4/))
    real::  Coef_II(6,4) 

    real, parameter::  Coef1_II(6,3) = reshape(([&
         -5./21, 9./14, 4.0,  -45./7, 27./7,  -5./6, &
         -5./42, 2./7,  4.0,  -40./7, 45./14, -2./3, &
         -2./63, 1./14, 8./3, -20./7, 10./7,  -5./18 &
         ]),&
         ([6,3]))
    real,parameter::  Distance1_II(4,3) = reshape((/&         
         -3,-2,3,4, &
         -4,-3,2,3, &
         -5,-4,1,2 &
         /), (/4,3/))

    real, parameter::  Coef2_II(6,2) = reshape(([&
         -1./7,  2./3, 5./3,   -2.0,  1.0,  -4./21, &
         -2./35, 1./5,  2.0,   -9./5, 4./5, -1./7   &
         ]),&
         ([6,2]))
    real,parameter::  Distance2_II(4,2) = reshape((/&         
         -2, -1, 2, 3, &
         -3, -2, 1, 2  &
         /), (/4,2/))

    real, parameter::  Coef_I(6) = ([1./20, -3./10, 3./4, 3./4, -3./10, 1./20])
    real, parameter::  Distance_I(4) = (/-2, -1, 1, 2/)

    integer:: iDir, jDir, kDir, iDir1, jDir1, kDir1, nDir1, iDir2, jDir2, kDir2
    integer:: jBegin, jEnd, iBegin, iEnd, kBegin, kEnd
    integer:: Di, Dj, Dk, i, j, k, i1, j1, k1
    integer:: Di1, Dj1, Dk1
    integer:: DiLevel, DiLevel1, Count, iVar, iCount1
    real:: Orig, CellValue_I(6)
    integer:: FineNeiIndex_II(3,3), IEdge_, IEdge1_, nEdge, iStage
    logical:: IsFineNei_I(3)
    integer, parameter:: Edge1_ = 1, Edge2_ = 2, Corner_ = 3
    character(len=*), parameter:: NameSub = 'correct_face_ghost_for_fine_block'
    !----------------------------------------------------------------------

    if(nI>7 .and. nJ>7 .and. (nDim ==2 .or. nK>7)) then
       ! If there are enough information, use 6 points interpolation, which
       ! is symmetric. 
       Coef_II = reshape(([&
            5./21, -15./14, 12./7, 3./7,  -3./7,  5./42, &
            5./14, -10./7,  12./7, 8./7, -15./14, 2./7, &
            2./7,  -15./14, 8./7,  12./7, -10./7, 5./14,&
            5./42, -3./7,   3./7,  12./7, -15./14,5./21 &
            ]),&
            ([6,4]))
    else
       ! Do not have enough infomation to do symmetric interpolation. 
       ! 5 points interpolation. 
       Coef_II = reshape(([&
            5./14,  -10./7,  2.0,    1./7,   -1./14,  0.0, &
            9./14,  -16./7, 12./5,  16./35,  -3./14,  0.0, &
            9./14,  -15./7,  2.0,    6./7,   -5./14,  0.0, &
            5./14,  -8./7,   1.0,    8./7,   -5./14,  0.0  &
            ]),&
            ([6,4]))
    endif

    if(nK == 1) then

       ! For example: 

       !--------- 
       !|   |    |
       !|8__|_9__|
       !|   |    |
       !|6  | 7  |
       ! -------------------
       !|        | 4  |  5 |
       !|   1    |____|____|
       !|        |    |    |
       !|        | 2  | 3  |
       !--------------------

       ! For block 4, the ghost cells -2 <= i <= 0, nJ - 3 <= j <= nJ
       ! will be corredted. 
       ! For block 7, the ghost cells nI-3 <= i <= nI, -2 <= j <= 0 
       ! will be corredted. 

       kDir = 0
       do iDir = -1, 1; do jDir = -1, 1
          if(iDir*jDir /=0) CYCLE
          DiLevel = DiLevelNei_IIIB(iDir,jDir,kDir,iBlock)
          if(DiLevel /= 1) CYCLE

          ! DiLevel = 1
          if(jDir == 0) then ! Resolution change in x-dir
             if(iDir ==1) then
                iBegin = nI+1; iEnd = nI+3; Di = 1
             else
                iBegin = 0; iEnd = -2; Di = -1
             endif

             do jDir1 = -1, 1, 2
                DiLevel1 = DiLevelNei_IIIB(iDir,jDir1,kDir,iBlock)
                if(DiLevel1 /=0) CYCLE

                if(jDir1 == -1) then
                   jBegin = 1; jEnd = 4; Dj = 1
                else
                   jBegin = nJ; jEnd = nJ - 3; Dj = -1
                endif

                k = 1
                
                do i = iBegin, iEnd, Di
                   do iVar = 1, nVar
                      CellValue_I(1) = Field_VG(iVar,i,jBegin-3*Dj,k)
                      CellValue_I(2) = Field_VG(iVar,i,jBegin-2*Dj,k)
                      CellValue_I(3) = Field_VG(iVar,i,jBegin-Dj,k)
                      CellValue_I(4) = Field_VG(iVar,i,jEnd+Dj,k)
                      CellValue_I(5) = Field_VG(iVar,i,jEnd+2*Dj,k)       
                      CellValue_I(6) = Field_VG(iVar,i,jEnd+3*Dj,k)

                      Count = 1
                      do j = jBegin, jEnd, Dj
                         Orig = sum(CellValue_I*Coef_II(:,Count))
                         Field_VG(iVar,i,j,k) = limit_interpolation(Orig,&
                              CellValue_I(2:5), Distance_II(:,Count))
                         Count = Count + 1                      
                      enddo
                   enddo ! iVar
                enddo ! i

             enddo ! jDir1
          elseif(iDir == 0) then
             if(jDir == 1) then
                jBegin = nJ+1; jEnd = nJ+3; Dj = 1
             else
                jBegin = 0; jEnd = -2; Dj = -1
             endif

             do iDir1 = -1, 1, 2
                DiLevel1 = DiLevelNei_IIIB(iDir1,jDir,kDir,iBlock)
                if(DiLevel1 /=0) CYCLE

                if(iDir1 == -1) then
                   iBegin = 1; iEnd = 4; Di = 1
                else
                   iBegin = nI; iEnd = nI - 3; Di = -1
                endif

                k = 1
                do j = jBegin, jEnd, Dj
                   do iVar = 1, nVar
                      CellValue_I(1) = Field_VG(iVar,iBegin-3*Di,j,k)
                      CellValue_I(2) = Field_VG(iVar,iBegin-2*Di,j,k)
                      CellValue_I(3) = Field_VG(iVar,iBegin-  Di,j,k)
                      CellValue_I(4) = Field_VG(iVar,iEnd  +  Di,j,k)
                      CellValue_I(5) = Field_VG(iVar,iEnd  +2*Di,j,k)
                      CellValue_I(6) = Field_VG(iVar,iEnd  +3*Di,j,k)

                      Count = 1
                      do i = iBegin, iEnd, Di
                         Orig = sum(CellValue_I*Coef_II(:,Count))
                         Field_VG(iVar,i,j,k) = limit_interpolation(Orig,&
                              CellValue_I(2:5),Distance_II(:,Count))
                         Count = Count + 1
                      enddo
                   enddo ! ivar
                enddo ! j
             enddo ! idir1
          endif

       enddo; enddo
    elseif(nK >1) then
       ! Example: 
       ! ___________________________________________
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |         7          |         8          |                        
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |____________________|____________________|
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |          5         |         6          |                        
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |____________________|____________________|
       !               TOP LAYER

       ! ___________________________________________
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |        3           |          4         |                        
       ! |                    |                    |                        
       ! |                    |                    |                        
       ! |                    |                    |        y                
       ! |____________________|____________________|        | 
       ! |                    |         |          |        |                
       ! |                    |   11    |   12     |        -----> x                     
       ! |                    |         |          |                        
       ! |         1          |_________|__________|                        
       ! |                    |         |          |                        
       ! |                    |   9     |   10     |                        
       ! |                    |         |          |                        
       ! |____________________|_________|__________|
       !                 BOTTOM LAYER
       ! 9-12 are the top layer of previous coares block. 

       ! To calculate the left face ghost cells of block 11, need to 
       ! consider the refinement of block 3, 5, and 7. 

       ! case 1: 3, 5, 7 are coarse, nothing need to do.
       ! case 2: only 3 or 5 is refined. 7 can be refined or not. 
       !         Similar with the 2D case above. 
       ! case 3: both 3 and 5 are refined. 7 can be refined or not. 
       !         Some face ghost cells can be interpolated in y 
       !         direction or z direction, and these cells use the 
       !         average of both direction. 
       ! case 4: only 7 is refined. The ghost cells: -2 .le. i .le. 0 .and. 
       !         nJ-3 .le. j .le. nJ .and. nK .le. k .le. nK . are not 5th
       !         order. They are calculated with 4th order accurate 
       !         interpolation in prolongation_high_order_for_face_ghost. 
       !         These cells may passed to block 9 and other block as 
       !         edge/corner ghost cells. But themselves, as the face ghost
       !         of block 11, they will be overwritten with 5th order value
       !         when iSendStage is 3. 

       do kDir = -1, 1; do jDir = -1, 1; do iDir = -1, 1
          ! Loop through 6 faces. 

          ! Check whether block 1 is coarse
          if(abs(kDir) + abs(jDir) + abs(iDir) /=1) CYCLE
          if(DiLevelNei_IIIB(iDir,jDir,kDir,iBlock) /=1) CYCLE

          ! Check whether block 3 and block 5 are refined. 
          IsFineNei_I = .false. 
          do kDir1 = -1, 1; do jDir1 = -1, 1; do iDir1 = -1, 1

             ! Only Check the neighbour blocks close to this coarse 
             ! face neighbour.
             if(iDir /=0 .and. iDir1 /= iDir) CYCLE
             if(jDir /=0 .and. jDir1 /= jDir) CYCLE
             if(kDir /=0 .and. kDir1 /= kDir) CYCLE
             if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlock) /= 0) CYCLE

             ! Only check edge neighbours. 
             nDir1 = abs(iDir1) + abs(jDir1) + abs(kDir1)
             if(nDir1 == 3) CYCLE

             ! At most 2 edges can be refined block. 
             if(.not. IsFineNei_I(Edge1_)) then                
                IEdge_ = Edge1_
             else
                IEdge_ = Edge2_
             endif
             if(IsFineNei_I(IEdge_)) call CON_stop(NameSub//': error1')
             IsFineNei_I(IEdge_) = .true.

             FineNeiIndex_II(IEdge_,1) = iDir1
             FineNeiIndex_II(IEdge_,2) = jDir1
             FineNeiIndex_II(IEdge_,3) = kDir1                
          enddo; enddo; enddo

          ! Check whether block 7 is refined. 
          if(.not. IsFineNei_I(Edge1_) .and. .not. IsFineNei_I(Edge2_)) then
             ! Only one corner neighbour can be refined now. 
             do kDir1 = -1, 1; do jDir1 = -1, 1; do iDir1 = -1, 1
                if(iDir /=0 .and. iDir1 /= iDir) CYCLE
                if(jDir /=0 .and. jDir1 /= jDir) CYCLE
                if(kDir /=0 .and. kDir1 /= kDir) CYCLE

                if(DiLevelNei_IIIB(iDir1,jDir1,kDir1,iBlock) /= 0) CYCLE
                nDir1 = abs(iDir1) + abs(jDir1) + abs(kDir1)

                if(nDir1 /= 3) CYCLE
                if(IsFineNei_I(Corner_)) call CON_stop(NameSub//': error2')
                IsFineNei_I(Corner_) = .true.
                FineNeiIndex_II(Corner_,1) = iDir1
                FineNeiIndex_II(Corner_,2) = jDir1
                FineNeiIndex_II(Corner_,3) = kDir1             
             enddo; enddo; enddo
          endif

          nEdge = 1 ! case 2
          if(IsFineNei_I(Edge1_) .and. IsFineNei_I(Edge2_)) nEdge = 2 !case 3
          IsCorrected_G = .false. 
          do iStage = 1, nEdge
             ! If nEdge is 2:
             ! iStage 1: calculate these cells can only interpolated in 
             !           one direction
             ! iStage 2: for these cells can be interpolated in two 
             !           directions, use the average. 
             do IEdge_ = 1, 2
                if(.not. IsFineNei_I(IEdge_)) CYCLE

                iDir1 = FineNeiIndex_II(IEdge_,1)
                jDir1 = FineNeiIndex_II(IEdge_,2)
                kDir1 = FineNeiIndex_II(IEdge_,3)

                if(nEdge == 2) then
                   IEdge1_ = mod(IEdge_,2) + 1
                   iDir2 = FineNeiIndex_II(IEdge1_,1)
                   jDir2 = FineNeiIndex_II(IEdge1_,2)
                   kDir2 = FineNeiIndex_II(IEdge1_,3)
                endif

                if(iDir /=0) then
                   if(iDir == 1) then
                      iBegin = nI + 1;  iEnd = nI + 3; Di = 1
                   elseif(iDir == -1) then
                      iBegin = 0; iEnd = -2; Di = -1
                   endif

                   if(jDir1 /=0) then
                      if(nEdge == 1) then
                         kBegin = 1; kEnd = nK; Dk = 1
                      else
                         if(iStage == 1) then
                            if(kDir2 ==1) then
                               kBegin = 1; kEnd = nK-4; Dk = 1
                            elseif(kDir2 == -1) then
                               kBegin = nK; kEnd = 5; Dk = -1
                            endif
                         else ! iStage == 2                               
                            if(kDir2 ==1) then
                               kBegin = nK - 3; kEnd = nK; Dk = 1
                            elseif(kDir2 == -1) then
                               kBegin = 4; kEnd = 1; Dk = -1
                            endif

                         endif

                      endif

                      if(jDir1 == 1) then
                         jBegin = nJ; jEnd = jBegin; Dj = -1
                      elseif(jDir1 == -1) then
                         jBegin = 1; jEnd = jBegin; Dj = 1
                      endif

                      Di1 = 0; Dj1 = Dj; Dk1 = 0

                   elseif(kDir1 /=0) then
                      if(nEdge == 1) then
                         jBegin = 1; jEnd = nJ; Dj = 1
                      else
                         if(iStage == 1) then
                            if(jDir2 ==1) then
                               jBegin = 1; jEnd = nJ-4; Dj = 1
                            elseif(jDir2 == -1) then
                               jBegin = nJ; jEnd = 5; Dj = -1
                            endif
                         else ! iStage = 2
                            if(jDir2 ==1) then
                               jBegin = nJ - 3; jEnd = nJ; Dj = 1
                            elseif(jDir2 == -1) then
                               jBegin = 4; jEnd = 1; Dj = -1
                            endif
                         endif
                      endif

                      if(kDir1 == 1) then
                         kBegin = nK; kEnd = kBegin; Dk = -1
                      elseif(kDir1 == -1) then
                         kBegin = 1; kEnd = kBegin; Dk = 1
                      endif

                      Di1 = 0; Dj1 = 0; Dk1 = Dk
                   endif

                elseif(jDir /=0) then
                   if(jDir == 1) then
                      jBegin = nJ + 1; jEnd = nJ + 3; Dj = 1
                   elseif(jDir == -1) then
                      jBegin = 0; jEnd = -2; Dj = -1
                   endif

                   if(iDir1 /=0) then
                      if(nEdge == 1) then
                         kBegin = 1; kEnd = nK; Dk = 1
                      else
                         if(iStage == 1) then
                            if(kDir2 ==1) then
                               kBegin = 1; kEnd = nK-4; Dk = 1
                            elseif(kDir2 == -1) then
                               kBegin = nK; kEnd = 5; Dk = -1
                            endif
                         else ! iStage == 2
                            if(kDir2 ==1) then
                               kBegin = nK - 3; kEnd = nK; Dk = 1
                            elseif(kDir2 == -1) then
                               kBegin = 4; kEnd = 1; Dk = -1
                            endif
                         endif ! iStage

                      endif

                      if(iDir1 == 1) then
                         iBegin = nI; iEnd = iBegin; Di = -1
                      elseif(iDir1 == -1) then
                         iBegin = 1; iEnd = iBegin; Di = 1
                      endif

                      Di1 = Di; Dj1 = 0; Dk1 = 0

                   elseif(kDir1 /=0) then
                      if(nEdge == 1) then
                         iBegin = 1; iEnd = nI; Di = 1
                      else
                         if(iStage == 1) then
                            if(iDir2 ==1) then
                               iBegin = 1; iEnd = nI-4; Di = 1
                            elseif(iDir2 == -1) then
                               iBegin = nI; iEnd = 5; Di = -1
                            endif
                         else
                            if(iDir2 ==1) then
                               iBegin = nI-3; iEnd = nI; Di = 1
                            elseif(iDir2 == -1) then
                               iBegin = 4; iEnd = 1; Di = -1
                            endif
                         endif
                      endif

                      if(kDir1 == 1) then
                         kBegin = nK; kEnd = kBegin; Dk = -1
                      elseif(kDir1 == -1) then
                         kBegin = 1; kEnd = kBegin; Dk = 1
                      endif

                      Di1 = 0; Dj1 = 0; Dk1 = Dk
                   endif

                elseif(kDir /=0) then
                   if(kDir == 1) then 
                      kBegin = nK + 1; kEnd = nK + 3; Dk = 1
                   elseif(kDir == -1) then
                      kBegin = 0; kEnd = -2; Dk = -1
                   endif


                   if(iDir1 /=0) then
                      if(nEdge == 1) then
                         jBegin = 1; jEnd = nJ; Dj = 1
                      else
                         if(iStage == 1) then
                            if(jDir2 ==1) then
                               jBegin = 1; jEnd = nJ-4; Dj = 1
                            elseif(jDir2 == -1) then
                               jBegin = nJ; jEnd = 5; Dj = -1
                            endif
                         else
                            if(jDir2 ==1) then
                               jBegin = nJ - 3; jEnd = nJ; Dj = 1
                            elseif(jDir2 == -1) then
                               jBegin = 4; jEnd = 1; Dj = -1
                            endif
                         endif

                      endif

                      if(iDir1 == 1) then
                         iBegin = nI; iEnd = iBegin; Di = -1
                      elseif(iDir1 == -1) then
                         iBegin = 1; iEnd = iBegin; Di = 1
                      endif

                      Di1 = Di; Dj1 = 0; Dk1 = 0

                   elseif(jDir1 /=0) then
                      if(nEdge == 1) then
                         iBegin = 1; iEnd = nI; Di = 1
                      else
                         if(iStage == 1) then
                            if(iDir2 ==1) then
                               iBegin = 1; iEnd = nI-4; Di = 1
                            elseif(iDir2 == -1) then
                               iBegin = nI; iEnd = 5; Di = -1
                            endif
                         else ! iStage == 2
                            if(iDir2 ==1) then
                               iBegin = nI-3; iEnd = nI; Di = 1
                            elseif(iDir2 == -1) then
                               iBegin = 4; iEnd = 1; Di = -1
                            endif
                         endif

                      endif

                      if(jDir1 == 1) then
                         jBegin = nJ; jEnd = jBegin; Dj = -1
                      elseif(jDir1 == -1) then
                         jBegin = 1; jEnd = jBegin; Dj = 1
                      endif
                      Di1 = 0; Dj1 = Dj; Dk1 = 0
                   endif
                endif

                do k=kBegin,kEnd,Dk; do j=jBegin,jEnd,Dj; do i=iBegin,iEnd,Di
                   do iVar = 1, nVar
                      CellValue_I(1) = Field_VG(iVar,&
                           i-3*Di1,j-3*Dj1,k-3*Dk1)
                      CellValue_I(2) = Field_VG(iVar,&
                           i-2*Di1,j-2*Dj1,k-2*Dk1)
                      CellValue_I(3) = Field_VG(iVar,&
                           i-  Di1,j-  Dj1,k-  Dk1)
                      CellValue_I(4) = Field_VG(iVar,&
                           i+4*Di1,j+4*Dj1,k+4*Dk1)
                      CellValue_I(5) = Field_VG(iVar,&
                           i+5*Di1,j+5*Dj1,k+5*Dk1)
                      CellValue_I(6) = Field_VG(iVar,&
                           i+6*Di1,j+6*Dj1,k+6*Dk1)

                      do Count = 1, 4
                         Orig = sum(CellValue_I*Coef_II(:,Count))

                         if(.not. IsCorrected_G(i+(Count-1)*Di1,&
                              j+(Count-1)*Dj1,k+(Count-1)*Dk1)) then
                            Field_VG(iVar,i+(Count-1)*Di1,&
                                 j+(Count-1)*Dj1,k+(Count-1)*Dk1)&
                                 = limit_interpolation(Orig,&
                                 CellValue_I(2:5),Distance_II(:,Count))

                            IsCorrected_G(i+(Count-1)*Di1,&
                                 j+(Count-1)*Dj1,k+(Count-1)*Dk1) = .true.
                         else
                            Field_VG(iVar,i+(Count-1)*Di1,&
                                 j+(Count-1)*Dj1,k+(Count-1)*Dk1)&
                                 = 0.5*(Field_VG(iVar,i+(Count-1)*Di1,&
                                 j+(Count-1)*Dj1,k+(Count-1)*Dk1)  &
                                 + & 
                                 limit_interpolation(Orig,&
                                 CellValue_I(2:5),Distance_II(:,Count)))
                         endif
                      enddo
                   enddo

                enddo; enddo; enddo

             enddo ! IEdge_
          enddo ! iStage

       enddo; enddo; enddo
    endif

  end subroutine correct_face_ghost_for_fine_block
  !======================================================================

  real function prolongation_high_order_amr(Cell_III)
    ! Calc 5th order refined cell for AMR.
    use BATL_size, ONLY: kRatio
    real, intent(in):: Cell_III(5,5,5)
    real:: Cell_I(5), Cell_II(5,5)
    integer:: i, j, k
    real:: Temp, Distance_I(4)=(/-7,-3,1,5/)
    real:: CellLimit_I(4)
    !----------------------------------------------------------------------

    if(kRatio == 2) then
       do j = 1, 5; do i = 1, 5
          Cell_II(i,j) = interpolate_in_coarse_blk_1d(&
               Cell_III(i,j,:), .true.)
       enddo; enddo
       do i = 1, 4
          CellLimit_I(i) = Cell_III(i,i,i)
       enddo
    else
       k = 1
       Cell_II = Cell_III(:,:,k)
       do i = 1, 4
          CellLimit_I(i) = Cell_III(i,i,1)
       enddo
    endif

    do i = 1, 5 ! Eliminate j dimension
       Cell_I(i) = interpolate_in_coarse_blk_1d(Cell_II(i,:), .true.)
    enddo
    Temp = interpolate_in_coarse_blk_1d(Cell_I, .true.)
    prolongation_high_order_amr = limit_interpolation(Temp,CellLimit_I,&
         Distance_I)
  end function prolongation_high_order_amr
  !======================================================================

  real function restriction_high_order_amr(Cell_III)
    ! Calc 6th order coarsened cell for AMR. 
    use BATL_size, ONLY: kRatio
    real, intent(in):: Cell_III(6,6,6)
    real:: Cell_I(6), Cell_II(6,6)
    integer:: i, j, k

    if(kRatio == 2) then 
       do j = 1, 6; do i = 1, 6
          Cell_II(i,j) = calc_face_value(Cell_III(i,j,:), .true.)
       enddo; enddo
    else
       k = 1
       Cell_II = Cell_III(:,:,k) 
    endif

    do i = 1, 6
       Cell_I(i) = calc_face_value(Cell_II(i,:), .true.)
    enddo
    restriction_high_order_amr = calc_face_value(Cell_I, .true.)
  end function restriction_high_order_amr

end module BATL_high_order