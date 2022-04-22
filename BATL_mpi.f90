!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_mpi

  use ModUtilities, ONLY: CON_stop
  use ModMpi

  implicit none

  SAVE

  private ! except

  public:: init_mpi       ! initialize this module
  public:: clean_mpi      ! finalize MPI
  public:: barrier_mpi    ! use an MPI barrier to synchronize processors

  integer, public:: iComm       ! MPI communicator for the group of processors
  integer, public:: nProc = 1   ! number of processors in this group
  integer, public:: iProc = 0   ! processor rank from 0 to nProc-1
  integer, public:: nThread = 1 ! number of threads per MPI process
  integer, public:: iThread = 0 ! thread rank from 0 to nThread-1
  integer, public:: nGpuDev=-1 ! total number of GPUs on this node
  integer, public:: iGpuDev=-1 ! the GPU the current MPI process is using
  !$acc declare create(iComm, nProc, iProc, nThread, iThread)

contains
  !============================================================================
#ifdef _OPENACC
  !  subroutine acc_error_handler(sMsg) bind(C)
  subroutine acc_error_handler() bind(C)
    ! use iso_c_binding, only: C_CHAR, C_NULL_CHAR
    use ModUtilities, only: CON_Stop
    ! character(kind=C_CHAR,len=1024), intent(in) :: sMsg
    ! integer i

    ! do i=1,len(sMsg)
    !    if (sMsg(i:i)==C_NULL_CHAR) exit
    ! end do

    !--------------------------------------------------------------------------
    write (*,*) 'DEBUG: Fortran OpenACC error handler is calling CON_Stop()'
    ! call CON_Stop("OpenACC error: " // sMsg(1:i-1))
    call CON_Stop("OpenACC error!")
    stop "UNREACHABLE CODE" ! hint for the compiler
  end subroutine acc_error_handler
  !============================================================================
#endif
#ifdef _OPENACC
  subroutine set_acc_error_handler()
    use iso_c_binding, only: C_FUNPTR, c_funloc
    interface
       subroutine acc_set_error_routine(FuncPtr) BIND(C)
         use iso_c_binding, only: C_FUNPTR
         type(C_FUNPTR), intent(in), value :: FuncPtr
       end subroutine acc_set_error_routine
    end interface
    type(C_FUNPTR) :: FuncPtr
    !--------------------------------------------------------------------------
    FuncPtr = c_funloc(acc_error_handler)
    call acc_set_error_routine(FuncPtr)
  end subroutine set_acc_error_handler
  !============================================================================
#endif
#ifdef _OPENACC
  subroutine init_gpu(iComm, iProc)

    use openacc

    integer, intent(in) :: iComm, iProc
    integer :: iLocalComm, iLocalProc, nLocalProc, iError

    character(len=*), parameter:: NameSub = 'init_gpu'
    !--------------------------------------------------------------------------

    ! Get the node-local (shared-memory capable) communicator
    call MPI_Comm_Split_Type(iComm, MPI_COMM_TYPE_SHARED, iProc, &
         MPI_INFO_NULL, iLocalComm, iError)

    ! Determine the number of local processes and the local rank
    call MPI_Comm_Size(iLocalComm, nLocalProc, iError)
    call MPI_Comm_Rank(iLocalComm, iLocalProc, iError)

    ! Determine the number of GPUs
    nGpuDev = acc_get_num_devices(ACC_DEVICE_NVIDIA)

    if (nGpuDev <= 0) call CON_stop(NameSub//': No GPUs detected on the node')

    iGpuDev = iLocalProc
    if (nLocalProc > nGpuDev) then ! we have more processes than GPUs
       if (iLocalProc==0) write (*,*) NameSub, " WARNING:", &
            ' iProc, nLocalProc > nGpu=', iProc, nLocalProc, nGpuDev
       iGpuDev = mod(iLocalProc, nGpuDev)
    end if

    ! set the device number we will be operating on
    !$acc set device_num(iGpuDev)
    !$acc init device_num(iGpuDev)

    call MPI_Comm_free(iLocalComm, iError)

    call set_acc_error_handler()

  end subroutine init_gpu
  !============================================================================
#endif
  subroutine init_mpi(iCommIn)

    ! Initialize iComm, nProc and iProc. If iCommIn is not present, set it
    ! to MPI_COMM_WORLD and also call MPI_init

    use omp_lib

    integer, optional, intent(in):: iCommIn
    integer :: iError

    ! OpenMP support levels
    integer, parameter :: lSupportRequired = MPI_THREAD_SINGLE
    integer            :: lSupportProvided = MPI_THREAD_SINGLE
    !--------------------------------------------------------------------------
    if(.not.present(iCommIn))then
       iComm = MPI_COMM_WORLD
       call MPI_init_thread(lSupportRequired, lSupportProvided, iError)
       ! Check the threading support level
       if(lSupportProvided < lSupportRequired) then
          if(iProc == 0) write(*,*) &
               "Warning:  This MPI implementation provides ",   &
               "insufficient threading support. Switching to pure MPI..."
          !$ call omp_set_num_threads(1)
       end if
       !$ nThread = omp_get_max_threads()
    else
       iComm = iCommIn
    end if
    call MPI_comm_rank(iComm, iProc, iError)
    call MPI_comm_size(iComm, nProc, iError)
#ifdef _OPENACC
    call init_gpu(iComm, iProc)
#endif
    !$acc update device(iComm, iProc, nProc)

  end subroutine init_mpi
  !============================================================================

  subroutine clean_mpi

    ! This should only be called if the whole application is finished

    integer :: iError
    !--------------------------------------------------------------------------
    call MPI_finalize(iError)

  end subroutine clean_mpi
  !============================================================================
  subroutine barrier_mpi

    use ModUtilities, ONLY: flush_unit
    use ModIoUnit,    ONLY: STDOUT_

    integer:: iError
    !--------------------------------------------------------------------------
    call flush_unit(STDOUT_)
    call MPI_barrier(iComm, iError)

  end subroutine barrier_mpi
  !============================================================================

end module BATL_mpi
!==============================================================================
