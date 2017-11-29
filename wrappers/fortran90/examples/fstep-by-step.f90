!
! @file fstep-by-step.f90
!
! Fortran 90 example using a matrix read with the spm driver.
!
! @copyright 2017-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
!                      Univ. Bordeaux. All rights reserved.
!
! @version 6.0.0
! @author Mathieu Faverge
! @date 2017-01-01
!
program fsimple
  use iso_c_binding
  use pastix_enums
  use spmf
  use pastixf
  ! use mpi_f08
  implicit none

  type(pastix_data_t),        pointer                    :: pastix_data
  type(pastix_order_t),       pointer                    :: order
  type(pastix_spm_t),         target                     :: spm
  type(pastix_spm_t),         pointer                    :: spm2
  integer(kind=pastix_int_t), target                     :: iparm(iparm_size)
  real(kind=c_double),        target                     :: dparm(dparm_size)
  integer(c_int)                                         :: info
  integer(kind=pastix_int_t)                             :: nrhs
  real(kind=c_double), dimension(:), allocatable, target :: x0, x, b
  type(c_ptr)                                            :: x0_ptr, x_ptr, b_ptr
  integer(kind=pastix_int_t), dimension(:), pointer      :: permtab
  integer                                                :: i, j, nfact, nsolv

  nfact = 2
  nsolv = 3

  !
  ! Initialize the problem
  !   1- The matrix
  call spmReadDriver( PastixDriverLaplacian, "d:10:10:10:2.", spm, 0, info )

  call spmCheckAndCorrect( spm, spm2 )
  if (.not. c_associated(c_loc(spm), c_loc(spm2))) then
     call spmExit( spm )
     spm = spm2
  end if

  call spmPrintInfo( spm, c_null_ptr )

  !   2- The right hand side
  nrhs = 1
  allocate(x0(spm%n))
  allocate(x(spm%n))
  allocate(b(spm%n))
  x0_ptr = c_loc(x0)
  x_ptr  = c_loc(x)
  b_ptr  = c_loc(b)

  ! 1- Initialize the parameters and the solver
  call pastixInitParam( iparm, dparm )
  call pastixInit( pastix_data, 0, iparm, dparm )

  ! 2- Perform ordering, symbolic factorization, and analyze steps
  call pastix_subtask_order( pastix_data, spm, null(), info )
  call pastix_subtask_symbfact( pastix_data, info )
  call pastix_subtask_reordering( pastix_data, info )
  call pastix_subtask_blend( pastix_data, info )

  ! If needed, get the generated ordering
  call pastixOrderGet( pastix_data, order )

  ! Convert the permtab to Fortran array
  call c_f_pointer( order%permtab, permtab, [order%vertnbr] )
  print *, permtab(1:10)

  ! 3- Factorize nfact times the matrix
  do i=0,nfact
     ! Perform the numerical factorization
     call pastix_subtask_spm2bcsc( pastix_data, spm, info )
     call pastix_subtask_bcsc2ctab( pastix_data, info )
     call pastix_subtask_sopalin( pastix_data, info )

     ! Perform nsolv solve steps
     do j=0,nsolv

        call spmGenRHS( PastixRhsRndX, nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, info )
        x = b

        ! 4- Solve the problem
        call pastix_task_solve( pastix_data, nrhs, x_ptr, spm%n, info )

        ! 5- Refine the solution
        call pastix_task_refine( pastix_data, x_ptr, nrhs, b_ptr, info )

        ! Check the solution
        call spmCheckAxb( nrhs, spm, x0_ptr, spm%n, b_ptr, spm%n, x_ptr, spm%n, info )

     end do
  end do

  ! 6- Destroy the C data structure
  call pastixFinalize( pastix_data )

  call spmExit( spm )
  deallocate(x0)
  deallocate(x)
  deallocate(b)

end program fsimple