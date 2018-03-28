!::::::::::::::::::::::::
! MADE BY GENIS LLEOPART
! PARALLEL VERSION BY ADRIÀ VICENS
!::::::::::::::::::::::::

! Contains the subroutine momentum, which calculate
! the total intertial momentum of our system

! Variables IN:
!   vel --> Velocities of the particles | DIMENSION(N,3), REAL(8)
!   myFirstPart --> Index of the first particle | integer
!   myLastPart --> Index of the first particle | integer
!   partialMomentum --> partial momentum of each processor | DIMENSION(3), REAL(8)

! Variables OUT:
!   totalMomentum --> total momentum | DIMENSION(3), REAL(8)
module moment_module
use mpi
implicit none
contains

subroutine momentum(myFirstPart, myLastPart, vel, totalMomentum)
implicit none
real(8), dimension(:,:), intent(in) :: vel
real(8), dimension(3), intent(out) :: totalMomentum
integer, intent(in) :: myFirstPart, myLastPart
real(8), dimension(3) :: partialMomentum
integer :: ierror, i
integer, parameter :: rMaster = 0

partialMomentum = 0.0d0

totalMomentum = 0.0d0
partialMomentum = sum(vel(myFirstPart:myLastPart, :), dim = 1)

call mpi_barrier( mpi_comm_world, ierror)
call mpi_reduce(partialMomentum, totalMomentum, 3, mpi_real8, mpi_sum, rMaster, mpi_comm_world, ierror)

end subroutine momentum


end module moment_module
