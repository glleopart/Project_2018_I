module rows_per_proc_module
use mpi
implicit none
contains
subroutine rows_per_proc(nPart,myFirstPart,myLastPart)
implicit none
integer, intent(in)                       :: nPart
integer, intent(out)                      :: myFirstPart, myLastPart
integer                                   :: numPart, rank, numProcs, ierror

call mpi_comm_rank(mpi_comm_world, rank, ierror)
call mpi_comm_size(mpi_comm_world, numProcs, ierror)

numPart = nPart/numProcs
myFirstPart = rank*numPart + 1
myLastPart = myFirstPart + numPart - 1
if (rank == numProcs - 1) myLastPart = nPart

end subroutine rows_per_proc

end module rows_per_proc_module
