module send_recv_module
use mpi
implicit none
contains

subroutine send_recv_array(dats, myFirstRow, myLastRow, rank, nRow, status, full_dats)
implicit none
integer, intent(in)                                     :: nRow, rank, status
integer                                                 :: ierror
integer, intent(in)                                     :: myFirstRow, myLastRow
real(8), dimension(:,:), intent(inout)                  :: dats
real(8), optional, dimension(:,:), intent(out)          :: full_dats
integer, parameter                                      :: rMaster = 0
integer                                                 :: i, nDats, numProcs
integer                                                 :: numRows, procFirstRow, procLastRow

call mpi_comm_size(mpi_comm_world, numProcs, ierror)
numRows = nRow/numProcs
if (rank == rMaster) then
        if (.not.present(full_dats)) then
                print*, "ERROR WHILE CALLING send_recv_array()"
                print*, "'full_dats' not given for the MASTER THREAD"
                print*, "EXITTING PROGRAM"
                call exit()
        end if

        do i = 0, numProcs - 1, 1
                procFirstRow = i*numRows + 1
                if (i /= numProcs - 1) procLastRow = procFirstRow + numRows - 1
                if (i == numProcs - 1) procLastRow = nRow
                if (i == rMaster) then
                        full_dats(procFirstRow:procLastRow,:) = dats(:,:)
                        cycle
                end if
                nDats = 3*(procLastRow - procFirstRow + 1)
                call mpi_recv(&
                & full_dats(procFirstRow:procLastRow,:),        &
                & nDats,                                        &
                & mpi_real8,                                    &
                & i,                                            &
                & mpi_any_tag,                                  &
                & mpi_comm_world,                               &
                & status,                                       &
                & ierror                                        &
                &)

                
        end do
else if (rank /= rMaster) then
        nDats = 3*(myLastRow - myFirstRow + 1)
        call mpi_send(&
        & dats,                 &
        & nDats,                &
        & mpi_real8,            &
        & rMaster,              &
        & 2001,                 &
        & mpi_comm_world,       &
        & ierror                &
        &)

end if
end subroutine send_recv_array
end module send_recv_module

