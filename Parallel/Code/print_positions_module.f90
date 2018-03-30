module print_positions_module
implicit none
contains

subroutine print_positions(un, nPart, pos, time)
implicit none
integer, intent(in)                             :: un, nPart
real(8), dimension(nPart,3), intent(in)         :: pos
real(8), intent(in)                             :: time
integer                                         :: i
character(25)                                   :: formOut, formOut2

formOut  = '(A3, 3F10.3)'
formOut2 = '(A18, F10.3, A7)'
write(un,'(I6)') nPart
write(un,formOut2) 'Simulation Time = ', time,'seconds'
do i = 1, nPart, 1
        write(un,formOut) 'C', pos(i,:)
end do

end subroutine print_positions

end module print_positions_module
