module vel_verlet_module
use mpi
use lj_module
use pbc_module
use send_recv_module
implicit none
contains

subroutine vel_verlet(time, dt, pos, vel, nPart, eps, sig, boxSize, cutOff, V, F, myFirstPart, myLastPart&
                &, rank, status)
implicit none
integer, intent(in)                                             :: nPart, rank, myFirstPart, myLastPart, status
real(8), dimension(nPart,3),intent(inout)                       :: vel, pos 
real(8), dimension(myFirstPart:myLastPart,3), intent(inout)     :: F 
real(8), intent(inout)                                          :: time
real(8), intent(out)                                            :: V
real(8), intent(in)                                             :: dt, eps, sig, boxSize, cutOff
real(8), dimension(myFirstPart:myLastPart,3)                    :: F_aux
real(8), dimension(3)                                           :: vec
integer                                                         :: i, j, k, ierror
integer, parameter                                              :: rMaster = 0

if ((time == 0).and.(rank == rmaster)) then 
        call LJ_pot(nPart, myFirstPart, myLastPart, pos, eps, sig, boxSize, cutOff, F_aux, V)
endif
if (time == 0) then
        do i = myFirstPart, myLastPart, 1
                vec(:) = pos(i,:) + vel(i,:)*dt + F_aux(i,:)*dt**3./2.
                call pbc(vec,boxsize)
                pos(i,:) = vec(:)
        end do
endif
if ((time /= 0).and.(rank == rmaster)) then
        F_aux(:,:) = F(:,:)
        call LJ_pot(nPart, myFirstPart, myLastPart, pos, eps, sig, boxSize, cutOff,F, V)
endif
if (time/=0) then 
        do i = myFirstPart,myLastPart, 1
                vel(i,:) = vel(i,:)*dt +(F_aux(i,:)+ F(i,:))*dt/2.
        end do
endif
call send_recv_array(vel(myFirstPart:myLastPart,:),myFirstPart,myLastPart,rank,nPart,status,vel) 
call send_recv_array(pos(myFirstPart:myLastPart,:),myFirstPart,myLastPart,rank,nPart,status,pos) 

end subroutine
end module 
