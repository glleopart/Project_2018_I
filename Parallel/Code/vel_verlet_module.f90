module vel_verlet_module
use mpi
use lj_module
use pbc_module
contains
    subroutine vel_verlet(time,dt, pos, vel, nPart, eps, sig, boxSize,cutOff,V,F_aux,myFirstPart,myLastPart)
    integer,intent(in)         :: nPart,rank,ierror,myFirstPart,myLastPart
    integer :: i,j,k
    integer,parameter :: rmaster = 0
    real(8),dimension(myFirstPart:myLastPart,3),intent(inout)        :: vel,pos,F 
    real(8),intent(in),dimension(myFirstPart:myLastPart,3)           :: F_aux 
    real(8),dimension(3)         :: vec
    real(8),intent(inout)        :: time
    real(8),intent(out)          :: V
    real(8),intent(in)           :: dt, eps,sig,boxsize,cutOff
    real(8),allocatable          :: full_dats
        if ((time == 0).and.(rank == rmaster)) then 
            call LJ_pot(nPart, myFirstPart, myLastPart, pos, eps, sig, boxSize, cutOff,F_aux, V)
        endif
        if (time == 0) then
            do i = myFirstPart , myLastPart, 1
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
        call send_recv_array(vel(myFirstPart:myLastPart,:),myFirstPart,myLastPart,rank,nRow,status,vel) 
        call send_recv_array(pos(myFirstPart:myLastPart,:),myFirstPart,myLastPart,rank,nRow,status,pos) 
        
            
        












    end subroutine
end module 
