module vel_verlet_module
use pbc_module
use lj_module
contains
subroutine velocity_verlet(time, dt, pos, vel, nPart, eps, sig, boxSize, cutOff, V, F, forcesTotT)
implicit none
integer, intent(in)                             :: nPart
real(8), intent(in)                             :: dt, eps, sig, boxSize, cutOff
real(8), intent(inout)                          :: time
real(8), intent(out)                            :: V
real(8), dimension(nPart,3), intent(inout)      :: pos, vel, F
real(4), intent(inout)                             :: forcesTotT
real(8), dimension(nPart,3)                     :: F_aux
integer                                         :: i, j, k
real(8), dimension(3)                           :: vec
real(4)                                         :: forcesIniT, forcesFinT

!Starting the simulation the forces are obtained calling LJ_pot
if (time == 0) then
         call cpu_time(forcesIniT)
         call LJ_pot(nPart, pos, eps, sig, boxSize, cutOff, F_aux, V)
         call cpu_time(forcesFinT)
         forcesTotT = forcesTotT + (forcesFinT - forcesIniT)
end if
!Once the simulation is running the forces enter as an input value F(t)
!Also F(t+dt) is the output of the subroutine
if (time /= 0) F_aux(:,:) = F(:,:)
time = time + dt

!calculation of Pos(t+dt)
do i = 1, nPart, 1
        vec(:) = pos(i,:) + vel(i,:)*dt + F_aux(i,:)*dt**3./2.
        call pbc(vec, boxSize)
        pos(i,:) = vec(:)
end do
!Calculation of forces at time t+dt
call cpu_time(forcesIniT)
call LJ_pot(nPart, pos, eps, sig, boxSize, cutOff, F, V)
call cpu_time(forcesFinT)
forcesTotT = forcesTotT + (forcesFinT - forcesIniT)

!calculation of vel(t+dt)
do i = 1, nPart, 1
        vel(i,:) = vel(i,:) + (F_aux(i,:) + F(i,:))*dt/2.
end do
end subroutine velocity_verlet

end module vel_verlet_module
