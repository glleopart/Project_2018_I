!------------------------
!MADE BY ORIOL PIQUÉ
!-------------------------

module lj_module
use mpi
use pbc_module
implicit none
contains

!This module contains the subroutine that calculates the resulting forces that appear because of the LJ potential interaction between the particles.


!Variables(in):

!Number of particles (nPart)
!Positions array (pos)
!Epsilon (eps)
!Sigma (sig)
!Size of the simulation box (boxSize)
!Cutoff distance (cutOff)


!Variables(out):

!Forces array (F)
!Potential (V)


subroutine LJ_pot(nPart, myFirstPart, myLastPart, pos, eps, sig, boxSize, cutOff, F, V)
implicit none
integer, intent(in)                                         :: nPart, myFirstPart, myLastPart
real(8), dimension(nPart, 3), intent(in)                    :: pos
real(8), intent(in)                                         :: eps, sig, boxSize, cutOff
real(8), dimension(myFirstPart:myLastPart,3), intent(out)   :: F
real(8), intent(out)                                        :: V
real(8)                                                     :: vPartial
real(8), dimension(3)                                       :: dist
real(8)                                                     :: rij, dV
integer                                                     :: i, j, k, ierror
integer, parameter                                          :: rMaster = 0

vPartial = 0.
V = 0.
F(:,:) = 0.
do i = myFirstPart, myLastPart, 1; do j = i + 1, nPart, 1
        dist(:) = pos(i,:) - pos(j,:)
        call pbc(dist, boxSize)
        rij = dsqrt(dot_product(dist,dist))
        if (rij < cutOff) then
                dist(:) = dist(:)/rij
                vPartial  = vPartial + 4.*eps*((sig/rij)**12. - (sig/rij)**6.)
                dV = 4*eps*(12.*sig**12./rij**13. - 6.*sig**6./rij**7)
                F(i,:) = F(i,:) + dV*dist(:)
        end if
end do; end do

call mpi_reduce(vPartial, V, 1, mpi_real8, mpi_sum, rMaster, mpi_comm_world, ierror)

end subroutine LJ_pot
end module lj_module
