module kinetic_energy_module
use mpi
contains
subroutine kinetic_energy(vel, KETotal, Tinst, myFirstPart, myLastPart, nPart)
implicit none
integer, intent(in)                             :: nPart, myFirstPart, myLastPart
real(8), dimension(nPart,3), intent(in)         :: vel
real(8), intent(out)                            :: KETotal, Tinst
real(8), dimension(3)                           :: vec
real(8)                                         :: modV, KE
integer                                         :: i, ierror
integer, parameter                              :: rMaster = 0

KE = 0.0D0
do i = myFirstPart, myLastPart, 1
        vec(:) = vel(i,:)
        modV = dsqrt(dot_product(vec, vec))
        KE = KE + modV**2.
end do


call mpi_reduce(KE, KETotal, 1, MPI_real8, MPI_SUM, rMaster, mpi_comm_world, ierror)
KETotal = KETotal*0.5
Tinst = 2.0*KETotal/(3.0*float(nPart))


end subroutine kinetic_energy
end module kinetic_energy_module
