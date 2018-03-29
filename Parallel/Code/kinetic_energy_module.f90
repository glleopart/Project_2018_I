module kinetic_energy_module
use mpi
implicit none
contains
subroutine kinetic_energy(vel, KETotal, Tinst, myFirstPart, myLastPart, nPart, rank)
implicit none
integer, intent(in)                             :: nPart, myFirstPart, myLastPart, rank
real(8), dimension(nPart,3), intent(in)         :: vel
real(8), intent(out)                            :: KETotal, Tinst
real(8), dimension(3)                           :: vec
real(8)                                         :: modV, KE
integer                                         :: i, ierror
integer, parameter                              :: rMaster = 0

KE = 0.0D0
ketotal = 0.0D0
do i = myFirstPart, myLastPart, 1
        vec(:) = vel(i,:)
        modV = dsqrt(dot_product(vec, vec))
        KE = KE + modV**2.
end do
if (rank.eq.rMaster) then
print *, "Sóc el màster i he acabat el loop de kinetic"
endif

print *, "Hola sóc el", rank, "Aquesta és la KE", KE, 'ketotal=', KETotal

call mpi_reduce(KE, KETotal, 1, MPI_real8, MPI_SUM, rMaster, mpi_comm_world, ierror)
print *, "He fet el reduce en kinetic energy", rank
if (rank == rMaster) then
        KETotal = KETotal*0.5
        Tinst = 2.0*KETotal/(3.0*float(nPart))
        print *, "Hola sóc el màster he fet el Tinst."
end if


end subroutine kinetic_energy
end module kinetic_energy_module
