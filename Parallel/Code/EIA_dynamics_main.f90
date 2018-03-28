program dynamics
!use the necessary modules
use moment_module
use positions_module
use velocities_module
use andersen_therm_module
use distort_module
use vel_verlet_module
use lj_module
use pbc_module
use read_data_module
use print_positions_module
use kinetic_energy_module
use print_data_module
use send_rec_module
use rows_per_proc_module

implicit none
character(25)                           :: fName, fff
integer                                 :: i, j, k, l, m, n
integer                                 :: fStat, un
real(4)                                 :: start, finish
real(8)                                 :: dt, boxSize, cutOff, T, density
integer                                 :: nPartDim, nPart, nSteps
real(8), allocatable, dimension(:,:)    :: pos, F, vel
real(8)                                 :: V, eps, sig, time, KE, Tinst
real(8), allocatable, dimension(:)      :: total_momentum
integer                                 :: seed, trjCount, thermCount
integer                                 :: initUn, finUn, trajUn, dataUn, velUn, paramUn
integer                                 :: ierror, rank, numProcs, status, numParts
integer                                 :: myFirstPart, myLastPart
integer, parameter                      :: rMaster = 0

call mpi_init(ierror)
call mpi_comm_rank(mpi_comm_world, rank, ierror)
call mpi_comm_size(mpi_comm_world, numProcs, ierror)

if (rank == rMaster) then 
        call cpu_time(start)
        call get_command_argument(1, fName, status=fStat)
        if (fStat /= 0) then
                print*, 'Any file given ---> Exitting program'
                call mpi_finalize(ierror)
                call exit()
        end if 
        un = 100
        open(unit=un, file=trim(fName), status='old') 
    
        ! Llegir les dades del fitxer de entrada i alocatar la variable que contindrà
        ! la posició de totes les partícules.
        call readData(un, dt, boxSize, cutOff, nPartDim, T, eps, sig, nSteps, density, seed)
        nPart = nPartDim**3
        if (mod(float(nPart),8.) /= 0.) then
                print*, 'Total number of particles must be'
                print*, 'divisible by 8 to fit a SC cell'
                print*, 'Exitting Program'
                call mpi_finalize(ierror)
                call exit()
        end if
end if
call mpi_bcast(dt, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(boxSize, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(cutOff, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(T, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(eps, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(sig, 1, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(nPart, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
call mpi_bcast(nSteps, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
call mpi_bcast(seed, 1, mpi_integer, rMaster, mpi_comm_world, ierror)
seed = seed + rank

call rows_per_proc(nPart, myFirstPart, myLastPart)
allocate(pos(nPart,3), F(nPart,3), vel(nPart,3), total_momentum(3))

if (rank == rMaster) then
        ! Generar la posició inicial de totes les partícules conforme a un cristall 
        ! d'una cel·la Simple Cubic. I distorsiona les posicions inicials per tal de que
        ! no tingui tanta energia reticular.
        call SC_init_conditions(nPart, pos, boxSize)
        call distort_geometry(nPart, pos, boxSize, seed)

        ! Genera una distribució de velocitats gaussiana amb l'algoritme de Box_Muller
        call IN_velocities(nPart, T, seed, vel)
    
        ! Obre els arxius on s'imprimiran els resultats de la simulació
        initUn = 101; finUn = 102; trajUn = 103; dataUn = 104; velUn = 105; paramUn = 106
        open(unit=initUn, file='initial.out')           ! Coord. Inicials
        open(unit=finUn , file='final.out')             ! Coord. Finals
        open(unit=trajUn, file='traj.xyz')              ! Trajectoria
        open(unit=dataUn, file='data.out')              ! T, Ken, V, temps...
        open(unit=velUn,  file='velocity.out')          ! T, Ken, V, temps...
        open(unit=paramUn,file='parameters.out')        ! Parameters needed for statistics
                                                ! (MB and RDF)
        write(paramUn,*) nSteps/100 - 1, nPart
        write(paramUn,*) boxSize
end if

call mpi_bcast(pos, 3*nPart, mpi_real8, rMaster, mpi_comm_world, ierror)
call mpi_bcast(vel, 3*nPart, mpi_real8, rMaster, mpi_comm_world, ierror)

if (rank == rMaster) then
        call print_positions(initUn, nPart, pos, time)
        call print_positions(trajUn, nPart, pos, 0.0D0)
        call print_positions(velUn,  nPart, vel, 0.0D0)
end if

time = 0.0D0; trjCount = 0; thermCount = 0
do i = 1, nSteps, 1
        call velocity_verlet(time, dt, pos, vel, nPart, eps, sig, boxSize, cutOff, V, F)
        if (mod(trjCount,100) == 0) then
                call kinetic_energy(vel, KE, Tinst, nPart)
                call momentum(nPart, vel, total_momentum)
                if (rank == rMaster) then
                        call print_positions(trajUn, nPart, pos, time)
                        call print_positions(velUn,  nPart, vel, time)
                        call print_data(time, V, KE, Tinst, total_momentum, dataUn)
                end if
                trjCount = 0
        end if
        if (mod(thermCount,10) == 0) then
                call andersen_thermo(dt, T, nPart, i, vel)
                thermCount = 0
        end if
        trjCount = trjCount + 1
        thermCount = thermCount + 1
end do
if (rank == rMaster) then
        call print_positions(finUn, nPart, pos, time)
        call cpu_time(finish)
        write(dataUn,*) "# CPU TIME: ", finish - start
        close(un); close(initUn); close(finUn); close(trajUn); close(dataUn)
end if
contains

end program dynamics
    
