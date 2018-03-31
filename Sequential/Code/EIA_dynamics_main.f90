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
real(4)                                 :: printIniT,  printFinT,  printTotT
real(4)                                 :: verletIniT, verletFinT, verletTotT, verletTime
real(4)                                 :: thermoIniT, thermoFinT, thermoTotT
real(4)                                 :: forcesIniT, forcesFinT, forcesTotT
real(4)                                 :: momentIniT, momentFinT, momentTotT
real(4)                                 :: kinetiIniT, kinetiFinT, kinetiTotT
    
call cpu_time(start)
call get_command_argument(1, fName, status=fStat)
if (fStat /= 0) then
        print*, 'Any file given ---> Exitting program'
        call exit()
end if 
un = 100
open(unit=un, file=trim(fName), status='old') 
    
! Llegir les dades del fitxer de entrada i alocatar la variable que contindrà
! la posició de totes les partícules.
call readData(un, dt, boxSize, cutOff, nPartDim, T, eps, sig, nSteps, density, seed)
nPart = nPartDim**3
allocate(pos(nPart,3), F(nPart,3), vel(nPart,3), total_momentum(3))
if (mod(float(nPart),8.) /= 0.) then
        print*, 'Total number of particles must be'
        print*, 'divisible by 8 to fit a SC cell'
        print*, 'Exitting Program'
        call exit()
end if
    
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
write(dataUn,*), '# time, V, KE, KE + V, mom, Tins'
write(paramUn,*) nSteps/100 - 1, nPart
write(paramUn,*) boxSize

call print_positions(initUn, nPart, pos, time)
! Començem la MD inicialitzant el temps a 0 i fent el numero de pasos
! especificats
time = 0.0D0; trjCount = 0; thermCount = 0
call print_positions(trajUn, nPart, pos, time)
call print_positions(velUn,  nPart, vel, time)

! Posem temps a zero
verletTotT = 0.0
thermoTotT = 0.0
printTotT  = 0.0
forcesTotT = 0.0
momentTotT = 0.0
kinetiTotT = 0.0


do i = 1, nSteps, 1
        call cpu_time(verletIniT)
        call velocity_verlet(time, dt, pos, vel, nPart, eps, sig, boxSize, cutOff, V, F,forcesTotT)
        call cpu_time(verletFinT)
        verletTotT = verletTotT + (verletFinT - verletIniT)
        if (mod(trjCount,100) == 0) then
                call cpu_time(kinetiIniT)
                call kinetic_energy(vel, KE, Tinst, nPart)
                call cpu_time(kinetiFinT)
                kinetiTotT = kinetiTotT + (kinetiFinT - kinetiIniT)
                
                call cpu_time(momentIniT)
                call momentum(nPart, vel, total_momentum)
                call cpu_time(momentFinT)
                momentTotT = momentTotT + (momentFinT - momentIniT)
                
                call cpu_time(printIniT)
                call print_positions(trajUn, nPart, pos, time) 
                call print_positions(velUn,  nPart, vel, time)
                call print_data(time, V, KE, Tinst, total_momentum, dataUn)
                call cpu_time(printFinT)
                printTotT = printTotT + (printFinT - printIniT)

                trjCount = 0
        end if
        if (mod(thermCount,10) == 0) then
                call cpu_time(thermoIniT)
                call andersen_thermo(dt, T, nPart, i, vel)
                call cpu_time(thermoFinT)
                thermoTotT = thermoTotT + (thermoFinT - thermoIniT)
                thermCount = 0
        end if
        trjCount = trjCount + 1
        thermCount = thermCount + 1
end do
call print_positions(finUn, nPart, pos, time)


call cpu_time(finish)
verletTime = verletTotT - forcesTotT
        write(dataUn,*) '# VERLET TOTAL_CPU_TIME  : ', verletTime
        write(dataUn,*) '# THERMO TOTAL_CPU_TIME  : ', thermoTotT
        write(dataUn,*) '# FORCES TOTAL_CPU_TIME  : ', forcesTotT
        write(dataUn,*) '# MOMENT TOTAL_CPU_TIME  : ', momentTotT
        write(dataUn,*) '# KINETIC TOTAL_CPU_TIME : ', kinetiTotT
        write(dataUn,*) '# PRINTING TIME : ', printTotT
        write(dataUn,*) "# TOTAL CPU TIME         : ", finish - start
        write(dataUn,*)
        write(dataUn,*)
        write(dataUn,*) '# RELATIVE TIME TO THE TOTAL CPU_TIME (IN %):'
        write(dataUn,*) '# VERLET   : ', 100*verletTime/(finish-start)
        write(dataUn,*) '# THERMO   : ', 100*thermoTotT/(finish-start)
        write(dataUn,*) '# FORCES   : ', 100*forcesTotT/(finish-start)
        write(dataUn,*) '# MOMENT   : ', 100*momentTotT/(finish-start)
        write(dataUn,*) '# KINETIC  : ', 100*kinetiTotT/(finish-start)
        write(dataUn,*) '# PRINTING : ', 100*printTotT/(finish-start)
close(un); close(initUn); close(finUn); close(trajUn); close(dataUn)
contains

end program dynamics
    
