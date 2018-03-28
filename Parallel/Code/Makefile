SHELL:=/bin/bash
path = $(shell pwd)/../
dir_name= test_run
new_dir = $(path)$(dir_name)
num_procs = 2
all_exec = compile-main compile-MB compile-rdf
to_compile = compile-pbc compile-lj compile-mom compile-init_pos compile-init_vel compile-thermo compile-verlet compile-pos compile-KE compile-read_data compile-print_data compile-distort compile-send_array
.PHONY: clean backup run-main run-MB run-rdf variables module-load

compile-all: module-load $(all_exec) 

module-load:
	$(shell module load openmpi/1.6.5_gcc-4.8.3)

compile-main: EIA_dynamics_main.f90 $(to_compile)
	mpif90 -g -o 	dynamics_mpi \
			pbc_module.o \
			distort_module.o \
		       	lj_module.o \
			moment_module.o \
			position_module.o \
			velocities_module.o \
			andersen_therm_module.o \
			vel_verlet_module.o \
			print_positions_module.o \
			kinetic_energy_module.o \
			read_data_module.o \
			print_data_module.o \
			EIA_dynamics_main.f90
			
compile-MB: check_MB.f90
	mpif90 -g -o 	check_MB_mpi \
			check_MB.f90

compile-rdf: rdf.f90
	mpif90 -g -o 	rdf_mpi \
			rdf.f90 \
			pbc_module.f90

compile-pbc: pbc_module.f90
	mpif90 -g -c 	pbc_module.f90

compile-lj: lj_module.f90
	mpif90 -g -c 	lj_module.f90

compile-mom: moment_module.f90
	mpif90 -g -c 	moment_module.f90

compile-init_pos: position_module.f90
	mpif90 -g -c 	position_module.f90

compile-init_vel: velocities_module.f90
	mpif90 -g -c 	velocities_module.f90

compile-thermo: andersen_therm_module.f90
	mpif90 -g -c 	andersen_therm_module.f90

compile-verlet: vel_verlet_module.f90
	mpif90 -g -c 	vel_verlet_module.f90

compile-pos: print_positions_module.f90
	mpif90 -g -c 	print_positions_module.f90

compile-KE: kinetic_energy_module.f90
	mpif90 -g -c 	kinetic_energy_module.f90

compile-read_data: read_data_module.f90
	mpif90 -g -c 	read_data_module.f90

compile-print_data: print_data_module.f90
	mpif90 -g -c 	print_data_module.f90

compile-distort: distort_module.f90
	mpif90 -g -c 	distort_module.f90

compile-send_array: send_rec_module.f90
	mpif90 -g -c 	send_rec_module.f90

compile-rows_per_proc: rows_per_proc_module.f90
	mpif90 -g -c 	rows_per_proc_module.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f dynamics_mpi rdf_mpi check_MB_mpi
	rm -f *.out traj.xyz

backup:
	@if [ -d 'backup' ]; then \
		rm -f ./backup/*; \
	else \
		mkdir ./backup; \
	fi 
	cp *.f90 ./backup

run-main: compile-main
	@if [ ! -d $(new_dir) ]; then \
		mkdir $(new_dir); \
		cp ./input.dat $(new_dir); \
		cp ./dynamics_mpi $(new_dir); \
		cd $(new_dir); \
		mpirun -np $(num_procs) dynamics_mpi input.dat; \
	else \
		echo 'THE DIRECTORY ALREADY EXISTS'; \
		echo 'PLEASE CHOOSE AN OTHER NAME' ; \
	fi 

run-MB: compile-MB
	@if [ -d $(new_dir) ]; then \
		cp ./check_MB $(new_dir); \
		cd $(new_dir); \
		mpirun -np $(num_procs) check_MB_mpi velocity.out parameters.out
	else \
		echo 'THE DIRECTORY DOES NOT EXIST'; \
		echo 'PLEASE TARGET A VALID ONE'; \
	fi

run-rdf: compile-rdf
	@if [ -d $(new_dir) ]; then \
		cp ./rdf $(new_dir); \
		cd $(new_dir); \
		mpirun -np $(num_procs) rdf_mpi traj.xyz parameters.out; \
	else \
		echo 'THE DIRECTORY DOES NOT EXIST'; \
		echo 'PLEASE TARGET A VALID ONE'; \
	fi

variables: 
	@echo '-------------------------------------------------'
	@echo 
	@echo "'dir_name' default is "$(dir_name)
	@echo 
	@echo 
	@echo "'num_procs' default is "$(num_procs)
	@echo 
	@echo '-------------------------------------------------'




























