path = $(shell pwd)/../
dir_name=test_run
new_dir = $(path)$(dir_name)
all_exec = compile-main compile-MB compile-rdf
.PHONY: clean backup run-main run-MB run-rdf

compile-all: $(all_exec)

compile-main: EIA_dynamics_main.f90
	gfortran 	pbc_module.f90 \
			distort_module.f90 \
		       	lj_module.f90 \
			moment_module.f90 \
			position_module.f90 \
			velocities_module.f90 \
			andersen_therm_module.f90 \
			vel_verlet_module.f90 \
			print_positions_module.f90 \
			kinetic_energy_module.f90 \
			read_data_module.f90 \
			print_data_module.f90 \
			EIA_dynamics_main.f90 \
			-o dynamics

compile-MB: check_MB.f90
	gfortran 	check_MB.f90 \
			-o check_MB

compile-rdf: rdf.f90
	gfortran	pbc_module.f90 \
			rdf.f90 \
			-o rdf

compile-pbc: pbc_module.f90
	gfortran -c 	pbc_module.f90

compile-lj: lj_module.f90
	gfortran -c 	lj_module.f90

compile-mom: moment_module.f90
	gfortran -c 	moment_module.f90

compile-init_pos: position_module.f90
	gfortran -c 	position_module.f90

compile-init_vel: velocities_module.f90
	gfortran -c 	velocities_module.f90

compile-thermo: andersen_therm_module.f90
	gfortran -c 	andersen_therm_module.f90

compile-verlet: vel_verlet_module.f90
	gfortran -c 	vel_verlet_module.f90

compile-pos: print_positions_module.f90
	gfortran -c 	print_positions_module.f90

compile-KE: kinetic_energy_module.f90
	gfortran -c 	kinetic_energy_module.f90

compile-read_data: read_data_module.f90
	gfortran -c 	read_data_module.f90

compile-print_data: print_data_module.f90
	gfortran -c 	print_data_module.f90

compile-distort: distort_module.f90
	gfortran -c 	distort_module.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f dynamics rdf check_MB
	rm -f *.out traj.xyz

backup:
	if [ -d 'backup' ]; then \
		rm -f ./backup/*; \
	else \
		mkdir ./backup; \
	fi 
	cp *.f90 ./backup

run-main: compile-main
	if [ ! -d $(new_dir) ]; then \
		mkdir $(new_dir); \
		cp ./input.dat $(new_dir); \
		cp ./dynamics $(new_dir); \
		cd $(new_dir); \
		nohup ./dynamics input.dat > nohup.err; \
	else \
		echo 'THE DIRECTORY ALREADY EXISTS'; \
		echo 'PLEASE CHOOSE AN OTHER NAME' ; \
	fi 

run-MB: compile-MB
	if [ -d $(new_dir) ]; then \
		cp ./check_MB $(new_dir); \
		cd $(new_dir); \
		./check_MB velocity.out parameters.out; \
	else \
		echo 'THE DIRECTORY DOES NOT EXIST'; \
		echo 'PLEASE TARGET A VALID ONE'; \
	fi

run-rdf: compile-rdf
	if [ -d $(new_dir) ]; then \
		cp ./rdf $(new_dir); \
		cd $(new_dir); \
		./rdf traj.xyz parameters.out; \
	else \
		echo 'THE DIRECTORY DOES NOT EXIST'; \
		echo 'PLEASE TARGET A VALID ONE'; \
	fi
plot-energy: 
	if [ -d $(new_dir) ]; then \
                cd $(new_dir); \
                gnuplot energy_all.gp; \
        else \
                echo 'THE DIRECTORY DOES NOT EXIST'; \
                echo 'PLEASE TARGET A VALID ONE'; \
        fi

plot-momentum: 
	if [ -d $(new_dir) ]; then \
                cd $(new_dir); \
                gnuplot momentum_all.gp; \
        else \
                echo 'THE DIRECTORY DOES NOT EXIST'; \
                echo 'PLEASE TARGET A VALID ONE'; \
        fi

plot-temperature: 
	if [ -d $(new_dir) ]; then \
                cd $(new_dir); \
                gnuplot temperature_all.gp; \
        else \
                echo 'THE DIRECTORY DOES NOT EXIST'; \
                echo 'PLEASE TARGET A VALID ONE'; \
        fi



help: 
	@echo '-------------------------------------------------'
	@echo Type "make" to compile everything
	@echo Type "make run-main" to run the MD simulation in ../test_dir and 2 threads by default.
	@echo '        'To change default values: "make run-main num_procs=# dir_name='your/dir/'"
	@echo Type "make run-rdf" to run the RDF computation program.
	@echo Type "make run-MB" to perform the checking the velocity distribution.
	@echo '        '\(This two commands must be executed AFTER running MD.
	@echo '        'But the Makefile already have this feature into account.
	@echo '        'Default options for "run-main" apply to "make run-rdf" and "make run-MB" too.\)
	@echo Type "make clean" to remove *.o, *.mod, *.xyz, *.out and the exec files.
	@echo Type "make backup" to create a backup directory with the *.f90 files.
	@echo Type "make help" to show the aviable options.
	@echo 
	@echo Type make plot-energy dir_name=toplot ,in order to plot Potential Energy, press enter to watch kinetic energy, press enter again 
	@echo to watch mechanical energy . 
	@echo 
	@echo 
	@echo Type : make plot-momentum dir_name=toplot , in order to watch momentum behaviour for i=4,12,2 particles per dimension  
	@echo 
	@echo Type : make plot-temperature dir_name=toplot , in order to watch temperature behaviour for i=4,12,2 particles per dimension  
	@echo 
	@echo "'dir_name' default is "$(dir_name)
	@echo 
	@echo 
	@echo "'num_procs' default is "$(num_procs)
	@echo 
	@echo '-------------------------------------------------'

























