# PIGFEM

========================================
INSTALLATION INSTRUCTIONS FOR PIGFEM
========================================

Developed by David Brandyberry (brandyb2@illinois.edu)
Last Update 8/8/2016

Inside this zip folder you should find 12 things:

	1. petsc-lite-3.7.3.tar.gz
	2. metis-5.1.0.tar.gz
	3. gsl-2.1.tar.gz
	4. include folder (PCIGFEM)
	5. src folder (PCIGFEM)
	6. obj folder (PCIFEM)
	7. Tests folder (PCIGFEM)
	8. Input folder (PCIGFEM)
	9. Output folder (PCIGFEM)
	10. Makefile (PCIGFEM)
	11. setup script
	12. This file

If you don't have all of these things something is wrong, please email the developer.

NOTE: If you don't have an internet connections, these instructions will not work as external packages are downloaded

Installation instructions:

	1. Run the setup script with the command "./setup"
		- This script will take ~30 minutes to run
		- This will decompress and compile PETSc, METIS, and GSL (GNU Scientific Library)
		- This default command build the optimized versions of all of these packages
		- To build the debug versions run the script "./setup debug"
		- This file also downloads a version of MPI that works well with PETSc. all of the executables are stored in the PETSc folder. the path to this folder is written to you ~/.bashrc file so can use this version of MPI from any directory
	
	2. After the script has finished running, close all open terminals
		- This is to reload your ~/.bashrc file

	3. Run the command "which mpiexec"
		- If the end of the output of this command is not "/petsc-3.7.3/arch-linux2-c-opt/bin" then you may have difficulties running PCIGFEM in parallel. Please read the running instruction 6 in detail.

	4. You should be good to run the code!




Running Instructions:

	1. You should have an exectuable in your current directory called Main. PCIGFEM is called from this file using command line parameters to control input
	
	2. Currently Main is only structured to read and parse input files containing "x,y,r" data that describes circular fibers. To run a file matching this descriptions continue with these instructions. Otherwise you're out of luck until the developer passes Quals.

	3. Place your input file in the Input folder. There are a few files already in the Input folder which follow the necessary input syntax

	4. Return to the folder with the Main executable

	5. To run the code in serial simply run the command "./Main your_file_name_here"

	6. How to run the code in parallel depends on your success in Installation step 3
		- If the output of Installation step 3 was as described, you should be able to run in parllel with the command "mpiexec -n number_of_procs_to_use ./Main your_file_name_here"
		- If the output of Installation step 3 was as not described, you should be able to run in parllel with the command "petsc-3.7.3/arch-linux2-c-opt/bin/mpiexec -n number_of_procs_to_use ./Main your_file_name_here"
		- NOTE: there are still some major issues running larger problems on multiple cores that will cause the program to fail. I'm still not sure what the issue is but I'm working on it. To get around this you should be able to just use a smaller number of cores and your solution will just take longer. Sorry about that...

	7. You should see output on the screen describing the nonlinear solver convergence history. Additionally, VTK output files will be output to the Output folder. The nonlinear history data will show up in Output/NL_history.##.vtk where ## is the load step number. Additionally, the final result should be output to a file called "Output/IGFEM_out_<#_of_fibers>.vtk" These files can be visualized with Paraview





Modification Instructions:
	If you wish to modify the problem being run (material parameters, mesh size, domain size, etc...) feel free to modify the file Tests/Main.cpp

	1. Open Tests/Main.cpp

	2. Parameters that you would generally wish to change are limited to the #define statements towards the begining of the file. Feel free to modify these.

	3. After modifications, return to the main directoy and run the instruction "make"
		- This recompiles Main.cpp and links it with everything else

	4. You should be good to run again

	5. If you break anything and can't fix it please email the developer
