/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "SolverLinear.h"
#include "Problem.h"
#include "Assembler.h"
#include "Mesh.h"
#include "DofObject.h"
#include "NodalData.h"
#include "InternalVars.h"
#include "BodyLoad.h"
#include "Utilities.h"
#include "SensitivitySolver.h"
#include <iostream>



PETScLinearKSPSolver::PETScLinearKSPSolver()
	: _phase(0)
{
}

PETScLinearKSPSolver::~PETScLinearKSPSolver()
{
}




// Solves the global system using PETSc.
// First assembles the matricies and vectors using assemble_linear_elasticity
// Then solves the partitioned system using the PETSc KSP context
PetscErrorCode PETScLinearKSPSolver::solve()
{
	//	Variables
	//
	//	K		-	Matrix that defines the linear system
	//	ksp		-	KSP context
	//	b, u	-	RHS, exact solution vectors

	KSP ksp;
	// AO ao;   				// Application Ordering object
	Vec temp1, temp2, Res_f;
	PetscErrorCode ierr;
	Utilities::timer timer;

	// Print to screen
	PIGFEMPrint("\nSolving the Linear Problem...\n\n\n");

	// Things I'll actually never need in the linear Problem. There are no ISVs in a linear problem
	std::vector<std::vector<std::vector<double> > > update_ISVs;
	_prob->get_internal_vars()->preallocate_internal_vars_object(update_ISVs);

	// Preallocate all memory for matricies and vectors
	ierr = VecDuplicate(_F_ext.f, &Res_f);CHKERRQ(ierr);
	ierr = VecDuplicate(_F_ext.f, &temp1);CHKERRQ(ierr);		  // Copy correct structure
	ierr = VecDuplicate(_F_ext.p, &temp2);CHKERRQ(ierr);
	
	// Actually assemble the linear system
	Assembler* assem = _prob->get_assembler();
	timer.start();
	ierr = assem->assemble_new_load_step(_U.p, _F_ext, 1.0);CHKERRQ(ierr);
	if (getNonzeroStart())
		{ierr = VecAXPY(_U.p, 1.0, _Up_initial);CHKERRQ(ierr);}
	ierr = assem->assemble(_K, _P_int, 1.0, update_ISVs, _prob->get_solution(), true, true);CHKERRQ(ierr);
	_prob->_assemble_time = timer.time_elapsed();

	// Create the KSP context and set operators
	ierr = initKSP(ksp, true);CHKERRQ(ierr);

	// Start the solve process
	timer.start();
	ierr = VecCopy(_F_ext.f, Res_f);CHKERRQ(ierr);
	ierr = MatMult(_K.fp, _U.p, temp1);CHKERRQ(ierr);    // Will store the product of Kfp*up in temp1
	ierr = VecAXPY(Res_f, -1.0, temp1);CHKERRQ(ierr);  // Stored Res_f = Ff - Kfp*Up. This will be te rhs of our solve
	ierr = VecAXPY(Res_f, -1.0, _P_int.f);CHKERRQ(ierr); // Store Res_f = Ff - Kfp*Up - Pf 
	
	// Actually Solve the system!!!!
	ierr = KSPSolve(ksp, Res_f, _U.f);CHKERRQ(ierr);  // Solve the system
	
	// Compute the reaction forces
	ierr = MatMult(_K.pf, _U.f, temp2);CHKERRQ(ierr);
	ierr = MatMult(_K.pp, _U.p, _F_ext.p);CHKERRQ(ierr);
	ierr = VecAXPY(_F_ext.p, 1.0, temp2);CHKERRQ(ierr);
	ierr = VecAXPY(_F_ext.p, 1.0, _P_int.p);CHKERRQ(ierr); // Stores Fp = Pp + Kfp*Uf + Kpp*Up
	_prob->_solve_time = timer.time_elapsed();
	
	// Stores the solution in the Problem object in an easier to access structure
	store_solution(_U, _prob->get_solution());
	store_solution(_F_ext, _prob->get_external_load());

	if (_prob->sensitivity())
	{
		timer.start();
		_prob->get_sensitivity_solver()->solveProblem();
		_prob->_sensitivity_time += timer.time_elapsed();
	}

	Output(_phase, 1.0);
	
	// Clean up local memory
	ierr = VecDestroy(&temp1);CHKERRQ(ierr);
	ierr = VecDestroy(&temp2);CHKERRQ(ierr);
	ierr = VecDestroy(&Res_f);CHKERRQ(ierr);
	ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

	_prob->get_solved() = true;
	_phase++;
	return ierr;
}

