/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "Problem.h"
#include "Assembler.h"
#include "BoundaryObject.h"
#include "DofObject.h"
#include "InternalVars.h"
#include "Mesh.h"
#include "NodalData.h"
#include "Solver.h"
#include "Writer.h"
#include "Writer_ASCII_VTK.h" // Default
#include "Writer_ASCII_VTK_Nodal.h"
#include "Writer_Binary_VTK.h"
#include "WriterSensitivity.h"
#include "WriterCohesiveFailureLengths.h"
#include "WriterPath.h"
#include "WriterHomogenization.h"
#include "WriterCohesiveOpenings.h"
#include "RefinerGeometricalInclusion.h"
#include "RefinerGeometricalInclusionTraversal.h"
#include "elem.h"
#include "material.h"
#include "BodyLoad.h"
#include "Utilities.h"
#include "SensitivityDirect.h"
#include "SensitivityAdjoint.h"
#include "SensitivityMaterialParameter.h"
#include "SensitivityShapeParameter.h"
#include "SensitivityFunction.h"
#include "Options.h"
#include "Inclusion.h"
#include "Integrator.h"
#include <iostream>


Problem::Problem()
	: _mesh(NULL), _boundary(NULL), _dofs(NULL), _assembler(NULL), _solver(NULL), _sensitivity_solver(NULL), _internal_vars(NULL), _solution(NULL), _external_load(NULL),
	_writer(NULL), _coh_writer(NULL), _sensitivity_writer(NULL), _coh_failure_writer(NULL), _path_writer(NULL), _coh_openings_writer(NULL), _homogenization_writer(NULL),
	_solved(false), _method(PENALTY), _subscale(false), _init(false), _setup(false), _sensitivity(false),
	_assemble_time(0.0), _solve_time(0.0), _misc_solve_time(0.0), _solver_write_time(0.0), _sensitivity_time(0.0), 
	_IGFEM_nodal_detection_time(0.0), _IGFEM_element_detection_time(0.0), _IGFEM_enrichments_time(0.0)
{
	// One of the first things I'll always is read in the options database
	_options = new Options;
	_boundary = new BoundaryObject;
	_dofs = new DofObject;
}

Problem::~Problem()
{
	delete _dofs;
	delete _assembler;
	delete _solver;
	delete _solution;
	delete _external_load;
	delete _internal_vars;
	delete _writer;
	delete _options;
	delete _boundary;
	if (_coh_writer != NULL)
		delete _coh_writer;
	if (_sensitivity_solver != NULL)
		delete _sensitivity_solver;
	if (_sensitivity_writer != NULL)
		delete _sensitivity_writer;
	if (_coh_failure_writer != NULL)
		delete _coh_failure_writer;
	if (_path_writer != NULL)
		delete _path_writer;
	if (_coh_openings_writer != NULL)
		delete _coh_openings_writer;
	if (_homogenization_writer != NULL)
		delete _homogenization_writer;
	for (id_type i=0; i<n_SensitivityFunctions(); ++i)
		delete _sensitivity_functions[i];
	for (id_type i=0; i<n_SensitivityParameters(); ++i)
		delete _parameters[i];
	_sensitivity_functions.clear();
	_parameters.clear();
}

Writer* Problem::get_cohesive_writer() const
{
	if (_options->getBoolOption("-OUTPUT_COH") && _mesh->is_cohesive())
		return _coh_writer;
	else
		err_message("Attempting to get the cohesive writer for a non-cohesive problem!");
}



void Problem::init()
{
	if(_mesh==NULL)
		err_message("Please attach a mesh to the problem before attempting to solve it.");

	// First thing's first, lets check all of the materials in the mesh to make sure they match the problem type (Will produce error if there's a mismatch)
	check_materials();

	// Add all of the enrichments here since this is the only place where I'm sure all inclusions are done being added
	{
		// (1)
		if (_mesh->IGFEM())
		{
			// Clear any old enrichments
			_mesh->clearEnrichments();

			Utilities::timer timer;

			// (1a)
			PIGFEMPrint("\n\tDetecting Nodes... ");
			timer.start();
			_mesh->analyze_nodes();
			_IGFEM_nodal_detection_time = timer.time_elapsed();
			PIGFEMPrint(_IGFEM_nodal_detection_time << " seconds");

			// (1b)
			// PIGFEMPrint("\n\tRefining the mesh...");
			// time.start();
			// RefinerGeometricalInclusion refiner(this);
			// refiner.refine();
			// _mesh->set_constraint_values();
			// _refine_time = timer.time_elapsed();
			// PIGFEMPrint(_refine_time << " seconds");


			// (1c)
			PIGFEMPrint("\n\tDetecting Elements... ");
			timer.start();
			_mesh->analyze_elements();
			_IGFEM_element_detection_time = timer.time_elapsed();
			PIGFEMPrint(_IGFEM_element_detection_time << " seconds");

			// (1d)
			PIGFEMPrint("\n\tAdding Enrichment Nodes... ");
			timer.start();
			_mesh->add_enrichments();
			_IGFEM_enrichments_time = timer.time_elapsed();
			PIGFEMPrint(_IGFEM_enrichments_time << " seconds");
		}
	}


	// Create and initialize all of the objects needed to solve the problem
	{
		// (2) Precalculate all shape function and gradient values on the current mesh configuration
		_mesh->store_shape_functions();

		// (5) Preallocate storage for the solution and loads
		if (_solution!=NULL)
			delete _solution;
		_solution = new NodalData(_mesh);
		_solution->preallocate_storage(nndof());
		if (_external_load != NULL)
			delete _external_load;
		_external_load = new NodalData(_mesh);
		_external_load->preallocate_storage(nndof());

		// Generate the assembler for this problem
		if (_assembler != NULL)
			delete _assembler;
		generate_assembler();
		_assembler->attach_problem(this);
		_assembler->storeBmats(); // Some classes of problems show improvement from pre-storing the shape function gradients in a certain arrangement

		// Generate the solver for this problem
		if (_solver != NULL)
			delete _solver;
		generate_solver();
		_solver->attach_problem(this);

		if (_internal_vars!=NULL)
			delete _internal_vars;
		_internal_vars = new InternalVars(this, false);
	}


	// Initialize the sensitivity formulation
	{
		if (n_SensitivityFunctions() > 0)
		{
			_sensitivity = true;
			if (_sensitivity_solver!=NULL)
				delete _sensitivity_solver;
			if (n_SensitivityParameters() > n_SensitivityFunctions())
			{
				if (_internal_vars->haveISVs())
				{
					if (_options->getBoolOption("-SENSITIVITY_FINAL_ONLY"))
						err_message("For problems with ISVs, sensitivities must be calculated at all load steps");
					_sensitivity_solver = new SensitivityDirect; // I don't know how to handle these with ana adjoint method yet
				}
				else
					_sensitivity_solver = new SensitivityAdjoint;
			}
			else
				_sensitivity_solver = new SensitivityDirect;

			_sensitivity_solver->attachProblem(this);

			// Add all of the parameters
			id_type id = 0;
			for (id_type i=0; i<_parameters.size(); ++i)
			{
				_parameters[i]->set_id(id);
				_sensitivity_solver->addParameter(_parameters[i]);
				id++;
			}

			// Add the functions
			for (id_type i=0; i<_sensitivity_functions.size(); ++i)
			{
				_sensitivity_solver->addFunction(_sensitivity_functions[i]);
			}

			_sensitivity_solver->init();
		}
	}


	// Create all of the Writers
	{
		// Creater the VTK Writer context
		if (_options->getBoolOption("-WRITE_VTK"))
		{
			if (_writer!=NULL)
				delete _writer;
			_writer = new Writer_ASCII_VTK(this); // default
			// _writer = new Writer_Binary_VTK(this); // Binary writer still doesn't work...
			_writer->store_mesh(); // Store this initial state of the mesh
			if (_options->getBoolOption("-OUTPUT_COH") && _mesh->is_cohesive())
			{
				if (_coh_writer!=NULL)
					delete _coh_writer;
				_coh_writer = new Writer_ASCII_VTK(this); // default
				_coh_writer->setParameter("COHESIVE", true);
				_coh_writer->store_mesh(); // Store this initial state of the mesh
			}
		}

		// If the user set the output cohesive failure lengths flag, create a writer
		if (_mesh->is_cohesive() && _options->getBoolOption("-OUTPUT_COH_FAIL"))
		{
			if (_coh_failure_writer!=NULL)
				delete _coh_failure_writer;
			_coh_failure_writer = new WriterCohesiveFailureLengths(this);
			_coh_failure_writer->setParameter("INDIVIDUAL", _options->getBoolOption("-OUTPUT_ALL_COH_FAIL"));
			_coh_failure_writer->store_mesh();
		}

		// If the user wants to output the solution along a path, build the writer here
		if (_options->getBoolOption("-WRITE_PATH"))
		{
			std::ifstream myfile;
			myfile.open( _options->getOption("-PATH_IN") );
			if (!myfile.good())
				err_message("Invalid input file name!");

			// Make the path writer
			if (_path_writer!=NULL)
				delete _path_writer;
			_path_writer = new WriterPath(this);
			_path_writer->store_mesh();

			// Store all points that are contained in the input file
			for (std::string line; std::getline(myfile, line); )
		    {
		        // Remove whitespace
		        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

		        // Parse into comma-separated values
		        std::vector<std::string> vals = Utilities::splitString(line, ",");
		        if (vals.size() != _mesh->dim())
		        	err_message("Path input file must contain points with the same dimension as the problem being run!");
		        else
		        {
		        	std::vector<double> point(vals.size());
		        	for (id_type d=0; d<vals.size(); ++d)
		        		point[d] = std::stod(vals[d]);
		        	_path_writer->setParameter("POINT", point);
		        }
		    }

		    _path_writer->store_mesh();
		}

		// If the use set the cohesive openings flag, make the writer here
		if (_mesh->is_cohesive() && _options->getBoolOption("-WRITE_COH_OPENINGS"))
		{
			if (_coh_openings_writer!=NULL)
				delete _coh_openings_writer;
			_coh_openings_writer = new WriterCohesiveOpenings(this);
		}

		// Create a homogenization writer
		if (_options->getBoolOption("-HOMOGENIZE"))
			_homogenization_writer = new WriterHomogenization(this);

		// Create the sensitivities writer
		if (sensitivity())
		{
			if (_sensitivity_writer!=NULL)
				delete _sensitivity_writer;
			_sensitivity_writer = new WriterSensitivity(this);
		}
	}

	_init = true;
}


void Problem::setup()
{
	if (_boundary == NULL || _dofs==NULL)
		err_message("Attempting to set up the finite element system before creating a mesh!");

	// Apply any BCs that were assciated with nodesets or sidesets now
	_boundary->apply_BCs();
		
	// Now that all nodes have been added, distribute the degrees of freedom
	_dofs->distribute_dofs(nndof());

	// Set up the finite element system inside the solver itself
	_solver->setup();

	// If the problem has sensitivities, set up the sensitivity solve
	if (sensitivity())
		_sensitivity_solver->setup();

	_setup = true;
}


























void Problem::attach_mesh(Mesh* mesh)
{
	_mesh = mesh;
	_dofs->attachProblem(this);
	_boundary->attachProblem(this);
}




// Check that all of the materials in the mesh are the correct classification for the given problem type
void Problem::check_materials()
{
	std::vector<Material*>& mats = _mesh->get_materials();

	for (id_type m=0; m<mats.size(); ++m)
		if (mats[m]->get_classification() != get_classification())
			err_message("Invalid material classification for given problem type.");
}







// Actually calls the solver to solve te problem at hand
void Problem::solve_problem()
{
	// Call the initialization routine
	if (!_init)
		err_message("Please initialize the problem before calling solve!");

	setup();

	// Actually solve the problem!
	Utilities::timer timer;
	timer.start();
	//_solver->purely_explicit_solve();
	_solver->solve();
	double total_solve_time = timer.time_elapsed();
	_misc_solve_time = total_solve_time - _assemble_time - _solve_time - _solver_write_time - _sensitivity_time;
	_solved = true;

	writeStatistics();
}








myFunctor<std::vector<double> >* Problem::getStrainKernel()
{
	return _assembler->getStrainKernel();
}
myFunctor<std::vector<double> >* Problem::getStressKernel()
{
	return _assembler->getStressKernel();
}

























// Computes the L2 Norm of the error of the solution with respect to the analytical function provided as input
double Problem::Compute_L2_Norm(std::vector<double> (*analytical_disp)(std::vector<double>))
{
	MyL2NormFunctor functor;
	functor._ptr = (void*) analytical_disp;
	Integrator<double> integrator;
	double e2 = integrator.integrate(this, &functor);
	return sqrt(e2);
}






// Get a pointer to the assembler's homogenization kernel
// Default ot integrate over volumetric portion of domain only
// Function is virtual to allow for additional functionality if desired though.
void Problem::homogenize()
{
	myFunctor<std::vector<double> >* functor = _assembler->getHomogenizationKernel();
	functor->_ptr = (void*) this;
	Integrator<std::vector<double> > integrator;
	_homogneizedValue = integrator.integrate(this, functor);

	VolumeFunctor vol;
	Integrator<double> volIntegrator;
	double totalVolume = volIntegrator.integrateSimple(this, &vol);
	_homogneizedValue = Utilities::scale(_homogneizedValue, 1./totalVolume);

	delete functor;
}

// Defaults to just being a wrapper for the assembler function
// Function is virtual to allow for additional functionality if desired though
std::vector<double> Problem::computeMacroscaleSolutionContribution(const std::vector<double>& coords)
{
	return _assembler->computeMacroscaleSolutionContribution(coords); 
}














































void Problem::interpolate_nodal_enrichment(NodalData& solution, NodalData& external_load)
{
	// Store the current normal solution
	for (id_type n=0; n<_mesh->n_local_nodes(); ++n)
		for (id_type d=0; d<nndof(); ++d)
		{
			solution.get_value_local(n, d) = _solution->get_value_local(n,d);
			external_load.get_value_local(n, d) = _external_load->get_value_local(n,d);
		}

	// Interpolate between the edges to get the solution value at the enrichment nodes
	for(Mesh::enrich_node_iterator it=_mesh->enrich_nodes_begin(), end=_mesh->enrich_nodes_end(); it!=end; ++it)
	{
		double interpolate = _mesh->getEnrichmentEdgeDistance( (*it)->get_id() );
		double scale_factor = _mesh->getEnrichementScaleFactorGlobal( (*it)->get_id() );
		id_type l_node = _mesh->global_to_local_node((*it)->get_id());
		std::vector<id_type> edge = _mesh->get_enriched_edge_local( l_node );
		id_type n1 = _mesh->global_to_local_node( edge[0] );
		id_type n2 = _mesh->global_to_local_node( edge[1] );

		for(id_type d=0; d<nndof(); ++d)
		{
			double sol1 = get_solution()->get_value_local(n1, d);
			double sol2 = get_solution()->get_value_local(n2, d);
			double load1 = get_external_load()->get_value_local(n1, d);
			double load2 = get_external_load()->get_value_local(n2, d);
			// solution.get_value_local(l_node, d) = sol1*(1.0-interpolate) + sol2*interpolate + (get_solution()->get_value_local(l_node, d) / scale_factor);
			// external_load.get_value_local(l_node, d) = load1*(1.0-interpolate) + load2*interpolate + (get_external_load()->get_value_local(l_node, d) / scale_factor);
			solution.get_value_local(l_node, d) = sol1*(1.0-interpolate) + sol2*interpolate + (get_solution()->get_value_local(l_node, d) * scale_factor);
			external_load.get_value_local(l_node, d) = load1*(1.0-interpolate) + load2*interpolate + (get_external_load()->get_value_local(l_node, d) * scale_factor);
		}
	}
}






























void Problem::writeStatistics()
{
	std::string filename = _options->getOption("-OUT_FOLDER") + _options->getOption("-STATS_FILE");
	if (_mesh->serial())
	{
		std::ofstream myfile;
		myfile.open(filename.c_str(), std::ofstream::out);
		if (!myfile.good())
			err_message("Error opening statistics output file!");
		myfile.precision(16);
		myfile << "PROBLEM DETAILS\n---------------------------------------------------\n";
		myfile << "Problem Type: " << problem_type_names[get_type()] << std::endl;
		myfile << "Number of Mesh Partitions: " << _mesh->n_ranks() << std::endl;
		myfile << "Number of Inclusions: " << _mesh->n_inclusions() << std::endl;
		myfile << "Number of Elements: " << _mesh->n_global_elem() << std::endl;
		myfile << "Number of Nodes: " << _mesh->n_global_nodes()+_mesh->n_global_enrich_nodes() << " (" << _mesh->n_global_enrich_nodes() << " enriched)" << std::endl;
		myfile << "Number of dofs: " << _dofs->n_global_dofs()  << " (" << _dofs->n_global_free_dofs() << " free)" << std::endl;

		myfile << "\nTIMING DETAILS\n---------------------------------------------------\n";
		myfile << "IGFEM nodal detection time: " << _IGFEM_nodal_detection_time << " seconds" << std::endl;
		myfile << "IGFEM element detection time: " << _IGFEM_element_detection_time << " seconds" << std::endl;
		myfile << "IGFEM enrichment addition time: " << _IGFEM_enrichments_time << " seconds" << std::endl;
		myfile << "\tIGFEM total time: " << _IGFEM_nodal_detection_time + _IGFEM_element_detection_time + _IGFEM_enrichments_time << " seconds" << std::endl;
		myfile << "Total assembly time: " << _assemble_time << " seconds" << std::endl;
		myfile << "Total linear solver time: " << _solve_time << " seconds" << std::endl;
		myfile << "Total file I/O time: " << _solver_write_time << " seconds" << std::endl;
		myfile << "Total miscellaneous solver time: " << _misc_solve_time << " seconds" << std::endl;
		if (sensitivity())
		{
			double solve = _sensitivity_solver->get_solve_time();
			double assemble = _sensitivity_solver->get_assembly_time();
			double subs = _sensitivity_solver->get_substitution_time();
			myfile << "Total sensitivity solver time: " << solve << " seconds" << std::endl;
			myfile << "Total sensitivity assembly time: " << assemble << " seconds" << std::endl;
			myfile << "Total sensitivity substitution time: " << subs << " seconds" << std::endl;
			myfile << "Total additional sensitivity time: " << _sensitivity_time << " seconds" << std::endl;
			myfile << "\tTotal Solver time: " << _assemble_time+_solve_time+_solver_write_time+_misc_solve_time+_sensitivity_time << " seconds" << std::endl;
		}
		else
			myfile << "\tTotal Solver time: " << _assemble_time+_solve_time+_solver_write_time+_misc_solve_time << " seconds" << std::endl;

		myfile.close();
	}
	else
	{
		double tot_assem_time, tot_solve_time, tot_misc_time, tot_write_time, tot_node_time, tot_elem_time, tot_enrich_time;
		MPI_Reduce(&_IGFEM_nodal_detection_time, &tot_node_time, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
		MPI_Reduce(&_IGFEM_element_detection_time, &tot_elem_time, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
		MPI_Reduce(&_IGFEM_enrichments_time, &tot_enrich_time, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
		MPI_Reduce(&_assemble_time, &tot_assem_time, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
		MPI_Reduce(&_solve_time, &tot_solve_time, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
		MPI_Reduce(&_solver_write_time, &tot_write_time, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
		MPI_Reduce(&_misc_solve_time, &tot_misc_time, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
		id_type n_ranks = _mesh->n_ranks();
		tot_node_time /= n_ranks;
		tot_elem_time /= n_ranks;
		tot_enrich_time /= n_ranks;
		tot_assem_time /= n_ranks;
		tot_solve_time /= n_ranks;
		tot_write_time /= n_ranks;
		tot_misc_time /= n_ranks;
		double tot_sens_solve(0.0), tot_sens_assem(0.0), tot_sens_subs(0.0), tot_sens(0.0);
		if (sensitivity())
		{
			double solve = _sensitivity_solver->get_solve_time();
			double assemble = _sensitivity_solver->get_assembly_time();
			double subs = _sensitivity_solver->get_substitution_time();
			MPI_Reduce(&solve, &tot_sens_solve, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
			MPI_Reduce(&assemble, &tot_sens_assem, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
			MPI_Reduce(&subs, &tot_sens_subs, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
			MPI_Reduce(&_sensitivity_time, &tot_sens, 1, MPI_DOUBLE, MPI_SUM, 0, _mesh->get_comm());
			tot_sens_solve /= n_ranks;
			tot_sens_assem /= n_ranks;
			tot_sens_subs /= n_ranks;
			tot_sens /= n_ranks;
		}

		if (_mesh->get_rank() == 0)
		{
			std::ofstream myfile;
			myfile.open(filename.c_str(), std::ofstream::out);
			if (!myfile.good())
				err_message("Error opening statistics output file!");
			myfile.precision(16);

			myfile << "PROBLEM DETAILS\n---------------------------------------------------\n";
			myfile << "Problem Type: " << problem_type_names[get_type()] << std::endl;
			myfile << "Number of Mesh Partitions: " << _mesh->n_ranks() << std::endl;
			myfile << "Number of Inclusions: " << _mesh->n_inclusions() << std::endl;
			myfile << "Number of Elements: " << _mesh->n_global_elem() << std::endl;
			myfile << "Number of Nodes: " << _mesh->n_global_nodes()+_mesh->n_global_enrich_nodes() << " (" << _mesh->n_global_enrich_nodes() << " enriched)" << std::endl;
			myfile << "Number of dofs: " << _dofs->n_global_dofs()  << " (" << _dofs->n_global_free_dofs() << " free)" << std::endl;

			myfile << "\nTIMING DETAILS\n---------------------------------------------------\n";
			myfile << "IGFEM nodal detection time (average): " << tot_node_time << " seconds" << std::endl;
			myfile << "IGFEM element detection time (average): " << tot_elem_time << " seconds" << std::endl;
			myfile << "IGFEM enrichment addition time (average): " << tot_enrich_time << " seconds" << std::endl;
			myfile << "\tIGFEM total time (average): " << tot_node_time + tot_elem_time + tot_enrich_time << " seconds" << std::endl;
			myfile << "Total assembly time (average): " << tot_assem_time << " seconds" << std::endl;
			myfile << "Total linear solver time (average): " << tot_solve_time << " seconds" << std::endl;
			myfile << "Total file I/O time (average): " << tot_write_time << " seconds" << std::endl;
			myfile << "Total miscellaneous solver time (average): " << tot_misc_time << " seconds" << std::endl;
			if (sensitivity())
			{
				myfile << "Total sensitivity solver time (average): " << tot_sens_solve << " seconds" << std::endl;
				myfile << "Total sensitivity assembly time (average): " << tot_sens_assem << " seconds" << std::endl;
				myfile << "Total sensitivity substitution time (average): " << tot_sens_subs << " seconds" << std::endl;
				myfile << "Total additional sensitivity time (average): " << tot_sens << " seconds" << std::endl;
				myfile << "\tTotal Solver time (average): " << tot_assem_time+tot_solve_time+tot_write_time+tot_misc_time+tot_sens << " seconds" << std::endl;
			}
			else
				myfile << "\tTotal Solver time (average): " << tot_assem_time+tot_solve_time+tot_write_time+tot_misc_time << " seconds" << std::endl;

			myfile.close();
		}
	}
}






void Problem::addSensitivityFunction(SensitivityFunction* func)
{
	_sensitivity_functions.push_back(func->allocate_and_copy());
}
// Functions to add parameters
void Problem::addMaterialSensitivityParameter(std::string mat_name, std::string param_name)
{
	// Find the correct material from the mesh
	int mat_num = _mesh->find_mat_number(mat_name);
	if (mat_num < 0)
		err_message("Invalid material name while addint a material sensitivity parameter");
	Material* mat = _mesh->get_material(mat_num);

	// Get the Sensitivity parameter directly from the material
	SensitivityParameter* param = mat->getSensitivityParameter(param_name);
	_parameters.push_back(param);
}
void Problem::addShapeSensitivityParameter(id_type inclusion_id, std::string param_name)
{
	Inclusion* inc = _mesh->get_inclusion(inclusion_id);

	SensitivityParameter* param = inc->getSensitivityParameter(param_name, get_classification());
	_parameters.push_back(param);
}
void Problem::addLoadSensitivityParameter(/* Not sure what parameters to place here yet */)
{
	err_message("Not totally sure how to do sensitivity wrt a load yet");
}