/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated September 2017

##################################################################################
*/
#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_
#include "common.h"
#include "DenseMatrix.h"
#include "Utilities.h"
#include "material.h"
#include "InternalVars.h"
#include "Mesh.h"
#include "elem.h"
#include "Problem.h"




// Define the Structs and Functors that will be used to actually call this integrator
// =========================================================================

struct elemState {
	DenseMatrix<double>* Uel;
	DenseMatrix<double>* X;
	std::vector<double>* N;
	std::vector<std::vector<double> >* dN;
	Material* Mat;
	Material::input_params* input;
	std::vector<double>* ISVs;
};

template<class T>
class myFunctor
{
	public:
		void* _ptr;	// General void pointer to store anything needed by the function (Use function pointer in L2 norm)
		myFunctor() {_ptr=NULL;};
		virtual T evaluate(elemState* state) = 0;
};


class MyL2NormFunctor : public myFunctor<double>
{
	public:
		typedef myFunctor<double> Base; // alias
		typedef std::vector<double> (*fptr)(std::vector<double>);
		
		virtual double evaluate(elemState* state)
		{
			// Cast the void pointer to a function pointer
			fptr func = (fptr) _ptr;
			
			// Get the current coordinates and local solution vector
			std::vector<double> gcoords = (*(state->N)) * (*(state->X));
			std::vector<double> sol = (*(state->N)) * (*(state->Uel));

			std::vector<double> analytical = func(gcoords);
			std::vector<double> diff = Utilities::minus(analytical, sol);
			return Utilities::dot(diff, diff);
		}
};


class VolumeFunctor : public myFunctor<double>
{
	public:
		virtual double evaluate(elemState* state)
		{
			return 1.0;
		}
};






// Define the integrator function itself
// ===================================================

template <class T>
class Integrator
{
	private:
		void AXPY(double alpha, T& X, T& Y)
		{
			Y += (alpha * X);
		}
		T reduceDomains(T& domainVal, MPI_Comm communicator)
		{
			err_message("Integrator not defined for parameter type");
		}

	public:

		T integrate(Problem* prob, myFunctor<T>* integrand);
		T integrateElem(Problem* prob, Elem* el, myFunctor<T>* integrand, int int_domain);

		// Simple calls only provide the position location (X and N), not the solution location to the kernel
		T integrateSimple(Problem* prob, myFunctor<T>* integrand);
		T integrateElemSimple(Problem* prob, Elem* el, myFunctor<T>* integrand, int int_domain);
};

















// Main Integrator functions
// ======================================================
template<class T>
T Integrator<T>::integrate(Problem* prob, myFunctor<T>* integrand)
{
	Mesh* mesh = prob->get_mesh();

	// initialize integrand
	T domainVal;

	// Loop over all of the local elements and compute the local contribution to the L2 norm
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		T elemVal = integrateElem(prob, (*it), integrand);
		AXPY(1.0, elemVal, domainVal);
	}
		
	// Now perform a global reduction to get the global value
	if(mesh->serial())
		return domainVal;
	else
		return reduceDomains(domainVal, mesh->get_comm());
}


template<class T>
T Integrator<T>::integrateElem(Problem* prob, Elem* el, myFunctor<T>* integrand, int int_domain=-1)
{
	Mesh* mesh = prob->get_mesh();
	T elemVal(0);
	id_type l_elem = mesh->global_to_local_elem(el->get_id());

	// Assemble the current solution vector
	DenseMatrix<double> elem_sol(0);
	ProblemUtilities::formElementalMatrix(elem_sol, el, prob->get_solution());

	// Get the element's node positions (Don't do anything with enrichment nodes because position isn't based off of them)
	DenseMatrix<double> X((el->n_nodes()+el->n_enrich_nodes()), mesh->dim());
	for (id_type n=0; n<el->n_nodes(); ++n)
		for (id_type d=0; d<mesh->dim(); ++d)
			X(n, d) = (*el)(n)(d);

	Material::input_params input;
	input.dim = mesh->dim();
	input.delta_t = 0.0;
	if (prob->get_classification() == STRUCTURAL)
	{
		std::vector<int> switch_dim = {1, 3, 6};
		input.strain.resize(switch_dim[input.dim - 1]); // I don't know if I need this here. Just a placeholder in case I don't wanna compute the actual strain
		if (prob->get_parameter("plane_strain") == 0.0)
			input.plane_strain = false;
		else
			input.plane_strain = true;
	}

	// Define the current Element state that I will pass to the function
	elemState currentState;
	currentState.Uel = &elem_sol;
	currentState.X = &X;
	currentState.input = &input;

	id_type curr_qp = 0;
	id_type start, stop;
	if (int_domain == -1)
	{
		start = 0;
		stop = el->n_integration_elem();
	}
	else
	{
		start = int_domain;
		stop = int_domain+1;
		for(id_type ie=0; ie<start; ++ie)
			curr_qp += el->get_integration_elem(ie)->n_q_points();
	}
	for(id_type ie=start; ie<stop; ++ie)
	{
		Elem* int_el = el->get_integration_elem(ie);
		Material* curr_mat = mesh->get_element_material_local(l_elem, ie);

		id_type nqp = int_el->n_q_points();
		for(id_type qp=0; qp<nqp; ++qp)
		{
		
			currentState.N = &mesh->get_shape(l_elem,curr_qp);
			currentState.dN = &mesh->get_shape_grad(l_elem,curr_qp);
			currentState.Mat = curr_mat;
			currentState.ISVs = &prob->get_internal_vars()->get_internal_vars_local(l_elem, curr_qp);

			// Compute the output
			T localVal = integrand->evaluate(&currentState);

			// Compute contibution to local L2 norm
			double& W = mesh->get_W(l_elem, curr_qp);
			double& J = mesh->get_J(l_elem, curr_qp);
			double wgt = J * W;
			
			// Add to the existing elemental value
			AXPY(wgt, localVal, elemVal);

			curr_qp++;
		}
	}

	return elemVal;
}




template<class T>
T Integrator<T>::integrateSimple(Problem* prob, myFunctor<T>* integrand)
{
	Mesh* mesh = prob->get_mesh();

	// initialize integrand
	T domainVal;

	// Loop over all of the local elements and compute the local contribution to the L2 norm
	for(Mesh::element_iterator it=mesh->active_elements_begin(), end=mesh->active_elements_end(); it!=end; ++it)
	{
		T elemVal = integrateElemSimple(prob, (*it), integrand, -1);
		AXPY(1.0, elemVal, domainVal);
	}
		
	// Now perform a global reduction to get the global value
	if(mesh->serial())
		return domainVal;
	else
		return reduceDomains(domainVal, mesh->get_comm());
}



template<class T>
T Integrator<T>::integrateElemSimple(Problem* prob, Elem* el, myFunctor<T>* integrand, int int_domain=-1)
{
	Mesh* mesh = prob->get_mesh();
	T elemVal(0);
	id_type l_elem = mesh->global_to_local_elem(el->get_id());

	// Get the element's node positions (Don't do anything with enrichment nodes because position isn't based off of them)
	DenseMatrix<double> X((el->n_nodes()+el->n_enrich_nodes()), mesh->dim());
	for (id_type n=0; n<el->n_nodes(); ++n)
		for (id_type d=0; d<mesh->dim(); ++d)
			X(n, d) = (*el)(n)(d);


	// Define the current Element state that I will pass to the function
	elemState currentState;
	currentState.X = &X;

	id_type curr_qp = 0;
	id_type start, stop;
	if (int_domain == -1)
	{
		start = 0;
		stop = el->n_integration_elem();
	}
	else
	{
		start = int_domain;
		stop = int_domain+1;
		for(id_type ie=0; ie<start; ++ie)
			curr_qp += el->get_integration_elem(ie)->n_q_points();
	}
	for(id_type ie=start; ie<stop; ++ie)
	{
		Elem* int_el = el->get_integration_elem(ie);

		id_type nqp = int_el->n_q_points();
		for(id_type qp=0; qp<nqp; ++qp)
		{
		
			currentState.N = &mesh->get_shape(l_elem,curr_qp);

			// Compute the output
			T localVal = integrand->evaluate(&currentState);

			// Compute contibution to local L2 norm
			double& W = mesh->get_W(l_elem, curr_qp);
			double& J = mesh->get_J(l_elem, curr_qp);
			double wgt = J * W;
			
			// Add to the existing elemental value
			AXPY(wgt, localVal, elemVal);

			curr_qp++;
		}
	}

	return elemVal;
}





























// Helper functions for various data types
// ============================================================

template<>
inline
void Integrator<std::vector<double> >::AXPY(double alpha, std::vector<double>& X, std::vector<double>& Y)
{
	if (Y.size() == 0)
		Y.resize(X.size());
	cblas_daxpy(Y.size(), alpha, X.data(), 1, Y.data(), 1);
}

template<>
inline
void Integrator<DenseMatrix<double> >::AXPY(double alpha, DenseMatrix<double>& X, DenseMatrix<double>& Y)
{
	if (Y.n_rows() == 0)
		Y.resize(X.n_rows(), X.n_cols());
	Y.AXPY(alpha, X);
}

template<>
inline
void Integrator<std::vector<float> >::AXPY(double alpha, std::vector<float>& X, std::vector<float>& Y)
{
	if (Y.size() == 0)
		Y.resize(X.size());
	cblas_saxpy(Y.size(), alpha, X.data(), 1, Y.data(), 1);
}

template<>
inline
void Integrator<DenseMatrix<float> >::AXPY(double alpha, DenseMatrix<float>& X, DenseMatrix<float>& Y)
{
	if (Y.n_rows() == 0)
		Y.resize(X.n_rows(), X.n_cols());
	Y.AXPY(alpha, X);
}

template<>
inline
void Integrator<std::vector<int> >::AXPY(double alpha, std::vector<int>& X, std::vector<int>& Y)
{
	if (Y.size() == 0)
		Y.resize(X.size());
	for (id_type i=0; i<Y.size(); ++i) // There's no BLAS AXPY routine for ints
		Y[i] = alpha * X[i] + Y[i];
}

template<>
inline
void Integrator<DenseMatrix<int> >::AXPY(double alpha, DenseMatrix<int>& X, DenseMatrix<int>& Y)
{
	if (Y.n_rows() == 0)
		Y.resize(X.n_rows(), X.n_cols());
	Y.AXPY(alpha, X);
}














template<>
inline
double Integrator<double>::reduceDomains(double& domainVal, MPI_Comm communicator)
{
	double globalVal;
	MPI_Allreduce(&domainVal, &globalVal, 1, MPI_DOUBLE, MPI_SUM, communicator);
	return globalVal;
}

template<>
inline
std::vector<double> Integrator<std::vector<double> >::reduceDomains(std::vector<double>& domainVal, MPI_Comm communicator)
{
	std::vector<double> globalVal(domainVal.size());
	MPI_Allreduce(domainVal.data(), globalVal.data(), domainVal.size(), MPI_DOUBLE, MPI_SUM, communicator);
	return globalVal;
}

template<>
inline
DenseMatrix<double> Integrator<DenseMatrix<double> >::reduceDomains(DenseMatrix<double>& domainVal, MPI_Comm communicator)
{
	id_type size = domainVal.n_rows() * domainVal.n_cols();
	DenseMatrix<double> globalVal(domainVal.n_rows(), domainVal.n_cols());
	MPI_Allreduce(domainVal.data(), globalVal.data(), size, MPI_DOUBLE, MPI_SUM, communicator);
	return globalVal;
}

template<>
inline
float Integrator<float>::reduceDomains(float& domainVal, MPI_Comm communicator)
{
	float globalVal;
	MPI_Allreduce(&domainVal, &globalVal, 1, MPI_FLOAT, MPI_SUM, communicator);
	return globalVal;
}

template<>
inline
std::vector<float> Integrator<std::vector<float> >::reduceDomains(std::vector<float>& domainVal, MPI_Comm communicator)
{
	std::vector<float> globalVal(domainVal.size());
	MPI_Allreduce(domainVal.data(), globalVal.data(), domainVal.size(), MPI_FLOAT, MPI_SUM, communicator);
	return globalVal;
}

template<>
inline
DenseMatrix<float> Integrator<DenseMatrix<float> >::reduceDomains(DenseMatrix<float>& domainVal, MPI_Comm communicator)
{
	id_type size = domainVal.n_rows() * domainVal.n_cols();
	DenseMatrix<float> globalVal(domainVal.n_rows(), domainVal.n_cols());
	MPI_Allreduce(domainVal.data(), globalVal.data(), size, MPI_FLOAT, MPI_SUM, communicator);
	return globalVal;
}

template<>
inline
int Integrator<int>::reduceDomains(int& domainVal, MPI_Comm communicator)
{
	int globalVal;
	MPI_Allreduce(&domainVal, &globalVal, 1, MPI_INT, MPI_SUM, communicator);
	return globalVal;
}

template<>
inline
std::vector<int> Integrator<std::vector<int> >::reduceDomains(std::vector<int>& domainVal, MPI_Comm communicator)
{
	std::vector<int> globalVal(domainVal.size());
	MPI_Allreduce(domainVal.data(), globalVal.data(), domainVal.size(), MPI_INT, MPI_SUM, communicator);
	return globalVal;
}

template<>
inline
DenseMatrix<int> Integrator<DenseMatrix<int> >::reduceDomains(DenseMatrix<int>& domainVal, MPI_Comm communicator)
{
	id_type size = domainVal.n_rows() * domainVal.n_cols();
	DenseMatrix<int> globalVal(domainVal.n_rows(), domainVal.n_cols());
	MPI_Allreduce(domainVal.data(), globalVal.data(), size, MPI_INT, MPI_SUM, communicator);
	return globalVal;
}



#endif