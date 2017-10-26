/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _SSSTRUCT_ASSEMBLER_H_
#define _SSSTRUCT_ASSEMBLER_H_
#include "Assembler.h"
#include "common.h"
#include "Utilities.h"













class AssemblerSmallStrainStructural : public Assembler
{
	protected:
		class MyStressFunctor : public myFunctor<std::vector<double> >
		{
			public:
				typedef myFunctor<std::vector<double> > Base; // alias
				
				virtual std::vector<double> evaluate(elemState* state)
				{
					// Cast the void pointer to a function pointer
					Problem* prob = (Problem*) _ptr;

					// Get the strain
					Material::input_params& input = *(state->input);
					ProblemUtilities::SmallStrainFast(input.strain, *(state->dN), *(state->Uel));

					// Add the current macroscopic strain to it
					if (prob->subscale())
					{
						std::vector<double> curr_val = prob->get_assembler()->getCurrentMacroscopicValue();
						if (input.strain.size() == curr_val.size())
							Utilities::_VecAXPY(1.0, curr_val, input.strain);
						else
							err_message("Invalid macroscopic value applied");
					}

					// Compute the constitutive relations
					Material::output_params* output = state->Mat->Constitutive(input);
					
					return output->stress;
				}
		};

		class MyStrainFunctor : public myFunctor<std::vector<double> >
		{
			public:
				typedef myFunctor<std::vector<double> > Base; // alias
				
				virtual std::vector<double> evaluate(elemState* state)
				{
					// Cast the void pointer to a function pointer
					Problem* prob = (Problem*) _ptr;

					// Get the strain
					Material::input_params& input = *(state->input);
					ProblemUtilities::SmallStrainFast(input.strain, *(state->dN), *(state->Uel));

					// Add the current macroscopic strain to it
					if (prob->subscale())
					{
						std::vector<double> curr_val = prob->get_assembler()->getCurrentMacroscopicValue();
						if (input.strain.size() == curr_val.size())
							Utilities::_VecAXPY(1.0, curr_val, input.strain);
						else
							err_message("Invalid macroscopic value applied");
					}
					
					return input.strain;
				}
		};

	public:
		virtual ~AssemblerSmallStrainStructural() {};

		virtual bool storedBmats();


		// Functions having to do with multiscale homogenization
		virtual myFunctor<std::vector<double> >* getHomogenizationKernel() {return new MyStressFunctor;};
		// virtual MyFunctor* getHomogenizationKernel() {return new MyStrainFunctor;};
		virtual std::vector<double> computeMacroscaleSolutionContribution(const std::vector<double>& coords);


		// Functions having to do with visualization output
		// Get stresses and strains. I want to think of a more general name for these but these will do for the moment
		virtual myFunctor<std::vector<double> >* getStrainKernel() {return new MyStrainFunctor;};
		virtual myFunctor<std::vector<double> >* getStressKernel() {return new MyStressFunctor;};


	protected:

		// Kernel function that actually does the computes the physics behind the specific problem
		// Computes the local contribution to the stiffness matrix and internal load vector
		// Also updates the internal variable storage in the input variable
		virtual void KernelVolumetric(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
									  const std::vector<double>& shape, const DenseMatrix<double>& B,
									  Material* mat, Material::input_params& input,
									  std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac);
		virtual void KernelVolumetric(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
									  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
									  Material* mat, Material::input_params& input,
									  std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac);
		virtual void KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
									 const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
									 Material* mat, Material::input_params& input,
									 std::vector<double>& coh_U_curr, bool assembleFunc, bool assembleJac);

		virtual void fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& shape_grad);
};


#endif