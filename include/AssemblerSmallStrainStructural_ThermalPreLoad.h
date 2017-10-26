/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated August 2017

##################################################################################
*/
#ifndef _SSSTRUCT_THERMALPRELOAD_ASSEMBLER_H_
#define _SSSTRUCT_THERMALPRELOAD_ASSEMBLER_H_
#include "Assembler.h"
#include "common.h"











class AssemblerSmallStrainStructural_ThermalPreLoad : public Assembler
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
		AssemblerSmallStrainStructural_ThermalPreLoad();
		virtual ~AssemblerSmallStrainStructural_ThermalPreLoad();

		virtual void clear();


		// Function to assemble the presribed displacements and free loads at the beginning of each load step
		// (Or at the beginning of the analysis for linear)
		// This is made virtual so that calling it can update the assembler state at every load step is desired
		virtual PetscErrorCode assemble_new_load_step(Vec& Up, pvector& F_ext, double current_time_frac);


		// Functions having to do with multiscale homogenization
		virtual myFunctor<std::vector<double> >* getHomogenizationKernel() {return new MyStressFunctor;};
		// virtual MyFunctor* getHomogenizationKernel() {return new MyStrainFunctor;};
		virtual std::vector<double> computeMacroscaleSolutionContribution(const std::vector<double>& coords);


		// Functions having to do with visualization output
		// Get stresses and strains. I want to think of a more general name for these but these will do for the moment
		virtual myFunctor<std::vector<double> >* getStrainKernel() {return new MyStrainFunctor;};
		virtual myFunctor<std::vector<double> >* getStressKernel() {return new MyStressFunctor;};

		
		/*
		 * Function used to set a general single parameter
		 */
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);

		virtual bool storedBmats();

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
	private:

		// The total change in the temperature (negative for cooling off)
		double _delta_T;
		double _curr_delta_T; // Determined by hw far into the simulation we are

		// quick boolean check
		bool _temp_set;
		bool _relaxation_complete;
};


#endif