/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "AssemblerThermal.h"
#include "Problem.h"
#include "Utilities.h"



bool AssemblerThermal::storedBmats()
{
	return true;
}

// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
void AssemblerThermal::KernelVolumetric(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
											  const std::vector<double>& shape, const DenseMatrix<double>& B,
											  Material* mat, Material::input_params& input,
											  std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	ProblemUtilities::B_U(input.strain, B, elem_U_curr);

	// Add the current macroscopic strain to it
	if (_prob->subscale())
	{
		if (_currentMacroscopicValue.size() == 1)
			input.temp += _currentMacroscopicValue[0];
		else
			err_message("Invalid macroscopic value applied");
	}

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::Bt_D_B(K_el, B, output->Dmat);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::Bt_F(P_el_int, B, output->stress);
}


void AssemblerThermal::KernelVolumetric(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
											  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
											  Material* mat, Material::input_params& input,
											  std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Just do the same thing as the B-matrix version, just assemble the b matrix here
	DenseMatrix<double> B;
	ProblemUtilities::Assemble_Thermal_B_Mat(B, shape_grad);

	// Compute the current strain
	ProblemUtilities::B_U(input.strain, B, elem_U_curr);

	// Add the current macroscopic strain to it
	if (_prob->subscale())
	{
		if (_currentMacroscopicValue.size() == 1)
			input.temp += _currentMacroscopicValue[0];
		else
			err_message("Invalid macroscopic value applied");
	}

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::Bt_D_B(K_el, B, output->Dmat);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::Bt_F(P_el_int, B, output->stress);
}




void AssemblerThermal::KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
													const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
													Material* mat, Material::input_params& input,
													std::vector<double>& coh_U_curr, bool assembleFunc, bool assembleJac)
{
	K_coh.clear();
	K_coh.resize(shape.size(), shape.size());
	P_int_coh.clear();
	P_int_coh.resize(shape.size());
}



// Select the appropriate B-matrix fill function
void AssemblerThermal::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& shape_grad)
{
	ProblemUtilities::Assemble_Thermal_B_Mat(B, shape_grad);
}