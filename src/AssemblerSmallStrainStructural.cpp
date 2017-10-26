/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "AssemblerSmallStrainStructural.h"
#include "material.h"
#include "Utilities.h"
#include "Problem.h"
#include <algorithm>




// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
//	B-matrix variant
void AssemblerSmallStrainStructural::KernelVolumetric(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
													  const std::vector<double>& shape, const DenseMatrix<double>& B,
													  Material* mat, Material::input_params& input,
													  std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	ProblemUtilities::B_U(input.strain, B, elem_U_curr);

	// Add the current macroscopic strain to it
	if (_prob->subscale())
	{
		if (input.strain.size() == _currentMacroscopicValue.size())
			Utilities::_VecAXPY(1.0, _currentMacroscopicValue, input.strain);
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


void AssemblerSmallStrainStructural::KernelVolumetric(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
													  const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
													  Material* mat, Material::input_params& input,
													  std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	ProblemUtilities::SmallStrainFast(input.strain, shape_grad, elem_U_curr);

	// Add the current macroscopic strain to it
	if (_prob->subscale())
	{
		if (input.strain.size() == _currentMacroscopicValue.size())
			Utilities::_VecAXPY(1.0, _currentMacroscopicValue, input.strain);
		else
			err_message("Invalid macroscopic value applied");
	}

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::SmallStrainKFast(K_el, shape_grad, output->Dmat);
	
	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::SmallStrainInternalForceFast(P_el_int, shape_grad, output->stress);
}





void AssemblerSmallStrainStructural::KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
													const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
													Material* mat, Material::input_params& input,
													std::vector<double>& coh_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current opening vector in the ttn coordinate system
	ProblemUtilities::cohesiveDeltaFast(input.delta, shape, coh_U_curr, rotation_matrix);

	// Compute the constitutive relations
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::cohesiveKFast(K_coh, shape, output->Dmat, rotation_matrix);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::cohesiveInternalForceFast(P_int_coh, shape, output->traction, rotation_matrix);	
}




bool AssemblerSmallStrainStructural::storedBmats()
{
	return true;
}

// Select the appropriate B-matrix fill function
void AssemblerSmallStrainStructural::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& shape_grad)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, shape_grad);
}





std::vector<double> AssemblerSmallStrainStructural::computeMacroscaleSolutionContribution(const std::vector<double>& coords)
{
	id_type dim = _prob->get_mesh()->dim();
	const std::vector<double>& origin = _prob->get_mesh()->origin();
	std::vector<double> coordsTranslated = Utilities::minus(coords, origin);

	if (dim == 1)
	{
		std::vector<double> ret(coords.size());
		ret[0] = (coordsTranslated[0]) * _currentMacroscopicValue[0];
		return ret;
	}

	// 2D
	else if (dim == 2)
	{
		std::vector<double> ret(coords.size());
		ret[0] = (coordsTranslated[0]) * _currentMacroscopicValue[0] + (coordsTranslated[1]) * _currentMacroscopicValue[2]; // Not 100% sure about the shear strains. Might need a 1/2
		ret[1] = (coordsTranslated[0]) * _currentMacroscopicValue[2] + (coordsTranslated[1]) * _currentMacroscopicValue[1];
		return ret;
	}

	// 3D
	else if (dim == 3)
	{
		std::vector<double> ret(coords.size());
		ret[0] = (coordsTranslated[0]) * _currentMacroscopicValue[0] + (coordsTranslated[1]) * _currentMacroscopicValue[5] + (coordsTranslated[2]) * _currentMacroscopicValue[4]; // Not 100% sure about the shear strains. Might need a 1/2
		ret[1] = (coordsTranslated[0]) * _currentMacroscopicValue[5] + (coordsTranslated[1]) * _currentMacroscopicValue[1] + (coordsTranslated[2]) * _currentMacroscopicValue[3];
		ret[2] = (coordsTranslated[0]) * _currentMacroscopicValue[4] + (coordsTranslated[1]) * _currentMacroscopicValue[3] + (coordsTranslated[2]) * _currentMacroscopicValue[2];
		return ret;
	}
	
	else
		err_message("Unknown dimension somehow");
}