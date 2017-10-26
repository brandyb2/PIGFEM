/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated August 2017

##################################################################################
*/
#include "AssemblerSmallStrainStructural_ThermalPreLoad.h"
#include "Problem.h"
#include "Utilities.h"
#include <algorithm>



AssemblerSmallStrainStructural_ThermalPreLoad::AssemblerSmallStrainStructural_ThermalPreLoad()
	: _delta_T(0.0), _curr_delta_T(0.0), _temp_set(false), _relaxation_complete(false)
{
}

AssemblerSmallStrainStructural_ThermalPreLoad::~AssemblerSmallStrainStructural_ThermalPreLoad()
{
}
void AssemblerSmallStrainStructural_ThermalPreLoad::clear()
{
	_delta_T = 0.0;
	_curr_delta_T = 0.0;
	_temp_set = false;
	_relaxation_complete = false;
	Assembler::clear();
}



// Kernel function that actually does the computes the physics behind the specific problem
// Computes the local contribution to the stiffness matrix and internal load vector
// Also updates the internal variable storage in the input variable
//	B-matrix variant
void AssemblerSmallStrainStructural_ThermalPreLoad::KernelVolumetric(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
													const std::vector<double>& shape, const DenseMatrix<double>& B,
													Material* mat, Material::input_params& input,
													std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	std::vector<double>& strain = input.strain;
	ProblemUtilities::B_U(strain, B, elem_U_curr);

	// Add the current macroscopic strain to it
	if (_prob->subscale() && _relaxation_complete)
	{
		if (input.strain.size() == _currentMacroscopicValue.size())
			Utilities::_VecAXPY(1.0, _currentMacroscopicValue, input.strain);
		else
			err_message("Invalid macroscopic value applied");
	}

	// Compute the constitutive relations
	input.temp_change = _curr_delta_T; // Inform the material that there is an applied temperature change
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::Bt_D_B(K_el, B, output->Dmat);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::Bt_F(P_el_int, B, output->stress);
}

void AssemblerSmallStrainStructural_ThermalPreLoad::KernelVolumetric(DenseMatrix<double>& K_el, std::vector<double>& P_el_int,
													const std::vector<double>& shape, const std::vector<std::vector<double> >& shape_grad,
													Material* mat, Material::input_params& input,
													std::vector<double>& elem_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current strain
	std::vector<double>& strain = input.strain;
	ProblemUtilities::SmallStrainFast(strain, shape_grad, elem_U_curr);

	// Add the current macroscopic strain to it
	if (_prob->subscale() && _relaxation_complete)
	{
		if (input.strain.size() == _currentMacroscopicValue.size())
			Utilities::_VecAXPY(1.0, _currentMacroscopicValue, input.strain);
		else
			err_message("Invalid macroscopic value applied");
	}

	// Compute the constitutive relations
	input.temp_change = _curr_delta_T; // Inform the material that there is an applied temperature change
	Material::output_params* output = mat->Constitutive(input);

	// Compute the stiffness matrix
	if (assembleJac)
		ProblemUtilities::SmallStrainKFast(K_el, shape_grad, output->Dmat);

	// Compute the internal force
	if (assembleFunc)
		ProblemUtilities::SmallStrainInternalForceFast(P_el_int, shape_grad, output->stress);
}




void AssemblerSmallStrainStructural_ThermalPreLoad::KernelCohesive(DenseMatrix<double>& K_coh, std::vector<double>& P_int_coh,
															 const std::vector<double>& shape, const DenseMatrix<double>& rotation_matrix,
															 Material* mat, Material::input_params& input,
															 std::vector<double>& coh_U_curr, bool assembleFunc, bool assembleJac)
{
	// Compute the current opening vector in the ntt coordinate system
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





bool AssemblerSmallStrainStructural_ThermalPreLoad::storedBmats()
{
	return true;
}

// Select the appropriate B-matrix fill function
void AssemblerSmallStrainStructural_ThermalPreLoad::fillBmat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& shape_grad)
{
	ProblemUtilities::Assemble_Small_Strain_B_Mat(B, shape_grad);
}












// FIXME: I DON'T ACTUALLY KNOW IF THIS IS CORRECT.
// HOW DO BOUNDARY CONDITIONS WORK IN A PERIODIC MEDIA FOR A MATERIAL THAT HAS BEEN PRESTRAINED???
// DO THE SOLUTION CONTRIBUTIONS NEED TO BE BASED ON THE DEFORMED CONFIGURATION HERE?
std::vector<double> AssemblerSmallStrainStructural_ThermalPreLoad::computeMacroscaleSolutionContribution(const std::vector<double>& coords)
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















// Function to assemble the presribed displacements and free loads at the beginning of each load step
// (Or at the beginning of the analysis for linear)
// This is made virtual so that calling it can update the assembler state at every load step is desired
PetscErrorCode AssemblerSmallStrainStructural_ThermalPreLoad::assemble_new_load_step(Vec& Up, pvector& F_ext, double current_time_frac)
{
	// Make sure I've assigned a strain value
	if (!_temp_set)
		err_message("Please set a temperature change value before attemping to solve a residual strain problem.");

	// Update the current temp here
	if (_relaxation_complete)
		_curr_delta_T = _delta_T;
	else
		_curr_delta_T = _delta_T * current_time_frac;

	// Call the base class function to do all of the same things
	return Assembler::assemble_new_load_step(Up, F_ext, current_time_frac);
}



/*
 * Function used to easily set a scalar-valued parameter
 */
void AssemblerSmallStrainStructural_ThermalPreLoad::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="TEMP" || name=="DELTA_T" || name=="DELTA_TEMP")
	{
		_delta_T = val;
		_curr_delta_T = 0.0;
		_temp_set = true;
	}
	else if (name=="COMPLETE" || name=="CONVERGED")
	{
		if (val==0.0) // stupid hack because virtual functions can't be templated
			_relaxation_complete = false;
		else
			_relaxation_complete = true;
	}
	else
		err_message(name << " is not a valid parameter name");
}
double AssemblerSmallStrainStructural_ThermalPreLoad::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="TEMP" || name=="DELTA_T" || name=="DELTA_TEMP")
	{
		if (_temp_set)
			return _delta_T;
		else
			err_message("Please set a temperature change value before attempting to get it.");
	}
	else if (name=="COMPLETE" || name=="CONVERGED")
	{
		if (_relaxation_complete) // stupid hack because virtual functions can't be templated
			return 1.0;
		else
			return 0.0;
	}
	else if (name=="CURRENT TEMP" || name=="CURRENT_TEMP")
		return _curr_delta_T;
	else
		err_message(name << " is not a valid parameter name");
}