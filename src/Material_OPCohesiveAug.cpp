/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated August 2017

##################################################################################
*/
#include "Material_OPCohesiveAug.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>


Material* OPCohesiveAugMaterial::allocate()
{
	return new OPCohesiveAugMaterial;
}


std::vector<std::string> OPCohesiveAugMaterial::internal_vars_name()
{
	std::vector<std::string> vec = {"maximum_opening"};
	return vec;
}

std::vector<bool> OPCohesiveAugMaterial::internal_vars_print()
{
	std::vector<bool> vec = {false};
	return vec;
}



Material::output_params* OPCohesiveAugMaterial::Constitutive(Material::input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Ortiz-Pandolfi Cohesive material.");

	// Get the internal variables and references to the output
	if((*input.internal_vars).size() != 1)
		err_message("Invalid number of internal variables for a Ortiz-Pandolfi Cohesive material.");
	double& delta_max = (*input.internal_vars)[0];
	double dmax = delta_max; // Just store a copy in case this interface is under compression
	DenseMatrix<double>& Dmat = _output.Dmat;
	std::vector<double>& traction_ttn = _output.traction;

	// Get the dimension of the problem and the input opening vector (int ttn coordinates hopefully)
	// ttn stands for tangential-tangential-normal
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;

	// Get the effective opening
	double delta_eff = getDeltaEff(delta_ttn);

	// Effective traction based on loading (Eq 28) or unloading (Eq 30)
	double trac_eff = 0;
	double Dtraction_Ddelta = 0;
	if(delta_eff>=delta_max || delta_eff<0.0) // Loading, continue along curve and update ISV
	{
		trac_eff = _sigma_c * (delta_eff/_delta_c) * exp(1.0-delta_eff/_delta_c);
		Dtraction_Ddelta = (_sigma_c/_delta_c) * (1.0-delta_eff/_delta_c) * exp(1.0-delta_eff/_delta_c);
		delta_max = delta_eff; // Update the internal state variable!
	}
	else	// Unloading. Unload to origin and don't update ISV
	{
		double max_traction = _sigma_c * (delta_max/_delta_c) * exp(1.0-delta_max/_delta_c);  // Traction at max opening
		trac_eff = max_traction * (delta_eff/delta_max);
		Dtraction_Ddelta = max_traction/delta_max;
	}

	// Actual traction vector
	std::vector<double> bracket_term(delta_ttn.size(), 0.0); // Term in brackets of Eq (43)
	traction_ttn.resize(delta_ttn.size());
	for(id_type i=0; i<delta_ttn.size(); ++i)
		bracket_term[i] = _beta2 * delta_ttn[i];
	bracket_term[dim-1] += (1.0-_beta2) * delta_ttn[dim-1]; // Add the second term of the bracketed expression (in the ttn coord system the dot product is just equal to the normal component of delta_ttn)
	for(id_type i=0; i<delta_ttn.size(); ++i)
		traction_ttn[i] = (trac_eff/delta_eff) * bracket_term[i];
	// If the normal component of opening vector is negative, overwrite the previous computation of the traction vector
	if (delta_ttn[dim-1] < 0)
		traction_ttn[dim-1] = (_sigma_c / _alpha) * exp(1.0 - dmax / _delta_c) * sinh(_alpha * delta_ttn[dim-1] / _delta_c);

	// Derivatives of delta_eff wrt delta_ttn
	std::vector<double> Ddeltaeff_Ddeltattn(dim);
	Ddeltaeff_Ddeltattn[dim-1] = delta_ttn[dim-1]/delta_eff;
	for(id_type i=0; i<(dim-1); i++)
		Ddeltaeff_Ddeltattn[i] = _beta2 * delta_ttn[i]/delta_eff;

	// Computation of the D matrix (Derivative of Eq 43 wrt delta_ttn)
	Dmat.resize(dim, dim);
	for(id_type i=0; i<dim; ++i)
		for(id_type k=0; k<dim; ++k)
			Dmat(i, k) = bracket_term[i] * (Dtraction_Ddelta*delta_eff - trac_eff)*(Ddeltaeff_Ddeltattn[k]/pow(delta_eff,2)) + 
						(trac_eff/delta_eff) * (_beta2*Utilities::kron(i,k)+(1.0-_beta2)*Utilities::kron(i,dim-1)*Utilities::kron(k,dim-1)); // NOTE: this is assuming that the normal vector is [0,0,1] (3D)
	// Again, if the normal component is actually negative, overwrite the compnents of the last row of the tangent matrix
	if (delta_ttn[dim-1] < 0)
	{
		for (id_type i=0; i<(dim-1); ++i)
			Dmat(dim-1, i) = 0.0;
		Dmat(dim-1, dim-1) = (_sigma_c / _delta_c) * exp(1.0 - dmax / _delta_c) * cosh(_alpha * delta_ttn[dim-1] / _delta_c);
	}

	// Return a pointer to the output structure
	return &_output;
}





Material::sensitivity_output_params* OPCohesiveAugMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// Get the internal variables and references to the output
	double& delta_max = (*input.internal_vars)[0];
	id_type param_id = input.parameter_id;
	double& ddelta_max_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars()];
	std::vector<double>& dtraction_ttn_dd = _sensitivity_output.dtraction_dd;

	// Get the effective opening
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;
	double delta_eff = getDeltaEff(delta_ttn);

	// Traction vector partial derivatives
	dtraction_ttn_dd.clear();
	dtraction_ttn_dd.resize(delta_ttn.size());

	// Used a lot so I'll just compute it here
	std::vector<double> bracket_term(delta_ttn.size(), delta_ttn[dim-1]); // Term in brackets of Eq (43)
	for(id_type i=0; i<(dim-1); ++i)
		bracket_term[i] = _beta2 * delta_ttn[i];

	// just store boolean here for ease of use
	bool loading = (delta_eff >= delta_max);

	// Partial derivative is only non-zero if the material I am taing the derivative wrt is the same as the current material I am in
	if (input.sensitivity_mat_name == _name) 
	{
		// Sigma_c sensitivity
		if (input.sensitivity_param_name == "SIGMA_C")
		{
			// 3b PARTIAL derivative
			double dtrac_eff_dsigmac;
			if (loading)
				dtrac_eff_dsigmac = delta_eff/_delta_c*exp(1-delta_eff/_delta_c);
			else
				dtrac_eff_dsigmac = delta_eff/_delta_c*exp(1-delta_max/_delta_c);

			for(id_type i=0; i<delta_ttn.size(); ++i)
				dtraction_ttn_dd[i] = (dtrac_eff_dsigmac/delta_eff) * bracket_term[i];
			if (delta_ttn[dim-1] < 0)
				dtraction_ttn_dd[dim-1] = (1.0 / _alpha) * exp(1.0 - delta_max / _delta_c) * sinh(_alpha * delta_ttn[dim-1] / _delta_c);
		}

		// Delta_c sensitivity
		else if (input.sensitivity_param_name == "DELTA_C")
		{
			// 3c PARTIAL derivative
			double dtrac_eff_ddeltac;
			if (loading)
				dtrac_eff_ddeltac = (_sigma_c*delta_eff/(_delta_c*_delta_c)*exp(1-delta_eff/_delta_c)) * (delta_eff/_delta_c - 1);
			else
				dtrac_eff_ddeltac = (_sigma_c*delta_eff/(_delta_c*_delta_c)*exp(1-delta_max/_delta_c)) * (delta_max/_delta_c - 1);

			for(id_type i=0; i<delta_ttn.size(); ++i)
				dtraction_ttn_dd[i] = (dtrac_eff_ddeltac/delta_eff) * bracket_term[i];
			if (delta_ttn[dim-1] < 0)
				dtraction_ttn_dd[dim-1] = (_sigma_c / _alpha) * exp(1.0 - delta_max / _delta_c) * ( sinh(_alpha * delta_ttn[dim-1] / _delta_c) * (delta_max / pow(_delta_c,2)) + 
																									cosh(_alpha * delta_ttn[dim-1] / _delta_c) * (-_alpha * delta_ttn[dim-1] / pow(_delta_c,2)) );
		}

		// Beta sensitivity
		else if (input.sensitivity_param_name == "BETA")
		{
			err_message("Ortiz-Pandolphi Beta sensitivity has not been derived yet!");
		}
	}

	// Total derivative wrt the internal state variable term is present regardless of whether or not I'm in the same material
	if (!loading)
	{
		double dtrac_eff_ddeltamax = (_sigma_c * delta_eff / _delta_c) * exp(1-delta_max/_delta_c) * (-1.0/_delta_c);
		for(id_type i=0; i<(dim-1); ++i)
			dtraction_ttn_dd[i] += ddelta_max_dd * (dtrac_eff_ddeltamax/delta_eff) * bracket_term[i];
		if (delta_ttn[dim-1] >= 0.0)
			dtraction_ttn_dd[dim-1] += ddelta_max_dd * (dtrac_eff_ddeltamax/delta_eff) * bracket_term[dim-1];
		else
			dtraction_ttn_dd[dim-1] += ddelta_max_dd * ((_sigma_c / _alpha) * exp(1.0 - delta_max / _delta_c) * sinh(_alpha * delta_ttn[dim-1] / _delta_c) * (-1.0 / _delta_c));
	}

	return &_sensitivity_output;
}



void OPCohesiveAugMaterial::updateSensitivityISVs(Material::sensitivity_input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// Get the internal variables and references to the output
	double& delta_max = (*input.internal_vars)[0];
	id_type param_id = input.parameter_id;
	double& ddelta_max_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars()];

	// Get the effective opening
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;
	double delta_eff = getDeltaEff(delta_ttn);

	// Loading, continue along curve and update ISV
	if(delta_eff >= delta_max)
	{
		std::vector<double>& ddelta_ttn_dd = input.ddelta_dd;

		// 2b
		double ddelta_eff_dd = 2.0*delta_ttn[dim-1]*ddelta_ttn_dd[dim-1];
		for (id_type i=0; i<(dim-1); ++i)
			ddelta_eff_dd += 2.0*_beta2*delta_ttn[i]*ddelta_ttn_dd[i];
		ddelta_eff_dd *= 1.0/(2.0*delta_eff);

		// Update the sensitivity internal state variables!
		ddelta_max_dd = ddelta_eff_dd;
	}
}