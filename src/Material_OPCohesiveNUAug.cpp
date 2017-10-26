/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated August 2017

##################################################################################
*/
#include "Material_OPCohesiveNUAug.h"
#include "SensitivityMaterialParameterStructural.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>


Material* OPCohesiveNUAugMaterial::allocate()
{
	return new OPCohesiveNUAugMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void OPCohesiveNUAugMaterial::copy(Material* other_mat)
{
	other_mat->set_name(_name);
	other_mat->set_id(_id);
	if(_s_set)
		other_mat->set_parameter("SIGMA", _sigma_c);
	if(_d_set)
		other_mat->set_parameter("DELTA", _delta_c);
	if(_b_set)
		other_mat->set_parameter("BETA", sqrt(_beta2));
}


OPCohesiveNUAugMaterial::OPCohesiveNUAugMaterial()
	: _sigma_c(0), _delta_c(0), _beta2(0), _alpha(50),
	  _s_set(false), _d_set(false), _b_set(false), ready(false)
{
}

double OPCohesiveNUAugMaterial::getDeltaEff(std::vector<double>& delta_ttn)
{
	// Replacing delta_ttn with a small number to prevent blowup at the beginning of loading
	double norm = 0;
	for(id_type i=0; i<delta_ttn.size(); ++i)
		norm += delta_ttn[i]*delta_ttn[i];
	if(norm < 1e-40) // Since this is norm squared
	{
		for (id_type i=0; i< delta_ttn.size(); ++i)
			delta_ttn[i] = 1e-20;
	}

	// Effective opening diplacement (Eq 23, Ortiz-Pandolfi (1999))
	double delta_eff = pow(delta_ttn[delta_ttn.size()-1], 2); // Normal contribution
	for(id_type i=0; i<(delta_ttn.size()-1); i++)			 // Tangential contributions
		delta_eff += _beta2*pow(delta_ttn[i], 2);
	delta_eff = sqrt(delta_eff);
	return delta_eff;
}


Material::output_params* OPCohesiveNUAugMaterial::Constitutive(Material::input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Ortiz-Pandolfi No Unloading Cohesive material.");

	// Get the internal variables and references to the output
	DenseMatrix<double>& Dmat = _output.Dmat;
	std::vector<double>& traction_ttn = _output.traction;

	// Get the dimension of the problem and the input opening vector (int ttn coordinates hopefully)
	// ttn stands for normal-tangential-tangential
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;

	// Get the effective opening
	double delta_eff = getDeltaEff(delta_ttn);

	// Effective traction based on loading (Eq 28) or unloading (Eq 30)
	double trac_eff = 0;
	double Dtraction_Ddelta = 0;
	trac_eff = _sigma_c * (delta_eff/_delta_c) * exp(1.0-delta_eff/_delta_c);
	Dtraction_Ddelta = (_sigma_c/_delta_c) * (1.0-delta_eff/_delta_c) * exp(1.0-delta_eff/_delta_c);

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
		traction_ttn[dim-1] = (_sigma_c * exp(1.0) / _alpha) * sinh(_alpha * delta_ttn[dim-1] / _delta_c);

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
		Dmat(dim-1, dim-1) = (_sigma_c * exp(1.0) / _delta_c) * cosh(_alpha * delta_ttn[dim-1] / _delta_c);
	}

	// Return a pointer to the output structure
	return &_output;
}





Material::sensitivity_output_params* OPCohesiveNUAugMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Ortiz-Pandolfi No Unloading Cohesive material.");

	// Get the internal variables and references to the output
	std::vector<double>& dtraction_ttn_dd = _sensitivity_output.dtraction_dd;

	// Get the dimension of the problem and the input opening vector (int ttn coordinates)
	// ttn stands for normal-tangential-tangential
	id_type dim = input.dim;
	std::vector<double>& delta_ttn = input.delta;

	// Traction vector partial derivatives
	dtraction_ttn_dd.clear();
	dtraction_ttn_dd.resize(delta_ttn.size());

	// Partial derivative is only non-zero if the material I am taing the derivative wrt is the same as the current material I am in
	if (input.sensitivity_mat_name == _name) 
	{
		// Get the effective opening
		double delta_eff = getDeltaEff(delta_ttn);

		std::vector<double> bracket_term(delta_ttn.size(), delta_ttn[dim-1]); // Term in brackets of Eq (43)
		for(id_type i=0; i<(dim-1); ++i)
			bracket_term[i] = _beta2 * delta_ttn[i];


		// Sigma_c sensitivity
		if (input.sensitivity_param_name == "SIGMA_C")
		{
			// 3b PARTIAL derivative
			double dtrac_eff_dsigmac = delta_eff/_delta_c*exp(1-delta_eff/_delta_c);

			for(id_type i=0; i<delta_ttn.size(); ++i)
				dtraction_ttn_dd[i] = (dtrac_eff_dsigmac/delta_eff) * bracket_term[i];
			if (delta_ttn[dim-1] < 0)
				dtraction_ttn_dd[dim-1] = (exp(1.0) / _alpha) * sinh(_alpha * delta_ttn[dim-1] / _delta_c);
		}

		// Delta_c sensitivity
		else if (input.sensitivity_param_name == "DELTA_C")
		{
			// 3c PARTIAL derivative
			double dtrac_eff_ddeltac = (_sigma_c*delta_eff/_delta_c/_delta_c*exp(1-delta_eff/_delta_c)) * (delta_eff/_delta_c - 1);

			for(id_type i=0; i<delta_ttn.size(); ++i)
				dtraction_ttn_dd[i] = (dtrac_eff_ddeltac/delta_eff) * bracket_term[i];
			if (delta_ttn[dim-1] < 0)
				dtraction_ttn_dd[dim-1] = (_sigma_c * exp(1.0) / _alpha) * cosh(_alpha * delta_ttn[dim-1] / _delta_c) * (-_alpha * delta_ttn[dim-1] / pow(_delta_c,2));
		}

		// Beta sensitivity
		else if (input.sensitivity_param_name == "BETA")
		{
			err_message("Ortiz-Pandolphi Beta sensitivity has not been derived yet!");
		}
	}



	// Return a pointer to the output structure
	return &_sensitivity_output;
}









SensitivityMaterialParameter* OPCohesiveNUAugMaterial::getSensitivityParameter(std::string param_name)
{
	SensitivityMaterialParameter* ret = new SensitivityMaterialParameterStructural;
	ret->set_mat_name(_name);
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	if(param_name=="SIGMA" || param_name=="SIGMA_C" || param_name=="SIGMA C" || param_name=="SIGMA_CRITICAL" || param_name=="SIGMA CRITICAL")
		ret->set_param_name("SIGMA_C");
	else if(param_name=="DELTA" || param_name=="DELTA_C" || param_name=="DELTA C" || param_name=="DELTA_CRITICAL" || param_name=="DELTA CRITICAL")
		ret->set_param_name("DELTA_C");
	else if(param_name=="BETA")
		ret->set_param_name("BETA");
	else
	{
		delete ret;
		err_message(param_name << " is not a valid parameter name");
	}

	return ret;
}












double OPCohesiveNUAugMaterial::getNormalizedCohesiveFailure(std::vector<double>& delta_ttn)
{
	double delta_eff = getDeltaEff(delta_ttn);
	return (delta_eff / _delta_c);
}
id_type OPCohesiveNUAugMaterial::get_cohesive_damage_zone(std::vector<double>& delta)
{
	double norm = getNormalizedCohesiveFailure( delta );

	if (norm <= 0)
		return 0;
	else if (norm <= 1)
		return 1;
	else if (norm <= 5.0) // Kinda arbitrary selection
		return 2;
	else
		return 3;
}




void OPCohesiveNUAugMaterial::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="SIGMA" || name=="SIGMA_C" || name=="SIGMA C" || name=="SIGMA_CRITICAL" || name=="SIGMA CRITICAL")
		set_sc(val);
	else if(name=="DELTA" || name=="DELTA_C" || name=="DELTA C" || name=="DELTA_CRITICAL" || name=="DELTA CRITICAL")
		set_dc(val);
	else if(name=="BETA")
		set_b(val);
	else if(name=="ALPHA")
		set_alpha(val);
	else
		err_message(name << " is not a valid parameter name");
}
double OPCohesiveNUAugMaterial::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="SIGMA" || name=="SIGMA_C" || name=="SIGMA C" || name=="SIGMA_CRITICAL" || name=="SIGMA CRITICAL")
		return get_sc();
	else if(name=="DELTA" || name=="DELTA_C" || name=="DELTA C" || name=="DELTA_CRITICAL" || name=="DELTA CRITICAL")
		return get_dc();
	else if(name=="BETA")
		return get_b();
	else if(name=="ALPHA")
		return _alpha;
	else
		err_message(name << " is not a valid parameter name");
}




void OPCohesiveNUAugMaterial::set_sc(double sc)
{
	_sigma_c = sc;
	_s_set = true;
	if(_s_set && _d_set && _b_set)
		ready = true;
}
double OPCohesiveNUAugMaterial::get_sc()
{
	if(_s_set)
		return _sigma_c;
	else
		err_message("Sigma critical could not be determined.");
}




void OPCohesiveNUAugMaterial::set_dc(double dc)
{
	_delta_c = dc;
	_d_set = true;
	if(_s_set && _d_set && _b_set)
		ready = true;
}
double OPCohesiveNUAugMaterial::get_dc()
{
	if(_d_set)
		return _delta_c;
	else
		err_message("Delta critical could not be determined.");
}




void OPCohesiveNUAugMaterial::set_b(double b)
{
	_beta2 = b*b;
	_b_set = true;
	if(_s_set && _d_set && _b_set)
		ready = true;
}
double OPCohesiveNUAugMaterial::get_b()
{
	if(_b_set)
		return sqrt(_beta2);
	else
		err_message("Beta could not be determined.");
}

void OPCohesiveNUAugMaterial::set_alpha(double a)
{
	if (a < 0.)
		err_message("Invalid value for alpha, must be positive!");

	_alpha = a;
}
