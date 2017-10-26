/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated October 2017

##################################################################################
*/
#include "Material_SmoothedTrilinearCohesiveNU.h"
#include "SensitivityMaterialParameterStructural.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>


Material* TrilinearSmoothedCohesiveNUMaterial::allocate()
{
	return new TrilinearSmoothedCohesiveNUMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void TrilinearSmoothedCohesiveNUMaterial::copy(Material* other_mat)
{
	other_mat->set_name(_name);
	other_mat->set_id(_id);
	if (_s_set)
		other_mat->set_parameter("SIGMA", _sigma_c);
	if (_d1_set)
		other_mat->set_parameter("DELTA1", _delta_c1);
	if (_d2_set)
		other_mat->set_parameter("DELTA2", _delta_c2);
	if (_d3_set)
		other_mat->set_parameter("DELTA3", _delta_c3);
	if (_G_set)
		other_mat->set_parameter("G", _G);

	other_mat->set_parameter("BETA", sqrt(_beta2));
	other_mat->set_parameter("ALPHA", sqrt(_alpha));
}


TrilinearSmoothedCohesiveNUMaterial::TrilinearSmoothedCohesiveNUMaterial()
	: _sigma_c(0.0), _delta_c1(0.0), _delta_c2(0.0), _delta_c3(0.0), _G(0.0), _beta2(1.0), _alpha(50.0),
	  _s_set(false), _d1_set(false), _d2_set(false), _d3_set(false), _G_set(false), ready(false)
{
}

double TrilinearSmoothedCohesiveNUMaterial::getDeltaEff(std::vector<double>& delta_ttn)
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


Material::output_params* TrilinearSmoothedCohesiveNUMaterial::Constitutive(Material::input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Trilinear No Unloading Cohesive material.");

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
	if (delta_eff >= 0 && delta_eff < _delta_c1)
	{
		double val = delta_eff / _delta_c1;
		trac_eff = _sigma_c * ( 2.0*val - pow(val, 2) );
		Dtraction_Ddelta = _sigma_c * ( (2.0/_delta_c1) * (1.0 - val) );
	}
	else if (delta_eff >= _delta_c1 && delta_eff<_delta_c2)
	{
		trac_eff = _sigma_c;
		Dtraction_Ddelta = 0.0;
	}
	else if (delta_eff >= _delta_c2 && delta_eff<_delta_c3)
	{
		double val = (delta_eff-_delta_c2) / (_delta_c3 - _delta_c2);
		trac_eff = _sigma_c * ( 2.0*pow(val, 3) - 3.0*pow(val, 2) + 1.0 );
		Dtraction_Ddelta = _sigma_c * ( 6.0 * (val/(_delta_c3-_delta_c2)) * (val-1.0) );
	}
	else
	{
		trac_eff = 0;
		Dtraction_Ddelta = 0;
	}

	// Actual traction vector
	std::vector<double> bracket_term(delta_ttn.size(), 0.0); // Term in brackets of Eq (43)
	traction_ttn.resize(delta_ttn.size());
	for(id_type i=0; i<delta_ttn.size(); ++i)
		bracket_term[i] = _beta2 * delta_ttn[i];
	bracket_term[dim-1] += (1.0-_beta2) * delta_ttn[dim-1]; // Add the second term of the bracketed expression (in the ttn coord system the dot product is just equal to the normal component of delta_ttn)
	for(id_type i=0; i<delta_ttn.size(); ++i)
		traction_ttn[i] = (trac_eff/delta_eff) * bracket_term[i];

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
		double c = 2.0*_sigma_c/_alpha;

		traction_ttn[dim-1] = c * sinh(_alpha * delta_ttn[dim-1] / _delta_c1);
		for (id_type i=0; i<(dim-1); ++i)
			Dmat(dim-1, i) = 0.0;
		Dmat(dim-1, dim-1) = c * (_alpha/_delta_c1) * cosh(_alpha * delta_ttn[dim-1] / _delta_c1);
	}

	// Return a pointer to the output structure
	return &_output;
}





Material::sensitivity_output_params* TrilinearSmoothedCohesiveNUMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	err_message("Sensitivity is not implemented for the bilinear cohesive model");
	// if(!ready)
	// 	err_message("Not all material parameters have been set for the Ortiz-Pandolfi No Unloading Cohesive material.");

	// // Get the internal variables and references to the output
	// std::vector<double>& dtraction_ttn_dd = _sensitivity_output.dtraction_dd;

	// // Get the dimension of the problem and the input opening vector (int ttn coordinates)
	// // ttn stands for normal-tangential-tangential
	// id_type dim = input.dim;
	// std::vector<double>& delta_ttn = input.delta;

	// // Traction vector partial derivatives
	// dtraction_ttn_dd.clear();
	// dtraction_ttn_dd.resize(delta_ttn.size());

	// // Partial derivative is only non-zero if the material I am taing the derivative wrt is the same as the current material I am in
	// if (input.sensitivity_mat_name == _name) 
	// {
	// 	// Get the effective opening
	// 	double delta_eff = getDeltaEff(delta_ttn);

	// 	std::vector<double> bracket_term(delta_ttn.size(), delta_ttn[dim-1]); // Term in brackets of Eq (43)
	// 	for(id_type i=0; i<(dim-1); ++i)
	// 		bracket_term[i] = _beta2 * delta_ttn[i];


	// 	// Sigma_c sensitivity
	// 	if (input.sensitivity_param_name == "SIGMA_C")
	// 	{
	// 		// 3b PARTIAL derivative
	// 		double dtrac_eff_dsigmac = delta_eff/_delta_c*exp(1-delta_eff/_delta_c);

	// 		for(id_type i=0; i<delta_ttn.size(); ++i)
	// 			dtraction_ttn_dd[i] = (dtrac_eff_dsigmac/delta_eff) * bracket_term[i];
	// 		if (delta_ttn[dim-1] < 0)
	// 			dtraction_ttn_dd[dim-1] = (exp(1.0) / _alpha) * sinh(_alpha * delta_ttn[dim-1] / _delta_c);
	// 	}

	// 	// Delta_c sensitivity
	// 	else if (input.sensitivity_param_name == "DELTA_C")
	// 	{
	// 		// 3c PARTIAL derivative
	// 		double dtrac_eff_ddeltac = (_sigma_c*delta_eff/_delta_c/_delta_c*exp(1-delta_eff/_delta_c)) * (delta_eff/_delta_c - 1);

	// 		for(id_type i=0; i<delta_ttn.size(); ++i)
	// 			dtraction_ttn_dd[i] = (dtrac_eff_ddeltac/delta_eff) * bracket_term[i];
	// 		if (delta_ttn[dim-1] < 0)
	// 			dtraction_ttn_dd[dim-1] = (_sigma_c * exp(1.0) / _alpha) * cosh(_alpha * delta_ttn[dim-1] / _delta_c) * (-_alpha * delta_ttn[dim-1] / pow(_delta_c,2));
	// 	}

	// 	// Beta sensitivity
	// 	else if (input.sensitivity_param_name == "BETA")
	// 	{
	// 		err_message("Ortiz-Pandolphi Beta sensitivity has not been derived yet!");
	// 	}
	// }



	// // Return a pointer to the output structure
	// return &_sensitivity_output;
}









SensitivityMaterialParameter* TrilinearSmoothedCohesiveNUMaterial::getSensitivityParameter(std::string param_name)
{
	err_message("Sensitivity is not implemented for the bilinear cohesive model");
	// SensitivityMaterialParameter* ret = new SensitivityMaterialParameterStructural;
	// ret->set_mat_name(_name);
	// std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	// if(param_name=="SIGMA" || param_name=="SIGMA_C" || param_name=="SIGMA C" || param_name=="SIGMA_CRITICAL" || param_name=="SIGMA CRITICAL")
	// 	ret->set_param_name("SIGMA_C");
	// else if(param_name=="DELTA" || param_name=="DELTA_C" || param_name=="DELTA C" || param_name=="DELTA_CRITICAL" || param_name=="DELTA CRITICAL")
	// 	ret->set_param_name("DELTA_C1");
	// else if(param_name=="BETA")
	// 	ret->set_param_name("BETA");
	// else
	// {
	// 	delete ret;
	// 	err_message(param_name << " is not a valid parameter name");
	// }

	// return ret;
}












double TrilinearSmoothedCohesiveNUMaterial::getNormalizedCohesiveFailure(std::vector<double>& delta_ttn)
{
	double delta_eff = getDeltaEff(delta_ttn);
	return (delta_eff / _delta_c1);
}
id_type TrilinearSmoothedCohesiveNUMaterial::get_cohesive_damage_zone(std::vector<double>& delta)
{
	double delta_eff = getDeltaEff(delta);
	if (delta[delta.size()-1] < 0)
		return 0;
	else if (delta_eff <= _delta_c1)
		return 1;
	else if (delta_eff <= _delta_c3)	// I should probably extrend this to an additions region. Oh well
		return 2;
	else
		return 3;
}










































void TrilinearSmoothedCohesiveNUMaterial::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="SIGMA" || name=="SIGMA_C" || name=="SIGMA C" || name=="SIGMA_CRITICAL" || name=="SIGMA CRITICAL")
		set_sc(val);
	else if (name=="DELTA1" || name=="DELTA_C1" || name=="DELTA C1" || name=="DELTA_CRITICAL" || name=="DELTA CRITICAL")
		set_dc1(val);
	else if (name=="DELTA2" || name=="DELTA_C2" || name=="DELTA C2")
		set_dc2(val);
	else if (name=="DELTA3" || name=="DELTA_C3" || name=="DELTA C3" || name=="DELTA_FINAL" || name=="DELTA FINAL")
		set_dc3(val);
	else if (name=="G" || name=="J" || name=="FRACTURE_ENERGY" || name=="FRACTURE ENERGY")
		set_G(val);
	else if (name=="BETA")
		set_b(val);
	else if (name=="ALPHA")
		set_alpha(val);
	else
		err_message(name << " is not a valid parameter name");
}
double TrilinearSmoothedCohesiveNUMaterial::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if (name=="SIGMA" || name=="SIGMA_C" || name=="SIGMA C" || name=="SIGMA_CRITICAL" || name=="SIGMA CRITICAL")
		return get_sc();
	else if (name=="DELTA1" || name=="DELTA_C1" || name=="DELTA C1" || name=="DELTA_CRITICAL" || name=="DELTA CRITICAL")
		return get_dc1();
	else if (name=="DELTA2" || name=="DELTA_C2" || name=="DELTA C2")
		return get_dc2();
	else if (name=="DELTA3" || name=="DELTA_C3" || name=="DELTA C3" || name=="DELTA_FINAL" || name=="DELTA FINAL")
		return get_dc3();
	else if (name=="G" || name=="J" || name=="FRACTURE_ENERGY" || name=="FRACTURE ENERGY")
		return get_G();
	else if (name=="BETA")
		return sqrt(_beta2);
	else if (name=="ALPHA")
		return _alpha;
	else
		err_message(name << " is not a valid parameter name");
}




void TrilinearSmoothedCohesiveNUMaterial::set_dc1(double dc)
{
	if (_d2_set && _d3_set && _s_set && _G_set)
		err_message("Can only specify 4 parameters for a trilinear cohesive model");
	else{
		_delta_c1 = dc;
		_d1_set = true;

		if (_delta_c1 <= 0.0)
			err_message("Incompatible parameter delta_c1 added to trilinear cohesive law");

		if (_d2_set && (dc > _delta_c2))
			err_message("Incompatible parameter delta_c1 added to trilinear cohesive law");
		if (_d3_set && (dc > _delta_c3))
			err_message("Incompatible parameter delta_c1 added to trilinear cohesive law");

		if (_d2_set && _d3_set && _s_set){
			ready = true;
			_G = get_G();
		}
		else if (_d2_set && _s_set && _G_set){
			ready = true;
			_delta_c3 = get_dc3();
			if (_delta_c3 < _delta_c2)
				err_message("Incompatible parameter delta_c1 added to trilinear cohesive law");
		}
		else if (_d3_set && _s_set && _G_set){
			ready = true;
			_delta_c2 = get_dc2();
			if (_delta_c2 < _delta_c1 || _delta_c2 > _delta_c3)
				err_message("Incompatible parameter delta_c1 added to trilinear cohesive law");
		}
		else if (_d2_set && _d3_set && _G_set){
			ready = true;
			_sigma_c = get_sc();
		}
	}
}
double TrilinearSmoothedCohesiveNUMaterial::get_dc1()
{
	if (_d2_set)
		return _delta_c2;
	else if (ready)     // The other 3 parameters are specified. Compute this one from them{
		return 3.0 * (_delta_c2/2.0 + _delta_c3/2.0 - _G/_sigma_c);
	else
		err_message("Peak opening displacement could not be determined.");
}




void TrilinearSmoothedCohesiveNUMaterial::set_dc2(double dc)
{
	if (_d1_set && _d3_set && _s_set && _G_set)
		err_message("Can only specify 4 parameters for a trilinear cohesive model");
	else{
		_delta_c2 = dc;
		_d2_set = true;

		if (_delta_c2 <= 0.0)
			err_message("Incompatible parameter delta_c2 added to trilinear cohesive law");

		if (_d1_set && (dc < _delta_c1))
			err_message("Incompatible parameter delta_c2 added to trilinear cohesive law");
		if (_d3_set && (dc > _delta_c3))
			err_message("Incompatible parameter delta_c2 added to trilinear cohesive law");

		if (_d1_set && _d3_set && _s_set){
			ready = true;
			_G = get_G();
		}
		else if (_d1_set && _s_set && _G_set){
			ready = true;
			_delta_c3 = get_dc3();
			if (_delta_c3 < _delta_c2)
				err_message("Incompatible parameter delta_c2 added to trilinear cohesive law");
		}
		else if (_d3_set && _s_set && _G_set){
			ready = true;
			_delta_c1 = get_dc1();
			if (_delta_c1 > _delta_c2 || _delta_c1 <= 0.0)
				err_message("Incompatible parameter delta_c2 added to trilinear cohesive law");
		}
		else if (_d1_set && _d3_set && _G_set){
			ready = true;
			_sigma_c = get_sc();
		}
	}
}
double TrilinearSmoothedCohesiveNUMaterial::get_dc2()
{
	if (_d2_set)
		return _delta_c2;
	else if (ready)     // The other 3 parameters are specified. Compute this one from them{
		return 2.0 * (_G/_sigma_c + _delta_c1/3.0 - _delta_c3/2.0);
	else
		err_message("Final plateau opening displacement could not be determined.");
}




void TrilinearSmoothedCohesiveNUMaterial::set_dc3(double dc)
{
	if (_d1_set && _d2_set && _s_set && _G_set)
		err_message("Can only specify 4 parameters for a trilinear cohesive model");
	else{
		_delta_c3 = dc;
		_d3_set = true;

		if (_delta_c3 <= 0.0)
			err_message("Incompatible parameter delta_c3 added to trilinear cohesive law");

		if (_d1_set && (dc < _delta_c1))
			err_message("Incompatible parameter delta_c3 added to trilinear cohesive law");
		if (_d2_set && (dc < _delta_c2))
			err_message("Incompatible parameter delta_c3 added to trilinear cohesive law");

		if (_d1_set && _d2_set && _s_set){
			ready = true;
			_G = get_G();
		}
		else if (_d1_set && _s_set && _G_set){
			ready = true;
			_delta_c2 = get_dc2();
			if (_delta_c2 < _delta_c1 || _delta_c2 > _delta_c3)
				err_message("Incompatible parameter delta_c3 added to trilinear cohesive law");
		}
		else if (_d2_set && _s_set && _G_set){
			ready = true;
			_delta_c1 = get_dc1();
			if (_delta_c1 > _delta_c2 || _delta_c1 <= 0.0)
				err_message("Incompatible parameter delta_c3 added to trilinear cohesive law");
		}
		else if (_d1_set && _d2_set && _G_set){
			ready = true;
			_sigma_c = get_sc();
		}
	}
}
double TrilinearSmoothedCohesiveNUMaterial::get_dc3()
{
	if (_d3_set)
		return _delta_c3;
	else if (ready)     // The other 3 parameters are specified. Compute this one from them{
		return 2.0 * (_G/_sigma_c + _delta_c1/3.0 - _delta_c2/2.0);
	else
		err_message("Final opening displacement could not be determined.");
}




void TrilinearSmoothedCohesiveNUMaterial::set_sc(double sc)
{
	if (_d1_set && _d2_set && _d3_set && _G_set)
		err_message("Can only specify 4 parameters for a trilinear cohesive model");
	else{
		_sigma_c = sc;
		_s_set = true;

		if (_sigma_c <= 0.0)
			err_message("Incompatible parameter sigma_c added to trilinear cohesive law");

		if (_d1_set && _d2_set && _d3_set){
			ready = true;
			_G = get_G();
		}
		else if (_d1_set && _d3_set && _G_set){
			ready = true;
			_delta_c2 = get_dc2();
			if (_delta_c2 < _delta_c1 || _delta_c2 > _delta_c3)
				err_message("Incompatible parameter sigma_c added to trilinear cohesive law");
		}
		else if (_d2_set && _d3_set && _G_set){
			ready = true;
			_delta_c1 = get_dc1();
			if (_delta_c1 > _delta_c2 || _delta_c1 <= 0.0)
				err_message("Incompatible parameter sigma_c added to trilinear cohesive law");
		}
		else if (_d1_set && _d2_set && _G_set){
			ready = true;
			_delta_c3 = get_dc3();
			if (_delta_c3 < _delta_c2)
				err_message("Incompatible parameter sigma_c added to trilinear cohesive law");
		}
	}
}
double TrilinearSmoothedCohesiveNUMaterial::get_sc()
{
	if (_s_set)
		return _sigma_c;
	else if (ready)     // The other 3 parameters are specified. Compute this one from them{
		return _G/(_delta_c2/2.0 + _delta_c3/2.0 - _delta_c1/3.0);
	else
		err_message("Peak traction could not be determined.");
}



void TrilinearSmoothedCohesiveNUMaterial::set_G(double G)
{
	if (_d1_set && _d2_set && _d3_set && _s_set)
		err_message("Can only specify 4 parameters for a trilinear cohesive model");
	else{
		_G = G;
		_G_set = true;

		if (_G <= 0.0)
			err_message("Incompatible parameter G added to trilinear cohesive law");

		if (_d1_set && _d2_set && _d3_set){
			ready = true;
			_sigma_c = get_sc();
		}
		else if (_d1_set && _d3_set && _s_set){
			ready = true;
			_delta_c2 = get_dc2();
			if (_delta_c2 < _delta_c1 || _delta_c2 > _delta_c3)
				err_message("Incompatible parameter G added to trilinear cohesive law");
		}
		else if (_d2_set && _d3_set && _s_set){
			ready = true;
			_delta_c1 = get_dc1();
			if (_delta_c1 > _delta_c2 || _delta_c1 <= 0.0)
				err_message("Incompatible parameter G added to trilinear cohesive law");
		}
		else if (_d1_set && _d2_set && _s_set){
			ready = true;
			_delta_c3 = get_dc3();
			if (_delta_c3 < _delta_c2)
				err_message("Incompatible parameter G added to trilinear cohesive law");
		}
	}	
}
double TrilinearSmoothedCohesiveNUMaterial::get_G()
{
	if (_G_set)
		return _G;
	else if (ready)     // The other 3 parameters are specified. Compute this one from them{
		return _sigma_c * (_delta_c2/2.0 + _delta_c3/2.0 - _delta_c1/3.0);
	else
		err_message("Fracture energy could not be determined.");
}




void TrilinearSmoothedCohesiveNUMaterial::set_b(double b)
{
	if (b > 0)
		_beta2 = b*b;
	else
		err_message("Beta must be a positive value");
}

void TrilinearSmoothedCohesiveNUMaterial::set_alpha(double a)
{
	if (a > 0)
		_alpha = a;
	else
		err_message("Alpha must be a positive value");
}
