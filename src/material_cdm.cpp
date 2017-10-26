/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#include "material_cdm.h"
#include "Utilities.h"
#include "SensitivityMaterialParameterStructural.h"
#include "Utilities.h"
#include <cmath>
#include <algorithm>
#include <iostream>



Material* ContinuumDamageModelMaterial::allocate()
{
	return new ContinuumDamageModelMaterial;
}
// This function assumes the current material is set up correctly. Which it should't be possible for it not to be
void ContinuumDamageModelMaterial::copy(Material* other_mat)
{
	other_mat->set_name(_name);
	other_mat->set_id(_id);
	if(E_set)
		other_mat->set_parameter("E", _E);
	if(nu_set)
		other_mat->set_parameter("nu", _nu);
	if(lambda_set)
		other_mat->set_parameter("lambda", _lambda);
	if(mu_set)
		other_mat->set_parameter("mu", _mu);
	if(K_set)
		other_mat->set_parameter("K", _K);
	if(P1_set)
		other_mat->set_parameter("P1", _P1);
	if(P1_set)
		other_mat->set_parameter("P2", _P2);
	if(P1_set)
		other_mat->set_parameter("Yin", _Yin);
	if(P1_set)
		other_mat->set_parameter("mu_visc", _mu_visc);
	if (_thermal_exp_set)
		other_mat->set_parameter("THERMAL_EXPANSION", _thermal_exp);
}


ContinuumDamageModelMaterial::ContinuumDamageModelMaterial()
	: n_set(0), E_set(false), nu_set(false), lambda_set(false), mu_set(false), K_set(false),
	  P1_set(false), P2_set(false), Yin_set(false), mu_visc_set(false), ready(false),
	  _thermal_exp(0.0), _thermal_exp_set(false)
{
}



std::vector<std::string> ContinuumDamageModelMaterial::internal_vars_name()
{
	std::vector<std::string> vec = {"damage", "damage_threshold"};
	return vec;
}

std::vector<bool> ContinuumDamageModelMaterial::internal_vars_print()
{
	std::vector<bool> vec = {true, false};
	return vec;
}



void ContinuumDamageModelMaterial::computeLinearDmat(Material::input_params& input)
{
	if(n_set<2)
		err_message("Please set two independant parameters for a linear elastic isotropic material.");

	id_type dim = input.dim;
	bool plane_strain = input.plane_strain;

	// Compute the constitutive matrix
	if(dim==1)
	{
		_linear_Dmat.resize(1,1);
		_linear_Dmat(0,0) = _E;
	}
	else if(dim==2)
	{
		_linear_Dmat.clear();
		_linear_Dmat.resize(3, 3);
		if(plane_strain)
		{
			double coef = _E / ((1.0+_nu) * (1.0-2.0*_nu));
			_linear_Dmat(0,0) = coef * (1.0-_nu);
			_linear_Dmat(0,1) = coef * (_nu);
			_linear_Dmat(1,0) = coef * (_nu);
			_linear_Dmat(1,1) = coef * (1.0-_nu);
			_linear_Dmat(2,2) = coef * (1.0-2.0*_nu)/2.0;
		}
		else
		{
			double coef = _E / (1.0-_nu*_nu);
			_linear_Dmat(0,0) = coef;
			_linear_Dmat(0,1) = coef * (_nu);
			_linear_Dmat(1,0) = coef * (_nu);
			_linear_Dmat(1,1) = coef;
			_linear_Dmat(2,2) = coef * (1.0-_nu)/2.0;	// Division by 2 is to account for the engineering strain we use
		}
	}
	else if(dim==3)
	{
		_linear_Dmat.clear();
		_linear_Dmat.resize(6, 6);
		double coef = _E / ((1.0+_nu) * (1.0-2.0*_nu));
		_linear_Dmat(0,0) = coef * (1.0-_nu);
		_linear_Dmat(0,1) = coef * (_nu);
		_linear_Dmat(0,2) = coef * (_nu);
		_linear_Dmat(1,0) = coef * (_nu);
		_linear_Dmat(1,1) = coef * (1.0-_nu);
		_linear_Dmat(1,2) = coef * (_nu);
		_linear_Dmat(2,0) = coef * (_nu);
		_linear_Dmat(2,1) = coef * (_nu);
		_linear_Dmat(2,2) = coef * (1.0-_nu);
		_linear_Dmat(3,3) = coef * (1.0-2.0*_nu)/2;
		_linear_Dmat(4,4) = coef * (1.0-2.0*_nu)/2;
		_linear_Dmat(5,5) = coef * (1.0-2.0*_nu)/2;
	}
	else
		err_message("Dimensions for the D matrix must be less than or equal to 3.");
}


void ContinuumDamageModelMaterial::apply_damage(DenseMatrix<double>& D, double damage)
{
	unsigned char nrows = D.n_rows();

	switch (nrows)
	{
		case 3: // Probably the most common case (2D)
			for (unsigned char i=0; i<2; ++i)
				for (unsigned char j=0; j<2; ++j)
					D(i,j) = (1.0 - damage) * _linear_Dmat(i, j);
			D(2,2) = (1.0 - damage) * _linear_Dmat(2,2);
			break;


		case 6:	// (3D)
			for (unsigned char i=0; i<3; ++i)
				for (unsigned char j=0; j<2; ++j)
					D(i,j) = (1.0 - damage) * _linear_Dmat(i, j);
			D(3,3) = (1.0 - damage) * _linear_Dmat(3,3);
			D(4,4) = (1.0 - damage) * _linear_Dmat(4,4);
			D(5,5) = (1.0 - damage) * _linear_Dmat(5,5);
			break;
	

		case 1: // (1D)
			D(0,0) = (1.0 - damage) * _linear_Dmat(0,0);
			break;
		

		default:
			err_message("Invalid number of rows in the apply damage function.");
	}
}
void ContinuumDamageModelMaterial::compute_stress(std::vector<double>& stress, const DenseMatrix<double>& D, const std::vector<double>& strain)
{
	unsigned char nrows = D.n_rows();

	switch (nrows)
	{
		case 3: // Probably the most common case (2D)
			stress[0] = D(0,0)*strain[0] + D(0,1)*strain[1];
			stress[1] = D(1,0)*strain[0] + D(1,1)*strain[1];
			stress[2] = D(2,2)*strain[2];
			break;


		case 6:	// (3D)
			stress[0] = D(0,0)*strain[0] + D(0,1)*strain[1] + D(0,2)*strain[2];
			stress[1] = D(1,0)*strain[0] + D(1,1)*strain[1] + D(1,2)*strain[2];
			stress[2] = D(2,0)*strain[0] + D(2,1)*strain[1] + D(2,2)*strain[2];
			stress[3] = D(3,3)*strain[3];
			stress[4] = D(4,4)*strain[4];
			stress[5] = D(5,5)*strain[5];
			break;
	

		case 1: // (1D)
			stress[0] = D(0,0) * strain[0];
			break;
		

		default:
			err_message("Invalid number of rows in the apply damage function.");
	}
}




Material::output_params* ContinuumDamageModelMaterial::Constitutive(Material::input_params& input)
{
	if(!ready)
		err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// Get the internal variables and references to the output
	double& curr_damage = (*input.internal_vars)[0];
	double& curr_damage_thresh = (*input.internal_vars)[1];
	double dt = input.delta_t;
	double dmu = dt*_mu_visc;
	DenseMatrix<double>& Dmat_bar = _output.Dmat;
	std::vector<double>& stress = _output.stress;

	// Get the dimension of the problem and compute the strain
	id_type dim = input.dim;
	std::vector<double>& strain = input.strain;
	if (_thermal_exp_set)
	{
		std::vector<double> thermal_strain(strain.size(), 0.0);
		double strain_val = -1.0 * input.temp_change * _thermal_exp; // Get the volumetric strain component from the material thermal expansion ratio
		switch (input.dim)
		{
			case 1:
				thermal_strain[0] = strain_val;
				break;
			case 2:
				std::fill(thermal_strain.begin(), thermal_strain.begin()+2, strain_val);
				break;
			case 3:
				std::fill(thermal_strain.begin(), thermal_strain.begin()+3, strain_val);
				break;
			default:
				err_message("Wrong dimension somehow...");
		}
		Utilities::_VecAXPY(1.0, thermal_strain, strain); // strain = 1.0*thermal_strain + strain
	}

	// Compute the linear D matrix for this material
	if(dim != _curr_dim)
	{
		computeLinearDmat(input);
		_curr_dim = dim;
		Dmat_bar = _linear_Dmat;
		stress.resize(3);
	}

	// Compute the Strain energy (0.5*\epsilon'*D*\epsilon) (Known as Y_bar)
	std::vector<double> dYbar_depsilon = strain*_linear_Dmat;
	double strain_energy = Utilities::dot(dYbar_depsilon, strain);
	strain_energy *= 0.5;

	// Check the criteria to see if there is new damage present
	// First Check: If the strain energy is greater than the initial threshold (Yin)
	// Second Check: If the damage function exceeds the curent damage threshold (internal variable)
	bool new_damage = false;
	double tmp = 0.0;
	double damage_func = 0.0;
	if(strain_energy >= _Yin)
	{
		tmp = (strain_energy - _Yin)/(_P1*_Yin);    // This is purely to make the code more readable. This value is used several times. Could probably some up with a better name for it I guess..
		damage_func = 1.0 - exp(-1.0*pow(tmp, _P2));
		if(damage_func > curr_damage_thresh)
			new_damage = true;
		else
			new_damage = false;
	}
	else
		new_damage = false;



	// If there is new damage present in the material, do all of the updating
	if(new_damage)
	{
		// Update the damage parameters
		curr_damage += (dmu/(1.0+dmu)) * (damage_func-curr_damage_thresh);
		curr_damage_thresh = (curr_damage_thresh+dmu*damage_func) / (1.0+dmu);

		// Compute the linear stress
		compute_stress(stress, _linear_Dmat, strain);
		//stress = _linear_Dmat*strain;

		// Compute the derivative of damage wrt the strain
		double dG_dYbar = (_P2/(_P1*_Yin)) * pow(tmp, _P2-1.0) * exp(-1.0*pow(tmp, _P2));
		double domega_dYbar = (dmu/(1.0+dmu)) * dG_dYbar;
		
		std::vector<double> domega_depsilon(dYbar_depsilon.size());
		for(id_type i=0; i<domega_depsilon.size(); ++i)      // I really wish there was a syntax to mutliply a std::vector by a scalar...
			domega_depsilon[i] = domega_dYbar * dYbar_depsilon[i];

		// Compute the local tangent material matrix (D_bar = -domega_depsilon*\sigma + (1-damage)*D)
		apply_damage(Dmat_bar, curr_damage);
		//Dmat_bar = (1.0-curr_damage)*_linear_Dmat;
		for(id_type i=0; i<domega_depsilon.size(); ++i)
			for(id_type j=0; j<stress.size(); ++j)
				Dmat_bar(i,j) -= domega_depsilon[i] * stress[j];
		
		// Update the stress to include the damage (this is actual sigma_bar)
		for(id_type i=0; i<stress.size(); ++i)
			stress[i] = (1.0-curr_damage)*stress[i];
	}
	
	// Otherwise there is no new damage present, simply use the current damage
	else
	{
		// Compute he current D matrix
		apply_damage(Dmat_bar, curr_damage);

		// Compute the stress using the current D matrix
		compute_stress(stress, Dmat_bar, strain);
/*
		// Compute the stress
		stress = _linear_Dmat*strain;
		for(id_type i=0; i<stress.size(); ++i)
			stress[i] = (1.0-curr_damage)*stress[i];

		// Compute the local tangent material matrix (D_bar =(1-damage)*D)
		Dmat_bar = (1.0-curr_damage)*_linear_Dmat;
*/
	}

	// Return a pointer to the output structure
	return &_output;
}






Material::sensitivity_output_params* ContinuumDamageModelMaterial::SensitivityConstitutive(Material::sensitivity_input_params& input)
{
	err_message("Sensitivity analysis is not implemented for a continuum damage material yet!");

	// Return a pointer to the output structure
	return &_sensitivity_output;
}












SensitivityMaterialParameter* ContinuumDamageModelMaterial::getSensitivityParameter(std::string param_name)
{
	SensitivityMaterialParameter* ret = new SensitivityMaterialParameterStructural;
	ret->set_mat_name(_name);
	std::transform(param_name.begin(), param_name.end(), param_name.begin(), ::toupper); // Capitilize the name
	
	err_message("Sensitivity analysis is not implemented for a continuum damae material yet!");

	return ret;
}


















void ContinuumDamageModelMaterial::updateSensitivityISVs(Material::sensitivity_input_params& input)
{
	// if (!ready)
	// 	err_message("Not all material parameters have been set for the Continuum Damage Material.");

	// // Get the internal variables and references to the output
	// double& curr_damage = (*input.internal_vars)[0];
	// double& curr_damage_thresh = (*input.internal_vars)[1];
	// double& ddamage_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars()];
	// double& dthresh_dd = (*input.internal_vars_sensitivity)[param_id * n_internal_vars() + 1];
	// double dt = input.delta_t;
	// double dmu = dt*_mu_visc;

	// // Get the dimension of the problem
	// id_type dim = input.dim;
	// std::vector<double>& strain = input.strain;

	// // Compute the linear D matrix for this material
	// if(dim != _curr_dim)
	// {
	// 	computeLinearDmat(input);
	// 	_curr_dim = dim;
	// 	Dmat_bar = _linear_Dmat;
	// 	stress.resize(3);
	// }

	// // Compute the Strain energy (0.5*\epsilon'*D*\epsilon) (Known as Y_bar)
	// std::vector<double> dYbar_depsilon = strain*_linear_Dmat;
	// double strain_energy = Utilities::dot(dYbar_depsilon, strain);
	// strain_energy *= 0.5;

	// // Check the criteria to see if there is new damage present
	// // First Check: If the strain energy is greater than the initial threshold (Yin)
	// // Second Check: If the damage function exceeds the curent damage threshold (internal variable)
	// bool new_damage = false;
	// double tmp = 0.0;
	// double damage_func = 0.0;
	// if(strain_energy >= _Yin)
	// {
	// 	tmp = (strain_energy - _Yin)/(_P1*_Yin);    // This is purely to make the code more readable. This value is used several times. Could probably some up with a better name for it I guess..
	// 	damage_func = 1.0 - exp(-1.0*pow(tmp, _P2));
	// 	if(damage_func > curr_damage_thresh)
	// 		new_damage = true;
	// 	else
	// 		new_damage = false;
	// }
	// else
	// 	new_damage = false;
}







void ContinuumDamageModelMaterial::set_parameter(std::string name, double val)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="E" || name=="YOUNG'S MODULUS" || name=="YOUNGS MODULUS")
		set_E(val);
	else if(name=="NU" || name=="POISSON'S RATIO" || name=="POISSONS RATIO")
		set_nu(val);
	else if(name=="LAMBDA" || name=="FIRST LAME CONSTANT")
		set_lambda(val);
	else if(name=="MU" || name=="G" || name=="SECOND LAME CONSTANT" || name=="SHEAR MODULUS")
		set_mu(val);
	else if(name=="K" || name=="BULK MODULUS")
		set_K(val);
	else if(name=="P1")
		set_P1(val);
	else if(name=="P2")
		set_P2(val);
	else if(name=="YIN" || name=="DAMAGETHRESHOLD" || name=="DAMAGE THRESHOLD")
		set_Yin(val);
	else if(name=="MU_VISC" || name=="MU_VISCOUS" || name=="MU VISC" || name=="MU VISCOUS")
		set_mu_visc(val);
	else if (name=="THERMAL_EXPANSION" || name=="THERMAL EXPANSION" || name=="THERMAL_SHINKAGE" || name=="THERMAL SHRINKAGE")
		set_thermal_exp(val);
	else
		err_message(name << " is not a valid parameter name");
}
double ContinuumDamageModelMaterial::get_parameter(std::string name)
{
	std::transform(name.begin(), name.end(), name.begin(), ::toupper); // Capitilize the name
	if(name=="E" || name=="YOUNG'S MODULUS" || name=="YOUNGS MODULUS")
		return get_E();
	else if(name=="NU" || name=="POISSON'S RATIO" || name=="POISSONS RATIO")
		return get_nu();
	else if(name=="LAMBDA" || name=="FIRST LAME CONSTANT")
		return get_lambda();
	else if(name=="MU" || name=="G" || name=="SECOND LAME CONSTANT" || name=="SHEAR MODULUS")
		return get_mu();
	else if(name=="K" || name=="BULK MODULUS")
		return get_K();
	else if(name=="P1")
		return get_P1();
	else if(name=="P2")
		return get_P2();
	else if(name=="YIN" || name=="DAMAGETHRESHOLD" || name=="DAMAGE THRESHOLD")
		return get_Yin();
	else if(name=="MU_VISC" || name=="MU_VISCOUS" || name=="MU VISC" || name=="MU VISCOUS")
		return get_mu_visc();
	else if (name=="THERMAL_EXPANSION" || name=="THERMAL EXPANSION" || name=="THERMAL_SHINKAGE" || name=="THERMAL SHRINKAGE")
		return get_thermal_exp();
	else
		err_message(name << " is not a valid parameter name");
}


void ContinuumDamageModelMaterial::set_E(double E)
{
	if(n_set>=2 && !E_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_E = E;
	if(!E_set) // setting value for the first time
	{
		n_set = n_set + 1;
		E_set = true;
	}
	if(n_set==2)
	{
		if(nu_set) // Parameters E and nu were set
		{
			_lambda = (_E*_nu)/((1.0+_nu)*(1.0-2.0*_nu));
			_mu = _E/(2.0*(1.0+_nu));
			_K = _E/(3.0*(1-2.0*_nu));
		}
		else if(lambda_set) // Parameters E and lambda were set
		{
			double R = sqrt(_E*_E + 9.0*_lambda*_lambda + 2.0*_E*_lambda);
			_nu = 2.0*_lambda/(_E+_lambda+R);
			_mu = (_E-3.0*_lambda+R)/4.0;
			_K = (_E+3.0*_lambda+R)/6.0;
		}
		else if(mu_set)
		{
			_lambda = _mu*(_E-2.0*_mu)/(3.0*_mu-_E);
			_nu = _E/(2.0*_mu) - 1.0;
			_K = _E*_mu/(3.0*(3.0*_mu-_E));
		}
		else if(K_set)
		{
			_lambda = 3.0*_K*(3.0*_K-_E)/(9.0*_K-_E);
			_nu = (3.0*_K-_E)/(6.0*_K);
			_mu = 3.0*_K*_E/(9.0*_K-_E);
		}
	}
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_E()
{
	if(n_set>=2 || E_set)
		return _E;
	else
		err_message("Young's modulus could not be determined.");
}




void ContinuumDamageModelMaterial::set_nu(double nu)
{
	if(n_set>=2 && !nu_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_nu = nu;
	if(!nu_set) // setting value for the first time
	{
		n_set = n_set + 1;
		nu_set = true;
	}
	if(n_set==2)
	{
		if(E_set) // Parameters E and nu were set
		{
			_lambda = (_E*_nu)/((1.0+_nu)*(1.0-2.0*_nu));
			_mu = _E/(2.0*(1.0+_nu));
			_K = _E/(3.0*(1-2.0*_nu));
		}
		else if(lambda_set) // Parameters nu and lambda were set
		{
			_E = _lambda*(1.0+_nu)*(1.0-2.0*_nu)/_nu;
			_mu = _lambda*(1.0-2.0*_nu)/(2.0*_nu);
			_K = _lambda*(1.0+_nu)/(3.0*_nu);
		}
		else if(mu_set)
		{
			_lambda = 2.0*_mu*_nu/(1.0-2.0*_nu);
			_E = 2.0*_mu*(1.0+_nu);
			_K = 2.0*_mu*(1.0+_nu)/(3.0*(1.0-2.0*_nu));
		}
		else if(K_set)
		{
			_lambda = 3.0*_K*_nu/(1+_nu);
			_E = 3.0*_K*(1.0-2.0*_nu);
			_mu = 3.0*_K*(1.0-2.0*_nu)/(2.0*(1.0+_nu));
		}
	}
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_nu()
{
	if(n_set>=2 || nu_set)
		return _nu;
	else
		err_message("Poisson's ratio could not be determined.");
}




void ContinuumDamageModelMaterial::set_lambda(double lambda)
{
	if(n_set>=2 && !lambda_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_lambda = lambda;
	if(!lambda_set) // setting value for the first time
	{
		n_set = n_set + 1;
		lambda_set = true;
	}
	if(n_set==2)
	{
		if(E_set) // Parameters E and lambda were set
		{
			double R = sqrt(_E*_E + 9.0*_lambda*_lambda + 2.0*_E*_lambda);
			_nu = 2.0*_lambda/(_E+_lambda+R);
			_mu = (_E-3.0*_lambda+R)/4.0;
			_K = (_E+3.0*_lambda+R)/6.0;
		}
		else if(nu_set) // Parameters nu and lambda were set
		{
			_E = _lambda*(1.0+_nu)*(1.0-2.0*_nu)/_nu;
			_mu = _lambda*(1.0-2.0*_nu)/(2.0*_nu);
			_K = _lambda*(1.0+_nu)/(3.0*_nu);
		}
		else if(mu_set)
		{
			_nu = _lambda/(2.0*(_lambda+_mu));
			_E = _mu*(3.0*_lambda+2.0*_mu)/(_lambda+_mu);
			_K = _lambda+2.0*_mu/3.0;
		}
		else if(K_set)
		{
			_nu = _lambda/(3.0*_K-_lambda);
			_E = 9.0*_K*(_K-_lambda)/(3.0*_K-_lambda);
			_mu = 3.0*(_K-_lambda)/2.0;
		}
	}
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_lambda()
{
	if(n_set>=2 || lambda_set)
		return _lambda;
	else
		err_message("First Lame constant could not be determined.");
}




void ContinuumDamageModelMaterial::set_mu(double mu)
{
	if(n_set>=2 && !mu_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_mu = mu;
	if(!mu_set) // setting value for the first time
	{
		n_set = n_set + 1;
		mu_set = true;
	}
	if(n_set==2)
	{
		if(E_set) // Parameters E and mu were set
		{
			_lambda = _mu*(_E-2.0*_mu)/(3.0*_mu-_E);
			_nu = _E/(2.0*_mu) - 1.0;
			_K = _E*_mu/(3.0*(3.0*_mu-_E));
		}
		else if(nu_set) // Parameters nu and lambda were set
		{
			_lambda = 2.0*_mu*_nu/(1.0-2.0*_nu);
			_E = 2.0*_mu*(1.0+_nu);
			_K = 2.0*_mu*(1.0+_nu)/(3.0*(1.0-2.0*_nu));
		}
		else if(lambda_set)
		{
			_nu = _lambda/(2.0*(_lambda+_mu));
			_E = _mu*(3.0*_lambda+2.0*_mu)/(_lambda+_mu);
			_K = _lambda+2.0*_mu/3.0;
		}
		else if(K_set)
		{
			_nu = (3.0*_K-2.0*_mu)/(2.0*(3.0*_K+_mu));
			_E = 9.0*_K*_mu/(3.0*_K+_mu);
			_lambda = _K-2.0*_mu/3.0;
		}
	}
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_mu()
{
	if(n_set>=2 || mu_set)
		return _mu;
	else
		err_message("Shear modulus could not be determined.");
}




void ContinuumDamageModelMaterial::set_K(double K)
{
	if(n_set>=2 && !K_set)
		err_message("Can only set 2 independant parameters in a linear elastic isotropic material.");
	
	_K = K;
	if(!K_set) // setting value for the first time
	{
		n_set = n_set + 1;
		K_set = true;
	}
	if(n_set==2)
	{
		if(E_set) // Parameters E and K were set
		{
			_lambda = 3.0*_K*(3.0*_K-_E)/(9.0*_K-_E);
			_nu = (3.0*_K-_E)/(6.0*_K);
			_mu = 3.0*_K*_E/(9.0*_K-_E);
		}
		else if(nu_set) // Parameters nu and K were set
		{
			_lambda = 3.0*_K*_nu/(1+_nu);
			_E = 3.0*_K*(1.0-2.0*_nu);
			_mu = 3.0*_K*(1.0-2.0*_nu)/(2.0*(1.0+_nu));
		}
		else if(lambda_set)
		{
			_nu = _lambda/(3.0*_K-_lambda);
			_E = 9.0*_K*(_K-_lambda)/(3.0*_K-_lambda);
			_mu = 3.0*(_K-_lambda)/2.0;
		}
		else if(mu_set)
		{
			_nu = (3.0*_K-2.0*_mu)/(2.0*(3.0*_K+_mu));
			_E = 9.0*_K*_mu/(3.0*_K+_mu);
			_lambda = _K-2.0*_mu/3.0;
		}
	}
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_K()
{
	if(n_set>=2 || K_set)
		return _K;
	else
		err_message("Bulk modulus could not be determined.");
}




void ContinuumDamageModelMaterial::set_P1(double P1)
{
	_P1 = P1;
	P1_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_P1()
{
	if(P1_set)
		return _P1;
	else
		err_message("P1 could not be determined.");
}




void ContinuumDamageModelMaterial::set_P2(double P2)
{
	_P2 = P2;
	P2_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_P2()
{
	if(P2_set)
		return _P2;
	else
		err_message("P2 could not be determined.");
}




void ContinuumDamageModelMaterial::set_Yin(double Yin)
{
	_Yin = Yin;
	Yin_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_Yin()
{
	if(Yin_set)
		return _Yin;
	else
		err_message("Yin could not be determined.");
}




void ContinuumDamageModelMaterial::set_mu_visc(double mu_visc)
{
	_mu_visc = mu_visc;
	mu_visc_set = true;
	if((n_set>=2) && P1_set && P2_set && Yin_set && mu_visc_set)
		ready = true;
}
double ContinuumDamageModelMaterial::get_mu_visc()
{
	if(mu_visc_set)
		return _mu_visc;
	else
		err_message("mu_visc could not be determined.");
}


void ContinuumDamageModelMaterial::set_thermal_exp(double val)
{
	_thermal_exp = val;
	_thermal_exp_set = true;
}
double ContinuumDamageModelMaterial::get_thermal_exp()
{
	if (_thermal_exp_set)
		return _thermal_exp;
	else
		err_message("Thermal expansion ratio could not be determined.");
}