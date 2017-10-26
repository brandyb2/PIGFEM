/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _CDMMATERIAL_H_
#define _CDMMATERIAL_H_

#include "material.h"

class ContinuumDamageModelMaterial : public Material
{
	protected:
		// Linear elastic parameters
		double _E;
		double _nu;
		double _lambda;
		double _mu;
		double _K;
		int n_set;
		bool E_set;
		bool nu_set;
		bool lambda_set;
		bool mu_set;
		bool K_set;
		
		// Damage parameters
		double _P1;
		double _P2;
		double _Yin;
		double _mu_visc;
		bool P1_set;
		bool P2_set;
		bool Yin_set;
		bool mu_visc_set;
		bool ready; // Stores whether or not this material is ready for calls to constitutive
		
		void set_E(double E);
		double get_E();
		void set_nu(double nu);
		double get_nu();
		void set_lambda(double lambda);
		double get_lambda();
		void set_mu(double mu);
		double get_mu();
		void set_K(double K);
		double get_K();
		
		void set_P1(double P1);
		double get_P1();
		void set_P2(double P2);
		double get_P2();
		void set_Yin(double YIn);
		double get_Yin();
		void set_mu_visc(double mu_visc);
		double get_mu_visc();

		// Thermal expansion parameters
		double _thermal_exp;
		bool _thermal_exp_set;
		void set_thermal_exp(double t);
		double get_thermal_exp();

		// The linear constitutive matrix
		DenseMatrix<double> _linear_Dmat;

		// Function to compute the linear D matrix for the given dimension
		virtual void computeLinearDmat(Material::input_params& input);

		void apply_damage(DenseMatrix<double>& D, double damage);
		void compute_stress(std::vector<double>& stress, const DenseMatrix<double>& D, const std::vector<double>& strain);

	public:
		ContinuumDamageModelMaterial();
		virtual ~ContinuumDamageModelMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return CONTINUUM_DAMAGE;};
		virtual classification get_classification() {return STRUCTURAL;};
		virtual bool linear() {return false;};
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);

		// Internal variable functions
		virtual id_type n_internal_vars() {return 2;};
		virtual std::vector<double> init_internal_vars() {return std::vector<double>(2,0.0);};
		virtual std::vector<double> max_internal_vars() {return std::vector<double>(2,1.0);};
		virtual std::vector<std::string> internal_vars_name();
		virtual std::vector<bool> internal_vars_print();
		
		// Annoying functions needed for communicating a material
		virtual Material* allocate();			// virtual function
		virtual void copy(Material* other_mat);	// virtual function

		// Functions related to material sensitivity
		virtual SensitivityMaterialParameter* getSensitivityParameter(std::string name);
		virtual Material::sensitivity_output_params* SensitivityConstitutive(Material::sensitivity_input_params& input);
		virtual void updateSensitivityISVs(Material::sensitivity_input_params& input);
};


#endif
