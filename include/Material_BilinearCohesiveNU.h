/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated October 2017

##################################################################################
*/
#ifndef _NU_BL_AUG_COH_MATERIAL_H_
#define _NU_BL_AUG_COH_MATERIAL_H_

#include "material.h"

/*
 * Cohesive law based on bilinear model
 * Tractions described by an "effective" opening which
 *	is a combination of normal and tangential components
 */
class BilinearCohesiveNUMaterial : public Material
{
	protected:
		// Linear elastic parameters
		double _sigma_c;		// Maximum traction
		double _delta_c1;		// opening at peak traction
		double _delta_c2;		// Opening at fully failed
		double _G;				// Fracture toughness (in units of sigma_c * delta_c)
		double _beta2;			// Parameter governing the mode 1-mode 2 coupling
		double _alpha;			// Parameter governing the contact condition
		bool _s_set;
		bool _d1_set;
		bool _d2_set;
		bool _G_set;
		bool _b_set;
		bool ready; // Stores whether or not this material is ready for calls to constitutive
		
		void set_sc(double sc);
		double get_sc();
		void set_dc1(double dc);
		double get_dc1();
		void set_dc2(double dc);
		double get_dc2();
		void set_G(double dc);
		double get_G();
		void set_b(double b);
		void set_alpha(double a);

		double getDeltaEff(std::vector<double>& delta_ttn);

	public:
		BilinearCohesiveNUMaterial();
		virtual ~BilinearCohesiveNUMaterial() {};

		// Main constitutive subroutine
		virtual Material::output_params* Constitutive(Material::input_params& input);

		// Helper functions
		virtual material_type get_type() {return BL_COHESIVE_NO_UNLOADING;};
		virtual classification get_classification() {return STRUCTURAL;};
		virtual bool linear() {return false;};
		virtual void set_parameter(std::string name, double val);
		virtual double get_parameter(std::string name);

		// Functions specific to cohesive materials
		virtual bool is_cohesive() {return true;};
		virtual id_type get_cohesive_damage_zone(std::vector<double>& delta);
		virtual double getNormalizedCohesiveFailure(std::vector<double>& delta_ntt);
		
		// Annoying functions needed for copying a material
		virtual Material* allocate();			// virtual function
		virtual void copy(Material* other_mat);	// virtual function

		// Functions related to material sensitivity
		virtual SensitivityMaterialParameter* getSensitivityParameter(std::string name);
		virtual Material::sensitivity_output_params* SensitivityConstitutive(Material::sensitivity_input_params& input); 
};

#endif
