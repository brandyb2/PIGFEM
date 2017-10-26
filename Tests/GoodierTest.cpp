#include <iostream>
#include <fstream>
#include <cstddef>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include "Mesh.h"
#include "elem.h"
#include "Options.h"
#include "Problem.h"
#include "ProblemMaxPrincipalStress.h"
#include "ProblemNonlinearStructural.h"
#include "ProblemLinearStructural_ThermalPreLoad.h"
#include "ProblemNonlinearStructural_ThermalPreLoad.h"
#include "ProblemLinearElasticity.h"
#include "ProblemLinearThermal.h"
#include "BoundaryObject.h"
#include "DofObject.h"
#include "Inclusion.h"
#include "InclusionPlane.h"
#include "InclusionEllipse.h"
#include "InclusionEllipsoid.h"
#include "InclusionPolyhedron.h"
#include "InclusionCircle.h"
#include "InclusionSphere.h"
#include "material.h"
#include "material_cdm.h"
#include "material_lei.h"
#include "material_leti.h"
#include "material_leimps.h"
#include "material_leipcd.h"
#include "material_lt.h"
#include "Material_OPCohesiveAug.h"
#include "Material_OPCohesiveNUAug.h"
#include "Material_OPCohesive.h"
#include "Material_OPCohesiveNU.h"
#include "Material_XNCohesiveNU.h"
#include "Material_XNCohesive.h"
#include "Writer.h"
#include "Solver.h"
#include "STL_Reader.h"
#include "SensitivityRightLoad.h"
#include "common.h"
#include "Utilities.h"
#include "mpi.h"
#include <algorithm>
#include <cctype>
#include <iomanip>
using namespace std;

#define pi 3.141592653589799323
#define mu1 1.13e9
//#define mu2 1.13e9
#define mu2 31.15e9
#define nu1 0.33
//#define nu2 0.33
#define nu2 0.22
#define Force 50e6
#define a 2.11
#define d_size 4
//#define z_size 5e-2
#define z_size 4
#define x_disp 1e-1

#define nelem 4
#define n_refines 6
#define dimension 3

#define goodier true
// The size of the square domain in the positive and negative x and y directions (centers at (0,0))
// ~/projects/petsc-3.6.2/arch-linux2-c-opt/bin/mpiexec


// From Goodier (1933) Plane strain
std::vector<double> TwoD_Goodier(std::vector<double> gcoords);

// From Goodier (1933)
std::vector<double> ThreeD_Goodier(std::vector<double> gcoords);

// Simple homegeneous bar
std::vector<double> Analytical_x_tension(std::vector<double> gcoords);


template <typename T>
void print_vec(const std::vector<T> vec);




int main (int argc, char* argv[])
{
	Mesh mesh(&argc, &argv);

	std::vector<int> n_elem(n_refines);
	std::vector<double> Norm_err(n_refines);
	for (int i=0; i<n_refines; ++i)
	{
		mesh.init();
		n_elem[i] = nelem * pow(2, i);

		// Generate a new mesh
		PIGFEMPrint("Generating the mesh...");
		if (dimension==2)
			mesh.generate_mesh(QUAD4, -d_size, d_size, n_elem[i], -d_size, d_size, n_elem[i]);
		else if (dimension==3)
			mesh.generate_mesh(TET4, -d_size, d_size, n_elem[i], -d_size, d_size, n_elem[i], -d_size, d_size, n_elem[i]);

		// Define the problem
		ProblemLinearElasticity problem;
		problem.attach_mesh(&mesh);

		problem.getOptions()->parseInput(argc, argv);

		// Add a material to the mesh
		LinearElasticIsotropicMaterial material;
		material.set_parameter("mu", mu1);
		material.set_parameter("nu", nu1);
		material.set_name("Epoxy");
		mesh.add_material(&material);
		mesh.set_material("Epoxy");



		// Inclusion Materials
		if (goodier)
		{
			LinearElasticIsotropicMaterial material2;
			material2.set_parameter("mu", mu2);
			material2.set_parameter("nu", nu2);
			material2.set_name("Glass");
			if(mesh.dim()==2)
			{
				Circle_Inclusion circle;
				std::vector<double> center = {0, 0};
				circle.set_vec_parameter("center", center);
				circle.set_parameter("R", a);
				circle.set_material(&material2);
				mesh.add_inclusion(&circle);
			}
			else if(mesh.dim()==3)
			{
				Sphere_Inclusion sphere;
				std::vector<double> center = {0, 0, 0};
				sphere.set_vec_parameter("center", center);
				sphere.set_parameter("R", a);
				sphere.set_material(&material2);
				mesh.add_inclusion(&sphere);
			}
			else
				err_message("Unknown mesh dimension");



			// Goodier Boundary Conditions
			std::vector<std::string> sets;
			if(mesh.dim()==2)
				sets = {"left", "right", "top", "bottom"};
			else if(mesh.dim()==3)
				sets = {"left", "right", "front", "back", "top", "bottom"};
			else
				err_message("Unknown mesh dimension");
			for(unsigned int s=0; s<sets.size(); ++s)
			{
				std::set<id_type> curr_set = mesh.get_nodeset(sets[s]);
				for(auto it=curr_set.begin(), end=curr_set.end(); it!=end; ++it)
				{
					Node* node = mesh.get_node_global(*it);
					std::vector<double> gcoords = node->get_coords();
					std::vector<double> disp;
					if(mesh.dim()==2)
						disp = TwoD_Goodier(gcoords);
					else if(mesh.dim()==3)
						disp = ThreeD_Goodier(gcoords);
					else
						err_message("Unknown mesh dimension");
					for(unsigned int d=0; d<mesh.dim(); ++d)
						problem.get_boundary()->set_dirichlet_bc(*it, d, disp[d]);
				}
			}
		}
		else
		{
			BoundaryObject* boundary = problem.get_boundary();
			boundary->set_dirichlet_bcs_from_nodeset("right", 0, x_disp);
			boundary->set_dirichlet_bcs_from_nodeset("left", 0, 0.0);
			for(unsigned int d=0; d<mesh.dim(); ++d)
				boundary->set_dirichlet_bc(0, d, 0);

			if (mesh.dim() == 3)
			{
				id_type top_left_node = (n_elem[i]+1) * n_elem[i];
				boundary->set_dirichlet_bc(top_left_node, 2, 0.0);
			}
		}

		
		
		// Solve the problem!
		problem.init();
		PIGFEMPrint("\n\nPROBLEM DETAILS\n---------------------------------------------------\n");
		PIGFEMPrint("Problem Type: " << problem_type_names[problem.get_type()] << std::endl);
		PIGFEMPrint("Number of Mesh Partitions: " << mesh.n_ranks() << std::endl);
		PIGFEMPrint("Number of Inclusions: " << mesh.n_inclusions() << std::endl);
		PIGFEMPrint("Number of Elements: " << problem.get_mesh()->n_global_elem() << std::endl);
		PIGFEMPrint("Number of Nodes: " << problem.get_mesh()->n_global_nodes()+problem.get_mesh()->n_global_enrich_nodes()  << " (" << problem.get_mesh()->n_global_enrich_nodes() << " enriched)" << std::endl);
		problem.solve_problem();

		// Compute the L2 Norm of the error
		PIGFEMPrint("\n\nComputing the L2 Norm of the error...");
		if (goodier)
		{
			if(mesh.dim()==2)
				Norm_err[i] = problem.Compute_L2_Norm(&TwoD_Goodier);
			else
				Norm_err[i] = problem.Compute_L2_Norm(&ThreeD_Goodier);
		}
		else
			Norm_err[i] = problem.Compute_L2_Norm(&Analytical_x_tension);

		PIGFEMPrint("\n\tL2 Norm of the error: " << Norm_err[i] << "\n\n");

		mesh.clear();
	}

	PIGFEMPrint("Num Elem: ");
	print_vec(n_elem);
	PIGFEMPrint("\nL2 Norm: ");
	print_vec(Norm_err);
	PIGFEMPrint("\n\n\n");
}





template <typename T>
void print_vec(const std::vector<T> vec)
{
	if(vec.size() >= 1)
	{
		PIGFEMPrint("{");
		for(unsigned int i=0; i<(vec.size()-1); ++i)
			PIGFEMPrint(vec[i] << ", ");
		PIGFEMPrint(vec[vec.size()-1] << "}");
	}
}




// From Goodier (1933) Plane strain
std::vector<double> TwoD_Goodier(std::vector<double> gcoords)
{
	double theta = atan2(gcoords[1], gcoords[0]);
	double r = sqrt(pow(gcoords[0], 2) + pow(gcoords[1], 2));
	double E1 = 2.0*mu1*(1.0+nu1);

	double A = a*a* (Force/(4.0*mu1)) * ((1.0-2.0*nu2)*mu1-(1.0-2.0*nu1)*mu2)/((1.0-2.0*nu2)*mu1+mu2);
	double B = a*a*a*a* (Force/(4.0*mu1)) * (mu1-mu2)/(mu1+(3.0-4.0*nu1)*mu2);
	double C = a*a* (Force/(2.0*mu1)) * (mu1-mu2)/(mu1+(3.0-4.0*nu1)*mu2);

	double u_r, u_t;
	if(r < a) // Inside the inclusion
	{
		double H = (1.0/(2.0*nu2*a*a)) * (Force/(4.0*mu2) - (mu1/mu2)*(3.0*B/(a*a*a*a)-C/(a*a)) - B/(a*a*a*a) - (C/(a*a))*(1.0-2.0*nu1) - (Force/(2.0*E1))*(1.0+nu1));
		double G = B/(a*a*a*a) + (C/(a*a))*(1.0-2.0*nu1) + (Force/(2.0*E1))*(1.0+nu1) + (2.0*nu2-3.0)*H*a*a;
		double F = (1.0/a)*(A/a + (-B/(a*a*a) + (2.0*C/a)*(1.0-nu1) - G*a - 2.0*nu2*H*a*a*a)*cos(2.0*theta) + (Force*a/(2.0*E1))*(1.0+nu1)*(1.0-2.0*nu1+cos(2.0*theta)));
		u_r = F*r + (G*r + 2.0*nu2*H*r*r*r)*cos(2.0*theta);
		u_t = -1.0*(G*r + (3.0-2.0*nu2)*H*r*r*r)*sin(2.0*theta);
	}
	else 	// Outside the inclusion
	{
		u_r = A/r + (-1.0*B/(r*r*r) + (2.0*C/r)*(1.0-nu1))*cos(2.0*theta) + (Force*r/(2.0*E1))*(1.0+nu1)*(1.0-2.0*nu1+cos(2.0*theta));
		u_t = -1.0*(B/(r*r*r) + (C/r)*(1.0-2.0*nu1))*sin(2.0*theta) - (Force*r/(2.0*E1))*(1.0+nu1)*sin(2.0*theta);
	}

	// Transform from cylindrical coordinates to cartesian
	std::vector<double> ret(2);
	ret[0] = u_r*cos(theta) - u_t*sin(theta);
	ret[1] = u_r*sin(theta) + u_t*cos(theta);
	return ret;
}

// From Goodier (1933)
std::vector<double> ThreeD_Goodier(std::vector<double> gcoords)
{
	double r = sqrt(pow(gcoords[0], 2) + pow(gcoords[1], 2) + pow(gcoords[2], 2));
	double r_th = sqrt(pow(gcoords[0], 2) + pow(gcoords[1], 2)); // Distance from z-axis parallel to x-y plane
	double phi = atan2(gcoords[1], gcoords[0]); // Angle of rotation about the z axis
	double theta = atan2(r_th, gcoords[2]);
	double E1 = 2.0*mu1*(1.0+nu1);

	double A = a*a*a*( (-1.0*Force/(8.0*mu1)) * ((mu1-mu2)/((7.0-5.0*nu1)*mu1+(8.0-10.0*nu1)*mu2)) * 
					   (((1.0-2.0*nu2)*(6.0-5.0*nu1)*2.0*mu1+(3.0+19.0*nu2-20.0*nu1*nu2)*mu2)/((1.0-2.0*nu2)*2.0*mu1+(1.0+nu2)*mu2)) +
					   (Force/(4.0*mu1)) * ((((1.0-nu1)*(1.0+nu2)/(1.0+nu1)-nu2)*mu2-(1.0-2.0*nu2)*mu1)/((1.0-2.0*nu2)*2.0*mu1+(1.0+nu2)*mu2)) );
	double B = a*a*a*a*a * (Force/(8.0*mu1)) * ((mu1-mu2)/((7.0-5.0*nu1)*mu1+(8.0-10.0*nu1)*mu2));
	double C = a*a*a * (Force/(8.0*mu1)) * ((5.0*(1.0-2.0*nu1)*(mu1-mu2))/((7.0-5.0*nu1)*mu1+(8.0-10.0*nu1)*mu2));

	double u_r, u_t;
	if(r < a) // Inside the inclusion
	{
		double ur_a = -1.0*A/(a*a) - 3.0*B/(a*a*a*a) + (((5.0-4.0*nu1)/(1.0-2.0*nu1))*(C/(a*a))-9.0*(B/(a*a*a*a)))*cos(2.0*theta) + (Force*a/(2.0*E1))*(1.0-nu1+(1.0+nu1)*cos(2.0*theta));
		double ut_a = -1.0*(2.0*C/(a*a) + 6.0*B/(a*a*a*a))*sin(2.0*theta) - (Force*a/(2.0*E1))*(1.0+nu1)*sin(2.0*theta);
		double sigrt_a = 2.0*mu1*(-1.0*(2.0*(1.0+nu1)/(1.0-2.0*nu1))*(C/(a*a*a))+24.0*B/(a*a*a*a*a))*sin(2.0*theta) - (Force/2.0)*sin(2.0*theta);
		double G = (-1.0/(6.0*nu2*a*a*a)) * (a*sigrt_a/(2.0*mu2)-ut_a)/sin(2.0*theta);
		double F = (1.0/3.0)*((-1.0/(2.0*mu2*sin(2.0*theta)))*sigrt_a - (7.0+2.0*nu2)*G*a*a);
		double H = (1.0/a)*(ur_a - F*a - 2.0*nu2*G*a*a*a - (3.0*F*a+6.0*nu2*G*a*a*a)*cos(2.0*theta));
		u_r = H*r + F*r + 2.0*nu2*G*r*r*r + (3.0*F*r+6.0*nu2*G*r*r*r)*cos(2.0*theta);
		u_t = -1.0*(3.0*F*r+(7.0-4.0*nu2)*G*r*r*r)*sin(2.0*theta);
	}
	else 	// Outside the inclusion
	{
		u_r = -1.0*A/(r*r) - 3.0*B/(r*r*r*r) + (((5.0-4.0*nu1)/(1.0-2.0*nu1))*(C/(r*r))-9.0*B/(r*r*r*r))*cos(2.0*theta) + (Force*r/(2.0*E1))*(1.0-nu1+(1.0+nu1)*cos(2.0*theta));
		u_t = -1.0*(2.0*C/(r*r) + 6.0*B/(r*r*r*r))*sin(2.0*theta) - (Force*r/(2.0*E1))*(1.0+nu1)*sin(2.0*theta);
	}

	// Transform from spherical coordinates to cartesian
	std::vector<double> ret(3);
	ret[0] = cos(phi)*(sin(theta)*u_r + cos(theta)*u_t);
	ret[1] = sin(phi)*(sin(theta)*u_r + cos(theta)*u_t);
	ret[2] = cos(theta)*u_r - sin(theta)*u_t;
	return ret;
}

std::vector<double> Analytical_x_tension(std::vector<double> gcoords)
{
	double eps_x = x_disp/(2.0*d_size);
	double eps_t = -nu1*eps_x;
	std::vector<double> u(gcoords.size());
	u[0] = eps_x*(gcoords[0]+d_size);

	if(gcoords.size()==2) // 2D
		u[1] = (-eps_t/(nu1-1.0))*(gcoords[1]+d_size); // Plane strain
	else if(gcoords.size()==3)
	{
		u[1] = (gcoords[1]+d_size)*eps_t;
		u[2] = (gcoords[2]+z_size)*eps_t;
	}
	else
		err_message("Unknown mesh dimension");

	return u;
}