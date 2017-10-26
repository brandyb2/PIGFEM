/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "Writer_Binary_VTK.h"
#include "Problem.h"
#include "Mesh.h"
#include "InternalVars.h"
#include "NodalData.h"
#include "Utilities.h"
#include "mpi.h"



/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
Writer_Binary_VTK::Writer_Binary_VTK(Problem* prob)
	: Writer_VTK(prob), _isBigEndian(false)
{
	_isBigEndian = Utilities::isBigEndian();
}
Writer_Binary_VTK::~Writer_Binary_VTK()
{}
Writer_Binary_VTK::Writer_Binary_VTK(const Writer_Binary_VTK& other)
{
	_prob = other.get_prob();
	store_mesh();
}




/*
 * Virtual function to open the stream in ascii or binary mode
 */
void Writer_Binary_VTK::open(std::ofstream& myfile, std::string& filename)
{
	myfile.open(filename.c_str(), std::ofstream::out | std::ofstream::binary);
	if (!myfile.good())
	{
		std::string out = "Error opening the output file " + filename;
		err_message( out.data() );
	}
}




/*
 * Writes the standard ASCII header for a legacy VTK file
 */
void Writer_Binary_VTK::write_header(std::ofstream& myfile, double curr_t)
{
	Mesh* mesh = _prob->get_mesh();
	std::stringstream out;
	out << "# vtk DataFile Version 3.0\n";
	if (mesh->IGFEM())
		out << "IGFEM " << mesh->dim() << "D simulation: t=" << curr_t << "\n";
	else
		out << "FEM " << mesh->dim() << "D simulation: t=" << curr_t << "\n";
	out << "BINARY\nDATASET UNSTRUCTURED_GRID";

	writeASCII(myfile, out.str());
}




/*
 * Writes all of the nodal points to the VTK format. For serial ouput simply pass in empty remote lists
 */
void Writer_Binary_VTK::write_nodes(std::ofstream& myfile,
								   std::vector<id_type>& local_to_VTK,
								   const std::vector<id_type>& id_list, const std::vector<float>& coords_list,
								   const std::vector<std::vector<id_type> >& remote_id_lists, const std::vector<std::vector<float> >& remote_coords_lists)
{
	Mesh* mesh = _prob->get_mesh();

	// Only print if this is rank 0
	if (mesh->get_rank()==0)
	{
		std::stringstream out_str;
		id_type n_ranks = mesh->n_ranks();

		// Set up the mapping vector from global ids to VTK numbering
		local_to_VTK.resize(mesh->n_global_nodes() + mesh->n_global_enrich_nodes());

		// Print out the nodes themselves
		id_type curr_node = 0;
		out_str << "POINTS " << (mesh->n_global_nodes()+mesh->n_global_enrich_nodes()) << " float";
		writeASCII(myfile, out_str.str());

		for(id_type n=0; n<id_list.size(); ++n)
		{
			id_type id = id_list[n];
			if(id >= ENRICH_START)
				id = id - ENRICH_START + mesh->n_global_nodes();
			local_to_VTK[id] = curr_node;
			curr_node++;
		}
		writeBinary(myfile, &coords_list[0], coords_list.size());
		for(id_type part=1; part<n_ranks; ++part)
		{
			for(id_type n=0; n<remote_id_lists[part-1].size(); ++n)
			{
				id_type id = remote_id_lists[part-1][n];
				if(id >= ENRICH_START)
					id = id - ENRICH_START + mesh->n_global_nodes();
				local_to_VTK[id] = curr_node;
				curr_node++;
			}
			writeBinary(myfile, &remote_coords_lists[part-1][0], remote_coords_lists[part-1].size());
		}
	}
}




/*
 * Writes all of the elements to the VTK format. For serial ouput simply pass in empty remote list
 */
void Writer_Binary_VTK::write_elements(std::ofstream& myfile,
									  const std::vector<id_type>& local_to_VTK,
									  const std::vector<id_type>& eptr, const std::vector<id_type>& eind, const std::vector<unsigned char>& elem_to_type,
									  const std::vector<std::vector<id_type> >& remote_eptr, const std::vector<std::vector<id_type> >& remote_eind,
									  const std::vector<std::vector<unsigned char> >& remote_etype)
{
	Mesh* mesh = _prob->get_mesh();

	// Only print if this is rank 0
	if (mesh->get_rank()==0)
	{
		id_type n_ranks = mesh->n_ranks();

		id_type n_elem = eptr.size() - 1;
		id_type n_points = eind.size();
		for(id_type part=1; part<mesh->n_ranks(); ++part)
		{
			n_elem += (remote_eptr[part-1].size() - 1);
			n_points += remote_eind[part-1].size();
		}

		// Print all the element-node table information
		std::stringstream out_str;
		out_str << "\nCELLS " << n_elem << " " << (n_elem+n_points);
		writeASCII(myfile, out_str.str());

		for(id_type e=0; e<(eptr.size()-1); ++e)
		{
			id_type start = eptr[e];
			id_type end = eptr[e+1];
			std::vector<unsigned int> out(end-start+1); out[0] = end-start;
			for (id_type i=start; i<end; ++i)
			{
				id_type id = eind[i];
				id = (id < ENRICH_START) ? id : (id - ENRICH_START + mesh->n_global_nodes());
				out[i-start+1] = local_to_VTK[id];
			}
			writeBinary(myfile, &out[0], out.size());
		}
		for(id_type part=1; part<n_ranks; ++part)
		{
			for(id_type e=0; e<(remote_eptr[part-1].size()-1); ++e)
			{
				id_type start = remote_eptr[part-1][e];
				id_type end = remote_eptr[part-1][e+1];
				std::vector<int> out(end-start+1); out[0] = end-start;
				for (id_type i=start; i<end; ++i)
				{
					id_type id = remote_eind[part-1][i];
					id = (id < ENRICH_START) ? id : (id - ENRICH_START + mesh->n_global_nodes());
					out[i-start+1] = (int)(local_to_VTK[id]);
				}
				writeBinary(myfile, &out[0], out.size());
			}
		}

		// Print out the cell types data
		std::stringstream out_str2;
		out_str2 << "\nCELL_TYPES " << n_elem;
		writeASCII(myfile, out_str2.str());

		std::vector<int> cell_out(elem_to_type.begin(), elem_to_type.end());
		writeBinary(myfile, &cell_out[0], cell_out.size());
		for(id_type part=1; part<n_ranks; ++part)
		{
			std::vector<int> cell_out2(remote_etype[part-1].begin(), remote_etype[part-1].end());
			writeBinary(myfile, &cell_out2[0], cell_out2.size());
		}
	}
}




/*
 * Writes all of the nodal data to the VTK format. For serial ouput simply pass in empty remote lists
 */
void Writer_Binary_VTK::write_nodal_data(std::ofstream& myfile,
										const std::vector<float>& sol_list, const std::vector<float>& rhs_list,
										const std::vector<std::vector<float> >& remote_sol_lists, const std::vector<std::vector<float> >& remote_rhs_lists)
{
	Mesh* mesh = _prob->get_mesh();

	if (mesh->get_rank() == 0)
	{
		id_type nndof = _prob->nndof();
		id_type n_ranks = mesh->n_ranks();

		// Print out solution header info information
		std::stringstream out_str;
		out_str << "\nPOINT_DATA " << (mesh->n_global_nodes() + mesh->n_global_enrich_nodes()) << "\n";
		if(_prob->get_classification()==STRUCTURAL)
			out_str << "VECTORS displacement float";
		else
		{
			out_str << "SCALARS solution float " << nndof << "\n"; // I don't actually know if this is right
			out_str << "LOOKUP_TABLE default";
		}
		writeASCII(myfile, out_str.str());

		// Print out the actual solution info
		writeBinary(myfile, &sol_list[0], sol_list.size());
		for(id_type part=1; part<n_ranks; ++part)
			writeBinary(myfile, &remote_sol_lists[part-1][0], remote_sol_lists[part-1].size());

		// Optionally do the external loading information
		if (_printRHS)
		{
			std::stringstream out_str2;
			// Print out external load header information
			if(_prob->get_classification()==STRUCTURAL)
				out_str2 << "\nVECTORS external_force float";
			else
			{
				out_str2 << "SCALARS external_load float " << nndof << "\n"; // I don't actually know if this is right
				out_str2 << "LOOKUP_TABLE default";
			}
			writeASCII(myfile, out_str2.str());

			// Print out the actual load info
			writeBinary(myfile, &rhs_list[0], rhs_list.size());
			for(id_type part=1; part<n_ranks; ++part)
				writeBinary(myfile, &remote_rhs_lists[part-1][0], remote_rhs_lists[part-1].size());
		}
	}
}




/*
 * Writes all of the element data to the VTK format. For serial ouput simply pass in empty remote lists
 */
void Writer_Binary_VTK::write_element_data(std::ofstream& myfile,
										  const std::vector<id_type>& elem_mats, const std::map<std::string, std::vector<float> >& internal_vars_to_avgs,
										  const std::vector<std::vector<float> >& avg_strain, const std::vector<std::vector<float> >& avg_stress,
										  const std::vector<std::vector<id_type> >& remote_mats, const std::vector<std::vector<std::vector<float> > >& remote_int_vars,
										  const std::vector<std::vector<std::vector<float> > >& remote_strain, const std::vector<std::vector<std::vector<float> > >& remote_stress)
{
	// Get the mesh
	Mesh* mesh = _prob->get_mesh();

	if (mesh->get_rank() == 0)
	{
		id_type n_ranks = mesh->n_ranks();
		id_type n_elem = elem_mats.size();
		for (id_type p=0; p<remote_mats.size(); ++p)
			n_elem += remote_mats[p].size();

		// Print out the element materials
		// Print out the cell types data
		std::stringstream out_str;
		out_str << "\nCELL_DATA " << n_elem << "\n";
		out_str << "SCALARS materials int\n";
		out_str << "LOOKUP_TABLE default";
		writeASCII(myfile, out_str.str());
		
		writeBinary(myfile, &elem_mats[0], elem_mats.size());
		for(id_type part=1; part<n_ranks; ++part)
			writeBinary(myfile, &remote_mats[part-1][0], remote_mats[part-1].size());

		// Print out partitioning information
		if (_printPartition)
		{
			std::stringstream out_str6;
			out_str6 << "\nSCALARS partitions int\n";
			out_str6 << "LOOKUP_TABLE default";
			writeASCII(myfile, out_str6.str());
			
			std::vector<int> out(mesh->n_local_active_elem(), 0);
			writeBinary(myfile, &out[0], out.size());
			for(id_type part=1; part<n_ranks; ++part)
			{
				std::vector<int> out2(remote_mats[part-1].size(), part);
				writeBinary(myfile, &out2[0], out2.size());
			}
		}

		// Print out the internal variable information
		id_type iv = 0;
		for (auto it=internal_vars_to_avgs.begin(), end=internal_vars_to_avgs.end(); it!=end; ++it)
		{
			std::stringstream out_str2;
			out_str2 << "\nSCALARS " << it->first << " float\n";
			out_str2 << "LOOKUP_TABLE default";
			writeASCII(myfile, out_str2.str());
			
			const std::vector<float>& elem_avgs = it->second; // Get a reference
			writeBinary(myfile, &elem_avgs[0], elem_avgs.size());
			for(id_type part=1; part<n_ranks; ++part)
				writeBinary(myfile, &remote_int_vars[part-1][iv][0], remote_int_vars[part-1][iv].size());

			iv++;
		}

		// Print out stresses and strains
		if (_prob->get_classification() == STRUCTURAL)
		{
			id_type voigt = avg_strain[0].size();
			std::vector<std::string> voigt_map;
			if (voigt == 1)
				voigt_map = {"11"};
			else if (voigt == 3)
				voigt_map = {"11", "22", "12"};
			else if (voigt == 6)
				voigt_map = {"11", "22", "33", "23", "13", "12"};
			else
				err_message("Invalid number of average strains/stresses");

			// Output strains
			if (_printStrain)
			{
				for (id_type v=0; v<voigt; ++v)
				{
					std::stringstream out_str3;
					out_str3 << "\nSCALARS e" << voigt_map[v] << " float\n";
					out_str3 << "LOOKUP_TABLE default";
					writeASCII(myfile, out_str3.str());
					
					std::vector<float> out(avg_strain.size());
					for (id_type e=0; e<avg_strain.size(); ++e)
						out[e] = avg_strain[e][v];
					writeBinary(myfile, &out[0], out.size());

					for (id_type part=1; part<n_ranks; ++part)
					{
						out.resize(remote_strain[part-1].size());
						for (id_type e=0; e<remote_strain[part-1].size(); ++e)
							out[e] = remote_strain[part-1][e][v];
						writeBinary(myfile, &out[0], out.size());
					}
				}
			}

			// Output stresses
			if (_printStress)
			{
				for (id_type v=0; v<voigt; ++v)
				{
					std::stringstream out_str4;
					out_str4 << "\nSCALARS s" << voigt_map[v] << " float\n";
					out_str4 << "LOOKUP_TABLE default";
					writeASCII(myfile, out_str4.str());

					std::vector<float> out(avg_stress.size());
					for (id_type e=0; e<avg_stress.size(); ++e)
						out[e] = avg_stress[e][v];
					writeBinary(myfile, &out[0], out.size());

					for (id_type part=1; part<n_ranks; ++part)
					{
						out.resize(remote_stress[part-1].size());
						for (id_type e=0; e<remote_stress[part-1].size(); ++e)
							out[e] = remote_stress[part-1][e][v];
						writeBinary(myfile, &out[0], out.size());
					}
				}

				// Output von Mises stresses
				std::stringstream out_str5;
				out_str5 << "\nSCALARS von_mises float\n";
				out_str5 << "LOOKUP_TABLE default";
				writeASCII(myfile, out_str5.str());

				std::vector<float> out(avg_stress.size());
				for (id_type e=0; e<avg_stress.size(); ++e)
					out[e] = Utilities::von_mises(avg_stress[e]);
				writeBinary(myfile, &out[0], out.size());

				for (id_type part=1; part<n_ranks; ++part)
				{
					out.resize(remote_stress[part-1].size());
					for (id_type e=0; e<remote_stress[part-1].size(); ++e)
						out[e] = Utilities::von_mises(remote_stress[part-1][e]);
					writeBinary(myfile, &out[0], out.size());
				}
			}
		}
	}
}