/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated March 2017

##################################################################################
*/
#include "WriterHomogenization.h"
#include "Problem.h"
#include "Mesh.h"
#include "Assembler.h"
#include "SensitivitySolver.h"

/*
 * The Big 3
 * Constructor need the Problem
 * Destructor doesn't do anything
 * Copy Constructor copies the problem pointer
 */
WriterHomogenization::WriterHomogenization(Problem* prob)
	: Writer(prob)
{}
WriterHomogenization::WriterHomogenization()
{}
WriterHomogenization::~WriterHomogenization()
{}
WriterHomogenization::WriterHomogenization(const WriterHomogenization& other)
{
	_prob = other.get_prob();
	store_mesh();
}



/*
 * Main function used to write problem to any inherited file type
 */
void WriterHomogenization::writeConsecutive(std::string filename, double curr_t, id_type step)
{
	if (_prob->get_mesh()->get_rank() == 0)
	{
		// Open file for writing or appending base on the load step
		std::ofstream myfile;
		if (step==0)
			myfile.open(filename.c_str(), std::ofstream::out);
		else
			myfile.open(filename.c_str(), std::ofstream::app);
		if (!myfile.good())
		{
			std::string out = "Error opening the output file " + filename;
			err_message( out.data() );
		}
		myfile.precision(16);

		writeFromStream(myfile, curr_t);

		myfile.close();
	}
}



/*
 * The actual function that will do all of the writing (maybe by calling other functions)
 */
void WriterHomogenization::writeFromStream(std::ofstream& myfile, double curr_t)
{
	myfile << curr_t << ",";

	std::vector<double> currentMacroVal = _prob->get_assembler()->getCurrentMacroscopicValue();
	std::vector<double> homogenizedVal = _prob->getHomogenizedValue();
	if (homogenizedVal.size() != 0)
	{
		for (id_type i=0; i<currentMacroVal.size(); ++i)
			myfile << currentMacroVal[i] << ",";
		for (id_type i=0; i<(homogenizedVal.size()-1); ++i)
			myfile << homogenizedVal[i] << ",";
		myfile << homogenizedVal[homogenizedVal.size()-1];
	}
	else
	{
		for (id_type i=0; i<(currentMacroVal.size()-1); ++i)
			myfile << currentMacroVal[i] << ",";
		myfile << currentMacroVal[currentMacroVal.size()-1];
	}
	myfile << "\n";
}