/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	All rights reserved
	Last Updated September 2017

##################################################################################
*/
#ifndef _Writer_HOMOGENIZATION_H_
#define _Writer_HOMOGENIZATION_H_

#include "Writer.h"
#include "Utilities.h"




/*
 * This class is used to write the Problem sensitivity data to
 *  a single file over the course of the simulation
 */
class WriterHomogenization : public Writer
{
	public:


		/*
		 * The Big 3
		 * Constructor need the Problem
		 * Destructor doesn't do anything
		 * Copy Constructor copies the problem pointer and stores the mesh
		 */
		WriterHomogenization(Problem* prob);
		WriterHomogenization();
		virtual ~WriterHomogenization();
		WriterHomogenization(const WriterHomogenization& other);


		/*
		 * Main function used to write problem to any inherited file type
		 */
		virtual void writeConsecutive(std::string filename, double curr_t, id_type step);


	protected:


		/*
         * The actual function that will do all of the writing (maybe by calling other functions)
         */
        virtual void writeFromStream(std::ofstream& myfile, double curr_t);

};





#endif