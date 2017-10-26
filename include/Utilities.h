/*
##################################################################################

	Written by David Brandyberry, University of Illinois at Urbana-Champaign
	Please contact (brandyb2@illinois.edu) before distributing
	Last Updated September 2016

##################################################################################
*/
#ifndef _UTILITIES_H_
#define _UTILITIES_H_



#include "DenseMatrix.h"
#include "common.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <gsl_math.h>
#include <gsl_eigen.h> // For computation of the principal stresses in 3D
#include "mpi.h"
#include <time.h>
#include <sys/time.h>

/*
 * This file is to have several namespaces that define functions that are needed but don't really fit into any one call
 * They may be used across multiple classes
 */

class Elem;
class NodalData;

namespace ProblemUtilities
{
	// Function to assemble the "B" matrix used in assembly of the small strain elasticity stiffness matrix.
	void Assemble_Small_Strain_B_Mat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x);

	// Function to assemble the "B" matrix used in assembly of the finite strain elasticity stiffness matrix.
	void Assemble_Finite_Strain_B_Mat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x);

	// Function to assemble the "B" matrix used in assembly of the thermal stiffness matrix.
	void Assemble_Thermal_B_Mat(DenseMatrix<double>& B, const std::vector<std::vector<double> >& grad_x);

	// Function to assemble the "B" matrix used in assembly of the elasticity stiffness matrix
	void Assemble_Cohesive_N_Mat(DenseMatrix<double>& N, const std::vector<double>& shapes, id_type dim);

	// This function essentially does Bt*D*B taking into account all of the zeros in the B matrix so its faster
	void SmallStrainKFast(DenseMatrix<double>& K, const std::vector<std::vector<double> >& dN, const DenseMatrix<double>& D);

	// This function essentially does B*U taking into account all of the zeros in the B matrix so its faster
	void SmallStrainFast(std::vector<double>& strain, const std::vector<std::vector<double> >& dN, const std::vector<double>& U);
	void SmallStrainFast(std::vector<double>& strain, const std::vector<std::vector<double> >& dN, const DenseMatrix<double>& U);

	// This function performs Bt*sigma taking into accoun the zero structure of Bt
	void SmallStrainInternalForceFast(std::vector<double>& P, const std::vector<std::vector<double> >& dN, const std::vector<double>& sigma);

	// This function does Bt*D*B  (B-matrix variant)
	void Bt_D_B(DenseMatrix<double>& K, const DenseMatrix<double>& B, const DenseMatrix<double>& D);

	// This function does B*U (B-matrix variant)
	void B_U(std::vector<double>& strain, const DenseMatrix<double>& B, const std::vector<double>& U);

	// This function performs Bt*sigma (B-matrix variant)
	void Bt_F(std::vector<double>& P, const DenseMatrix<double>& B, const std::vector<double>& sigma);

	// This function computes the opening accross a cohesive surface from the shape functions, solution vector, and rotatino matrix
	void cohesiveDeltaFast(std::vector<double>& delta, const std::vector<double>& N, const std::vector<double>& U, const DenseMatrix<double>& R);

	// This function compute the stiffness matrix contribution from a cohesive opening
	void cohesiveKFast(DenseMatrix<double>& K, const std::vector<double>& N, const DenseMatrix<double>& D, const DenseMatrix<double>& R);

	// This function compute the internal force contribution contribution from a cohesive opening
	void cohesiveInternalForceFast(std::vector<double>& P, const std::vector<double>& N, const std::vector<double>& trac, const DenseMatrix<double>& R);

	// This function forms the elemental solution vector
	void formElementalVector(std::vector<double>& vec, Elem* el, NodalData* data);
	void formElementalMatrix(DenseMatrix<double>& elem_mat, Elem* el, NodalData* data);
}




namespace Utilities
{
	class timer {
		public:
			struct timespec start_time;

			inline
			void start()
			{
				clock_gettime(CLOCK_MONOTONIC, &start_time);
			}

			inline
			double time_elapsed()
			{
				struct timespec end;
				clock_gettime(CLOCK_MONOTONIC, &end);
				return (end.tv_sec - start_time.tv_sec) + (end.tv_nsec - start_time.tv_nsec)*(1e-9);
			}
	};

	// constrains the given angle to be within -pi to pi
	double constrain_angle(double x);

	// Kroeneker Delta
	int kron(int i, int j);

	// Converts stress or strain to the principal values (2D requires plane strain boolean, defaults to false)
	std::vector<double> principal_stress_strain(const std::vector<double>& vec/*, bool plane_strain = false*/);

	double von_mises(const std::vector<double>& stress);
	float von_mises(const std::vector<float>& stress);

	// General wrapper to recieve a message of unknown size using MPI
	template <typename T>
	void RecieveUnknown(std::vector<T>& recv_vec, int source, int tag, MPI_Datatype type, MPI_Comm comm)
	{
		MPI_Status status;
		int count;
		MPI_Probe(source, tag, comm, &status);
		MPI_Get_count(&status, type, &count);
		recv_vec.resize(count);
		MPI_Recv(recv_vec.data(), count, type, source, tag, comm, &status);
	}

	// Some basic vector operations
	double dot(const std::vector<double>& v1, const std::vector<double>& v2);
	double perp(const std::vector<double>& v1, const std::vector<double>& v2);
	std::vector<double> cross(const std::vector<double>& v1, const std::vector<double>& v2);
	std::vector<double> plus(const std::vector<double>& v1, const std::vector<double>& v2);
	std::vector<double> minus(const std::vector<double>& v1, const std::vector<double>& v2);
	std::vector<double> VecAXPY(double alpha, const std::vector<double>& x, const std::vector<double>& y);
	std::vector<double> scale(const std::vector<double>& vec, double alpha);
	std::vector<float> scale(const std::vector<float>& vec, float alpha);
	void _VecAXPY(double alpha, std::vector<double>& x, std::vector<double>& y);
	// void _scale(std::vector<double>& vec, double alpha);
	double _dot(std::vector<double> v1, std::vector<double> v2);
	double _perp(std::vector<double> v1, std::vector<double> v2);
	std::vector<double> _cross(std::vector<double> v1, std::vector<double> v2);
	std::vector<double> _plus(std::vector<double> v1, std::vector<double> v2);
	std::vector<double> _minus(std::vector<double> v1, std::vector<double> v2);
	std::vector<double> _scale(std::vector<double> vec, double alpha);

	double det(DenseMatrix<double>& J);
	DenseMatrix<double> inverse(DenseMatrix<double>& J);

	// Finds the intersection of Line starting a P0 in the direction dS and line atsrting in the direction V0 in the direction e
	// Returns the number of multiples of dS that occur bitween Po and the intersection
	double LineLineIntersection(const std::vector<double>& P0, const std::vector<double>& dS,
								 const std::vector<double>& V0, const std::vector<double>& e);

	// Find the intersection of a line starting at P0 in the dirction dS and the plane with normal n and point V0
	// Returns the number of multiples of dS that occur bitween Po and the intersection
	double LinePlaneIntersection(const std::vector<double>& P0, const std::vector<double>& dS,
								  const std::vector<double>& V0, const std::vector<double>& n);

	// Returns whether or not the given point P is inside of the projected triangle defined by a, b, c
	bool pointInTriangleSS(const std::vector<double>& P, const std::vector<double>& a,
						  const std::vector<double>& b, const std::vector<double>& c);
	bool pointInTriangleBary(const std::vector<double>& P, const std::vector<double>& a,
							const std::vector<double>& b, const std::vector<double>& c);
	// Helper function fo pointInTriangle. Deermines is P1 and P2 are in the same half plane defined by AB
	bool sameSide(const std::vector<double>& P1, const std::vector<double>& P2,
				 const std::vector<double>& a, const std::vector<double>& b);

	// Splits a string into pieces using the delimiter specified
	std::vector<std::string> splitString(std::string obj, std::string del);




	inline
	float dtos(double val)
	{
		return (float) val;
	}

	inline
	std::vector<float> dtos(std::vector<double> val)
	{
		std::vector<float> ret(val.size());
		for (id_type i=0; i<val.size(); ++i)
			ret[i] = (float) val[i];
		return ret;
	}
	


	/*
	 * Functions to do with converting between endianness formats
	 */
	// Checks to see if he system this is being run on stores data in big Engian format
	// Important for binary VTK output
	int isBigEndian(void);

	template <class T>
	T swap32( T f )
	{
		union
		{
			T f;
			unsigned char b[4];
		} dat1, dat2;
		dat1.f = f;
		dat2.b[0] = dat1.b[3];
		dat2.b[1] = dat1.b[2];
		dat2.b[2] = dat1.b[1];
		dat2.b[3] = dat1.b[0];
		return dat2.f;
	}

	template <class T>
	T swap64( T f )
	{
		union
		{
			T f;
			unsigned char b[8];
		} dat1, dat2;
		dat1.f = f;
		dat2.b[0] = dat1.b[7];
		dat2.b[1] = dat1.b[6];
		dat2.b[2] = dat1.b[5];
		dat2.b[3] = dat1.b[4];
		dat2.b[4] = dat1.b[3];
		dat2.b[5] = dat1.b[2];
		dat2.b[6] = dat1.b[1];
		dat2.b[7] = dat1.b[0];
		return dat2.f;
	}

	template <class T>
	T byteSwap(T data)
	{
		err_message("byteSwap is not defined for this data type");
	}
	template <>
	inline
	int byteSwap(int data)
	{
		return swap32(data);
	}
	template <>
	inline
	unsigned int byteSwap(unsigned int data)
	{
		return swap32(data);
	}
	template <>
	inline
	float byteSwap(float data)
	{
		return swap32(data);
	}
	template <>
	inline
	double byteSwap(double data)
	{
		return swap64(data);
	}
	template <>
	inline
	char byteSwap(char data)
	{
		return data; // only 1 byte
	}
	template <>
	inline
	unsigned char byteSwap(unsigned char data)
	{
		return data; // only 1 byte
	}


}



#endif