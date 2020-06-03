#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>
#include <fstream>
#include <string>
#include <cassert>

using namespace std;

#include "types.h"
#include "utility.h"
#include "mkl.h"



#ifndef NVECTOR_H
#define NVECTOR_H



// inline void SparsMatrixProduct(const int n, const int m, const int l,
//   const char *matdescra, const complex<double> * value,
//   const int  *rowindex , const int *colindex, const int nonzero,
//   const complex<double> * B, complex<double> * C)
// {

// 	mkl_zcoomm ( "N" , &n , &l , &m , &mkl_alpha , matdescra , value , rowindex ,colindex ,&nonzero , B , &l , &mkl_beta , C ,
// 	&l );
// }





/********************************************************/
/*    template nVector vide:: seul les spécialisation   */    
/*                 sont inmplimentées                   */
/********************************************************/


template <typename T, int d> class nVector
{

};

/********************************************************/
/*           spécialisation pour le cas 1d              */    
/*                   complex<double>                    */
/********************************************************/
template <> class nVector<complex<double>,1>
{
public:
	//nVector(const int );
	nVector(const int, ComputeMatrixExponential = NoExpM );
	nVector(nVector<complex<double>,1> const & );
	nVector<complex<double>,1> & operator= (nVector<complex<double>,1> const & );
	complex<double> & operator() (const int );
	complex<double> const & operator() (const int ) const;
	complex<double> *  ptrTo ();
	complex<double> *  ptrTo () const;
	void swap(nVector<complex<double>,1> & );
	int  Size() const;
	void CopyTo(nVector <complex<double>, 1> const &);
	void AddTo(nVector <complex<double>, 1> const & );
	void AddTo(nVector <complex<double>, 1> const &, nVector <complex<double>, 1> const & );
	void Mul(const complex<double> );
	void Mul(const double );
	void Mul(nVector <complex<double>, 1> const & );
	double Norm() const;
	void Exp();
	void Exp(const double);
	void Exp(const complex<double> );
	void Exp(nVector <complex<double>, 1> const &);
	void Exp(const complex<double>, nVector <complex<double>, 1> const &);
	void MatrixForm_HMatrixExp(const complex<double>, nVector <complex<double>, 1>  &);
	void MatrixForm_CopyTo(nVector <complex<double>, 1> const &, MatrixPart );
	void MatrixForm_AddTo(nVector <complex<double>, 1> const &, MatrixPart );
	void MatrixForm_Initialize(complex<double>, complex<double>, MatrixPart);
	void Write(ofstream  &ofstr);
	void Write(const string filename);
	void Read(ifstream & ifstr);
	void Read(const string filename);
	~nVector();
private:
	complex<double> * data_;
	int size_1_;
	int size_;

	ComputeMatrixExponential CmpExpM_;
	int sqrt_size_;
	complex<double> * eve_;
	double * eva_;
	complex<double> * tmpv_;
};

/********************************************************/
/*           spécialisation pour le cas 1d              */    
/*                   double                             */
/********************************************************/
template <> class nVector<double,1>
{
public:
	nVector(const int );
	nVector(nVector<double,1> const & );
	nVector<double,1> & operator= (nVector<double,1> const & );
	double & operator() (const int );
	double const & operator() (const int ) const;
	double *  ptrTo ();
	double*  ptrTo () const;
	void swap(nVector<double,1> & );
	int  Size() const;
	void CopyTo(nVector <double, 1> const &);
	// void AddTo(nVector <cdouble, 1> const & );
	// void AddTo(nVector <double, 1> const &, nVector <complex<double>, 1> const & );
	// void Mul(const double );
	// void Mul(nVector double, 1> const & );
	// double Norm();
	// void Exp();
	// void Exp(const double );
	// void Exp(nVector <cdouble, 1> const &);
	// void Exp(const double, nVector <double, 1> const &);
	~nVector();
private:
	double * data_;
	int size_1_;
	int size_;
};


/********************************************************/
/*           spécialisation pour le cas 2d              */    
/*                   complex<double>                    */
/********************************************************/

template <> class nVector<complex<double>,2>
{
public:
	nVector(const int, const int );
	nVector(nVector <complex<double>, 2> const & );
	nVector(const int, nVector <complex<double>, 1> const & );
	nVector <complex<double>, 2> & operator= (nVector <complex<double>, 2> const & );
	complex<double> & operator() (const int, const int );
	complex<double> const & operator() (const int, const int  ) const;
	void swap(nVector <complex<double>, 2> & );
	int  Size( ) const;
	int  Size(int) const;
	void CopyTo(nVector <complex<double>, 2> const &);
	void Mul(const double );
	void Mul(const complex<double> );
	void Mul(nVector <complex<double>, 2> const & );
	void MulByConj(nVector <complex<double>, 2> const & );
	void Exp();
	void Exp(const double);
	void Exp(const complex<double> );
	double Norm() const;
	double Norm(int ) const;
	void Write(ofstream & ofstr);
	void Write(const string filename);
	void Read(ifstream & ifstr);
	void Read(const string filename);
	~nVector();
private:
	complex<double> * data_;
	int size_1_;
	int size_2_;
	int size_;
};





#endif
