#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>
#include <fstream>
#include <string>
#include <iomanip>
#include <cassert>

using namespace std;

#include "mathematical_constants.h"
#include "physical_constants.h"
#include "types.h"
#include "utility.h"
#include "nvector.h"

#include "mkl.h"

const complex<double> mkl_alpha{1.f, 0.f};
const complex<double> mkl_beta {0.f, 0.f};

/*********** nVector<complex<double>,1>******************/

nVector<complex<double>,1>::nVector(const int size_1, ComputeMatrixExponential CmpExpM)
{
	assert(size_1 > 0);
	size_1_ = size_1;
	size_   = size_1;
	data_   = nullptr;
	data_   = (complex<double> *)mkl_calloc( size_t(size_), sizeof( complex<double> ), 64);

	CmpExpM_ = CmpExpM;
	sqrt_size_ = 0;
	eve_  = nullptr;
	eva_  = nullptr;
	tmpv_ = nullptr;

	if(CmpExpM_ == ExpM)
	{
		sqrt_size_ = static_cast<int>(sqrt(size_)); 
		assert(size_ == sqrt_size_ * sqrt_size_);
		eve_   = (complex<double> *)mkl_calloc( size_t(size_), sizeof( complex<double> ), 64);
		eva_   = (double *)mkl_calloc( size_t(sqrt_size_), sizeof( double ), 64);
		tmpv_  = (complex<double> *)mkl_calloc( size_t(sqrt_size_), sizeof( complex<double> ), 64);
	}
}
// nVector<complex<double>,1>::nVector(nVector <complex<double>, 1> const & rhs)
// {
// 	size_1_ = rhs.size_1_;
// 	size_   = rhs.size_;
// 	data_   = nullptr;
// 	data_   = (complex<double> *)mkl_calloc( size_t(size_), sizeof( complex<double>), 64);
// 	CopyTo(rhs);
// }

nVector<complex<double>,1>::nVector(nVector <complex<double>, 1> const & rhs)
{
	size_1_ = rhs.size_1_;
	size_   = rhs.size_;
	data_   = nullptr;
	data_   = (complex<double> *)mkl_calloc( size_t(size_), sizeof( complex<double>), 64);
	CopyTo(rhs);

	CmpExpM_   = rhs.CmpExpM_;
	sqrt_size_ = rhs.sqrt_size_;
	eve_  = nullptr;
	eva_  = nullptr;
	tmpv_ = nullptr;
	if(CmpExpM_ == ExpM)
	{
		eve_   = (complex<double> *)mkl_calloc( size_t(size_), sizeof( complex<double> ), 64);
		eva_   = (double *)mkl_calloc( size_t(sqrt_size_), sizeof( double ), 64);
		tmpv_  = (complex<double> *)mkl_calloc( size_t(sqrt_size_), sizeof( complex<double> ), 64);
	}
}

nVector <complex<double>, 1> & nVector<complex<double>,1>::operator= (nVector <complex<double>, 1> const & rhs)
{
	nVector <complex<double>, 1> other(rhs);
	swap(other);
	return (*this);
}

complex<double> & nVector<complex<double>,1>::operator() (const int i)
{
	assert(i >= 0); 
	assert(i < size_1_);
	return data_[i];
}

complex<double> const & nVector<complex<double>,1>::operator() (const int i) const
{
	assert(i >= 0); 
	assert(i < size_1_);
	return data_[i];
}

complex<double> *  nVector<complex<double>,1>::ptrTo()
{
	return data_;
}
complex<double> *  nVector<complex<double>,1>::ptrTo () const
{
	return data_;
}

void nVector<complex<double>,1>::swap(nVector <complex<double>, 1> & other)
{
	std::swap(size_1_,other.size_1_);
	std::swap(size_,other.size_);
	std::swap(data_,other.data_);

	std::swap(CmpExpM_,other.CmpExpM_);
	std::swap(sqrt_size_,other.sqrt_size_);
	std::swap(eve_,other.eve_);
	std::swap(eva_,other.eva_);
	std::swap(tmpv_,other.tmpv_);
}
nVector<complex<double>,1>::~nVector()
{
	mkl_free(data_);
	if(CmpExpM_ == ExpM)
	{
		mkl_free(eve_);
		mkl_free(eva_);
		mkl_free(tmpv_);
	}
}
int nVector<complex<double>,1>::Size() const
{
	return size_;
}

void nVector<complex<double>,1>::CopyTo(nVector <complex<double>, 1> const & v)
{
	assert(size_ == v.size_);
	cblas_zcopy(size_, v.data_, 1, data_, 1);
}

void nVector<complex<double>,1>::MatrixForm_CopyTo(nVector <complex<double>, 1> const & M, MatrixPart  part)
{
	int sqrt_size = static_cast<int>(sqrt(size_));
	assert(size_ == M.size_);
	assert(size_ == sqrt_size * sqrt_size);

	LAPACKE_zlacpy ( LAPACK_ROW_MAJOR, part , sqrt_size, sqrt_size, M.data_, sqrt_size, data_, sqrt_size);
	
}

void nVector<complex<double>,1>::MatrixForm_AddTo(nVector <complex<double>, 1> const & M, MatrixPart  part)
{
	int sqrt_size = static_cast<int>(sqrt(size_));
	assert(size_ == M.size_);
	assert(size_ == sqrt_size * sqrt_size);

	switch(part)
	{
		case UpperTriangularPart :
			for (int m = 0; m < sqrt_size; ++m)
			{
				for (int mp = m; mp < sqrt_size; ++mp)
				{
					data_[ToRowMajorIndex(m, mp, sqrt_size )] += M(ToRowMajorIndex(m, mp, sqrt_size ));
				}
			}
		    break;

		case LowerTriangularPart :
		    for (int m = 0; m < sqrt_size; ++m)
		    {
		    	for (int mp = m; mp < sqrt_size; ++mp)
		    	{
		    		data_[ToRowMajorIndex(mp, m, sqrt_size )] += M(ToRowMajorIndex(mp, m, sqrt_size ));
		    	}
		    }
		    break;  
		case All:
			AddTo(M);   
	}
}



void nVector<complex<double>,1>::MatrixForm_Initialize(complex<double> diag, complex<double> elsewhere, MatrixPart  part)
{
	int sqrt_size = static_cast<int>(sqrt(size_));
	assert(size_ == sqrt_size * sqrt_size);

	LAPACKE_zlaset (LAPACK_ROW_MAJOR , part, sqrt_size, sqrt_size, elsewhere , diag , data_ , sqrt_size);
}








void nVector<complex<double>,1>::AddTo(nVector <complex<double>, 1> const & v)
{
	assert(size_ == v.size_);
	vzAdd( size_, v.data_, data_, data_);
}

void nVector<complex<double>,1>::AddTo(nVector <complex<double>, 1> const & v1, nVector <complex<double>, 1> const & v2)
{
	assert(size_ == v1.size_);
	assert(size_ == v2.size_);
	vzAdd( size_, v1.data_, v2.data_, data_);
}


void nVector<complex<double>,1>::Mul(const complex<double> alpha)
{
	cblas_zscal(size_, &alpha, data_, 1);
}

void nVector<complex<double>,1>::Mul(const double alpha)
{
	cblas_zdscal(size_, alpha, data_, 1);
}

void nVector<complex<double>,1>::Mul(nVector <complex<double>, 1> const & v)
{
	assert(size_ == v.size_);
	vzMul( size_, v.data_, data_, data_ );
}

double nVector<complex<double>,1>::Norm() const
{
	return cblas_dznrm2( size_ ,data_, 1);
}

void nVector<complex<double>,1>::Exp()
{
	vzExp(size_, data_, data_); 
}
void nVector<complex<double>,1>::Exp(const double alpha)
{
	Mul(alpha);
	vzExp(size_, data_, data_); 
}

void nVector<complex<double>,1>::Exp(const complex<double> alpha)
{
	Mul(alpha);
	vzExp(size_, data_, data_); 
}

void nVector<complex<double>,1>::Exp( nVector <complex<double>, 1> const & v)
{
	assert(size_ == v.size_);
	vzExp(size_, v.data_, data_); 
}

void nVector<complex<double>,1>::Exp(const complex<double> alpha, nVector <complex<double>, 1> const & v)
{
	assert(size_ == v.size_);
	CopyTo(v);
	Mul(alpha);
	vzExp(size_, data_, data_); 
}

void nVector<complex<double>,1>::MatrixForm_HMatrixExp(const complex<double> alpha, nVector <complex<double>, 1>  & v)
{
	assert(CmpExpM_ == ExpM);
	assert(size_ == v.size_ * v.size_);

	cblas_zcopy(size_, data_, 1, eve_, 1);
	cblas_zcopy(v.size_, v.data_, 1, tmpv_, 1);
	LAPACKE_zheev ( LAPACK_ROW_MAJOR , 'V' , 'U' ,  sqrt_size_ , eve_ , sqrt_size_ , eva_);

	cblas_zgemv ( CblasRowMajor, CblasConjTrans, v.size_, v.size_, &mkl_alpha, eve_, v.size_, 
		tmpv_, 1, &mkl_beta, v.data_, 1);
	
	for (int i = 0; i < v.size_; ++i)
	{
		tmpv_[i] = exp(alpha * eva_[i]) *  v(i);
	}
	cblas_zgemv ( CblasRowMajor, CblasNoTrans, v.size_, v.size_, &mkl_alpha, eve_, v.size_, 
	 					tmpv_, 1, &mkl_beta, v.data_, 1);
}

void nVector<complex<double>,1>::Write(ofstream & ofstr)
{
	ofstr.precision(15);
	ofstr<<scientific;
	ofstr<<showpos;
	for (int i = 0; i < size_; ++i)
	{
		ofstr<< real(data_[i])<<" ";
		ofstr<< imag(data_[i])<<" ";
	}
	ofstr<<endl;
}
void nVector<complex<double>,1>::Write(const string filename)
{
	ofstream ofstr;
	ofstr.open(filename, ios::out);
	Write(ofstr);
	ofstr.close();
}

void nVector<complex<double>,1>::Read(ifstream &ifstr)
{
	double re = 0.0;
	double im = 0.0;
	for (int i = 0; i < size_; ++i)
	{
		ifstr >> re;
		ifstr >> im;
		data_[i ] = re + I * im;
	}
}


void nVector<complex<double>,1>::Read(const string filename)
{
	ifstream ifstr;
	ifstr.open(filename, ios::in);
	Read(ifstr);
	ifstr.close();
}


/*********** nVector<double,1>******************/

nVector<double,1>::nVector(const int size_1)
{
	assert(size_1 > 0);
	size_1_ = size_1;
	size_   = size_1;
	data_   = nullptr;
	data_   = (double *)mkl_calloc( size_t(size_), sizeof( double), 64);
}
nVector<double,1>::nVector(nVector <double, 1> const & rhs)
{
	size_1_ = rhs.size_1_;
	size_   = rhs.size_;
	data_   = nullptr;
	data_   = (double *)mkl_calloc( size_t(size_), sizeof( double), 64);
	CopyTo(rhs);
}

nVector<double,1> & nVector<double,1>::operator= (nVector <double, 1> const & rhs)
{
	nVector <double, 1> other(rhs);
	swap(other);
	return (*this);
}

double & nVector<double,1>::operator() (const int i)
{
	assert(i >= 0); 
	assert(i < size_1_);
	return data_[i];
}

double const & nVector<double,1>::operator() (const int i) const
{
	assert(i >= 0); 
	assert(i < size_1_);
	return data_[i];
}

double *  nVector<double,1>::ptrTo()
{
	return data_;
}
double *  nVector<double,1>::ptrTo () const
{
	return data_;
}

void nVector<double,1>::swap(nVector <double, 1> & other)
{
	std::swap(size_1_,other.size_1_);
	std::swap(size_,other.size_);
	std::swap(data_,other.data_);
}
nVector<double,1>::~nVector()
{
	mkl_free(data_);
}
int nVector<double,1>::Size() const
{
	return size_;
}

void nVector<double,1>::CopyTo(nVector <double, 1> const & v)
{
	assert(size_ == v.size_);
	cblas_dcopy(size_, v.data_, 1, data_, 1);
}




/********* 2222222222 ***********/

nVector<complex<double>,2>::nVector(const int size_1, const int size_2)
{
	size_1_ = size_1;
	size_2_ = size_2;
	size_   = size_1 * size_2;
	data_   = nullptr;
	data_   = (complex<double> *)mkl_calloc( size_t(size_), sizeof( complex<double> ), 64);	
}

nVector<complex<double>,2>::nVector(nVector <complex<double>, 2> const & rhs)
{
	size_1_ = rhs.size_1_;
	size_2_ = rhs.size_2_;
	size_   = rhs.size_;
	data_   = nullptr;
	data_   = (complex<double> *)mkl_calloc( size_t(size_), sizeof( complex<double>), 64);
	cblas_zcopy(size_, rhs.data_, 1, data_, 1);
}

nVector<complex<double>,2>::nVector(const int size_1, nVector <complex<double>, 1> const & rhs)
{
	size_1_ = size_1;
	size_2_ = rhs.Size();
	size_   = size_1 * rhs.Size();
	data_   = nullptr;
	data_   = (complex<double> *)mkl_calloc( size_t(size_), sizeof( complex<double>), 64);
	for (int i = 0; i < size_1_; ++i)
	{
		cblas_zcopy(size_2_, &rhs(0), 1, &data_[i * size_2_], 1);
	}
	
}


nVector <complex<double>, 2> & nVector<complex<double>,2>::operator= (nVector <complex<double>,2> const & rhs)
{
	nVector <complex<double>, 2> other(rhs);
	swap(other);
	return (*this);
}
complex<double> & nVector<complex<double>,2>::operator() (const int i, const int j)
{
	assert(i >= 0); 
	assert(j >= 0);
	assert(i < size_1_);
	assert(j < size_2_);
	return data_[i * size_2_ + j];
}

complex<double> const & nVector<complex<double>,2>::operator() (const int i, const int j) const
{
	assert(i >= 0); 
	assert(j >= 0);
	assert(i < size_1_);
	assert(j < size_2_);
	return data_[i * size_2_ + j];
}

void nVector<complex<double>,2>::swap(nVector <complex<double>, 2> & other)
{
	std::swap(size_1_,other.size_1_);
	std::swap(size_2_,other.size_2_);
	std::swap(size_,other.size_);
	std::swap(data_,other.data_);
}
nVector<complex<double>,2>::~nVector()
{
	mkl_free(data_);
}

int nVector<complex<double>,2>::Size(const int dim) const
{ 
	assert((dim == 1) || (dim == 2) );
	if( dim == 1) {return size_1_;}
	return size_2_;
}

int nVector<complex<double>,2>::Size() const
{ 
	return size_;
}

void nVector<complex<double>,2>::CopyTo(nVector <complex<double>, 2> const & v)
{
	assert(size_ == v.size_);
	cblas_zcopy(size_, v.data_, 1, data_, 1);
}

void nVector<complex<double>,2>::Mul(const double alpha)
{
	cblas_zdscal(size_, alpha, data_, 1);
	//cout<<"It's meeeeee!\n";
}

void nVector<complex<double>,2>::Mul( const complex<double> alpha)
{
	cblas_zscal(size_, &alpha, data_, 1);
}

void nVector<complex<double>,2>::Mul(nVector <complex<double>, 2> const & v)
{
	assert(size_ == v.size_);
	vzMul( size_, v.data_, data_, data_ );
}

void nVector<complex<double>,2>::MulByConj(nVector <complex<double>, 2> const & v)
{
	assert(size_ == v.size_);
	vzMulByConj( size_, data_, v.data_, data_ );
}

void nVector<complex<double>,2>::Exp()
{
	vzExp(size_, data_, data_); 
}

void nVector<complex<double>,2>::Exp(const double alpha)
{
	Mul(alpha);
	vzExp(size_, data_, data_); 
}

void nVector<complex<double>,2>::Exp(const complex<double> alpha)
{
	Mul(alpha);
	vzExp(size_, data_, data_); 
}


double nVector<complex<double>,2>::Norm() const
{
	return cblas_dznrm2( size_ ,data_, 1);
}

double nVector<complex<double>,2>::Norm(int i) const
{
	assert(i >= 0);
	assert(i < size_1_);
	return cblas_dznrm2( size_2_ ,&data_[i * size_2_], 1);
}

void nVector<complex<double>,2>::Write(ofstream &ofstr)
{
	ofstr.precision(15);
	ofstr<<scientific;
	ofstr<<showpos;
	for (int i = 0; i < size_1_; ++i)
	{
		for (int j = 0; j < size_2_; ++j)
		{
			ofstr<< real(data_[i * size_2_ + j ])<<" ";
			ofstr<< imag(data_[i * size_2_ + j ])<<" ";
		}
		ofstr<<endl;
	}
}
void nVector<complex<double>,2>::Write(const string filename)
{
	ofstream ofstr;
	ofstr.open(filename, ios::out);
	Write(ofstr);
	ofstr.close();
}

void nVector<complex<double>,2>::Read(ifstream &ifstr)
{
	double re = 0.0;
	double im = 0.0;
	for (int i = 0; i < size_1_; ++i)
	{
		for (int j = 0; j < size_2_; ++j)
		{
			ifstr >> re;
			ifstr >> im;
			data_[i * size_2_ + j ] = re + I * im;
		}
	}
}


void nVector<complex<double>,2>::Read(const string filename)
{
	ifstream ifstr;
	ifstr.open(filename, ios::in);
	Read(ifstr);
	ifstr.close();
}



