#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <cassert>

using namespace std;
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "types.h"
#include "utility.h"
#include "spin.h"

#include "mkl.h"

SpinMatrix::SpinMatrix(const int s, Axis r):
	nps(2 * s + 1)
{
	S = nullptr;
	S = (complex<double> *)mkl_calloc( size_t(nps * nps), sizeof( complex<double> ), 64);
	complex<double> scal{0,0};
	complex<double> * Sp = nullptr;
	complex<double> * Sm = nullptr;
	Sp = (complex<double> *)mkl_calloc( size_t(nps * nps), sizeof( complex<double> ), 64);
	Sm = (complex<double> *)mkl_calloc( size_t(nps * nps), sizeof( complex<double> ), 64);
	for (int m = -s; m < s ; m++)
	{
		Sp[ToRowMajorIndex(ToSpinIndex(s,m +1), ToSpinIndex(s,m), nps)] =  sqrt((s - m)*(s + 1 + m));
	}

	for (int m = -s + 1; m <= s ; m++)
	{
		Sm[ToRowMajorIndex(ToSpinIndex(s,m -1), ToSpinIndex(s,m), nps)] =  sqrt((s + m)*(s + 1 - m));
	}

	switch(r)
	{
		case x :
			scal = 0.5;
			vzAdd(nps * nps, Sp, Sm, S);
			cblas_zscal(nps * nps, &scal, S, 1);
		    break;
		case y:
		 	scal = -0.5 * I;
		 	vzSub(nps * nps, Sp, Sm, S);
		 	cblas_zscal(nps * nps, &scal, S, 1);
		 	break;
		case z :
		 	for (int m = s; m >= - s; --m)
		 	{
		 		S[ToRowMajorIndex(ToSpinIndex(s,m),ToSpinIndex(s,m),nps)] = m;
		 	}
		    break;
	}

	mkl_free(Sp);
	mkl_free(Sm);
}

complex<double> const & SpinMatrix::operator() (const int i) const
{
	assert(i >= 0); 
	assert(i < nps * nps);
	return S[i];
}

complex<double>  & SpinMatrix::operator() (const int i) 
{
	assert(i >= 0); 
	assert(i < nps * nps);
	return S[i];
}

SpinMatrix::~SpinMatrix(){}


SparsSpinMatrix::SparsSpinMatrix(const int s, Axis r):
	nps_(2 * s + 1)
{

	SpinMatrix S{s, r};
	complex<double> nullcmplx{0.0,0.0};
	
	matrix_descr Sx_descreption
	{
		SPARSE_MATRIX_TYPE_SYMMETRIC,
		SPARSE_FILL_MODE_UPPER,
		SPARSE_DIAG_NON_UNIT
	};

	matrix_descr Sy_descreption
	{
		SPARSE_MATRIX_TYPE_HERMITIAN,
		SPARSE_FILL_MODE_UPPER,
		SPARSE_DIAG_NON_UNIT
	};

	matrix_descr Sz_descreption
	{
		SPARSE_MATRIX_TYPE_DIAGONAL,
		SPARSE_FILL_MODE_UPPER,
		SPARSE_DIAG_NON_UNIT
	};
	

	
	switch(r)
	{
		case x :
			descreption = Sx_descreption;
			nonzero_ = nps_ - 1;
			break;
		case y:
		 	descreption = Sy_descreption;
		 	nonzero_ = nps_ - 1;
		 	break;
		case z :
			descreption = Sz_descreption;
			nonzero_ = nps_;
		    break;
	}
	//nonzero_ = nps_ - 1;
	value_     = nullptr;
	rowindex_  = nullptr;
	colindex_  = nullptr;
	value_     = (complex<double> *)mkl_calloc( size_t(nonzero_), sizeof( complex<double> ), 64);
	rowindex_  = (int *)mkl_calloc( size_t(nonzero_), sizeof( int ), 64);
	colindex_  = (int *)mkl_calloc( size_t(nonzero_), sizeof( int ), 64);

	int i = 0;
	if(r == x || r == y)
	for (int m = 0; m < nps_; ++m)
	{
		for (int mp = m; mp < nps_; ++mp)
		{
			if(S(ToRowMajorIndex(m,mp,nps_)) != nullcmplx)
			{
				value_[i]    = S(ToRowMajorIndex(m,mp,nps_));
				rowindex_[i] = m;
				colindex_[i] = mp;
				i++;
			}
		}
	}
	
	if(r == z)
	for (int m = 0; m < nps_; ++m)
	{
		value_[m]    = S(ToRowMajorIndex(m,m,nps_));
		rowindex_[m] = m;
		colindex_[m] = m;	
	}


	const int expected_calls = 100000;
	mkl_sparse_z_create_coo ( &spars, SPARSE_INDEX_BASE_ZERO, nps_, 
		nps_, nonzero_, rowindex_, colindex_, value_);

	 mkl_sparse_set_mv_hint ( spars , SPARSE_OPERATION_NON_TRANSPOSE,
		 descreption , expected_calls );
	mkl_sparse_set_memory_hint ( spars , SPARSE_MEMORY_AGGRESSIVE );
	mkl_sparse_optimize ( spars );

}

SparsSpinMatrix::~SparsSpinMatrix()
{
	mkl_free(rowindex_);
	mkl_free(colindex_);
	mkl_free(value_);
	mkl_sparse_destroy (spars);
}



