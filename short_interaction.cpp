#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <cassert>
#include <valarray>

using namespace std;

#include "mathematical_constants.h"
#include "physical_constants.h"
#include "types.h"
#include "utility.h"
#include "nvector.h"
#include "short_interaction.h"

#include "mkl.h"
#include "gsl/gsl_sf_coupling.h"

UgTensor::UgTensor(const int SPIN, const double alpha):
	size_1_((2 * SPIN + 1) * (2 * SPIN + 1)),
	size_2_((2 * SPIN + 1) * (2 * SPIN + 1)),
	descreption
	{
		SPARSE_MATRIX_TYPE_SYMMETRIC,
		SPARSE_FILL_MODE_UPPER,
		SPARSE_DIAG_NON_UNIT
	}
	
{
	const int nps = 2 * SPIN + 1;
	nVector<complex<double>, 2> Ug{size_1_, size_2_};
	
	const double ac[] = {a0,0.0,
						 a2,0.0,
						 a4,0.0,
						 a6};

	int m1;
	int m2;
	int m1p;
	int m2p;
	int M;
	int Fmax = 0;
	double sum_F = 0.0f;
	double TMP = 0.0f;
	double sum_M = 0.0f;
	complex<double> nullcmplx{0.0,0.0};
	nonzero_upper_ = 0;
	switch(SPIN)
	{
		case SPIN_1 :
			Fmax = 2;
		    break;
		 
		case SPIN_2 :
			Fmax = 4;
		    break;

		case SPIN_3 :
			Fmax = 6;
		    break;
	}



	for (int i_m1 = 0; i_m1 < nps; ++i_m1)
	{
		for (int j_m2 = 0; j_m2 < nps; ++j_m2)
		{
			for (int i_m1p = 0; i_m1p < nps; ++i_m1p)
			{
				for (int j_m2p = 0; j_m2p < nps; ++j_m2p)
				{
					
					sum_F = 0.0;
					for (int F = 0; F <= Fmax; F = F + 2)
					{

						sum_M = 0.0;
						for (int Mp = 0; Mp < 2 * F + 1; ++Mp)
						{
							m1  = i_m1  - SPIN;
							m2  = j_m2  - SPIN;
							m1p = i_m1p - SPIN;
							m2p = j_m2p - SPIN;
							M   = Mp   - F;
					
							TMP  = gsl_sf_coupling_3j (2 * SPIN, 2 * SPIN, 2 * F, 2 * m1, 2 * m2, 2* (- M));
							TMP *= pow(-1, M) * sqrt(1.0 + 2.0 * F ) ;
							TMP *= gsl_sf_coupling_3j (2 * SPIN , 2 * SPIN, 2 * F , 2 * m1p, 2 * m2p, 2 * (- M));
							TMP *= pow(-1, M) * sqrt(1.0 + 2.0 * F ) ;
							sum_M += TMP; 
						}
						sum_F += ac[F] * sum_M;
					}
					Ug(ToRowMajorIndex(i_m1, i_m1p, nps), ToRowMajorIndex(j_m2, j_m2p, nps)) =  sum_F;
				}
			}
		}
	}



	for (int i = 0; i < size_1_; ++i)
	{
		for (int j = i; j < size_2_; ++j)
		{
			if(Ug(i,j) != nullcmplx) 
			{
				nonzero_upper_++;
			}
		}
	}
	

	value_upper_     = nullptr;
	rowindex_upper_  = nullptr;
	colindex_upper_  = nullptr;
	value_upper_     = (complex<double> *)mkl_calloc( size_t(nonzero_upper_), sizeof( complex<double> ), 64);
	rowindex_upper_  = (int *)mkl_calloc( size_t(nonzero_upper_), sizeof( int ), 64);
	colindex_upper_  = (int *)mkl_calloc( size_t(nonzero_upper_), sizeof( int ), 64);
	
	int l = 0;
	for (int i = 0; i < size_1_; ++i)
	{
		for (int j = i; j < size_2_; ++j)
		{
			if(Ug(i,j) != nullcmplx) 
			{
				value_upper_[l] = Ug(i,j);
				rowindex_upper_[l] = i;
				colindex_upper_[l] = j;
				l++;
			}
		}
	}
	cblas_zdscal(nonzero_upper_, alpha, value_upper_, 1);

	

	const int expected_calls = 100000;
	mkl_sparse_z_create_coo ( &spars, SPARSE_INDEX_BASE_ZERO, size_1_, 
		size_2_, nonzero_upper_, rowindex_upper_, colindex_upper_, value_upper_);

	mkl_sparse_set_mv_hint ( spars , SPARSE_OPERATION_NON_TRANSPOSE,
		 descreption , expected_calls );
	mkl_sparse_set_memory_hint ( spars , SPARSE_MEMORY_AGGRESSIVE );
	mkl_sparse_optimize ( spars );

}
UgTensor::~UgTensor()
{
 	mkl_free(value_upper_);
 	mkl_free(rowindex_upper_);
 	mkl_free(colindex_upper_);
 	mkl_sparse_destroy (spars);
}


