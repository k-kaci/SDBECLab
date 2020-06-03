#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>
#include <cassert>

using namespace std;
#include "types.h"
#include "utility.h"
#include "nvector.h"
#include "fft.h"

#include "mkl.h"
FastFourierTransform::FastFourierTransform(
	SpaceDimension dimension,
	const valarray<int> npd,
	int input_distance,
	int number_of_transform)
{
	long   dimension_;
	long * npd_ = nullptr;
	
	switch(dimension)
	{
		case DIMENSION_1 :
		dimension_ = 1;
		npd_ = new long[1];
		npd_[x] = npd[x];
		DftiCreateDescriptor( &FFT_SPACE, DFTI_DOUBLE, DFTI_COMPLEX, dimension_, npd_[0]);
			
			break;
		case DIMENSION_2 :
		dimension_ = 2;
		npd_ = new long[2];
		npd_[x] = npd[x];
		npd_[y] = npd[y];
		DftiCreateDescriptor( &FFT_SPACE, DFTI_DOUBLE, DFTI_COMPLEX, dimension_, npd_);
			
		    break;
		case DIMENSION_3 : 
		dimension_ = 3;
		npd_ = new long[3];
		npd_[x] = npd[x];
		npd_[y] = npd[y];
		npd_[z] = npd[z];
		DftiCreateDescriptor( &FFT_SPACE, DFTI_DOUBLE, DFTI_COMPLEX, dimension_, npd_);  
	}

	DftiSetValue( FFT_SPACE, DFTI_FORWARD_SCALE, 1.0  / sqrt(input_distance));
	DftiSetValue( FFT_SPACE, DFTI_BACKWARD_SCALE, 1.0 / sqrt(input_distance));
	DftiSetValue( FFT_SPACE, DFTI_NUMBER_OF_TRANSFORMS, number_of_transform);
	DftiSetValue( FFT_SPACE, DFTI_INPUT_DISTANCE, input_distance);
	DftiCommitDescriptor( FFT_SPACE );
	delete[] npd_;
}

FastFourierTransform::~FastFourierTransform()
{
	DftiFreeDescriptor(&FFT_SPACE);
}