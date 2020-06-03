#include "types.h"
#include "utility.h"
#include "nvector.h"
#include "mkl.h"
#ifndef FFT_H
#define FFT_H
class FastFourierTransform : private NonCopyable 
{
public:
	FastFourierTransform(
		SpaceDimension dimension,
		const valarray<int> npd,
		int input_distance,
		int number_of_transform = 1);
	~FastFourierTransform();
	inline void Fourier(complex<double> * data)
	{
		DftiComputeBackward(FFT_SPACE, data);
	}
	inline void InverseFourier(complex<double> * data)
	{
		DftiComputeForward(FFT_SPACE, data);
	}
private:
	DFTI_DESCRIPTOR_HANDLE FFT_SPACE;
};

#endif


