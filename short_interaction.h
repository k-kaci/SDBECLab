
#include "physical_constants.h"
#include "types.h"
#include "utility.h"
#include "mkl.h"

#ifndef SHORT_INTERACTION_H
#define SHORT_INTERACTION_H
class UgTensor : private NonCopyable
{
public:
	UgTensor(const int, const double alpha);
	~UgTensor();
private:
	const int size_1_;
	const int size_2_;
	int   nonzero_upper_;
	int * rowindex_upper_;
	int * colindex_upper_;
	complex<double> * value_upper_;
public:
	sparse_matrix_t spars;
	matrix_descr descreption;

};

#endif
