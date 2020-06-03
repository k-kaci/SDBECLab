#include "utility.h"
#include "types.h"
#include "mkl.h"

#ifndef SPIN_H
#define SPIN_H
class SpinMatrix : private NonCopyable
{
public:
	SpinMatrix(const int s,  Axis r);
	complex<double> & operator() (const int);
	complex<double> const & operator() (const int i) const;
	~SpinMatrix();
private:
	complex<double> * S;
	const int nps;
};

class SparsSpinMatrix : private NonCopyable
{
public:
	SparsSpinMatrix(const int s,  Axis r);
	~SparsSpinMatrix();
private:
	const int nps_;
	int   nonzero_;
	int * rowindex_;
	int * colindex_;
	complex<double> * value_;
public:
	sparse_matrix_t spars;
	matrix_descr descreption;
};

#endif
