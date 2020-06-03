#ifndef TYPES_H
#define TYPES_H
//#define NDEBUG
//#define MKL_DIRECT_CALL 
#define USE_NCURSES
#define COMPLEX       complex<double>
#define MKL_Complex16 complex<double>
enum SpaceDimension 
{
	DIMENSION_1 = 1,
	DIMENSION_2 = 2,
	DIMENSION_3 = 3
};

enum SpinDimension
{
	SPIN_1 = 1,
	SPIN_2 = 2,
	SPIN_3 = 3
};


enum DynamcEvolution 
{ 
	RealTimeEvolution,
  	ImaginaryTimeEvolution
};

enum ComputeGroundState 
{ 
	DoImaginaryTimeEvolution,
	UseDataFile,
	UseDefaultGaussian
};

enum Axis
{
	x,
	y,
	z
};

enum DipolarInteractionKernel 
{ 
	DipolarInteraction,
	SpinConservingDipolarInteraction,
	NoDipolarInteraction
};

enum DoIt 
{ 
	No  = 0,
	Yes = 1
};




class ZeemanEffect
{
public:
	ZeemanEffect(double p_ = 0.0, double theta_ = 0.0, double phi_ = 0.0,
	double q_ = 0.0,double p_inf_ = 0.0,double tho_ = 0.0,double t0_ = 0.0)
	{
	    p           = p_;
		theta       = theta_;
		phi         = phi_;
		q           = q_;
		p_inf       = p_inf_;
		tho         = tho_;
		t0          = t0_;
	}
	~ZeemanEffect(){};
public:
	double p     ;
  	double theta;
  	double phi;
  	double q;
  	double p_inf     ;
  	double tho     ;
  	double t0     ;
};


class BfieldGradient
{
public:
	BfieldGradient(double gradient_x_ = 0.0, 
				   double gradient_y_ = 0.0,
				   double gradient_z_ = 0.0)
	{
		gradient_x = gradient_x_;
		gradient_y = gradient_y_;
		gradient_z = gradient_z_;
	}
	~BfieldGradient(){};
public:
  	double gradient_x;
  	double gradient_y;
  	double gradient_z;
};


class EulerAngles
{
public:
	EulerAngles(double alpha_ = 0.0, double beta_ = 0.0, double gamma_ = 0.0):
		alpha(alpha_),
		beta(beta_),
		gamma(gamma_){}
	~EulerAngles(){};
public:
  	const double alpha;
  	const double beta;
  	const double gamma;
};



class SpinRotation
{
public:
	SpinRotation(Axis axe_ = y, double angle_ = 0.0 )
	{
		axe   = axe_;
		angle = angle_;
	}
	~SpinRotation(){};
public:
	Axis axe;
  	double angle;	
};

class LossAtoms
{
public:
	LossAtoms(double zeeman_p = 0.0, DoIt instruction = No):
		p(zeeman_p),
		do_it(instruction)
	{	
	}
	~LossAtoms(){};
public:
	const double p;
	const DoIt do_it;
	
};


enum ComputeMatrixExponential 
{
	ExpM,
	NoExpM,
};

enum MatrixPart 
{
	UpperTriangularPart = 'U',
	LowerTriangularPart = 'L',
	All = 'A'
};

enum EnergyType 
{ 
	kinetic,
	potential,
	zeemanP,
	zeemanQ,
	contact,
  	dipolar,
  	total,
};


// struct FilesNames
// {
// 	string psi_m         = "./dynamic_results/psi_m.txt";
// 	string spin_density  = "./dynamic_results/spin_density.txt";
// 	string spin_dynamics = "./dynamic_results/spin_dynamics.txt";
// 	string energy        = "./dynamic_results/energy.txt";
// 	string dipolar_field = "./dynamic_results/dipolar_field.txt";
// };

#endif
