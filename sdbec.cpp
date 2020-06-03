#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>
#include <unistd.h>
#include <cassert>
#include <curses.h>


using namespace std;

#include "mathematical_constants.h"
#include "physical_constants.h"
#include "types.h"
#include "utility.h"
#include "nvector.h"
#include "fft.h"
#include "spin.h"
#include "short_interaction.h"
#include "sdbec.h"
#include "mkl.h"
#include <omp.h>


#ifdef USE_NCURSES
extern string msg_welcome;
extern string msg_in_progress;
extern string msg_Waiting;
extern string msg_completed;
extern WINDOW * simulation_status_win;
#endif



void GetRotation3D( nVector<complex<double>, 1> const & vector, nVector<complex<double>, 1>  & rotate_vector, EulerAngles euler )
{
	rotate_vector(0)  = vector(0) * (cos(euler.alpha) * cos(euler.beta) * cos(euler.gamma) - sin(euler.alpha) * sin(euler.gamma));
	rotate_vector(0) += vector(1) * (cos(euler.beta)  * cos(euler.gamma) * sin(euler.alpha) + cos(euler.alpha) * sin(euler.gamma));
	rotate_vector(0) += vector(2) * (-cos(euler.gamma) * sin(euler.beta));


	rotate_vector(1)  = vector(0) * (-cos(euler.gamma) * sin(euler.alpha) - cos(euler.alpha)  * cos(euler.beta) * sin(euler.gamma));
	rotate_vector(1) += vector(1) * (cos(euler.alpha) * cos(euler.gamma) - cos(euler.beta) * sin(euler.alpha)  * sin(euler.gamma));
	rotate_vector(1) += vector(2) * (sin(euler.beta) * sin(euler.gamma));

	rotate_vector(2)  = vector(0) * (cos(euler.alpha) * sin(euler.beta)); 
	rotate_vector(2) += vector(1) * (sin(euler.alpha) * sin(euler.beta)); 
	rotate_vector(2) += vector(2) * (cos(euler.beta)); 

}

sdbec::sdbec(SpinDimension spin, SpaceDimension dimension,
	const valarray<int> ngrid, const double delta, 
	const valarray<double> frequencies,EulerAngles trap_rotation_, const int atoms,
	const int ms):
	SPIN(spin),
	DIMENSION(dimension),
	nps(2 * SPIN + 1),
	npr(ngrid[ngrid.size() - 1]),
	npd(ngrid[slice(0,DIMENSION, 1)]),
	lgd(DIMENSION),
	dr(delta),

	nba(atoms),
	fr(frequencies),
	f0(fr.min()),
	trap_rotation(trap_rotation_),
	
	r_grid{DIMENSION,npd.max()},
 	k_grid{DIMENSION,npd.max()},

	psi_m_initial{nps, npr},
	psi_m_ground_state{nps, npr}
{ 
	

	for (int i = 0; i < DIMENSION; ++i)
	{
		lgd[i] = npd[i] * dr;
	}
	Grid();
	InitializePsiM(psi_m_initial, ms);
	InitializePsiM(psi_m_ground_state, ms);
	

	switch(DIMENSION)
	{
		case DIMENSION_1 :
			cout<<"sdbec::sdbec() \n<<"; 
			cout<<" ->> le cas 1D n'est pas encore disponible ...\n";
			assert(DIMENSION == DIMENSION_3);
			break;
		case DIMENSION_2 :
			cout<<"sdbec::sdbec() \n<<"; 
			cout<<" ->> le cas 2D n'est pas encore disponible ...\n";
			assert(DIMENSION == DIMENSION_3);
			l0 = sqrt(hBar / (mass * 2.0 * PI * f0));
			normalize_interaction = sqrt(8.0 * PI) / l0;
			cdd = normalize_interaction * add;
		    break;
		case DIMENSION_3 :
			l0 = sqrt(hBar / (mass * 2.0 * PI * f0));
			normalize_interaction = 4.0 * PI / l0;
			cdd =  normalize_interaction * add;
			break; 
	}

}
sdbec::~sdbec()
{
}

void sdbec::Grid()
{
	for (int d = 0; d < DIMENSION; ++d)
	{
		for (int i = 0; i < npd[d]; ++i)
		{
			r_grid(d,i) =  - (lgd[d] / 2.0) + double(i) * dr;
			if (i >=  npd[d]/2)
			{
				k_grid(d,i)  = 2.0 * PI * double(i -  npd[d]) / lgd[d] ;
				continue;
			}
			k_grid(d,i) = 2.0 * PI *  double(i) / lgd[d] ;
		}
	}
}

void sdbec::InitializePsiM(nVector<complex<double>, 2> & psi_m, const int m)
{
	complex<double> tmp_psi{0.0,0.0};
	nVector<complex<double>,1> r_vector{DIMENSION};
	nVector<complex<double>,1> r_vector_rotate{DIMENSION};



	switch(DIMENSION)
	{
		case DIMENSION_1 :
		   	for (int i = 0; i < npd[x]; ++i)
		   	{
		   		tmp_psi  = pow(r_grid(x,i), 2.0);
		   		tmp_psi *= -0.25;
		   		psi_m(ToSpinIndex(SPIN, m), i) = exp(tmp_psi);
		   	}
		    break;

		case DIMENSION_2 :
		    for (int i = 0; i < npd[x]; ++i)
		    {
		    	for (int j = 0; j < npd[y]; ++j)
		    	{
		    		tmp_psi  = pow(r_grid(x,i), 2.0);
		    		tmp_psi += pow(r_grid(y,j), 2.0);
		    		tmp_psi *= -0.25;
		    		psi_m(ToSpinIndex(SPIN, m), ToRowMajorIndex(i,j, npd[y])) = exp(tmp_psi);
		    	}
		    }
		    break;

		case DIMENSION_3 :
		    for (int i = 0; i < npd[x]; ++i)
		    {
		    	for (int j = 0; j < npd[y]; ++j)
		    	{
		    		for (int l = 0; l < npd[z]; ++l)
		    		{

		    			r_vector(x) = r_grid(x,i);
		    			r_vector(y) = r_grid(y,j);
		    			r_vector(z) = r_grid(z,l);
		    			
		    			GetRotation3D( r_vector, r_vector_rotate, trap_rotation );

		    			tmp_psi  =  pow( (fr[x] / f0) * r_vector_rotate(x), 2) ;
					    tmp_psi +=  pow( (fr[y] / f0) * r_vector_rotate(y), 2) ;
					    tmp_psi +=  pow( (fr[z] / f0) * r_vector_rotate(z), 2) ;
					    tmp_psi *= -1;
		    			psi_m(ToSpinIndex(SPIN, m), ToRowMajorIndex(i,j, npd[y], l, npd[z])) = exp(tmp_psi);
		    		}
		    	}
		    }
		    break;
	}



	NormalisePsiM(psi_m);	
}

void sdbec::InitializeKineticOperator(nVector<complex<double>, 1> & kinetic_operator)
{

	complex<double> tmp_hp{0.0,0.0};

	switch(DIMENSION)
	{
		case DIMENSION_1 :
			cout<<"sdbec::InitializeKineticOperator() \n<<"; 
			cout<<" ->> le cas 1D n'est pas encore disponible ...\n";
			assert(DIMENSION == DIMENSION_3);
		    break;

		case DIMENSION_2 :
			cout<<"sdbec::InitializeKineticOperator() \n<<"; 
			cout<<" ->> le cas 2D n'est pas encore disponible ...\n";
			assert(DIMENSION == DIMENSION_3);
			for (int i = 0; i < npd[x]; ++i)
			{
				for (int j = 0; j < npd[y]; ++j)
				{
					tmp_hp =  cos((2.0 * PI * i) / npd[x]);
					tmp_hp += cos((2.0 * PI * j) / npd[y]);
					tmp_hp += -2.0;
					tmp_hp *= -1.0 / (dr * dr);
					kinetic_operator(ToRowMajorIndex(i,j, npd[y])) = tmp_hp;
				}
			}
		    break;

		case DIMENSION_3 :
			for (int i = 0; i < npd[x]; ++i)
			{
				for (int j = 0; j < npd[y]; ++j)
				{
					for (int l = 0; l < npd[z]; ++l)
					{
						tmp_hp  = cos((2.0 * PI * i) / npd[x]);
						tmp_hp += cos((2.0 * PI * j )/ npd[y]);
						tmp_hp += cos((2.0 * PI * l )/ npd[z]);
						tmp_hp += -3.0;
						tmp_hp *= -1.0 / (dr * dr);
						kinetic_operator(ToRowMajorIndex(i,j, npd[y], l, npd[z])) = tmp_hp;
					}
				}
			}
		    break;      
	}
}

void sdbec::InitializeMomentumOperator(nVector<complex<double>, 2> & momentum_operator)
{
	complex<double> tmp_momentum_operator_x{0.0,0.0};
	complex<double> tmp_momentum_operator_y{0.0,0.0};
	complex<double> tmp_momentum_operator_z{0.0,0.0};

	switch(DIMENSION)
	{
		case DIMENSION_1 :
			cout<<"sdbec::initialize_MomentumOperator() \n<<"; 
			cout<<" ->> le cas 1D n'est pas encore disponible ...\n";
			assert(DIMENSION == DIMENSION_3);
		    break;

		case DIMENSION_2 :
			cout<<"sdbec::initialize_MomentumOperator() \n<<"; 
			cout<<" ->> le cas 2D n'est pas encore disponible ...\n";
			assert(DIMENSION == DIMENSION_3);
		    break;

		case DIMENSION_3 :
			for (int i = 0; i < npd[x]; ++i)
			{
				for (int j = 0; j < npd[y]; ++j)
				{
					for (int l = 0; l < npd[z]; ++l)
					{
						tmp_momentum_operator_x  = -(1.0/dr) * sin((2.0 * PI * i) / npd[x]);
						tmp_momentum_operator_y  = -(1.0/dr) * sin((2.0 * PI * j) / npd[y]);
						tmp_momentum_operator_z  = -(1.0/dr) * sin((2.0 * PI * l) / npd[z]);

						momentum_operator(x,ToRowMajorIndex(i,j, npd[y], l, npd[z]))  = tmp_momentum_operator_x;
						momentum_operator(y,ToRowMajorIndex(i,j, npd[y], l, npd[z]))  = tmp_momentum_operator_y;
						momentum_operator(z,ToRowMajorIndex(i,j, npd[y], l, npd[z]))  = tmp_momentum_operator_z;

					}
				}
			}
		    break;      
	}

}

void sdbec::InitializeHarmonicTrap(nVector<complex<double>, 1> & harmonic_trap)
{

	complex<double> tmp_ht{0.0,0.0};
	nVector<complex<double>,1> r_vector{DIMENSION};
	nVector<complex<double>,1> r_vector_rotate{DIMENSION};
	

	switch(DIMENSION)
	{
		case DIMENSION_1 :
			cout<<"sdbec::InitializeHarmonicTrap() \n<<"; 
			cout<<" ->> le cas 1D n'est pas encore disponible ...\n";
			assert(DIMENSION == DIMENSION_3);
			for (int i = 0; i < npd[x]; ++i)
			{
				tmp_ht  = pow( (fr[x] / f0) * r_grid(x,i) , 2.0);
				tmp_ht *= 0.5;
				harmonic_trap(i) = tmp_ht;		
			}
		    break;

		case DIMENSION_2 :
			cout<<"sdbec::InitializeHarmonicTrap() \n<<"; 
			cout<<" ->> le cas 2D n'est pas encore disponible ...\n";
			assert(DIMENSION == DIMENSION_3);
			for (int i = 0; i < npd[x]; ++i)
			{
				for (int j = 0; j < npd[y]; ++j)
				{
					tmp_ht  = pow( (fr[x] / f0) * r_grid(x,i) , 2.0);
					tmp_ht += pow( (fr[y] / f0) * r_grid(y,j) , 2.0);
					tmp_ht *= 0.5;
					harmonic_trap(ToRowMajorIndex(i,j, npd[y])) = tmp_ht;	
				}
			}
		    break;

		case DIMENSION_3 :
			for (int i = 0; i < npd[x]; ++i)
			{
				for (int j = 0; j < npd[y]; ++j)
				{
					for (int l = 0; l < npd[z]; ++l)
					{
						r_vector(x) = r_grid(x,i);
		    			r_vector(y) = r_grid(y,j);
		    			r_vector(z) = r_grid(z,l);
		    			
		    			GetRotation3D( r_vector, r_vector_rotate, trap_rotation );
					
					    tmp_ht  =  pow( (fr[x] / f0) * r_vector_rotate(x), 2) ;
					    tmp_ht +=  pow( (fr[y] / f0) * r_vector_rotate(y), 2) ;
					    tmp_ht +=  pow( (fr[z] / f0) * r_vector_rotate(z), 2) ;

						harmonic_trap(ToRowMajorIndex(i,j, npd[y], l, npd[z])) = 0.5 * tmp_ht;
					}
				}
			}
		    break;      
	}
}

void sdbec::InitializeLinearZeeman(nVector<complex<double>, 1> & linear_zeeman, const ZeemanEffect Zeeman)
{
	SpinMatrix Sx{SPIN, x};
	SpinMatrix Sy{SPIN, y};
	SpinMatrix Sz{SPIN, z};
	complex<double> tmp_linear_zeeman{0.0,0.0};
	for (int m = 0; m < nps; ++m)
	{
		for (int mp = 0; mp < nps; ++mp)
		{
			tmp_linear_zeeman  = sin(Zeeman.theta) * cos(Zeeman.phi) * Sx(ToRowMajorIndex(m, mp, nps));
			tmp_linear_zeeman += sin(Zeeman.theta) * sin(Zeeman.phi) * Sy(ToRowMajorIndex(m, mp, nps));
			tmp_linear_zeeman += cos(Zeeman.theta) * Sz(ToRowMajorIndex(m, mp, nps));
			linear_zeeman(ToRowMajorIndex(m, mp, nps)) = Zeeman.p * tmp_linear_zeeman; 
		}
	}
}
void sdbec::InitializeIBfieldGradient(nVector<complex<double>, 2> & magnetic_field_gradient, const BfieldGradient BGradient)
{

	complex<double> tmp_magnetic_field_gradient{0.0,0.0};
	
	nVector<complex<double>,1> tmp_bgrad_x{nps * nps};
	nVector<complex<double>,1> tmp_bgrad_y{nps * nps};
	nVector<complex<double>,1> tmp_bgrad_z{nps * nps};
	
	ZeemanEffect zeeman_x{BGradient.gradient_x * l0 / f0};
	ZeemanEffect zeeman_y{BGradient.gradient_y * l0 / f0};
	ZeemanEffect zeeman_z{BGradient.gradient_z * l0 / f0};

	InitializeLinearZeeman(tmp_bgrad_x,zeeman_x);
	InitializeLinearZeeman(tmp_bgrad_y,zeeman_y);
	InitializeLinearZeeman(tmp_bgrad_z,zeeman_z);

	nVector<complex<double>,1> r_vector{DIMENSION};
	nVector<complex<double>,1> r_vector_rotate{DIMENSION};
	
	for (int i = 0; i < npd[x]; ++i)
	{
		for (int j = 0; j < npd[y]; ++j)
		{
			for (int l = 0; l < npd[z]; ++l)
			{
				r_vector(x) = r_grid(x,i) + 0.0 * lgd[x]/2.0;
				r_vector(y) = r_grid(y,j) + 0.0 * lgd[y]/2.0;
				r_vector(z) = r_grid(z,l) + 0.0 * lgd[z]/2.0;
				    			
				GetRotation3D( r_vector, r_vector_rotate, trap_rotation );

				for (int m = 0; m < nps; ++m)
				{
					for (int mp = 0; mp < nps; ++mp)
					{
						tmp_magnetic_field_gradient  = tmp_bgrad_x(ToRowMajorIndex(m, mp, nps)) * r_vector_rotate(x) ;
						tmp_magnetic_field_gradient += tmp_bgrad_y(ToRowMajorIndex(m, mp, nps)) * r_vector_rotate(y) ;
						tmp_magnetic_field_gradient += tmp_bgrad_z(ToRowMajorIndex(m, mp, nps)) * r_vector_rotate(z) ;
						magnetic_field_gradient(ToRowMajorIndex(i,j, npd[y], l, npd[z]),ToRowMajorIndex(m, mp, nps)) = tmp_magnetic_field_gradient;
					}
					
				}
			}
		}
	}
}

void sdbec::InitializeQuadraticZeeman(nVector<complex<double>, 1> & quadratic_zeeman, const ZeemanEffect Zeeman)
{
	SpinMatrix Sz{SPIN, z};
	quadratic_zeeman.Mul(double(0.0));
	for (int m = 0; m < nps; ++m)
	{
		quadratic_zeeman(ToRowMajorIndex(m, m, nps)) = Zeeman.q * Sz(ToRowMajorIndex(m,m,nps)) * Sz(ToRowMajorIndex(m,m,nps)); 
	}
}
void sdbec::InitializeDipolarKernel(nVector<complex<double>, 2> & dipolar_kernel, DipolarInteractionKernel DipolarKernel)
{
	
	complex<double> k{0.0,0.0};
	complex<double> k2{0.0,0.0};
	complex<double> A{0.0,0.0};
	complex<double> R{lgd[0]/2.0,0.0};
	complex<double> cmplxnull{0.0,0.0};

	dipolar_kernel.Mul(0.0);

	for (int i = 0; i < npd[x]; ++i)
	{
		for (int j = 0; j <  npd[y]; ++j)
		{
			for (int l = 0; l < npd[z]; ++l)
			{
				k2  = k_grid(x,i) *  k_grid(x,i);
				k2 += k_grid(y,j) *  k_grid(y,j);
				k2 += k_grid(z,l) *  k_grid(z,l);
				k   = sqrt(k2);
				
				if(k != cmplxnull)
				{
					A  = 1.0/3.0;
					A += cos(k * R) / pow(k * R, 2);
					A -= sin(k * R) / pow(k * R, 3);

					switch(DipolarKernel)
					{
						case DipolarInteraction :
							dipolar_kernel(ToRowMajorIndex(x,x, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (1.0 - 3.0 * k_grid(x,i) * k_grid(x,i) / k2);
							dipolar_kernel(ToRowMajorIndex(y,y, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (1.0 - 3.0 * k_grid(y,j) * k_grid(y,j) / k2);
							dipolar_kernel(ToRowMajorIndex(z,z, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (1.0 - 3.0 * k_grid(z,l) * k_grid(z,l) / k2);

							dipolar_kernel(ToRowMajorIndex(x,y, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (0.0 - 3.0 * k_grid(x,i) * k_grid(y,j)/ k2 );
							dipolar_kernel(ToRowMajorIndex(x,z, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (0.0 - 3.0 * k_grid(x,i) * k_grid(z,l)/ k2 );

							dipolar_kernel(ToRowMajorIndex(y,x, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (0.0 - 3.0 * k_grid(y,j) * k_grid(x,i)/ k2 );
							dipolar_kernel(ToRowMajorIndex(y,z, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (0.0 - 3.0 * k_grid(y,j) * k_grid(z,l)/ k2 );

							dipolar_kernel(ToRowMajorIndex(z,x, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (0.0 - 3.0 * k_grid(z,l) * k_grid(x,i)/ k2 );
							dipolar_kernel(ToRowMajorIndex(z,y, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  A * (- 4.0 * PI) * (0.0 - 3.0 * k_grid(z,l) * k_grid(y,j)/ k2 );
							break;
						
						case SpinConservingDipolarInteraction :
							dipolar_kernel(ToRowMajorIndex(x,x, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  (-1.0/2.0) * A * (- 4.0 * PI) * (1.0 - 3.0 * k_grid(z,l) * k_grid(z,l) / k2);
							dipolar_kernel(ToRowMajorIndex(y,y, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  (-1.0/2.0) * A * (- 4.0 * PI) * (1.0 - 3.0 * k_grid(z,l) * k_grid(z,l) / k2);
							dipolar_kernel(ToRowMajorIndex(z,z, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =               A * (- 4.0 * PI) * (1.0 - 3.0 * k_grid(z,l) * k_grid(z,l) / k2);
						
						case NoDipolarInteraction :
							 //Do Nothing
						    break;   
					}

					
				}
				else
				{
					
					dipolar_kernel(ToRowMajorIndex(x,x, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;
					dipolar_kernel(ToRowMajorIndex(y,y, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;
					dipolar_kernel(ToRowMajorIndex(z,z, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;

					dipolar_kernel(ToRowMajorIndex(x,y, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;
					dipolar_kernel(ToRowMajorIndex(x,z, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;

					dipolar_kernel(ToRowMajorIndex(y,x, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;
					dipolar_kernel(ToRowMajorIndex(y,z, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;

					dipolar_kernel(ToRowMajorIndex(z,x, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;
					dipolar_kernel(ToRowMajorIndex(z,y, 3), ToRowMajorIndex(i,j, npd[y], l, npd[z])) =  0.0;
				}

			}
		}
	}

}

void sdbec::NormalisePsiM(nVector<complex<double>, 2> & psi_m)
{
	psi_m.Mul(sqrt(nba / GetNormPsiM(psi_m)));
}

double sdbec::GetNormPsiM(nVector<complex<double>, 2> const & psi_m) const
{
	return pow(psi_m.Norm(), 2) * pow(dr, DIMENSION);
}


void sdbec::GetNumberOfAtoms(nVector<complex<double>, 2> const & psi_m, nVector<complex<double>, 1> & nba_m)
{
	assert(nba_m.Size() == nps + 1);
	nba_m(nps) = 0.0;
	for (int m = 0; m < psi_m.Size(1); ++m)
	{
		nba_m(m) = pow(psi_m.Norm(m), 2) * pow(dr, DIMENSION);
		nba_m(nps) += nba_m(m);
	}
}

void sdbec::GetOrbitalAngularMomentum(nVector<complex<double>, 2> const & psi_m,
	nVector<complex<double>, 2> const & fft_psi_m,
	nVector<complex<double>, 1>  & orbital_angular_momentum)
{
	assert(orbital_angular_momentum.Size() == nps + 1);

	nVector<complex<double>, 2> momentum_operator{3,npr};
	InitializeMomentumOperator(momentum_operator);
	nVector<complex<double>, 2> px_psi_m{nps,npr};
	nVector<complex<double>, 2> py_psi_m{nps,npr};
	nVector<complex<double>, 2> lz_psi_m{nps,npr};
	
	
	complex<double>  y_px_psi_m{0.0,0.0};
	complex<double>  x_py_psi_m{0.0,0.0};
	

	FastFourierTransform FFT_SPACE{DIMENSION, npd, npr, nps};

	for (int m = 0; m < nps; ++m)
	{
		vzMul(npr, &momentum_operator(x,0),&fft_psi_m(m,0),&px_psi_m(m,0));
		vzMul(npr, &momentum_operator(y,0),&fft_psi_m(m,0),&py_psi_m(m,0));
	}

	FFT_SPACE.InverseFourier(&px_psi_m(0,0));
	FFT_SPACE.InverseFourier(&py_psi_m(0,0));
	
	orbital_angular_momentum(nps) = 0.0;
	for (int m = 0; m < nps; ++m)
	{
		for (int i = 0; i < npd[x]; ++i)
		{
			for (int j = 0; j < npd[y]; ++j)
			{
				for (int l = 0; l < npd[z]; ++l)
				{
					y_px_psi_m = px_psi_m(m,ToRowMajorIndex(i,j, npd[y], l, npd[z])) * r_grid(y,j);
					x_py_psi_m = py_psi_m(m,ToRowMajorIndex(i,j, npd[y], l, npd[z])) * r_grid(x,i);
					lz_psi_m(m,ToRowMajorIndex(i,j, npd[y], l, npd[z])) = x_py_psi_m - y_px_psi_m;
				}
			}
		}

		 cblas_zdotc_sub (npr , &psi_m(m,0) , 1 , &lz_psi_m(m,0) , 1 , &orbital_angular_momentum(m));
		 orbital_angular_momentum(m)   *= pow(dr, DIMENSION);
		 orbital_angular_momentum(nps) += orbital_angular_momentum(m);
	}

}

void sdbec::GetSpinAngularMomentum(nVector<complex<double>, 2> const & psi_m,
		nVector<complex<double>, 1> & spin_angular_momentum)
{

	assert(spin_angular_momentum.Size() == nps + 1);
	spin_angular_momentum(nps) = 0.0;
	for (int m = 0; m < psi_m.Size(1); ++m)
	{
		spin_angular_momentum(m) = pow(psi_m.Norm(m), 2) * pow(dr, DIMENSION) * (SPIN - m);
		spin_angular_momentum(nps) += spin_angular_momentum(m);
	}

}

void sdbec::DoSpinRotation(nVector<complex<double>, 2>  & psi_m, const SpinRotation spin_rotation)
{

	SpinMatrix Sx{SPIN, x};
	SpinMatrix Sy{SPIN, y};
	SpinMatrix Sz{SPIN, z};
	nVector<complex<double>, 1> psi_x{nps};
	nVector<complex<double>, 1> J_ROT{nps * nps, ExpM};
	
	switch(spin_rotation.axe)
	{
		case x :
			for (int i = 0; i < nps * nps; ++i)
			{
				J_ROT(i) = Sx(i);
			}
			 break;
		case y :
			for (int i = 0; i < nps * nps; ++i)
			{
				J_ROT(i) = Sy(i);
			}
		    break;
		case z :
			for (int i = 0; i < nps * nps; ++i)
			{
				J_ROT(i) = Sz(i);
			}
		    break;    
	}

	for (int i = 0; i < npr; ++i)
	{
		for (int m = 0; m < nps; ++m)
		{
			psi_x(m) = psi_m(m,i);
		}

		J_ROT.MatrixForm_HMatrixExp(- I * spin_rotation.angle , psi_x);

		for (int m = 0; m < nps; ++m)
		{
			psi_m(m,i) = psi_x(m);
		}	
	}
}

void sdbec::InitializeLossMatrix(nVector<complex<double>, 1> & loss_matrix,
	LossAtoms loss_atoms)
{

	double cdd_phys = (4.0 * PI * hBar * hBar / mass) * add;
	assert(loss_matrix.Size() == nps * nps);
	loss_matrix.Mul(0.0);
	if(loss_atoms.do_it == Yes)
	{
		loss_matrix.Read("./loss_matrix.txt");
		double tmp;
		tmp = pow(cdd_phys,2.0) * pow(mass/2.0,2.0) / pow(hBar,4.0);
		tmp *= 4.0 * sqrt(loss_atoms.p * hBar * 2.0 * PI / mass);
		tmp  /=  pow(l0,3.0) * (2.0 * PI * f0);
		loss_matrix.Mul(tmp);

	}


}
void sdbec::GroundState(
	const ZeemanEffect Zeeman,
	DipolarInteractionKernel DipolarKernel,
	const double dt,
	const double stop, FilesNames files_names)
{
	SpinRotation pulse;
	LossAtoms loss_atoms;
	PrintParametres print_parametres;
	BfieldGradient BGradient;
	Dynamic(pulse,loss_atoms,Zeeman, BGradient, DipolarKernel, dt, stop, files_names, print_parametres, ImaginaryTimeEvolution);
}
void sdbec::GroundState(FilesNames files_names)
{
	
	psi_m_ground_state.Read(files_names.file_psi_m);

}



void sdbec::Dynamic( const SpinRotation spin_rotation,
	const LossAtoms loss_atoms,
	const ZeemanEffect Zeeman,
	const BfieldGradient BGradient,
	DipolarInteractionKernel DipolarKernel,
	const double dt,
	const double stop,
	FilesNames files_names,
	PrintParametres print_parametres,
	DynamcEvolution DYNAMIC)
{

	nVector<complex<double>,2> psi_m{nps, npr};	
	nVector<complex<double>,1> kinetic_operator{npr};
	nVector<complex<double>,1> harmonic_trap{npr};
	nVector<complex<double>,1> linear_zeeman{nps * nps};
	nVector<complex<double>,2> magnetic_field_gradient{npr, nps * nps};
	nVector<complex<double>,1> quadratic_zeeman{nps * nps};
	
	
	nVector<complex<double>,1> Hr{nps * nps, ExpM};
	nVector<complex<double>,1> Hg{nps * nps};
	nVector<complex<double>,1> Hd{nps * nps};
	UgTensor Ug{SPIN, normalize_interaction};
	nVector<complex<double>,2> dipolar_kernel{9, npr};
	nVector<complex<double>,1> loss_matrix{nps * nps};

	nVector<complex<double>,2> spin_density{3, npr};
	nVector<complex<double>,2> dipolar_field{3, npr};

	nVector<complex<double>,1> energy{total + 1};
	nVector<complex<double>,1> numbre_of_atoms{nps + 1};


	InitializeKineticOperator(kinetic_operator);
	InitializeHarmonicTrap(harmonic_trap);
	InitializeLinearZeeman(linear_zeeman,Zeeman);
	InitializeIBfieldGradient(magnetic_field_gradient, BGradient);
	InitializeQuadraticZeeman(quadratic_zeeman, Zeeman);
	InitializeDipolarKernel(dipolar_kernel, DipolarKernel);
	InitializeLossMatrix(loss_matrix, loss_atoms);



	SpinMatrix Sx{SPIN, x};
	SpinMatrix Sy{SPIN, y};
	SpinMatrix Sz{SPIN, z};

	SparsSpinMatrix sparsSx{SPIN, x};
	SparsSpinMatrix sparsSy{SPIN, y};
	SparsSpinMatrix sparsSz{SPIN, z};

	nVector<complex<double>,2> expHp(nps, kinetic_operator);
	nVector<complex<double>,2> expHp_saved(nps, kinetic_operator);
	nVector<complex<double>,2> fft_psi_m{nps, npr};
	nVector<complex<double>,1> psi_x{nps};
	nVector<complex<double>,1> psi_density_x{nps};
	nVector<complex<double>,1> V{nps * nps};
	nVector<complex<double>,1> W{npr};
	nVector<complex<double>,2> Wx{nps, npr};
	nVector<complex<double>,2> Wy{nps, npr};
	nVector<complex<double>,2> Wz{nps, npr};
	nVector<complex<double>, 1> tmp_energy{total + 1};


	FastFourierTransform FFT_SPACE_FOR_PSI{DIMENSION, npd, npr,nps};
	FastFourierTransform FFT_SPACE_FOR_spin_density{DIMENSION, npd, npr, 3};
	

	complex<double> Idt{0.0,0.0};
	bool CONTUNUE = true;
	const int TIME_MAX = int (stop / ((dt/ (2.0 * PI * f0)) * 1e3))  + 1;
 




	ofstream out_spin_dynamics;
	ofstream out_magnetic_field;
	ofstream out_psi_m;
	ofstream out_psi_m_ground_state;

	switch(DYNAMIC)
	{
		case RealTimeEvolution :
			Idt = - I * dt;
			expHp.Exp(Idt * 0.5 );
			psi_m = psi_m_ground_state;

			DoSpinRotation(psi_m, spin_rotation);

			out_spin_dynamics.open(files_names.file_spin_dynamics, ios::out);
			out_magnetic_field.open(files_names.file_magnetic_field, ios::out);
			out_psi_m.open(files_names.file_psi_m, ios::out);

			
			kinetic_operator.Write(files_names.file_kinetic_operator);
			harmonic_trap.Write(files_names.file_harmonic_trap);
			linear_zeeman.Write(files_names.file_linear_zeeman);
			quadratic_zeeman.Write(files_names.file_quadratic_zeeman);
			dipolar_kernel.Write(files_names.file_dipolar_kernel);



			 break;
		case ImaginaryTimeEvolution :
			Idt = - 1.0 * dt;
			expHp.Exp(Idt * 0.5 );
			//psi_m = psi_m_initial;
			psi_m.Read(files_names.file_psi_m);


			kinetic_operator.Write(files_names.file_kinetic_operator);
			harmonic_trap.Write(files_names.file_harmonic_trap);
			linear_zeeman.Write(files_names.file_linear_zeeman);
			quadratic_zeeman.Write(files_names.file_quadratic_zeeman);
			dipolar_kernel.Write(files_names.file_dipolar_kernel);
		    break;    
	}


	
	double mu    = 0.0;
	double muold = 0.0;
	double dmu   = 1.0;
	int    TIME  = 0;

	
	ZeemanEffect zeeman_dynamic;
	double delta_time;
	double times = 0.0;
	nVector<complex<double>,1> magnetic_field{1};
	


	print_parametres.n_print_spin_dynamics  = 0;
	print_parametres.n_print_other_quantity = 0;

	
	const MKL_Complex16 mkl_alpha = 1.0;
	const MKL_Complex16 mkl_beta  = 0.0;


	
	do 
	{



		switch(DYNAMIC)
		{
			case RealTimeEvolution :
				#ifdef USE_NCURSES
				mvwprintw(simulation_status_win,3, 3,"Dynamic ............... %s  (time = %2.2lf ms) ",msg_in_progress.c_str(), TIME * dt * 1e3 / (2.0 * PI * f0) );
				wrefresh(simulation_status_win); 
				#endif
				GetNumberOfAtoms(psi_m, numbre_of_atoms);
				times = TIME * dt * 1e3 / (2.0 * PI * f0);
				if(times <= Zeeman.t0)
				{
					zeeman_dynamic.p = (Zeeman.p - Zeeman.p_inf) * exp(- pow(times/Zeeman.tho,3.0) / (1.0 + times/Zeeman.tho ))  + Zeeman.p_inf;
				}
				
				if (times > Zeeman.t0 &&  times <= 2.0 * Zeeman.t0)
				{
					delta_time = Zeeman.t0 - (times - Zeeman.t0);
					zeeman_dynamic.p = (Zeeman.p - Zeeman.p_inf) * exp(- pow(delta_time/Zeeman.tho,3.0) / (1.0 + delta_time/Zeeman.tho ))  + Zeeman.p_inf;
				}
				
				
				InitializeLinearZeeman(linear_zeeman,zeeman_dynamic.p);
				magnetic_field(0) = zeeman_dynamic.p;
				if (TIME % print_parametres.f_print_spin_dynamics == 0)
				{
					numbre_of_atoms.Write(out_spin_dynamics);
					magnetic_field.Write(out_magnetic_field);
					print_parametres.n_print_spin_dynamics ++;

					
				}
				if (TIME % print_parametres.f_print_other_quantity == 0)
				{
					psi_m.Write(out_psi_m);
					print_parametres.n_print_other_quantity ++;
				}
				break;
			case ImaginaryTimeEvolution :
				#ifdef USE_NCURSES
				mvwprintw(simulation_status_win,2, 3,"Ground State .......... %s (eps = %1.2lE > %1.2lE)", msg_in_progress.c_str(), dmu, stop);
				wrefresh(simulation_status_win); 
				#endif
			    break;   
		}
	
		
			
		



		

		{
			FFT_SPACE_FOR_PSI.Fourier(&psi_m(0,0));
			psi_m.Mul(expHp);
			FFT_SPACE_FOR_PSI.InverseFourier(&psi_m(0,0));

			spin_density.Mul(0.0);
			mkl_sparse_z_mm ( SPARSE_OPERATION_NON_TRANSPOSE , 1.0 , sparsSx.spars , sparsSx.descreption ,
				SPARSE_LAYOUT_ROW_MAJOR , &psi_m(0,0) , npr , npr , 0.0 , &Wx(0,0), npr );
			mkl_sparse_z_mm ( SPARSE_OPERATION_NON_TRANSPOSE , 1.0 , sparsSy.spars , sparsSy.descreption ,
				SPARSE_LAYOUT_ROW_MAJOR , &psi_m(0,0) , npr , npr , 0.0 , &Wy(0,0) , npr );
			mkl_sparse_z_mm ( SPARSE_OPERATION_NON_TRANSPOSE , 1.0 , sparsSz.spars , sparsSz.descreption ,
				SPARSE_LAYOUT_ROW_MAJOR , &psi_m(0,0) , npr , npr , 0.0 , &Wz(0,0) , npr );
			Wx.MulByConj(psi_m);
			Wy.MulByConj(psi_m);
			Wz.MulByConj(psi_m);
			for (int m = 0; m < nps; ++m)
			{
				vzAdd( npr, &Wx(m,0), &spin_density(x,0), &spin_density(x,0));
				vzAdd( npr, &Wy(m,0), &spin_density(y,0), &spin_density(y,0));
				vzAdd( npr, &Wz(m,0), &spin_density(z,0), &spin_density(z,0));
			}


			




			FFT_SPACE_FOR_spin_density.Fourier(&spin_density(0,0));
			dipolar_field.Mul(0.0);
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					vzMul (npr , &dipolar_kernel(ToRowMajorIndex(i,j,3),0) , &spin_density(j,0), &W(0));
					vzAdd( npr, &W(0), &dipolar_field(i,0), &dipolar_field(i,0));
				}
			}
			FFT_SPACE_FOR_spin_density.InverseFourier(&dipolar_field(0,0));

			

			
			for (int i = 0; i < npr; ++i)
			{
				
				for (int m = 0; m < nps; ++m)
				{
					psi_x(m) = psi_m(m,i);
					psi_density_x(m) = conj(psi_m(m,i)) * psi_m(m,i);;
				}
				
				
				for (int m2 = 0; m2 < nps; ++m2)
				{
					for (int m2p = 0; m2p < nps; ++m2p)
					{
						V(ToRowMajorIndex( m2, m2p, nps)) =  conj(psi_x(m2)) * psi_x(m2p);
					}
				}
				
				mkl_sparse_z_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, Ug.spars, Ug.descreption,
				 &V(0), 0.0, &Hg(0));

				Hd.Mul(0.0);
				for (int m = 0; m < nps-1; ++m)
				{
					Hd(ToRowMajorIndex(m,m + 1,nps))  = cdd * real(dipolar_field(x,i)) * Sx(ToRowMajorIndex(m,m + 1,nps));
					Hd(ToRowMajorIndex(m,m + 1,nps)) += cdd * real(dipolar_field(y,i)) * Sy(ToRowMajorIndex(m,m + 1,nps));
					Hd(ToRowMajorIndex(m,m,nps))      = cdd * real(dipolar_field(z,i)) * Sz(ToRowMajorIndex(m,m,nps));
				}
				Hd(ToRowMajorIndex(nps-1,nps-1,nps))  = cdd * real(dipolar_field(z,i)) * Sz(ToRowMajorIndex(nps-1,nps-1,nps));

				


				Hr.MatrixForm_Initialize(harmonic_trap(i), double(0), UpperTriangularPart);
				Hr.MatrixForm_AddTo(linear_zeeman, UpperTriangularPart );
				Hr.MatrixForm_AddTo(quadratic_zeeman, UpperTriangularPart );
				Hr.MatrixForm_AddTo(Hg, UpperTriangularPart );
				Hr.MatrixForm_AddTo(Hd, UpperTriangularPart );


				

				for (int m = 0; m < nps; ++m)
				{
					for (int mp = m; mp < nps; ++mp)
					{
						Hr(ToRowMajorIndex(m, mp, nps )) += magnetic_field_gradient(i,ToRowMajorIndex(m, mp, nps ));
					}
				}

	

				
				Hr.MatrixForm_HMatrixExp(Idt, psi_x);


				 
				for (int m = 0; m < nps; ++m)
				{
					psi_m(m,i) = psi_x(m);
				}

				if (loss_atoms.do_it == Yes)
				{
					cblas_zgemv ( CblasRowMajor, CblasNoTrans, nps, nps, &mkl_alpha, &loss_matrix(0), nps, 
					&psi_density_x(0), 1, &mkl_beta, &psi_x(0), 1);
					
					for (int m = 0; m < nps; ++m)
					{
						psi_m(m,i) *= exp(- I * Idt * psi_x(m));
					}	
				}	
			}
			
			FFT_SPACE_FOR_PSI.Fourier(&psi_m(0,0));
			psi_m.Mul(expHp);
			FFT_SPACE_FOR_PSI.InverseFourier(&psi_m(0,0));
		}
		
		

		switch(DYNAMIC)
		{
			case RealTimeEvolution :
				CONTUNUE = (TIME < TIME_MAX);
				break;
			case ImaginaryTimeEvolution :
				mu =  -( 1.0 / (2.0 * dt)) * log( GetNormPsiM(psi_m) / nba);
				dmu = fabs((mu - muold)/mu);
				muold = mu;
				NormalisePsiM(psi_m);
				CONTUNUE = (dmu > stop);
			    break;   
		}





		
		TIME++;

	}while(CONTUNUE);

	
	
	if(DYNAMIC == ImaginaryTimeEvolution)
	{
		psi_m_ground_state = psi_m;
		psi_m_ground_state.Write(files_names.file_psi_m);
		print_parametres.n_print_other_quantity ++;
	}

	if(DYNAMIC == RealTimeEvolution)
	{

		out_spin_dynamics.close();
		out_magnetic_field.close();
		out_psi_m.close();
	}

	

}

void sdbec::PrintResults(FilesNames files_names, PrintParametres print_parametres)
{

	nVector<complex<double>, 2> psi_m{nps, npr};
	
	nVector<complex<double>, 1> kinetic_operator{npr};
	nVector<complex<double>, 1> harmonic_trap{npr};
	nVector<complex<double>,1> linear_zeeman{nps * nps};
	nVector<complex<double>,1> quadratic_zeeman{nps * nps};
	nVector<complex<double>,1> Hg{nps * nps};
	nVector<complex<double>,1> Hd{nps * nps};
	UgTensor Ug{SPIN, normalize_interaction};
	nVector<complex<double>, 2> dipolar_kernel{9, npr};

	nVector<complex<double>, 2> spin_density{3, npr};
	nVector<complex<double>, 2> dipolar_field{3, npr};
	nVector<complex<double>,1> orbital_angular_momentum_z{nps +1};
	nVector<complex<double>,1> spin_angular_momentum_z{nps +1};

	nVector<complex<double>, 1> energy{total + 1};

	nVector<complex<double>,1> local_spin_length{npr};
	nVector<complex<double>,1> mean_spin_length{1};
	nVector<complex<double>,1> magnetization{4};


	SparsSpinMatrix sparsSx{SPIN, x};
	SparsSpinMatrix sparsSy{SPIN, y};
	SparsSpinMatrix sparsSz{SPIN, z};


	nVector<complex<double>,2> fft_psi_m{nps, npr};
	nVector<complex<double>,1> V{nps * nps};
	nVector<complex<double>,1> W{npr};
	nVector<complex<double>,2> Wx{nps, npr};
	nVector<complex<double>,2> Wy{nps, npr};
	nVector<complex<double>,2> Wz{nps, npr};
	nVector<complex<double>,1> tmp_energy{total + 1};





	FastFourierTransform FFT_SPACE_FOR_PSI{DIMENSION, npd, npr,nps};
	FastFourierTransform FFT_SPACE_FOR_spin_density{DIMENSION, npd, npr, 3};


	kinetic_operator.Read(files_names.file_kinetic_operator);
	harmonic_trap.Read(files_names.file_harmonic_trap);
	linear_zeeman.Read(files_names.file_linear_zeeman);
	quadratic_zeeman.Read(files_names.file_quadratic_zeeman);
	dipolar_kernel.Read(files_names.file_dipolar_kernel);

	

	ifstream in_psi_m;
	ofstream out_energy;
	ofstream out_spin_density;
	ofstream out_dipolar_field;
	ofstream out_orbtal_angular_momentum_z;
	ofstream out_spin_angular_momentum_z;
	ofstream out_local_spin_length;
	ofstream out_mean_spin_length;
	ofstream out_magnetization;
	in_psi_m.open(files_names.file_psi_m, ios::in);
	out_energy.open(files_names.file_energy, ios::out);
	out_spin_density.open(files_names.file_spin_density, ios::out);
	out_dipolar_field.open(files_names.file_dipolar_field, ios::out);
	out_orbtal_angular_momentum_z.open(files_names.file_orbital_angular_momentum, ios::out);
	out_spin_angular_momentum_z.open(files_names.file_spin_angular_momentum, ios::out);
	out_local_spin_length.open(files_names.file_local_spin_length, ios::out);
	out_mean_spin_length.open(files_names.file_mean_spin_length, ios::out);
	out_magnetization.open(files_names.file_magnetization, ios::out);

	

	for (int TIME = 1; TIME <= print_parametres.n_print_other_quantity; ++TIME)
	{
		psi_m.Read(in_psi_m);
		{


			spin_density.Mul(0.0);
			mkl_sparse_z_mm ( SPARSE_OPERATION_NON_TRANSPOSE , 1.0 , sparsSx.spars , sparsSx.descreption ,
				SPARSE_LAYOUT_ROW_MAJOR , &psi_m(0,0) , npr , npr , 0.0 , &Wx(0,0), npr );
			mkl_sparse_z_mm ( SPARSE_OPERATION_NON_TRANSPOSE , 1.0 , sparsSy.spars , sparsSy.descreption ,
				SPARSE_LAYOUT_ROW_MAJOR , &psi_m(0,0) , npr , npr , 0.0 , &Wy(0,0) , npr );
			mkl_sparse_z_mm ( SPARSE_OPERATION_NON_TRANSPOSE , 1.0 , sparsSz.spars , sparsSz.descreption ,
				SPARSE_LAYOUT_ROW_MAJOR , &psi_m(0,0) , npr , npr , 0.0 , &Wz(0,0) , npr );
			Wx.MulByConj(psi_m);
			Wy.MulByConj(psi_m);
			Wz.MulByConj(psi_m);
			for (int m = 0; m < nps; ++m)
			{
				vzAdd( npr, &Wx(m,0), &spin_density(x,0), &spin_density(x,0));
				vzAdd( npr, &Wy(m,0), &spin_density(y,0), &spin_density(y,0));
				vzAdd( npr, &Wz(m,0), &spin_density(z,0), &spin_density(z,0));
			}

			spin_density.Write(out_spin_density);


			local_spin_length.Mul(double(0.0));
			mean_spin_length.Mul(double(0.0));
			magnetization.Mul(double(0.0));
			for (int i = 0; i < npr; ++i)
			{
				for (int r = x; r <= z; ++r)
				{
					local_spin_length(i) += spin_density(r,i) * spin_density(r,i);
					magnetization(r)     +=  spin_density(r,i) * pow(dr, DIMENSION);
				}
				local_spin_length(i)  = sqrt(local_spin_length(i));
				mean_spin_length(0)  += local_spin_length(i) * pow(dr, DIMENSION);
			}
			for (int r = x; r <= z; ++r)
				{
					magnetization(3)     +=  magnetization(r) * magnetization(r);
				}
			magnetization(3)     =  sqrt( magnetization(3) );
			
			local_spin_length.Write(out_local_spin_length);
			mean_spin_length.Write(out_mean_spin_length);
			magnetization.Write(out_magnetization);
			


			FFT_SPACE_FOR_spin_density.Fourier(&spin_density(0,0));
			dipolar_field.Mul(0.0);
			for (int i = x; i <= z; ++i)
			{
				for (int j = x; j <= z; ++j)
				{
					vzMul (npr , &dipolar_kernel(ToRowMajorIndex(i,j,3),0) , &spin_density(j,0), &W(0));
					vzAdd( npr, &W(0), &dipolar_field(i,0), &dipolar_field(i,0));
				}
			}
			FFT_SPACE_FOR_spin_density.InverseFourier(&spin_density(0,0));
			FFT_SPACE_FOR_spin_density.InverseFourier(&dipolar_field(0,0));
			dipolar_field.Write(out_dipolar_field);



			fft_psi_m = psi_m;
			FFT_SPACE_FOR_PSI.Fourier(&fft_psi_m(0,0));

			
			
			GetOrbitalAngularMomentum(psi_m, fft_psi_m, orbital_angular_momentum_z);
			GetSpinAngularMomentum(psi_m, spin_angular_momentum_z);				

			orbital_angular_momentum_z.Write(out_orbtal_angular_momentum_z);
			spin_angular_momentum_z.Write(out_spin_angular_momentum_z);


			energy.Mul(double(0.0));
			
			
			for (int i = 0; i < npr; ++i)
			{	
			
				

				


				for (int m = 0; m < nps; ++m)
				{
					for (int mp = 0; mp < nps; ++mp)
					{
						V(ToRowMajorIndex( m, mp, nps)) =  conj(psi_m(m,i)) * psi_m(mp,i);
					}
				}
				
				mkl_sparse_z_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, Ug.spars, Ug.descreption,
				 &V(0), 0.0, &Hg(0));
				
				tmp_energy.Mul(double(0.0));
				for (int m = 0; m < nps; ++m)
				{
					tmp_energy(kinetic)   += kinetic_operator(i) *  conj(fft_psi_m(m,i)) * fft_psi_m(m,i);
					tmp_energy(potential) += harmonic_trap(i) *  conj(psi_m(m,i)) * psi_m(m,i);
				}
				cblas_zdotu_sub ( nps * nps, &linear_zeeman(0), 1, &V(0), 1, &tmp_energy(zeemanP));
				cblas_zdotu_sub ( nps * nps, &quadratic_zeeman(0), 1, &V(0), 1, &tmp_energy(zeemanQ));
				cblas_zdotu_sub ( nps * nps, &Hg(0), 1, &V(0), 1, &tmp_energy(contact));
				
				for (int r = x; r < z; ++r)
				{
					tmp_energy(dipolar)   += cdd * spin_density(r,i) * dipolar_field(r,i);
					
				}

				energy(kinetic)   += tmp_energy(kinetic)   * pow(dr, DIMENSION);
				energy(potential) += tmp_energy(potential) * pow(dr, DIMENSION);
				energy(zeemanP)   += tmp_energy(zeemanP)   * pow(dr, DIMENSION);
				energy(zeemanQ)   += tmp_energy(zeemanQ)   * pow(dr, DIMENSION);
				energy(contact)   += tmp_energy(contact)   * pow(dr, DIMENSION) * complex<double>(0.5);
				energy(dipolar)   += tmp_energy(dipolar)   * pow(dr, DIMENSION) * complex<double>(0.5);
				
				energy(total) = 0.0;
				for (int enrj = 0; enrj < total; ++ enrj)
				{
					energy(total) += energy(enrj);
				}
				
			}
			energy.Write(out_energy);
			
		}
		
	}

	
	in_psi_m.close();
	out_energy.close();
	out_spin_density.close();
	out_dipolar_field.close();
	out_orbtal_angular_momentum_z.close();
	out_spin_angular_momentum_z.close();
	out_local_spin_length.close();
	out_mean_spin_length.close();
	out_magnetization.close();
}
