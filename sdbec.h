
#include "mathematical_constants.h"
#include "physical_constants.h"
#include "types.h"
#include "utility.h"
#include "nvector.h"
#include "short_interaction.h"
#include "ncurses_parameters.h"



#ifndef SDBEC_H
#define SDBEC_H

class sdbec : private NonCopyable
{
public:
	sdbec(SpinDimension spin, SpaceDimension dimension,
	const valarray<int> ngrid, const double delta, 
	const valarray<double> frequencies, EulerAngles trap_rotation_, const int atoms,
	const int ms);
	~sdbec();
	void Grid();
	void InitializePsiM(nVector<complex<double>, 2> & psi_m, const int m);
	void InitializeKineticOperator(nVector<complex<double>, 1> & kinetic_operator);
	void InitializeMomentumOperator(nVector<complex<double>, 2> & momentum_operator);
	void InitializeHarmonicTrap(nVector<complex<double>, 1>  & harmonic_trap);
	void InitializeLinearZeeman(nVector<complex<double>, 1> & linear_zeeman, const ZeemanEffect Zeeman);
	void InitializeIBfieldGradient(nVector<complex<double>, 2> & magnetic_field_gradient,
	 	const BfieldGradient BGradient);
	void InitializeQuadraticZeeman(nVector<complex<double>, 1> & quadratic_zeeman, const ZeemanEffect Zeeman);
	void InitializeDipolarKernel(nVector<complex<double>, 2> & Qd, DipolarInteractionKernel DipolarKernel);
	void InitializeLossMatrix(nVector<complex<double>, 1> & loss_matrix,
	LossAtoms loss_atoms);
	void NormalisePsiM(nVector<complex<double>, 2> & psi_m);
	double GetNormPsiM(nVector<complex<double>, 2> const & psi_m) const;
	void GetNumberOfAtoms(nVector<complex<double>, 2> const & psi_m, nVector<complex<double>, 1> & nba_m);
	void GetOrbitalAngularMomentum(nVector<complex<double>, 2> const & psi_m,
		nVector<complex<double>, 2> const & fft_psi_m,
		nVector<complex<double>, 1>  & orbital_angular_momentum);
	void GetSpinAngularMomentum(nVector<complex<double>, 2> const & psi_m,
			nVector<complex<double>, 1> & spin_angular_momentum);
	void DoSpinRotation(nVector<complex<double>, 2>  & psi_m, const SpinRotation spin_rotation);
	void Dynamic(const SpinRotation spin_rotation,
		const LossAtoms loss_atoms, 
		const ZeemanEffect Zeeman,
		const BfieldGradient BGradient,
		DipolarInteractionKernel DipolarKernel,
		const double dt,
		const double stop, FilesNames files_names, PrintParametres print_parametres, DynamcEvolution DYNAMIC = RealTimeEvolution);
	void GroundState(
		const ZeemanEffect Zeeman,
		DipolarInteractionKernel DipolarKernel,
		const double dt,
		const double stop, FilesNames files_names);
	void GroundState(FilesNames files_names);
	void PrintResults(FilesNames files_names, PrintParametres print_parametres);
	void EverythingIsOK();
private:

	SpinDimension SPIN;
	SpaceDimension DIMENSION;

	const int   nps;
	const int   npr;
	const valarray<int> npd;
	valarray<double> lgd;
	const double dr;
	

	const int nba;
	const valarray<double> fr;
	const double f0;
	EulerAngles trap_rotation;
	
	double l0;
	double normalize_interaction;
	double cdd;
	
	nVector<complex<double>, 2> r_grid;
	nVector<complex<double>, 2> k_grid;

	
	nVector<complex<double>, 2> psi_m_initial;
	nVector<complex<double>, 2> psi_m_ground_state;
	




};



#endif