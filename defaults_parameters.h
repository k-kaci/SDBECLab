#include "types.h"
#include "mathematical_constants.h"
#ifndef DEFAULTS_PARAMETRES_H
#define DEFAULTS_PARAMETRES_H
/*********************************************************************/
/*                        System parameters                          */
/*********************************************************************/
SpinDimension   _SPIN_      = SPIN_3;
SpaceDimension  _DIMENSION_ = DIMENSION_3;
const int _NUMBER_OF_ATOMS_ = 800000;
const int _ms_ = -3;

const double _fx_ = 410.0;  /* en Hz*/ // ZEEMAN
const double _fy_ = 410.0;  /* en Hz*/ // VERTICAL
const double _fz_ = 410.0;  /* en Hz*/ // QUNATIF
const double _alpha_ = 0.0;  // (1) rotating Alpha around the current z axis
const double _beta_  = 0.0;  // (2) then Beta around the current y axis
const double _gamma_ = 0.0;  // (3) and then Gamma around the current z axis.
  
/*********************************************************************/
/*                          Grid parameters                          */
/*********************************************************************/
const int _nx_ = 64;
const int _ny_ = 64;
const int _nz_ = 64;

const double _lg_ = 25.0;   /* lg = max(lx,ly,lz)*/


const double _B_DD_       = 2.40779; //(h_bar omega)
const double _B_INF_      = 0.00 * _B_DD_;
const double _THO_        = 1.0;



const double _B0_         = _B_DD_ * 4.0; 
const double _T0_         = _THO_  * 4.0;
const double _T_MAX_      = _T0_   * 2.0 + 2.0; 



/*********************************************************************/
/*              Imaginary time evolution parameters                  */
/*********************************************************************/
const double _ITE_eps_ = 1e-10;
const double _ITE_dt_  = 0.5 * 1e-2;


const double _ITE_linear_zeeman_p_     = _B0_;
const double _ITE_linear_zeeman_theta_ = 0.0;
const double _ITE_linear_zeeman_phi_   = 0.0;
const double _ITE_quadratic_zeeman_q_  = 0.0;






ComputeGroundState  _GroundState_ = UseDataFile;
/*DoImaginaryTimeEvolution
UseDataFile
UseDefaultGaussian*/
const string _ITE_directory_  = "./ground_state";

/*********************************************************************/
/*                  Real time evolution parameters                   */
/*********************************************************************/

const double _RTE_TIME_MAX_ = _T_MAX_;           /*millisecond*/
const double _RTE_dt_       = 0.5 * 1e-3;

const double _RTE_linear_zeeman_p_     = _B0_;
const double _RTE_linear_zeeman_theta_ = 0.0;
const double _RTE_linear_zeeman_phi_   = 0.0;
const double _RTE_quadratic_zeeman_q_  = 0.0;

const double _RTE_linear_zeeman_p_inf_  = _B_INF_;
const double _RTE_linear_zeeman_tho_    = _THO_;
const double _RTE_linear_zeeman_t0_     = _T0_;






const double _RET_gradient_x_ =  0.0 ; // gradient suivant les axes du piege
const double _RET_gradient_y_ =  0.0 ; // en MHz/m
const double _RET_gradient_z_ =  0.0 ;

const double _RTE_SPIN_ROTATION_ANGLE_ = 0.0 * PI; 
Axis         _RTE_SPIN_ROTATION_AXIS_  = y; /*(x,y,z) = (0,1,2)*/

DoIt          _RTE_do_dipolar_relaxiation_loss_ = No;
const double  _RTE_do_dipolar_relaxiation_b0    = 0.0; /* en KHz*/

DipolarInteractionKernel _RTE_dipolar_kernel_ = DipolarInteraction;
/*DipolarInteraction,
SpinConservingDipolarInteraction,
NoDipolarInteraction*/

const string _RTE_directory_0_  = "systeme";
const string _RTE_directory_1_  = "results_1";
const string _RTE_directory_2_  = "results_2";

const int _RET_print_frequency_spin_dynamics_  = 10;
const int _RET_print_frequency_other_quantity_ = 1030;


#endif
