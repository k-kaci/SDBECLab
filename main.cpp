#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cassert>
#include <curses.h>

using namespace std;

#include "types.h"
#include "defaults_parameters.h"
#include "ncurses_parameters.h"
#include "utility.h"
#include "nvector.h"
#include "sdbec.h"
#include "spin.h"
#include "mkl.h"

#include <omp.h>



#ifdef USE_NCURSES
string msg_welcome     ="Spinor-dipolar BEC Simulation";
string msg_in_progress = "In progress";
string msg_waiting     = "Waiting";
string msg_completed   = "Completed.";
WINDOW * simulation_status_win = nullptr;
#endif

	
string DoubleToStr(double n) 
{
    stringstream result;
    result << n;
    return result.str();
}



int main()
{


	/*********************************************************************/
	/*                        System parameters                          */
	/*********************************************************************/
	SpinDimension  SPIN      = _SPIN_;
	SpaceDimension DIMENSION = _DIMENSION_;
	const valarray<double> frequencies = {_fx_, _fy_,_fz_};
	EulerAngles _trap_rotation_{_alpha_,_beta_,_gamma_};
	const int atoms = _NUMBER_OF_ATOMS_;
	const int ms = _ms_;
	/*********************************************************************/
	/*                          Grid parameters                          */
	/*********************************************************************/
	valarray<int>  ngrid(DIMENSION + 1);
	ngrid[x] = _nx_;
	ngrid[y] = _ny_;
	ngrid[z] = _nz_;
	ngrid[3] = ngrid[x] * ngrid[y] * ngrid[z];
	const double lg = _lg_;
	const double dr =  lg / ngrid[x];
	/*********************************************************************/
	/*              Imaginary time evolution parameters                  */
	/*********************************************************************/
	const double ITE_dt  = _ITE_dt_;
	const double ITE_eps =_ITE_eps_; 

	ZeemanEffect ITE_Zeeman{
		_ITE_linear_zeeman_p_,
	 	_ITE_linear_zeeman_theta_,
	 	_ITE_linear_zeeman_phi_,
	 	_ITE_quadratic_zeeman_q_};
	
	FilesNames files_names_ground_state{
		_ITE_directory_, _ITE_directory_};

	PrintParametres ITE_print_parametres;
	

	/*********************************************************************/
	/*                  Real time evolution parameters                   */
	/*********************************************************************/
	const double RTE_TIME_MAX = _RTE_TIME_MAX_;
	const double RTE_dt  = _RTE_dt_;
	
	ZeemanEffect RET_Zeeman{
		_RTE_linear_zeeman_p_,
	 	_RTE_linear_zeeman_theta_,
	 	_RTE_linear_zeeman_phi_,
	 	_RTE_quadratic_zeeman_q_,
	 	_RTE_linear_zeeman_p_inf_,
	 	_RTE_linear_zeeman_tho_,
	 	_RTE_linear_zeeman_t0_};
	
	BfieldGradient RET_BGradient{
		_RET_gradient_x_* 1e6,
		_RET_gradient_y_* 1e6,
		_RET_gradient_z_* 1e6};
	
	SpinRotation RTE_SPIN_ROTATION{
		_RTE_SPIN_ROTATION_AXIS_,
		_RTE_SPIN_ROTATION_ANGLE_};
	
	DipolarInteractionKernel RTE_dipolar_kernel = _RTE_dipolar_kernel_;
	
	FilesNames files_names_dynamics{
		_RTE_directory_0_,_RTE_directory_1_,_RTE_directory_2_};
	
	PrintParametres print_parametres{
		_RET_print_frequency_spin_dynamics_,
		_RET_print_frequency_other_quantity_};

	#ifdef USE_NCURSES
	initscr();
	
	mvprintw(0,(COLS - msg_welcome.size())/2,"%s",msg_welcome.c_str());
	
	mvprintw(get_row() + 1,1," -System parameters:");
	mvprintw(get_row() + 1,2,"   Spin ..................... %d", SPIN);
	mvprintw(get_row() + 1,2,"   Dimension ................ %d", DIMENSION);
	mvprintw(get_row() + 1,2,"   Number of atoms .......... %d", atoms);
	mvprintw(get_row() + 1,2,"   Harmonic trap ............ (%d,",int(frequencies[x]));
	mvprintw(get_row(),get_col(),"%d,%d) Hz",int(frequencies[y]), int(frequencies[z]));
	

	mvprintw(get_row() + 2,1," -Grid parameters:");
	mvprintw(get_row() + 1,2,"   (nx,ny,nz) .............. (%d,",ngrid[x]);
	mvprintw(get_row(),get_col(),"%d,%d)",ngrid[y], ngrid[z]);
	mvprintw(get_row() + 1,2,"   dr ...................... %1.2f", dr);
	

	mvprintw(get_row() + 2,1," -Imaginary time evolution parameters:");
	mvprintw(get_row() + 1,2,"  Zeeman effect -->");
	mvprintw(get_row() + 1,4,"  p ........................ %3.2lf",ITE_Zeeman.p);
	mvprintw(get_row() + 1,4,"  theta .................... %3.2lf",ITE_Zeeman.theta);
	mvprintw(get_row() + 1,4,"  phi ...................... %3.2lf",ITE_Zeeman.phi);
	mvprintw(get_row() + 1,4,"  q ........................ %3.2lf",ITE_Zeeman.q);	

	mvprintw(get_row() + 2,1," -Real time evolution parameters:");
	mvprintw(get_row() + 1,2,"  Zeeman effect -->");
	mvprintw(get_row() + 1,4,"  p ........................ %3.2lf",RET_Zeeman.p);
	mvprintw(get_row() + 1,4,"  theta .................... %3.2lf",RET_Zeeman.theta);
	mvprintw(get_row() + 1,4,"  phi ...................... %3.2lf",RET_Zeeman.phi);
	mvprintw(get_row() + 1,4,"  q ........................ %3.2lf",RET_Zeeman.q);	
	mvprintw(get_row() + 1,2,"  Spin rotation (pulse) -->");
	mvprintw(get_row() + 1,4,"  angle .................... %3.2lf PI",RTE_SPIN_ROTATION.angle/ PI );
	mvprintw(get_row() + 1,4,"  axis  .................... %3.2lf",RTE_SPIN_ROTATION.axe);
	refresh();

	simulation_status_win = newwin(4, COLS, get_row() + 1, 0);

	mvwprintw(simulation_status_win,1, 1,"-Simulation status");
	mvwprintw(simulation_status_win,2, 3,"Ground State .......... %s",msg_waiting.c_str() );
	mvwprintw(simulation_status_win,3, 3,"Dynamic ............... %s",msg_waiting.c_str());
	wrefresh(simulation_status_win);
	#endif

	
	


	// /*********** SIMULATION ***********/
	sdbec simulation{SPIN, DIMENSION, ngrid, dr, frequencies, _trap_rotation_, atoms, ms};
	
	#ifdef USE_NCURSES
	mvwprintw(simulation_status_win,2, 3,"Ground State .......... %s",msg_in_progress.c_str() );
	#endif

	switch(_GroundState_)
	{
		case DoImaginaryTimeEvolution :
			simulation.GroundState(ITE_Zeeman, DipolarInteraction,
			   	ITE_dt, ITE_eps, files_names_ground_state);
			simulation.GroundState(ITE_Zeeman, DipolarInteraction,
			   	ITE_dt, ITE_eps, files_names_ground_state);
			ITE_print_parametres.n_print_other_quantity = 1;
			simulation.PrintResults(files_names_ground_state, ITE_print_parametres);
			break;
		case UseDataFile :
			simulation.GroundState(files_names_ground_state);
			break;
		case UseDefaultGaussian:
			break;   
	}
	#ifdef USE_NCURSES
	wmove(simulation_status_win,2, 0);
	wclrtoeol(simulation_status_win);
	wmove(simulation_status_win,3, 0);
	wclrtoeol(simulation_status_win);
	
	mvwprintw(simulation_status_win,2, 3,"Ground State .......... %s",msg_completed.c_str() );
	mvwprintw(simulation_status_win,3, 3,"Dynamic ............... %s",msg_in_progress.c_str());
	wrefresh(simulation_status_win);
	#endif
	LossAtoms loss_atoms{_RTE_do_dipolar_relaxiation_b0 * 1e3, _RTE_do_dipolar_relaxiation_loss_};
	simulation.Dynamic(RTE_SPIN_ROTATION, loss_atoms, RET_Zeeman, RET_BGradient, RTE_dipolar_kernel, 
		RTE_dt, RTE_TIME_MAX, files_names_dynamics, print_parametres);
	
	print_parametres.n_print_other_quantity = _RTE_TIME_MAX_ * 5.0; 
	simulation.PrintResults(files_names_dynamics, print_parametres);
	
	#ifdef USE_NCURSES
	wmove(simulation_status_win,3, 1);
	wclrtoeol(simulation_status_win);
	mvwprintw(simulation_status_win,3, 3,"Dynamic ............... %s",msg_completed.c_str());
	wrefresh(simulation_status_win);
	mvprintw(LINES -1,0,"This screen has %d rows and %d columns\n",LINES,COLS);
	
	getch();               
	endwin();
	#endif
	

	return 0;
}
