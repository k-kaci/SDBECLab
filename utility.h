#ifndef UTILITY_H
#define UTILITY_H

inline int ToRowMajorIndex(const int i, const int j, const int nj, const int l, const int nl)
{
	return (i * nj +j)* nl + l;
}

inline int ToRowMajorIndex(const int i, const int j, const int nj)
{
	return (i * nj) + j;
}

inline int ToSpinIndex(const int s, const int m)
{
	return s - m;
}








class NonCopyable
{
public:
	NonCopyable(NonCopyable const &) = delete;
	NonCopyable & operator = (NonCopyable const &) = delete;
protected:
	NonCopyable(){}
	~NonCopyable(){}
};

/*
class FilesNames
{
public:
	FilesNames(const string directory):
	file_psi_m                    (directory + "/psi_m"         + ".txt"),
	file_spin_density             (directory + "/spin_density"  + ".txt"),
	file_spin_dynamics            (directory + "/spin_dynamics" + ".txt"),
	file_spin_angular_momentum    (directory + "/spin_angular_momentum" + ".txt"),
	file_orbital_angular_momentum (directory + "/orbital_angular_momentum" + ".txt"),
	file_energy                   (directory + "/energy"        + ".txt"),
	file_dipolar_field            (directory + "/dipolar_field" + ".txt"),
	file_local_spin_length        (directory + "/local_spin_length" + ".txt"),
	file_mean_spin_length         (directory + "/mean_spin_length" + ".txt"),
	file_magnetization            (directory + "/magnetization" + ".txt"),

	


	file_kinetic_operator  (directory + "/kinetic_operator"        + ".txt"),
	file_harmonic_trap     (directory + "/harmonic_trap"           + ".txt"),
	file_linear_zeeman     (directory + "/linear_zeeman"           + ".txt"),
	file_quadratic_zeeman  (directory + "/quadratic_zeeman"        + ".txt"),
	file_dipolar_kernel    (directory + "/dipolar_kernel"          + ".txt")
	{
	}
	FilesNames(const string directory, const string suffix):
	file_psi_m                    (directory + "/psi_m_"         + suffix + ".txt"),
	file_spin_density             (directory + "/spin_density_"  + suffix + ".txt"),
	file_spin_dynamics            (directory + "/spin_dynamics_" + suffix + ".txt"),
	file_spin_angular_momentum    (directory + "/spin_angular_momentum" + suffix + ".txt"),
	file_orbital_angular_momentum (directory + "/orbital_angular_momentum" + suffix + ".txt"),
	file_energy                   (directory + "/energy_"        + suffix + ".txt"),
	file_dipolar_field            (directory + "/dipolar_field_" + suffix + ".txt"),
	file_local_spin_length        (directory + "/local_spin_length_" + suffix + ".txt"),
	file_mean_spin_length         (directory + "/mean_spin_length_" + suffix + ".txt"),
	file_magnetization            (directory + "/magnetization _" + suffix + ".txt"),

	
	file_kinetic_operator  (directory + "/kinetic_operator"  + suffix + ".txt"),
	file_harmonic_trap     (directory + "/harmonic_trap"     + suffix + ".txt"),
	file_linear_zeeman     (directory + "/linear_zeeman"     + suffix + ".txt"),
	file_quadratic_zeeman  (directory + "/quadratic_zeeman"  + suffix + ".txt"),
	file_dipolar_kernel    (directory   + "/dipolar_kernel"  + suffix + ".txt")
	{
	}
	~FilesNames(){};
public:
	const string file_psi_m;        
	const string file_spin_density; 
	const string file_spin_dynamics;
	const string file_spin_angular_momentum;
	const string file_orbital_angular_momentum;
	const string file_energy;        
	const string file_dipolar_field;
	const string file_local_spin_length;
	const string file_mean_spin_length;
	const string file_magnetization;



	const string file_kinetic_operator;
	const string file_harmonic_trap;
	const string file_linear_zeeman;  
	const string file_quadratic_zeeman;
	const string file_dipolar_kernel;          
};

*/


class FilesNames
{
public:
	FilesNames(const string directory_0, const string directory_1 , const string directory_2):
	file_psi_m                    (directory_2 + "/psi_m"         + ".txt"),
	file_spin_density             (directory_2 + "/spin_density"  + ".txt"),
	file_spin_dynamics            (directory_1 + "/spin_dynamics" + ".txt"),
	file_spin_angular_momentum    (directory_1 + "/spin_angular_momentum" + ".txt"),
	file_orbital_angular_momentum (directory_1 + "/orbital_angular_momentum" + ".txt"),
	file_energy                   (directory_1 + "/energy"        + ".txt"),
	file_dipolar_field            (directory_2 + "/dipolar_field" + ".txt"),
	file_local_spin_length        (directory_2 + "/local_spin_length" + ".txt"),
	file_mean_spin_length         (directory_1 + "/mean_spin_length" + ".txt"),
	file_magnetization            (directory_1 + "/magnetization" + ".txt"),
	file_magnetic_field           (directory_1 + "/magnetic_field" + ".txt"),


	


	file_kinetic_operator  (directory_0 + "/kinetic_operator"        + ".txt"),
	file_harmonic_trap     (directory_0 + "/harmonic_trap"           + ".txt"),
	file_linear_zeeman     (directory_0 + "/linear_zeeman"           + ".txt"),
	file_quadratic_zeeman  (directory_0 + "/quadratic_zeeman"        + ".txt"),
	file_dipolar_kernel    (directory_0 + "/dipolar_kernel"          + ".txt")
	{
	}
	FilesNames(const string directory_0, const string directory_1):
	file_psi_m                    (directory_1 + "/psi_m"         + ".txt"),
	file_spin_density             (directory_1 + "/spin_density"  + ".txt"),
	file_spin_dynamics            (directory_1 + "/spin_dynamics" + ".txt"),
	file_spin_angular_momentum    (directory_1 + "/spin_angular_momentum" + ".txt"),
	file_orbital_angular_momentum (directory_1 + "/orbital_angular_momentum" + ".txt"),
	file_energy                   (directory_1 + "/energy"        + ".txt"),
	file_dipolar_field            (directory_1 + "/dipolar_field" + ".txt"),
	file_local_spin_length        (directory_1 + "/local_spin_length" + ".txt"),
	file_mean_spin_length         (directory_1 + "/mean_spin_length" + ".txt"),
	file_magnetization            (directory_1 + "/magnetization" + ".txt"),
	file_magnetic_field           (directory_1 + "/magnetic_field" + ".txt"),

	


	file_kinetic_operator  (directory_0 + "/kinetic_operator"        + ".txt"),
	file_harmonic_trap     (directory_0 + "/harmonic_trap"           + ".txt"),
	file_linear_zeeman     (directory_0 + "/linear_zeeman"           + ".txt"),
	file_quadratic_zeeman  (directory_0 + "/quadratic_zeeman"        + ".txt"),
	file_dipolar_kernel    (directory_0 + "/dipolar_kernel"          + ".txt")
	{
	}
	// FilesNames(const string directory, const string suffix):
	// file_psi_m                    (directory + "/psi_m_"         + suffix + ".txt"),
	// file_spin_density             (directory + "/spin_density_"  + suffix + ".txt"),
	// file_spin_dynamics            (directory + "/spin_dynamics_" + suffix + ".txt"),
	// file_spin_angular_momentum    (directory + "/spin_angular_momentum" + suffix + ".txt"),
	// file_orbital_angular_momentum (directory + "/orbital_angular_momentum" + suffix + ".txt"),
	// file_energy                   (directory + "/energy_"        + suffix + ".txt"),
	// file_dipolar_field            (directory + "/dipolar_field_" + suffix + ".txt"),
	// file_local_spin_length        (directory + "/local_spin_length_" + suffix + ".txt"),
	// file_mean_spin_length         (directory + "/mean_spin_length_" + suffix + ".txt"),
	// file_magnetization            (directory + "/magnetization _" + suffix + ".txt"),

	
	// file_kinetic_operator  (directory + "/kinetic_operator"  + suffix + ".txt"),
	// file_harmonic_trap     (directory + "/harmonic_trap"     + suffix + ".txt"),
	// file_linear_zeeman     (directory + "/linear_zeeman"     + suffix + ".txt"),
	// file_quadratic_zeeman  (directory + "/quadratic_zeeman"  + suffix + ".txt"),
	// file_dipolar_kernel    (directory + "/dipolar_kernel"  + suffix + ".txt")
	// {
	// }
	~FilesNames(){};
public:
	
	const string file_psi_m;        
	const string file_spin_density; 
	const string file_spin_dynamics;
	const string file_spin_angular_momentum;
	const string file_orbital_angular_momentum;
	const string file_energy;        
	const string file_dipolar_field;
	const string file_local_spin_length;
	const string file_mean_spin_length;
	const string file_magnetization;
	const string file_magnetic_field;



	const string file_kinetic_operator;
	const string file_harmonic_trap;
	const string file_linear_zeeman;  
	const string file_quadratic_zeeman;
	const string file_dipolar_kernel;          
};

class PrintParametres
{
public:
	PrintParametres(int print_frequency_spin_dynamics = 0,
		int print_frequency_other_quantity = 0):
		f_print_spin_dynamics(print_frequency_spin_dynamics),
		f_print_other_quantity(print_frequency_other_quantity)
	{
		n_print_spin_dynamics  = 0;
		n_print_other_quantity =0;
	}
	~PrintParametres(){};
public:
	const int f_print_spin_dynamics;
	const int f_print_other_quantity;
	int n_print_spin_dynamics;
	int n_print_other_quantity;
};

int get_row();
int get_col();

#endif
