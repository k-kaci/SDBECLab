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
#include "utility.h"
#include "nvector.h"

#include "mkl.h"

int get_row() 
{
	int r[2];
	getyx(stdscr, r[0], r[1]);
	return r[0];
}
int get_col() 
{
	int c[2];
	getyx(stdscr, c[0], c[1]);
	return c[1];
}



	

