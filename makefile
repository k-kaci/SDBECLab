# fichier makefile pour  

#choix du compilateur 
CC = icc

#option de compilation
CFLAGS= -g -W -Wall -std=c++11 -Wextra -pedantic  

#option de bibliotheques
MKLROOT = /softs/intel/compilers_and_libraries_2016.1.150/linux/mkl
MKL_LIBS  = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core 
MKL_LIBS += -lmkl_intel_thread -qopenmp -lpthread -lm -ldl 
MKL_LIBS_OPT = -qopenmp -I${MKLROOT}/include 


GSLROOT      = /usr/local/compiled/gsl-2.1
GSL_LIBS     = -L${GSLROOT}/lib/ -lgsl -lgslcblas -lm
GSL_LIBS_OPT = -I${GSLROOT}/include
LIBS=   -lm -ldl -lcurses 


# commande de netoyage! 
RM=/bin/rm -f

# sources-------------Editer les 2 lignes ci-dessous-------
MAIN= main.cpp            # <<  ici le fichier de programme principal
SUB=  nvector.cpp utility.cpp fft.cpp short_interaction.cpp  spin.cpp sdbec.cpp # << indiquer ici les fichiers de procedures
#----------------------------------------------------------

#nom des header des procedure
HEAD = $(SUB:.cpp=.h)

#nom des fichier objets 
OBJ_MAIN = $(MAIN:.cpp=.o)
OBJ_SUB = $(SUB:.cpp=.o)

#nom de l'exucutable à générer
EXEC = $(MAIN:.cpp=)
#EXEC = a.out_ok

$(EXEC): $(OBJ_MAIN) $(OBJ_SUB)
	$(CC) $(OBJ_MAIN) $(OBJ_SUB)  $(MKL_LIBS)  $(GSL_LIBS) $(LIBS)   -o $@
# $@ :: nom de la cible

$(OBJ_MAIN): $(MAIN) $(HEAD) types.h mathematical_constants.h physical_constants.h  utility.h 
	$(CC) $(CFLAGS)  $(MKL_LIBS_OPT)  $(GSL_LIBS_OPT)  -c  $< -o $@
# $< :: nom de la premiere dépendance

# Construction des fichier objet à partir des fichier .c
%.o: %.cpp %.h types.h mathematical_constants.h physical_constants.h utility.h 
	$(CC) $(CFLAGS)  $(MKL_LIBS_OPT) $(GSL_LIBS_OPT) -c  $< -o $@
# %.o: %.c :: construire un .o à partir d'un .c

# Notoyage complet
clean: 
	$(RM) *.out *.o $(EXEC)







