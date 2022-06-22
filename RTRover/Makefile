#************************************************************
#Nom ......... : Makefile
#Role ........ : compilation and library generation
#Auteur ...... : D. Laurichesse (CNES), A. Privat (CS-SI)
#Version ..... : V1.3 2/15/2016
#Licence ..... : see file LICENSE
# 
#CNES authorize you, free of charge, to circulate and
#distribute for no charge, for purposes other than commercial,
#the source and/or object code of COMPOSITE SOFTWARE on any
#present and future support.
# 
#************************************************************
CC=g++
CFLAGS=-O3 -DRTRover_STATIC_LIB

#****************

# Objects
FICLIB = libRTRover.a
FICOBJ = \
	rtrover.o\
	rtrover_broadcast.o\
	rtrover_date.o\
	rtrover_filter.o\
	rtrover_givens_filter.o\
	rtrover_kalman_filter.o\
	rtrover_map.o\
	rtrover_model.o\
	rtrover_ppp.o\
	rtrover_rover.o\
	rtrover_smoothed_measurement.o\
	rtrover_utility.o

FICH = \
	rtrover.h\
	rtrover_filter_modeled_measurement.h\
	rtrover_kalman_filter.h\
	rtrover_rover.h\
	rtrover_broadcast.h\
	rtrover_frequency.h\
	rtrover_map.h\
	rtrover_smoothed_measurement.h\
	rtrover_date.h\
	rtrover_givens_filter.h\
	rtrover_model.h\
	rtrover_utility.h\
	rtrover_filter.h\
	rtrover_interface.h\
	rtrover_ppp.h\
	rtrover_vector3D.h

#################
# Compilation

all : $(FICLIB)

$(FICLIB) : $(FICOBJ)
	@rm -f $@ 
	@echo "Creating libRTRover.a..."
	@ar -ru $@ $(FICOBJ)

%.o : %.cpp $(FICH)
	@rm -f $*.o
	@$(CC) -c $(CFLAGS) $<
    
clean :
	rm -f *.o $(FICLIB)
