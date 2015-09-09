#if using the GNU C++ compiler, uncomment the following:
CXX = g++
#	Note that -fsched-spec-load may be problematic.
#	Try using only -O3 if there are compilation errors.
CXXFLAGS = -O3 -fsched-spec-load 

#	The following may be helpful, but -fschedule-insns
#	has been problematic on recent version of g++.
#	Try including -fschedule-insns for best performance:
#CXXFLAGS = -O3 -fsched-spec-load -fschedule-insns
#	following is helpful for debugging
#CXXFLAGS = -O0 -ggdb

#for SGI and Sun, uncomment the following:
#CXX = CC
#CXXFLAGS = -O3

#if the Intel C compiler is available, use:
#CXX = icc
#CXXFLAGS = -O3

#for xlc compiler, e.g. on AIX, use:
#CXX = xlC
#CXXFLAGS = -O3

.SUFFIXES:	.cpp
.cpp.o:
	${CXX} ${CXXFLAGS} -c -o $@ $<

LINK=${CXX} ${CXXFLAGS} -o ../exe/OligoWalk

OBJ_FILES=algorithm.o rna_library.o structure.o globals.o platform.o \
		  siRNAfilter.o thermo.o random.o stochastic.o intermolecular.o alltrace.o pclass.o oligowalk.o

# rules
#####################
all: 
	@echo 'Type make OligoWalk to make executable file in exe/'

OligoWalk: ${OBJ_FILES}
	${LINK} ${OBJ_FILES}
	@echo "------"
	@echo "Installation complete."
	@echo "------"
	@echo 'A perl script, siRNAWalk.pl, is avaible to run OligoWalk and predict siRNA automatically'
	@echo 'Prerequisite: libsvm. (copy svm-predict and svm-scale to directory svm/)'


# Object Files 
##############################################
oligowalk.o: oligowalk.cpp intermolecular.h rna_library.h algorithm.h structure.h thermo.h pclass.h alltrace.h stochastic.h siRNAfilter.h interface.h platform.h defines.h globals.h stdafx.h 

intermolecular.o: intermolecular.cpp intermolecular.h siRNAfilter.h pclass.h alltrace.h stochastic.h defines.h interface.h platform.h rna_library.h structure.h algorithm.h stdafx.h

siRNAfilter.o: siRNAfilter.cpp siRNAfilter.h algorithm.h rna_library.h structure.h interface.h stdafx.h

pclass.o: pclass.cpp pclass.h rna_library.h algorithm.h structure.h interface.h defines.h

stochastic.o: stochastic.cpp stochastic.h random.h pclass.h stdafx.h platform.h structure.h rna_library.h algorithm.h interface.h defines.h pclass.h stdafx.h

alltrace.o: alltrace.cpp alltrace.h defines.h pclass.h rna_library.h algorithm.h structure.h interface.h platform.h stdafx.h

algorithm.o:  algorithm.cpp algorithm.h defines.h structure.h platform.h stdafx.h interface.h rna_library.h stdafx.h

rna_library.o: rna_library.cpp  rna_library.h platform.h structure.h stdafx.h

structure.o: structure.cpp structure.h platform.h defines.h stdafx.h 

globals.o: globals.cpp globals.h 

platform.o: platform.cpp platform.h 

thermo.o: thermo.cpp thermo.h algorithm.h defines.h interface.h rna_library.h structure.h 

random.o: random.cpp random.h

##############################################

#clean the files
clean:
	rm -f *.o

realclean: clean
	rm -f ../exe/OligoWalk
