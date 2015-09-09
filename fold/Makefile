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

LINK=${CXX} ${CXXFLAGS} -o $@


COMMON_OBJ_FILES=fold.o algorithm.o arrayclass.o  dotarray.o \
	forceclass.o rna_library.o stackclass.o stackstruct.o \
	structure.o TProgressDialog.o 
MCOMMON_OBJ_FILES=mmfold.o algorithm.o arrayclass.o  dotarray.o \
	forceclass.o rna_library.o stackclass.o stackstruct.o \
	structure.o TProgressDialog.o 


all:
	@echo 'Use `make mfold` to make executable file: mfold'

mfold: ${COMMON_OBJ_FILES} 
	${LINK} ${COMMON_OBJ_FILES} 

mmfold: ${MCOMMON_OBJ_FILES} 
	${LINK} ${MCOMMON_OBJ_FILES} 


# Common Object Files (COMMON_OBJ_FILES)

fold.o: fold.cpp algorithm.h arrayclass.h defines.h  \
		forceclass.h platform.h rna_library.h \
		structure.h TProgressDialog.h

mmfold.o: mmfold.cpp algorithm.h arrayclass.h defines.h  \
		forceclass.h platform.h rna_library.h \
		structure.h TProgressDialog.h


algorithm.o: algorithm.cpp algorithm.h arrayclass.h defines.h	\
	dotarray.h forceclass.h platform.h rna_library.h \
	stackclass.h stackstruct.h structure.h TProgressDialog.h

arrayclass.o: arrayclass.cpp arrayclass.h defines.h

dotarray.o: dotarray.cpp dotarray.h defines.h

forceclass.o: forceclass.cpp forceclass.h

rna_library.o: rna_library.cpp rna_library.h defines.h platform.h structure.h

stackclass.o: stackclass.cpp stackclass.h defines.h

stackstruct.o: stackstruct.cpp stackstruct.h

structure.o: structure.cpp structure.h defines.h platform.h

TProgressDialog.o: TProgressDialog.cpp TProgressDialog.h

#clean files
clean:
	rm -f *.o

realclean: clean
	rm -f *.o; rm -f mfold 
