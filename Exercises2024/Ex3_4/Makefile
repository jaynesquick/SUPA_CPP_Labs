CC=g++ #Name of compiler
FLAGS=-std=c++20 -w #Compiler flags (the s makes it silent)
TARGET=RunFiniteFunctions #Executable name
OBJECTS=RunFiniteFunctions.o FiniteFunctions.o InputFuncs.o
LIBS=-I ../../GNUplot/ -lboost_iostreams

#First target in Makefile is default
${TARGET}:${OBJECTS} #Make target from objects
	@echo "Linking..."
	${CC} ${FLAGS} ${OBJECTS} ${LIBS} -o ${TARGET}
	@make clean

RunFiniteFunctions.o : RunFiniteFunctions.cxx FiniteFunctions.h 
	${CC} ${FLAGS} ${LIBS} -c RunFiniteFunctions.cxx 

FiniteFunctions.o : FiniteFunctions.cxx FiniteFunctions.h
	${CC} ${FLAGS} ${LIBS} -c FiniteFunctions.cxx

InputFuncs.o : InputFuncs.cxx InputFuncs.h
	${CC} ${FLAGS} ${LIBS} -c InputFuncs.cxx

clean: #No targets just run shell command to rm object files and emacs backups
	@rm -f *.o *~

cleantarget: #Delete the exectuable
	@rm -f ${TARGET}
