GITHUBDIR:=/Users/Jordan/Documents/GitHub
ROOTDIR:=${GITHUBDIR}/hdp_mixture
SONLIBROOTDIR:=${GITHUBDIR}/sonLib
SONLIBDIR:=${SONLIBROOTDIR}/lib
IMPLDIR:=${ROOTDIR}/impl
INCDIR:=${ROOTDIR}/inc
BINDIR:=${ROOTDIR}/bin
LIBDIR:=${ROOTDIR}/lib
TESTDIR:=${ROOTDIR}/test
OBJDIR:=${ROOTDIR}/obj
EXTERNDIR:=${ROOTDIR}/external
CFLAGS:=-std=c99 -Wall -Wextra -Wpedantic
INC:=-I${INCDIR} -I${SONLIBDIR} -I${EXTERNDIR}
LINK:=-L${SONLIBDIR} -L${OBJDIR} -lm
SHELL:=/bin/bash

all: ${OBJDIR}/CuTest.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/rnglib.o ${OBJDIR}/ranlib.o ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/hdp.o ${OBJDIR}/nanopore_hdp.o ${BINDIR}/main ${BINDIR}/utils_tests ${BINDIR}/nanopore_hdp_tests ${BINDIR}/serialization_tests



${OBJDIR}/nanopore_hdp.o: ${OBJDIR}/hdp.o
	${CC} ${CFLAGS} -c ${IMPLDIR}/nanopore_hdp.c ${INC} -o ${OBJDIR}/nanopore_hdp.o

${OBJDIR}/hdp.o: ${SONLIBDIR}/sonLib.a ${OBJDIR}/rnglib.o ${OBJDIR}/ranlib.o ${OBJDIR}/hdp_math_utils.o
	${CC} ${CFLAGS} -c ${IMPLDIR}/hdp.c ${INC} -o ${OBJDIR}/hdp.o
	
${OBJDIR}/hdp_math_utils.o: ${SONLIBDIR}/sonLib.a
	${CC} ${CFLAGS} -c ${IMPLDIR}/hdp_math_utils.c ${INC} -o ${OBJDIR}/hdp_math_utils.o 

${OBJDIR}/ranlib.o: ${OBJDIR}/rnglib.o
	${CC} ${CFLAGS} -c ${EXTERNDIR}/ranlib.c -I${EXTERNDIR} -o ${OBJDIR}/ranlib.o

${OBJDIR}/rnglib.o:
	${CC} ${CFLAGS} -c ${EXTERNDIR}/rnglib.c -I${EXTERNDIR} -o ${OBJDIR}/rnglib.o

${OBJDIR}/CuTest.o:
	${CC} ${CFLAGS} -c ${EXTERNDIR}/CuTest.c -I${EXTERNDIR} -o ${OBJDIR}/CuTest.o

${SONLIBDIR}/sonLib.a:
	cd ${SONLIBROOTDIR}
	make
	cd ${ROOTDIR}


	
${BINDIR}/main: ${TESTDIR}/main.c ${OBJDIR}/hdp.o ${INCDIR}/hdp.h
	${CC} ${CFLAGS} -c ${TESTDIR}/main.c ${INC} -o ${OBJDIR}/main.o
	${CC} ${CFLAGS} ${OBJDIR}/main.o  ${OBJDIR}/hdp.o ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/ranlib.o ${OBJDIR}/rnglib.o ${SONLIBDIR}/sonLib.a ${INC} ${LINK} -lm -o ${BINDIR}/main
	chmod +x ${BINDIR}/main
	rm ${OBJDIR}/main.o

${BINDIR}/utils_tests: ${TESTDIR}/utils_tests.c ${OBJDIR}/hdp_math_utils.o ${INCDIR}/hdp_math_utils.h ${OBJDIR}/CuTest.o
	${CC} ${CFLAGS} -c ${TESTDIR}/utils_tests.c ${INC} -o ${OBJDIR}/utils_tests.o
	${CC} ${CFLAGS} ${OBJDIR}/utils_tests.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${INC} ${LINK} -o ${BINDIR}/utils_tests
	chmod +x ${BINDIR}/utils_tests
	rm ${OBJDIR}/utils_tests.o
	
${BINDIR}/nanopore_hdp_tests: ${TESTDIR}/nanopore_hdp_tests.c ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/CuTest.o
	${CC} ${CFLAGS} -c ${TESTDIR}/nanopore_hdp_tests.c ${INC} -o ${OBJDIR}/nanopore_hdp_tests.o
	${CC} ${CFLAGS} ${OBJDIR}/nanopore_hdp_tests.o ${OBJDIR}/nanopore_hdp.o ${OBJDIR}/hdp.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/ranlib.o ${OBJDIR}/rnglib.o ${OBJDIR}/CuTest.o ${INC} ${LINK} -o ${BINDIR}/nanopore_hdp_tests
	chmod +x ${BINDIR}/nanopore_hdp_tests
	rm ${OBJDIR}/nanopore_hdp_tests.o

${BINDIR}/serialization_tests: ${TESTDIR}/serialization_tests.c ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/hdp.o ${OBJDIR}/CuTest.o
	${CC} ${CFLAGS} -c ${TESTDIR}/serialization_tests.c ${INC} -o ${OBJDIR}/serialization_tests.o
	${CC} ${CFLAGS} ${OBJDIR}/serialization_tests.o ${OBJDIR}/hdp.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/ranlib.o ${OBJDIR}/rnglib.o ${OBJDIR}/CuTest.o ${INC} ${LINK} -o ${BINDIR}/serialization_tests
	chmod +x ${BINDIR}/serialization_tests
	rm ${OBJDIR}/serialization_tests.o	

clean:
	@if [ $$(find bin -type f | wc -l) -gt 0 ]; \
	then { \
		echo "The following will be deleted:"; \
		echo "------------------------------"; \
		find $(BINDIR) $(LIBDIR) $(OBJDIR) -type f; \
		echo "------------------------------"; \
		read -p "Continue (y/n)? " -n 1 -r CONTINUE; \
		echo; \
	}; \
	else echo "No files to delete."; \
	fi; \
	\
	if [[ $$CONTINUE =~ ^[Yy]$$ ]]; \
	then find $(BINDIR) $(LIBDIR) $(OBJDIR) -type f -delete; \
	else echo "Aborted"; \
	fi;