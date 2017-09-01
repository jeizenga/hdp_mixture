GITHUBDIR:=..
ROOTDIR:=${GITHUBDIR}/hdp_mixture
IMPLDIR:=${ROOTDIR}/impl
INCDIR:=${ROOTDIR}/inc
BINDIR:=${ROOTDIR}/bin
LIBDIR:=${ROOTDIR}/lib
TESTDIR:=${ROOTDIR}/test
OBJDIR:=${ROOTDIR}/obj
EXTERNDIR:=${ROOTDIR}/external
SONLIBROOTDIR:=${EXTERNDIR}/sonLib
SONLIBDIR:=${SONLIBROOTDIR}/lib
CFLAGS:=-std=c99 -Wall -Wextra -Wpedantic -fopenmp
INC:=-I${INCDIR} -I${SONLIBDIR} -I${EXTERNDIR}
LINK:=-L${SONLIBDIR} -L${OBJDIR} -lm
SHELL:=/bin/bash
CC:=gcc

.PHONY: clean .pre-build

all: .pre-build ${OBJDIR}/CuTest.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/rnglib.o ${OBJDIR}/ranlib.o ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/hdp.o ${OBJDIR}/nanopore_hdp.o ${BINDIR}/main ${BINDIR}/utils_tests ${BINDIR}/nanopore_hdp_tests ${BINDIR}/serialization_tests ${BINDIR}/distribution_comparison_tests ${BINDIR}/base_params_tests ${BINDIR}/parallel_tests ${BINDIR}/distr_script ${BINDIR}/distr_comparison_script

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
	if [ ! -d ${SONLIBROOTDIR} ]; then cd ${EXTERNDIR} && git clone https://github.com/benedictpaten/sonLib.git; fi
	cd ${SONLIBROOTDIR} && make
	
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

${BINDIR}/serialization_tests: ${TESTDIR}/serialization_tests.c ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/hdp.o ${OBJDIR}/CuTest.o ${OBJDIR}/nanopore_hdp.o
	${CC} ${CFLAGS} -c ${TESTDIR}/serialization_tests.c ${INC} -o ${OBJDIR}/serialization_tests.o
	${CC} ${CFLAGS} ${OBJDIR}/serialization_tests.o ${OBJDIR}/hdp.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/ranlib.o ${OBJDIR}/rnglib.o ${OBJDIR}/nanopore_hdp.o ${OBJDIR}/CuTest.o ${INC} ${LINK} -o ${BINDIR}/serialization_tests
	chmod +x ${BINDIR}/serialization_tests
	rm ${OBJDIR}/serialization_tests.o	

${BINDIR}/distribution_comparison_tests: ${TESTDIR}/distribution_comparison_tests.c ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/hdp.o ${OBJDIR}/CuTest.o ${OBJDIR}/nanopore_hdp.o
	${CC} ${CFLAGS} -c ${TESTDIR}/distribution_comparison_tests.c ${INC} -o ${OBJDIR}/distribution_comparison_tests.o
	${CC} ${CFLAGS} ${OBJDIR}/distribution_comparison_tests.o ${OBJDIR}/hdp.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/ranlib.o ${OBJDIR}/rnglib.o ${OBJDIR}/nanopore_hdp.o ${OBJDIR}/CuTest.o ${INC} ${LINK} -o ${BINDIR}/distribution_comparison_tests
	chmod +x ${BINDIR}/distribution_comparison_tests
	rm ${OBJDIR}/distribution_comparison_tests.o	
	
${BINDIR}/base_params_tests: ${TESTDIR}/base_params_tests.c ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/CuTest.o ${SONLIBDIR}/sonLib.a
	${CC} ${CFLAGS} -c ${TESTDIR}/base_params_tests.c ${INC} -o ${OBJDIR}/base_params_tests.o
	${CC} ${CFLAGS} ${OBJDIR}/base_params_tests.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/CuTest.o ${INC} ${LINK} -o ${BINDIR}/base_params_tests
	chmod +x ${BINDIR}/base_params_tests
	rm ${OBJDIR}/base_params_tests.o	

${BINDIR}/parallel_tests: ${TESTDIR}/parallel_tests.c ${OBJDIR}/hdp_math_utils.o ${OBJDIR}/CuTest.o ${SONLIBDIR}/sonLib.a
	${CC} ${CFLAGS} -c ${TESTDIR}/parallel_tests.c ${INC} -o ${OBJDIR}/parallel_tests.o
	${CC} ${CFLAGS} ${OBJDIR}/parallel_tests.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/CuTest.o ${INC} ${LINK} -o ${BINDIR}/parallel_tests
	chmod +x ${BINDIR}/parallel_tests
	rm ${OBJDIR}/parallel_tests.o	

${BINDIR}/distr_script: ${TESTDIR}/distr_script.c ${OBJDIR}/nanopore_hdp.o ${OBJDIR}/hdp_math_utils.o
	${CC} ${CFLAGS} -c ${TESTDIR}/distr_script.c ${INC} -o ${OBJDIR}/distr_script.o
	${CC} ${CFLAGS} ${OBJDIR}/distr_script.o ${OBJDIR}/nanopore_hdp.o ${OBJDIR}/hdp.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/ranlib.o ${OBJDIR}/rnglib.o ${INC} ${LINK} -o ${BINDIR}/distr_script
	chmod +x ${BINDIR}/distr_script
	rm ${OBJDIR}/distr_script.o	

${BINDIR}/distr_comparison_script: ${TESTDIR}/distr_comparison_script.c ${OBJDIR}/nanopore_hdp.o
	${CC} ${CFLAGS} -c ${TESTDIR}/distr_comparison_script.c ${INC} -o ${OBJDIR}/distr_comparison_script.o
	${CC} ${CFLAGS} ${OBJDIR}/distr_comparison_script.o ${OBJDIR}/nanopore_hdp.o ${OBJDIR}/hdp.o ${OBJDIR}/hdp_math_utils.o ${SONLIBDIR}/sonLib.a ${OBJDIR}/ranlib.o ${OBJDIR}/rnglib.o ${INC} ${LINK} -o ${BINDIR}/distr_comparison_script
	chmod +x ${BINDIR}/distr_comparison_script
	rm ${OBJDIR}/distr_comparison_script.o	

clean:
	@if [ $$(find bin -type f | wc -l) -gt 0 ]; \
	then { \
		echo "The following hdp_mixture files will be deleted:"; \
		echo "------------------------------"; \
		find $(BINDIR) $(OBJDIR) -type f; \
		echo "------------------------------"; \
		read -p "Continue (y/n)? " -n 1 -r CONTINUE; \
		echo; \
	}; \
	else echo "No files to delete."; \
	fi; \
	\
	if [[ $$CONTINUE =~ ^[Yy]$$ ]]; \
	then find $(BINDIR) $(OBJDIR) -type f -delete; \
	else echo "Aborted"; \
	fi;
	if [ -d ${SONLIBROOTDIR} ]; then cd ${SONLIBROOTDIR} && make clean; fi

.pre-build:
	if [ ! -d $(BINDIR) ]; then mkdir -p $(BINDIR); fi
	if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR); fi
	if [ `gcc --version | grep '\(4\.9\.[1-9]\d*\)\|\([5-9]\d*\.\d\{1,\}\(\.\d\{1,\}\)\{0,1\}\)' | wc -l` -eq 0 ]; then echo "WARNING: hdp_mixture requires gcc version 4.9.1 or above"; fi