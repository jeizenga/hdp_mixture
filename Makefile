GITHUBDIR:=/Users/Jordan/Documents/GitHub
ROOTDIR:=${GITHUBDIR}/hdp_mixture
SONLIBROOTDIR:=${GITHUBDIR}/sonLib
SONLIBDIR:=${SONLIBROOTDIR}/lib
IMPLDIR:=${ROOTDIR}/impl
INCDIR:=${ROOTDIR}/inc
BINDIR:=${ROOTDIR}/bin
LIBDIR:=${ROOTDIR}/lib
TESTDIR:=${ROOTDIR}/test
CFLAGS:=-Wall -Wextra -Wpedantic
INC:=-I${INCDIR} -I${SONLIBDIR}
SHELL:=/bin/bash

all: ${SONLIBDIR}/sonLib.a ${LIBDIR}/hdplib.a ${BINDIR}/main ${BINDIR}/utils_tests

${SONLIBDIR}/sonLib.a:
	cd ${SONLIBROOTDIR}
	make
	cd ${ROOTDIR}

${LIBDIR}/hdplib.a: ${IMPLDIR}/hdp.c ${IMPLDIR}/hdp_math_utils.c ${IMPLDIR}/ranlib.c ${IMPLDIR}/rnglib.c ${INCDIR}/hdp.h ${INCDIR}/hdp_math_utils.h ${INCDIR}/ranlib.h ${INCDIR}/rnglib.h ${SONLIBDIR}/sonLib.a
	${CC} ${CFLAGS} -c ${IMPLDIR}/*.c ${SONLIBDIR}/sonLib.a ${INC} -lm
	ar rcs ${LIBDIR}/hdplib.a ./*.o
	rm ./*.o

${BINDIR}/main: ${TESTDIR}/main.c ${LIBDIR}/hdplib.a
	${CC} ${CFLAGS} -c ${TESTDIR}/main.c ${INC} -o ${BINDIR}/main
	chmod +x ${BINDIR}/main

${BINDIR}/utils_tests: ${TESTDIR}/main.c ${LIBDIR}/hdplib.a
	${CC} ${CFLAGS} -c ${TESTDIR}/utils_tests.c ${INC} -o ${BINDIR}/utils_tests
	chmod +x ${BINDIR}/utils_tests

clean:
	@if [ $$(find bin -type f | wc -l) -gt 0 ]; \
	then { \
		echo "The following will be deleted:"; \
		echo "------------------------------"; \
		find $(BINDIR) $(LIBDIR) -type f; \
		echo "------------------------------"; \
		read -p "Continue (y/n)? " -n 1 -r CONTINUE; \
		echo; \
	}; \
	else echo "No files to delete."; \
	fi; \
	\
	if [[ $$CONTINUE =~ ^[Yy]$$ ]]; \
	then find $(BINDIR) $(LIBDIR) -type f -delete; \
	else echo "Aborted"; \
	fi;