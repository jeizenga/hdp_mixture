BIN=/Users/Jordan/Documents/GitHub/hdp_mixture/bin/
SONLIBC=/Users/Jordan/Documents/GitHub/sonLib/C/impl/
SONLIBH=/Users/Jordan/Documents/GitHub/sonLib/C/inc/
CUTEST=/Users/Jordan/Documents/GitHub/sonLib/externalTools/cutest/

all: sonLibList sonLibSet utils hdp main
	
main:
	gcc -c -Wall ${CUTEST}cutest.c ${SONLIBC}sonLibList.c ${SONLIBC}sonLibSet.c  hdp_math_utils.c  hdp.c main.c -I${SONLIBC} -I${SONLIBH} -I${CUTEST}

hdp:
	gcc -c -Wall ${CUTEST}cutest.c hdp.c -I${SONLIBH} -lm -o hdp.o

utils:
	gcc -c -Wall hdp_math_utils.c -lm -o hdp_math_utils.o

sonLibList:
	gcc -c -Wall ${SONLIBC}sonLibList.c -I${SONLIBC} -I${SONLIBH} -I${CUTEST} -o sonLibList.o

sonLibSet:
	gcc -c -Wall ${SONLIBC}sonLibSet.c -I${SONLIBC} -I${SONLIBH} -I${CUTEST} -o sonLibSet.o

clean:
	rm -f main