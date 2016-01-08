GITHUB=/Users/Jordan/Documents/GitHub
ROOTDIR=${GITHUB}/hdp_mixture
SONLIBROOT=${GITHUB}/sonLib
SONLIB=${SONLIBROOT}/lib/

test: 
	gcc -Wall -c impl/*.c -Iinc -I${SONLIB} -lm
	ar rcs bin/hdplib.a ./*.o

archive:
	${ar} 

main:
	gcc -c -Wall ${CUTEST}cutest.c ${SONLIBC}sonLibList.c ${SONLIBC}sonLibSet.c  hdp_math_utils.c  hdp.c main.c -I${SONLIBC} -I${SONLIBH} -I${CUTEST}

hdp:
	gcc -c -Wall ${CUTEST}cutest.c hdp.c -I${SONLIBH} -lm -o hdp.o

utils:
	gcc -Wall -c impl/hdp_math_utils.c -Iinc -lm -o bin/hdp_math_utils.o

sonLibList:
	gcc -c -Wall ${SONLIBC}sonLibList.c -I${SONLIBC} -I${SONLIBH} -I${CUTEST} -o sonLibList.o

sonLibSet:
	gcc -c -Wall ${SONLIBC}sonLibSet.c -I${SONLIBC} -I${SONLIBH} -I${CUTEST} -o sonLibSet.o

clean:
	rm ./*.o