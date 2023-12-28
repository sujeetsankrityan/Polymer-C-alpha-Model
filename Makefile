#FLAGS = -I. -O3 -Wall
FLAGS =  -I. -O3 -g -Wall
nonnat: nonnat.c rmsd1.o ran3n.o defs.h 
	gcc $(FLAGS) -o nonnat nonnat.c rmsd1.o ran3n.o -lm

rmsd1.o:	rmsd1.c
		gcc -O3 -Wall -I. -c rmsd1.c

ran3n.o:	ran3n.c
		gcc -O3 -Wall -I. -c ran3n.c



clean:
	rm -f *.o
	rm nonnat -f
	rm -f movie.pdb
	rm -f _rt
