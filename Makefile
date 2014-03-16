all:phirun

phirun: phi_main1.o Mersenne_Twister.o randGraph.o wolf.o
	gcc -o dmd.linux -lm phi_main1.o Mersenne_Twister.o randGraph.o wolf.o

phi_main1.o: phi_main1.c phi_lib.h
	gcc -c -g phi_main1.c

Mersenne_Twister.o: Mersenne_Twister.c phi_lib.h 
	gcc -c -g Mersenne_Twister.c

randGraph.o: randGraph.c phi_lib.h
	gcc -c -g randGraph.c

wolf.o: wolf.c phi_lib.h
	gcc -c -g wolf.c

clean:
	rm -rf *o phirun
