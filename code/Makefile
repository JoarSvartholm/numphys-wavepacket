#
# Makefile for program "wavepacket"
#
# C.M. Dion, A. Hashemloo, and G. Rahali
# April 2013
#

CC = gcc
LOADLIBES =  -L/usr/local/lib -lgsl -lgslcblas -lm
CFLAGS = -g -I/usr/local/include

CFLAGS = -O3 -std=c99

LDLIBS = -L/usr/local/lib -lfftw3 -lm

freepart:	wavepacket.o free_particle.o
	$(CC) -o $@ $^ $(LDLIBS)

potwell:	wavepacket.o pot_well.o
	$(CC) -o $@ $^ $(LDLIBS)

barrier:	wavepacket.o pot_barrier.o
	$(CC) -o $@ $^ $(LDLIBS)

diffScheme:	wavediff.o diff_pot_barrier.o
	$(CC) -o $@ $^ 

%.o:%.c
	$(CC) -c -o $@ $<


.PHONY: clean
clean:
	rm -f *.o iontrap3D co
