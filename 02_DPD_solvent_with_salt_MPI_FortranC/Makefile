CC = mpicc
FC = mpif77
FFLAGS = -O3 -fallow-argument-mismatch -Wall -o

INCLUDE1 = home/server/fftw/nfft2
INCLUDE2 = home/server/fftw/fftw-3.2alhpa2
INCLUDE3 = -lnfft -lfftw3
INCLUDE4 = -L/home/server/fftw/fftw-3.2alhpa2/lib

solvent_with_salts:
	enuf_nfft.o $(FC) $(FFLAGS) \
	solvent_with_salts solvent_with_salts.f \
	enuf_nfft.o $(INCLUDE3) \
	-L/$(INCLUDE1)/lib -L/$(INCLUDE2)/lib $(INCLUDE4)

enuf_nfft.o: 
	$(CC) $(FFLAGS) enuf_nfft.o -c enuf_nfft.c \
	-I/$(INCLUDE1)/include -I/$(INCLUDE2)/include

.PHONY: clean
clean:
	rm *.o solvent_with_salts
