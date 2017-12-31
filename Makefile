include config.mk

all:
	cd src; $(MAKE)
	cd test; $(MAKE)

fftw:
	cd fftw-3.3.3;                                                        \
	./configure CC=$(CC) MPICC=$(CC) $(FFTW_CONFIG_FLAGS);                \
	 make; make install

clean:
	rm -f *~
	cd src; $(MAKE) clean
	cd test; $(MAKE) clean
