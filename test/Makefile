include ../config.mk

all: 
	cd rv; $(MAKE)
	cd herm2d; $(MAKE)
	cd herm3d; $(MAKE)
	cd input; $(MAKE)
	cd kaimal; $(MAKE)
	cd simiu; $(MAKE)
	cd isotropic; $(MAKE)
	cd hyp2f1; $(MAKE)
	cd eigen; $(MAKE)

clean:
	rm -f main *.o *~
	cd rv; $(MAKE) clean
	cd herm2d; $(MAKE) clean
	cd herm3d; $(MAKE) clean
	cd input; $(MAKE) clean
	cd kaimal; $(MAKE) clean
	cd simiu; $(MAKE) clean
	cd isotropic; $(MAKE) clean
	cd hyp2f1; $(MAKE) clean
	cd eigen; $(MAKE) clean

