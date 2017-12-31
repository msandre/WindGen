include ../../config.mk

all: isotropic

isotropic: isotropic.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o isotropic isotropic.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f isotropic *.res *.png *.o *~
