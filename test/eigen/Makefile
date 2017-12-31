include ../../config.mk

all: eigen

eigen: eigen.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o eigen eigen.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f eigen *.o *~
