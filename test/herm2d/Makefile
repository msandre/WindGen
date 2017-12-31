include ../../config.mk

all: herm2d

herm2d: herm2d.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o herm2d herm2d.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f herm2d *.o *~
