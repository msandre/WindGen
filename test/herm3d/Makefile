include ../../config.mk

all: herm3d

herm3d: herm3d.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o herm3d herm3d.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f herm3d *.o *~
