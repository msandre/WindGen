include ../config.mk

all: hdf5_demo hdf5_demo2D

hdf5_demo: hdf5_demo.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o hdf5_demo hdf5_demo.c $(CFLAGS) $(LDFLAGS)

hdf5_demo2D: hdf5_demo2D.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o hdf5_demo2D hdf5_demo2D.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f hdf5_demo *.h5 *.msh *.res *.vv *.o *~
	rm -f hdf5_demo2D *.h5 *.msh *.res *.vv *.o *~
