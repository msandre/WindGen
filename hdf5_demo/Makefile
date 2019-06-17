include ../config.mk

default: hdf5_demo

hdf5_demo: hdf5_demo.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o hdf5_demo hdf5_demo.c $(CFLAGS) $(LDFLAGS)

clean:
	@rm -f hdf5_demo *.h5 *.msh *.res *.vv *.o *~ 2>/dev/null || true

.PHONY: default clean
