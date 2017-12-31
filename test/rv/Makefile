include ../../config.mk

all: rv

rv: rv.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o rv rv.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f rv *.o *~ *.png *.res
