include ../../config.mk

all: kaimal

kaimal: kaimal.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o kaimal kaimal.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f kaimal *.res *.png *.o *~
