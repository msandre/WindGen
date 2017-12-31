include ../../config.mk

all: input

input: input.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o input input.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f input *.o *~
