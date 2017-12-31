include ../../config.mk

all: hyp2f1

hyp2f1: hyp2f1.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o hyp2f1 hyp2f1.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f hyp2f1 *.o *~
