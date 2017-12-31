include ../../config.mk

all: simiu

simiu: simiu.c $(WFS_ROOT)/src/libwfs.a
	$(CC) -o simiu simiu.c $(CFLAGS) $(LDFLAGS)

clean:
	rm -f simiu *.res *.png *.o *~
