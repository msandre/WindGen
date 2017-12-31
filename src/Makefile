include ../config.mk

CDEPS = wfs.c input.c dft.c rv.c herm.c conv.c
OBJS = $(CDEPS:.c=.o)

all: libwfs.a

libwfs.a: $(OBJS)
	rm -f libwfs.a
	ar rcs libwfs.a $(OBJS)

clean:
	rm -f *.o *.a *~
