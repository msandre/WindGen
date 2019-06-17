include ../config.mk

CDEPS = wfs.c input.c dft.c rv.c herm.c conv.c range.c
OBJS = $(CDEPS:.c=.o)

default: libwfs.a

libwfs.a: $(OBJS) ../config.mk
	rm -f libwfs.a
	ar rcs libwfs.a $(OBJS)

clean:
	@rm -f *.o *.a *~ 2>/dev/null || true

.PHONY: default clean
