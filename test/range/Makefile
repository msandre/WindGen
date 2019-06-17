include ../../config.mk

all: range ../../config.mk

range: range.c $(WFS_ROOT)/src/libwfs.a
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	@rm -f range *.o *~ 2>/dev/null || true

.PHONY: all clean
