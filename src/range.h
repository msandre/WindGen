#include <stdio.h>
#include <stddef.h>

typedef struct _range range_t;

struct _range
{
  ptrdiff_t start;
  ptrdiff_t size;
};

range_t range_create(ptrdiff_t start, ptrdiff_t size);

range_t range_intersect(range_t r1, range_t r2);

range_t range_shift(range_t r, ptrdiff_t offset);

void range_copy(double * src, range_t src_range, double * dest, range_t dest_range, range_t local_range);
