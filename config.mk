WFS_ROOT              := /home/bkeith/src/ExaQUte/WindGen

# ----------------------------------------------------------------------
#				NO MPI
# ----------------------------------------------------------------------
# CC                    := gcc 
# FFTW_CONFIG_FLAGS     := --prefix=$(WFS_ROOT)/fftw-3.3.3/
# #CFLAGS                := -O2 -ffast-math -freciprocal-math -I $(WFS_ROOT)/src -I $(WFS_ROOT)/fftw-3.3.3/include
# CFLAGS                := -g -fsanitize=undefined -O0 -I $(WFS_ROOT)/src -I $(WFS_ROOT)/fftw-3.3.3/include
# LDFLAGS               := -L $(WFS_ROOT)/src -L $(WFS_ROOT)/fftw-3.3.3/lib
# LDFLAGS               += -lwfs -lfftw3

# ----------------------------------------------------------------------
#				MPI
# ----------------------------------------------------------------------
CC                    := mpicc 
FFTW_CONFIG_FLAGS     := --enable-mpi --prefix=$(WFS_ROOT)/fftw-3.3.3/
CFLAGS                := -DHAVE_MPI -O2 -ffast-math -freciprocal-math -I $(WFS_ROOT)/src -I $(WFS_ROOT)/fftw-3.3.3/include
LDFLAGS               := -L $(WFS_ROOT)/src -L $(WFS_ROOT)/fftw-3.3.3/lib
LDFLAGS               += -lwfs -lfftw3_mpi -lfftw3

# ----------------------------------------------------------------------
#				HDF5
# ----------------------------------------------------------------------
HDF5_ROOT              := /software/user/HDF5/HDF5
ZLIB_ROOT              := /usr/lib64/
SZIP_ROOT              := /software/user/HDF5/szip
CFLAGS                 += -DHAVE_HDF5 -I $/usr/include/openmpi-x86_64/
LDFLAGS                += /usr/lib64/openmpi/lib/libhdf5.so /usr/lib64/libz.so /usr/lib64/libsz.so
# ----------------------------------------------------------------------

LDFLAGS                += -lm
