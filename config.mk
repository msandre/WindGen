WFS_ROOT              := <path_to_repository>/WindGen

# ----------------------------------------------------------------------
#                               NO MPI
# ----------------------------------------------------------------------
CC                    := gcc-11
FFTW_CONFIG_FLAGS     := --prefix=$(WFS_ROOT)/fftw-3.3.10/
CFLAGS                := -O2 -ffast-math -freciprocal-math -I $(WFS_ROOT)/src -I $(WFS_ROOT)/fftw-3.3.10/include
LDFLAGS               := -L $(WFS_ROOT)/src -L $(WFS_ROOT)/fftw-3.3.10/lib
LDFLAGS               += -lwfs -lfftw3

# ----------------------------------------------------------------------
#				MPI
# ----------------------------------------------------------------------
#CC                    := mpicc 
#FFTW_CONFIG_FLAGS     := --enable-mpi --prefix=$(WFS_ROOT)/fftw-3.3.10/
#CFLAGS                := -DHAVE_MPI -O2 -ffast-math -freciprocal-math -I $(WFS_ROOT)/src -I $(WFS_ROOT)/fftw-3.3.10/include
#LDFLAGS               := -L $(WFS_ROOT)/src -L $(WFS_ROOT)/fftw-3.3.10/lib
#LDFLAGS               += -lwfs -lfftw3_mpi -lfftw3

# ----------------------------------------------------------------------
#                               HDF5
# ----------------------------------------------------------------------
HDF5_ROOT              := /usr/lib/x86_64-linux-gnu
ZLIB_ROOT              := /usr/lib/x86_64-linux-gnu
SZIP_ROOT              := /usr/lib/x86_64-linux-gnu
CFLAGS                 += -I/usr/include/hdf5/serial -DHAVE_HDF5
LDFLAGS                += $(HDF5_ROOT)/hdf5/serial/libhdf5.a
# ----------------------------------------------------------------------

LDFLAGS                += -pthread -lsz -lz -ldl -lm

