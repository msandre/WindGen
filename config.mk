WFS_ROOT              := /software/user/windgen

# ----------------------------------------------------------------------
#				NO MPI
# ----------------------------------------------------------------------
CC                    := gcc 
FFTW_CONFIG_FLAGS     := --prefix=$(WFS_ROOT)/fftw-3.3.3/
CFLAGS                := -O2 -ffast-math -freciprocal-math -I $(WFS_ROOT)/src -I $(WFS_ROOT)/fftw-3.3.3/include
LDFLAGS               := -L $(WFS_ROOT)/src -L $(WFS_ROOT)/fftw-3.3.3/lib
LDFLAGS               += -lwfs -lfftw3

# ----------------------------------------------------------------------
#				MPI
# ----------------------------------------------------------------------
#CC                    := mpicc 
#FFTW_CONFIG_FLAGS     := --enable-mpi --prefix=$(WFS_ROOT)/fftw-3.3.3/
#CFLAGS                := -DHAVE_MPI -O0 -I $(WFS_ROOT)/src -I $(WFS_ROOT)/fftw-3.3.3/include
#LDFLAGS               := -L $(WFS_ROOT)/src -L $(WFS_ROOT)/fftw-3.3.3/lib
#LDFLAGS               += -lwfs -lfftw3_mpi -lfftw3

# ----------------------------------------------------------------------
#				HDF5
# ----------------------------------------------------------------------
#HDF5_ROOT              := /software/user/HDF5/HDF5
#ZLIB_ROOT              := /software/user/HDF5/zlib
#SZIP_ROOT              := /software/user/HDF5/szip
#CFLAGS                 += -DHAVE_HDF5 -I $(HDF5_ROOT)/include
#LDFLAGS                += $(HDF5_ROOT)/lib/libhdf5.a $(ZLIB_ROOT)/lib/libz.a $(SZIP_ROOT)/lib/libsz.a
# ----------------------------------------------------------------------

LDFLAGS                += -lm
