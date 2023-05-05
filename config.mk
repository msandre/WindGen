WFS_ROOT              := /<path_to_software>/WindGen

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
HDF5_ROOT              := /<path_to_software>/HDF5/HDF5
ZLIB_ROOT              := /<path_to_software>/HDF5/zlib
SZIP_ROOT              := /<path_to_software>/HDF5/szip
CFLAGS                 += -DHAVE_HDF5 -I $(HDF5_ROOT)/include
LDFLAGS                += $(HDF5_ROOT)/lib/libhdf5.a $(ZLIB_ROOT)/lib/libz.a $(SZIP_ROOT)/lib/libsz.a
# ----------------------------------------------------------------------

LDFLAGS                += -pthread -lsz -lz -ldl -lm

