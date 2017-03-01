
# Custom CMake module for finding "GaussianLib" files
# Written by Lukas Hermanns on 24/08/2015

# Find library

find_path(GaussLib_INCLUDE_DIR Gauss/Gauss.h)

# Setup package handle

include(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(
	GaussLib
	DEFAULT_MSG
	GaussLib_INCLUDE_DIR
)
