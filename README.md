#############################################################################
DataTree structure of data. ROOT (v5 or v6) is needed to run this program.
#############################################################################
INSTALLATION:
cd <DataTree directory>/

Make build directory:
mkdir build

cd build

Link all source files via cmake:
cmake ../src/

Install soft via make:
make

Shared library libDataTree.so is ready to use.

USAGE:

Open ROOT:
root -l

Load shared library:
gSystem->Load("libDataTree.so")
#############################################################################
