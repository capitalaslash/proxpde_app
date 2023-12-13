#! /bin/bash

cmake \
  --no-warn-unused-cli \
  -S . \
  -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE:STRING=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE \
  -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/gcc-10 \
  -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/g++-10 \
  -DProXPDE_DIR:PATH=~/software/proxpde/lib/cmake/ProXPDE \
  $*

