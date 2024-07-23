#!/bin/bash
set -xe

# Install gfortran and cmake
yum install -y gcc-gfortran lapack-devel cmake

# Install NLopt
curl -LO https://github.com/stevengj/nlopt/archive/refs/tags/v2.7.1.tar.gz
tar -xzf v2.7.1.tar.gz
cd nlopt-2.7.1
cmake . && make && make install
cd ..
rm -rf nlopt-2.7.1 v2.7.1.tar.gz
