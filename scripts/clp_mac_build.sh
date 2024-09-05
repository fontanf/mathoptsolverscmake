#!/usr/bin/env bash

##############################################################
# Build COIN-OR Clp from sources for different architectures #
##############################################################

clp_version=1.17.10
archs=(arm64 x86_64)
libs=(Clp ClpSolver CoinUtils Osi OsiClp OsiCommonTests)
build_dir=build

# Create the build dir
mkdir -p ${build_dir}

# Get coinbrew
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew -P ${build_dir} -O ${build_dir}/coinbrew
chmod u+x ${build_dir}/coinbrew

# Fetch Clp
(cd ${build_dir} && ./coinbrew fetch Clp@${clp_version} --no-third-party)

# Build Clp for each architecture
for arch in ${archs[*]}
do
  (cd ${build_dir} && env /usr/bin/arch -${arch} ./coinbrew build Clp --no-prompt --reconfigure --build-dir=build/${arch} --prefix=dist/${arch} --disable-shared --disable-zlib --disable-bzlib --disable-readline --without-cholmod --without-amd --without-sample --without-lapack --no-third-party LT_LDFLAGS=-all-static --tests=none)
done

# Create universal lib dir
mkdir -p ${build_dir}/dist/universal/lib

# Copy includes
cp -r ${build_dir}/dist/${archs[0]}/include ${build_dir}/dist/universal/

# Combine each lib in single fat binary lib
for lib in ${libs[*]}
do
  files=""; for arch in ${archs[*]}; do files="${files} ${build_dir}/dist/${arch}/lib/lib${lib}.a"; done
  lipo ${files} -create -output ${build_dir}/dist/universal/lib/lib${lib}.a
done

# Create an archive with all static libs
tar -czf ${build_dir}/clp.tar.gz -C ${build_dir}/dist .
