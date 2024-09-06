#!/usr/bin/env bash

##############################################################
# Build COIN-OR Cbc from sources for different architectures #
##############################################################

project_name=Cbc
project_version=2.10.12

archs=(arm64 x86_64)
libs=(Cbc Cgl Clp ClpSolver CoinUtils Osi OsiCbc OsiClp OsiCommonTests)

. scripts/_mac_build.sh

build_coin_project
