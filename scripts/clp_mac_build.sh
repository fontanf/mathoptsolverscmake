#!/usr/bin/env bash

##############################################################
# Build COIN-OR Clp from sources for different architectures #
##############################################################

project_name=Clp
project_version=1.17.10

archs=(arm64 x86_64)
libs=(Clp ClpSolver CoinUtils Osi OsiClp OsiCommonTests)

. scripts/_mac_build.sh

build_coin_project
