#!/usr/bin/env bash

##############################################################
# Build COIN-OR Osi from sources for different architectures #
##############################################################

project_name=Osi
project_version=0.108.11

archs=(arm64 x86_64)
libs=(CoinUtils Osi OsiCommonTests)

. scripts/_mac_build.sh

build_coin_project
