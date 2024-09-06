#!/usr/bin/env bash

####################################################################
# Build COIN-OR CoinUtils from sources for different architectures #
####################################################################

project_name=CoinUtils
project_version=2.11.12

archs=(arm64 x86_64)
libs=(CoinUtils)

. scripts/_mac_build.sh

build_coin_project
