#!/bin/sh
# Template to set environment variables for the
# non-regression testing framework.

export VORPALINE_BUILD_CONFIG="Release"
export VORPALINE_SOURCE_DIR="/Users/lihongbo/Desktop/code/LpCVT_recon"
export VORPALINE_BUILD_DIR="/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release"
export VORPALINE_BIN_DIR="/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/bin"
export VORPALINE_LIB_DIR="/Users/lihongbo/Desktop/code/LpCVT_recon/cmake-build-release/lib"
export VORPATEST_ROOT_DIR="/Users/lihongbo/Desktop/code/LpCVT_recon/tests"
export DATADIR="/Users/lihongbo/Desktop/code/LpCVT_recon/tests/data"

args=
while [ -n "$1" ]; do
    case "$1" in
        --with-*=*)
            var=`echo "$1" | sed 's/--with-\([^=]*\)=\(.*\)$/VORPALINE_WITH_\U\1\E=\2/'`
            export "$var"
            shift
            ;;
        --with-*)
            var=`echo "$1" | sed 's/--with-\(.*\)$/VORPALINE_WITH_\U\1=1/'`
            export "$var"
            shift
            ;;
        *)
            args="$args $1"
            shift;
            ;;
    esac
done

