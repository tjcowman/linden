#!/bin/sh

# Parse the input parameters.
Build=;
Test=;

for arg in $@
do
    case $arg in
    --build)
        Build=1;;
    --test)
        Test=1;;
    esac
done


SCRIPT_DIR="$(dirname "$0")"
OUT_DIR=${SCRIPT_DIR}/out
BUILD_DIR=${OUT_DIR}/build

# Build the project.
if [ -n "$Build" ]
then
    if [ ! -d ${OUT_DIR} ]
    then
        mkdir ${OUT_DIR}
    fi

    if [ ! -d ${BUILD_DIR} ]
    then
        mkdir ${BUILD_DIR}
    fi

    # Configure the build.
    cmake -B ${BUILD_DIR}

    # Build the project
    cmake --build ${BUILD_DIR}

    # Install the project
    cmake --install ${BUILD_DIR} --prefix "${OUT_DIR}/install"
fi

# Run the test cases.
if [ -n "$Test" ]
then
    ${OUT_DIR}/install/test/linden-tests
fi