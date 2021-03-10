#!/usr/bin/env bash

SCRIPT_NAME=$1

if [ -z ${SCRIPT_NAME} ]; then
  echo "[ERROR]: No script name passed"
  exit 1
fi

OUTPUTNAME=${SCRIPT_NAME/C/exe}

g++ ${SCRIPT_NAME} -o ${OUTPUTNAME} -Wl,-rpath,$(pwd) WrappedProb3++.*.so $(root-config --cflags --libs)

