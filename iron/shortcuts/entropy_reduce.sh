#!/bin/bash


if ! type "Rscript" > /dev/null; then
  # install foobar here
  echo "Rscript is required but it's not installed.  Please install it and make sure its in your path."
  exit 1;
fi

pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd -P`
popd > /dev/null
Rscript $SCRIPTPATH/../utilities/entropy_reduce.R "$@"
