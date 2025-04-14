#!/bin/bash

# A small command line tool to comment out all passed requirement lines in a requirements file

set -e

if [ ! "$1" == "" ]; then
  REQUIREMENTS_FILE=$1
else
  # none given -> auto-detect
  if [ -f "requirements.txt" ]; then
    REQUIREMENTS_FILE="requirements.txt"
  elif [ -f "requirements/requirements.txt" ]; then
    REQUIREMENTS_FILE="requirements/requirements.txt"
  else
    echo "No requirements file found"
    exit 1
  fi
fi

LOOSE_REQUIREMENTS_FILE=${REQUIREMENTS_FILE/.txt/_loose.txt}

echo using $REQUIREMENTS_FILE $LOOSE_REQUIREMENTS_FILE

for a in "$@"; do
  sed -i "s/$a/### $a/" $REQUIREMENTS_FILE

  if [ -f $LOOSE_REQUIREMENTS_FILE ]; then
      sed -i "s/$a/### $a/" $LOOSE_REQUIREMENTS_FILE
  fi
done

echo $REQUIREMENTS_FILE:
cat  $REQUIREMENTS_FILE
