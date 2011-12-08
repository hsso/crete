#! /bin/bash
# Executable shell script to configure Ratran for Linux

export OS=`uname`

export RATRAN=`pwd`
export PATH=$PATH:$RATRAN/bin
export RATRANRUN=$RATRAN/run

mkdir -p $RATRANRUN

[ -e amc/Makefile ] && rm amc/Makefile
[ -e sky/Makefile ] && rm sky/Makefile

export fort=gfortran

$fort -g molec/readmol.f -o bin/readmol.lnx
ln -s amc.make.lnx amc/Makefile
ln -s sky.make.lnx sky/Makefile
echo Set up for $OS with compiler $fort
