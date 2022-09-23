#!/bin/bash

rm -rf *.lst

cd ./src
#############################################################################

####################################
echo
echo "Enter into AFTAN ... "
cd AFTAN
make cleanall
make
make install
make clean
cd ../
####################################

####################################
echo
echo "Enter into ANCC ... "
cd ANCC
make cleanall
make
make install
make clean
cd ../
####################################

####################################
echo
echo "Enter into TF_PWS ... "
cd TF_PWS
make cleanall
make
make install
make clean
cd ../
####################################

#############################################################################

####################################
echo
echo "Compiling done ... "
