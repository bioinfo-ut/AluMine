#!/bin/bash
# Starts simulations with 2 individuals (1 man and 1 woman)

a="/smallea/tmp/maido/sim"
b="/bigea/test_directory/maido/sim"

if (( $1%2 != 0 )); then
	dir=$a$1
	cd $dir
	pwd
	perl ~/autom4.pl
	bash ~/gtester_for_simulations.sh
else
	dir=$b$1
	cd $dir
	pwd
	perl ~/autom4.pl 
	bash ~/gtester_for_simulations.sh
fi
