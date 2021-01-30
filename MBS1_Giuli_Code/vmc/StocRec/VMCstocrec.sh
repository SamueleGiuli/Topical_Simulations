#!/bin/bash

rm graph.dat
rm pass.txt
touch graph.dat
touch pass.txt

# Starting alpha and beta
echo -e "1.7 \n 0.0 " >> pass.txt
for ia in {0..500}
do
    echo "step $ia$"
    rm data.dat
    rm pass2.txt
    
#                 t2    V      L      Nsteps  NstepsEachMeas
    echo -e "0.0 \n 5.0 \n 30 \n 1000000 \n 10" >> pass.txt
    ./vmc.x < pass.txt
    rm pass.txt
    touch pass2.txt
    touch pass.txt
#	       BinDim   BinToSkip   Delta
    echo -e "1000    \n 2      \n 0.02" >> pass2.txt
    (./vmc_sr_read.x < pass2.txt )
    cat variable.dat >> pass.txt
    rm variable.dat
	
done
    
	
