#!/bin/bash


rm graph.dat
rm pass.txt
touch graph.dat
touch pass.txt



#Starting alpha and beta
echo -e "0.65 \n 0.0 \n" >> pass.txt
for ia in {0..1000}
do
    rm data.dat
    rm pass2.txt
#   Param:  t2/t   V/t   L    Nsteps    FreqSampl
    echo -e "0.0 \n 1.0 \n 30 \n 1000000 \n 10" >> pass.txt
    ./vmc.x < pass.txt
    rm pass.txt
    touch pass2.txt
    touch pass.txt
#           BinLen   BinToSkip   Delta  
    echo -e "1000    \n 2        \n 0.02" >> pass2.txt
    (./vmc_sd_read.x < pass2.txt )
    cat variable.dat >> pass.txt
    rm variable.dat
    
done
