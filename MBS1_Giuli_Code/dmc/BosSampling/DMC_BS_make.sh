#!/bin/bash

rm graph.dat
touch graph.dat
rm pass.txt
rm Process.dat

echo "Making DMC_BS simulation..."
#              alfa  t2/t   V/t    L  X_start  Nsteps
echo -e      "1.68 \n 0.0 \n 5.0 \n 30 \n 2 \n 5000000 \n " >> dmc_bs_input.dat

./dmc_bs.x < dmc_bs_input.dat
rm dmc_bs_input.dat

echo "Analizing data..."
for (( ip = 10; ip<501; ip=ip+10 )); do
    touch pass.txt
    #        Nstep   BinDim  BinSkip  p
    echo "making p= ${ip}"
    echo -e "5000000  100000    1    ${ip}" >> pass.txt
    (./dmc_bs_read.x < pass.txt)
    cat p_sim.dat >> graph.dat
    rm pass.txt
    rm p_sim.dat
done
