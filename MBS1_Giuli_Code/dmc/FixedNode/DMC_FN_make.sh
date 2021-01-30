#!/bin/bash

rm graph.dat
touch graph.dat
rm pass.txt
rm Process.dat

echo "Making DMC_FN simulation..."
#              alfa  t2/t   V/t    L  X_start  Nsteps
echo -e     " 1.68 \n 0.0 \n 5.0 \n 30 \n 2 \n 10000000 \n  " >> dmc_fn_input.dat
./dmc_fn.x < dmc_fn_input.dat
rm dmc_fn_input.dat


echo "Analizing data..."
for (( ip = 20; ip<41; ip=ip+20 )); do
    touch pass.txt
    #        Nstep   BinDim  BinSkip  p
    echo "making p= ${ip}"
    echo -e "1000000  50000    1    ${ip}" >> pass.txt
    (./dmc_fn_read.x < pass.txt) >> graph.dat
    rm pass.txt
done
