#!/bin/bash

L=(6 7 8 9 10 11 12)

cat Data_All.dat

for i in {0..6}
do
    rm pass.txt
    echo -e "${L[i]} \n 0 " >> pass.txt

    ./Exact.x < pass.txt
     cat outputs.dat >> Data_All.dat
     rm outputs.dat
     
done
    
