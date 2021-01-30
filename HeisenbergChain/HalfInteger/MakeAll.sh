#!/bin/bash

L=(12 14 16 18 20 22 24 26)
Nlanc=(10 20 30 40 50 50 50 50)

cat Data_All.dat

for i in {0..8}
do
    rm pass.txt
    echo -e "${L[i]} \n 0 \n ${Nlanc[i]}" >> pass.txt

    ./Half.x < pass.txt
     cat outputs.dat >> Data_All.dat
     rm outputs.dat
     
done
    
