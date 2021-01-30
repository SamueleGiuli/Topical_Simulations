#!/bin/bash

L=(6 7 8 9 10 11 12 13 14)
Nlanc=(10 20 30 40 40 50 50 50 50)

cat Data_All.dat

for i in {0..8}
do
    rm pass.txt
    echo -e "${L[i]} \n 0 \n ${Nlanc[i]}" >> pass.txt

    ./prova.x < pass.txt
     cat outputs.dat >> Data_All.dat
     rm outputs.dat
     
done
    
