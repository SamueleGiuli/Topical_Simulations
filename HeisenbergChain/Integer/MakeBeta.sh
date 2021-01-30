#!/bin/bash

L=(8 9 10 11 12 13 14)
Nlanc=(30 35 40 50 50 50 50)
Beta=( 1.0 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5)

cat Data_All.dat

for i in {0..6}
do
    for b in {0..0}
    do
	
	rm pass.txt
	echo -e "${L[i]} \n 0 \n ${Nlanc[i]} \n ${Beta[b]} \n " >> pass.txt
	
	./Beta.x < pass.txt
	cat outputs.dat >> Data_All_b.dat
	rm outputs.dat
	
	
    done
done
mv Data_All_b.dat Data_All_b1.dat
    
