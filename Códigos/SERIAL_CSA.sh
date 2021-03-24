#!/bin/bash

#SBATCH --partition=cluster
#SBATCH --job-name=serial_csa
#SBATCH --output=serial_csa.out
#SBATCH --error=serial_csa.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=0-15:00
#SBATCH --hint=compute_bound
#SBATCH --exclusive

#Carrega os módulos do sistema
rm serial_csa.txt
module load compilers/gnu/7.3
module load compilers/intel/2018.2.199
eval $loadIntelLibs

#Compila o código
g++ -g -Wall -o serial_csa CSA_Problem1.cpp SERIAL_CSA.cpp

tentativas=5 #Quantas vezes o código será executado

for function in 2001 2003 2006 #função utilizada
do
	for cores in 2 4 8 16 32 #número de cores utilizados
	do
		for size in 5 10 50 100 #tamanho do problema			
		do   									
			echo -e "\n$function\t$cores\t$size\t\t\c" >> "serial_csa.txt" 
				
			for tentativa in $(seq $tentativas) #Cria uma vetor de 1 a "tentativas"				
			do
				echo -e `./serial_csa $cores $size $function`						
			done			
		done
	done	
done

exit 
