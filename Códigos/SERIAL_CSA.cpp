/*
 *	Simulated Annealing Acoplado
 *	
 *  COMPILE: g++ -g -Wall -o serial_csa CSA_Problem1.cpp SERIAL_CSA.cpp
 *  
 *	RUN: ./serial_csa  número de otimizadores dimensão do problema  function number
 *
 *	Leandro S. Ferreira
 *	Lucas F. Lucena
 *
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include "CSA_Problem1.h"
#include <iostream>
#include <time.h>
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif


int main(int argc, char* argv[]) 
{
	int       k                   = 0;                                                  //Iterador
	double    tmp_gen             = 100.0;                                             //Temperatura de geração inicial
	double    tmp_ac              = 100.0;                                            //Temperatura de aceitação inicial
    int       m_threads           = strtol(argv[1], NULL, 10);                       //Número total de otimizadores
   	int       n_size              = strtol(argv[2], NULL, 10);                      //Dimensão inicial
   	int       func_num            = strtol(argv[3], NULL, 10);                     //CSA_EvalCost -> Function Number
   	double    variancia           = 0.00;                                        //Variância
   	double    max_custo_corrente  = 0.0;                                       //Custo mais alto entre os otimizadores
   	double    gamma               = 0.0;                                      //Termo de acoplamento gamma
   	double    somatorio_A         = 0.0;                                    //Guarda o valor do somatoria das probabilidades de aceitação
	struct    timeval start, stop;                                         //Relógio
	double    custo_sol_melhor    = 0.0;
	double    m_threads_d         = static_cast<double>(m_threads);                     //Número de threads em double
	double    variancia_d         = 0.99 * ((m_threads_d - 1)/pow(m_threads_d,2.0));  //Variância desejada

	gettimeofday(&start, NULL); //Start the clock

	//Gera as matrizes de solução, os vetores de custo e das probabilidades de avaliação
	double sol_corrente[m_threads][n_size];		//Solução corrente
    double sol_nova [m_threads][n_size];	   //Solução nova	
    double custo_sol_corrente[m_threads];
    double custo_sol_nova[m_threads];
    double A[m_threads];
    
    //Numero Aleatório
	double num_aleatorio = 0.0;

	//Gera soluções iniciais - Xi
   	for (int i = 0; i < m_threads; ++i)
   	{
   		for (int j = 0; j < n_size; ++j)
   		{
   			srand(i*time(NULL));
   			num_aleatorio      = ((double) rand() / ((double)(RAND_MAX))); //Gera um número entre 0 e 1
			sol_corrente[i][j] = (2.0 * num_aleatorio) - 1.0;
   		}
	}

	//Avalia o custo da solução corrente - X
	for (int i = 0; i < m_threads; ++i)
   	{
   		custo_sol_corrente[i] = CSA_EvalCost(sol_corrente[i], n_size, func_num);
			
		//Avalia o melhor custo e o custo máximo entre os otimizadores
		if (i == 0)
		{
			custo_sol_melhor   = custo_sol_corrente[i];
			max_custo_corrente = custo_sol_corrente[i];
		}
		else
		{
			if (custo_sol_melhor > custo_sol_corrente[i])
			{
				custo_sol_melhor = custo_sol_corrente[i];
			}
			else
			{
				if (custo_sol_corrente[i] > max_custo_corrente)
				{
					max_custo_corrente = custo_sol_corrente[i];
				}
			}
		}
	}

	//Avalia GAMMA
	for (int i = 0; i < m_threads; ++i)
	{
		gamma += exp((custo_sol_corrente[i] - max_custo_corrente)/tmp_ac);
	}

	//Calcula a probabilidade de aceitação
	for (int i = 0; i < m_threads; ++i)
	{
		double e = exp((custo_sol_corrente[i] - max_custo_corrente)/tmp_ac); 
		A[i] = e/gamma;
	}

	//Iniciar o loop principal - critério de parada
	for (k = m_threads; k < 1000000;)
	{			
		//Gera as novas soluções - y
   		for (int i = 0; i < m_threads; i++)
   		{
   			for (int j = 0; j < n_size; ++j)
   			{
   				num_aleatorio  = ((double) rand() / ((double)(RAND_MAX))); //Gera um número entre 0 e 1
				sol_nova[i][j] = fmod(sol_corrente[i][j] + (tmp_gen*(tan(PI * (num_aleatorio - 0.5)))), 1.0);
   			}
		}
		
		//Avalia o custo da solução nova - y
		for (int i = 0; i < m_threads; ++i)
   		{
   			custo_sol_nova[i] = CSA_EvalCost(sol_nova[i], n_size, func_num);
   			k++;
   		}
			
		//Joga a moeda
		for (int i = 0; i < m_threads; ++i)
   		{
   			num_aleatorio  = ((double) rand() / ((double)(RAND_MAX)));
   			
   			if (custo_sol_corrente[i] > custo_sol_nova[i] or A[i] > num_aleatorio)
   			{
   				for (int j = 0; j < n_size; ++j)
   				{
  					sol_corrente[i][j] = sol_nova[i][j]; //Atualiza a solução corrente com o valor da solução nova	
   				}

   				custo_sol_corrente[i] = custo_sol_nova[i]; //Atualiza o custo da solução corrente com o custo da solução nova
   				
   				//Verifica se a nova solucao corrente eh a melhor
   				if (custo_sol_corrente[i] < custo_sol_melhor)
				{
					custo_sol_melhor = custo_sol_corrente[i];
				}
				
				//Avalia o melhor custo e o custo máximo entre os otimizadores
				if (i == 0)
				{
					max_custo_corrente = custo_sol_corrente[i];
				}
				else
				{
					if (custo_sol_corrente[i] > max_custo_corrente)
					{
						max_custo_corrente = custo_sol_corrente[i];
					}
				}
   			}
   		}
				
		//Reavalia GAMMA
		gamma = 0.00;
		for (int i = 0; i < m_threads; ++i)
		{
			gamma += exp((custo_sol_corrente[i] - max_custo_corrente)/tmp_ac);
		}

		//Recalcula a probabilidade de aceitação
		for (int i = 0; i < m_threads; ++i)
		{
			double e = exp((custo_sol_corrente[i] - max_custo_corrente)/tmp_ac); 
			A[i] = e/gamma;
		}

		//Somatório de A
		somatorio_A = 0.0;
		for (int i = 0; i < m_threads; ++i)
		{
			somatorio_A += pow(A[i], 2.0);
		}

		variancia = ((1/(m_threads_d))*(somatorio_A)) - (1/(pow(m_threads_d, 2.0)));

		//Atualiza as temperaturas de aceitação e geração
		if (variancia < variancia_d)
		{
			tmp_ac = 0.95 * tmp_ac;
		}
		else
		{
			if (variancia >= variancia_d)
			{
				tmp_ac = 1.05 * tmp_ac;
			}
		}
		tmp_gen = 0.99992 * tmp_gen;			
	}

	gettimeofday(&stop, NULL);

	//Abre o arquivo txt
   	char arquivo[] = "serial_csa.txt";
    FILE *fp;     
    fp = fopen(arquivo, "a");                 
       
    //Salva o tempo de execução na tabela
    if (fp == NULL) exit(1);
    else 
    {
		fprintf(fp, "%.3f\t", (stop.tv_sec + stop.tv_usec*1e-6)-(start.tv_sec + start.tv_usec*1e-6));
        fclose(fp);     
    }

	//Imprime qual é a melhor solução e o código da função de teste utilizada
	printf("\n%d OTIMIZADORES, DIMENSÃO %d e UTILIZANDO A FUNÇÃO #%d\n", m_threads, n_size, func_num);
	printf("\n MELHOR CUSTO --> %1.2e\n", custo_sol_melhor);

	return 0; 
}  
/* main */
