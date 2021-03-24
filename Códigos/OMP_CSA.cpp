/*
 *	Simulated Annealing Acoplado em OpenMP
 *	
 *	Codigo 4 - OpenMP utilizando #pragma omp for
 *
 *  COMPILE: g++ -g -Wall -fopenmp -o omp_csa CSA_Problem1.cpp OMP_CSA.cpp
 *  
 *	RUN: ./omp_csa  número de threads  dimensão do problema  function number
 *
 *	Leandro S. Ferreira
 *	Lucas F. Lucena
 *
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h> 
#include <math.h>
#include <pthread.h>
#include <string.h>
#include "CSA_Problem1.h"
#include <iostream>
using namespace std;

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

int main(int argc, char* argv[]) 
{
	int       k                   = 1;                                                  //Iterador
	double    tmp_gen             = 100.0;                                             //Temperatura de geração inicial
	double    tmp_ac              = 100.0;                                            //Temperatura de aceitação inicial
    int       m_threads           = strtol(argv[1], NULL, 10);                       //Número total de otimizadores
   	int       n_size              = strtol(argv[2], NULL, 10);                      //Dimensão inicial
   	int       func_num            = strtol(argv[3], NULL, 10);                     //CSA_EvalCost -> Function Number
   	double    variancia_d         = 0.99 * ((m_threads - 1)/pow(m_threads,2.0));  //Variância desejada
   	double    variancia           = 0.00;                                        //Variância
   	double    custos_correntes [m_threads];                                     //Vetor com os custos de cada otimizador
   	double    max_custo_corrente  = 0.0;                                       //Custo mais alto entre os otimizadores
   	double    gamma               = 0.0;                                      //Termo de acoplamento gamma
   	double    best_sol [m_threads];                                          //Vetor que guarda os melhores custos_ de cada thread
   	double    best                = 0.0;                                    //Melhor solução geral
   	double    media               = 0.0;								   //Média dos custos
   	double    probabilidades_A [m_threads];                               //Vetor que guarda o A de cada otimizador
   	double    somatorio_A         = 0.0;                                 //Guarda o valor do somatoria das probabilidades de aceitação
	struct    timeval start, stop;                                      //Relógio
	int       modular             = 1000000 % m_threads;               //Calcula se 1M é divisivel pelo numero de threads
	int       parada              = (1000000/(m_threads)) + modular;  //Calcula a condição de parada do laço principal
	double    m_threads_d         = static_cast<double>(m_threads);  //Número de threads em double


	gettimeofday(&start, NULL); //Start the clock

	#pragma omp parallel num_threads(m_threads) default(none) shared(m_threads_d, tmp_gen, tmp_ac, m_threads, n_size, func_num, variancia_d, variancia, custos_correntes, max_custo_corrente, gamma, best_sol, probabilidades_A, somatorio_A, stop, modular, parada) private(k)
	{
		int    my_rank = omp_get_thread_num();                           //Ranking de cada otimizador 
		double sol_corrente [n_size];                                   //Solução Corrente
		double sol_nova [n_size];                                      //Solução Nova 
		double A = 0.0;                                               //Probabilidade de aceitação
		double custo_sol_corrente, custo_sol_nova, custo_sol_melhor; //Custos das soluções

		//Semente
		struct drand48_data buffer;
  		srand48_r(my_rank*time(NULL),(&buffer)); //Gera semente

		double num_aleatorio = 0.0;

		//Gera soluções iniciais - Xi
   		for (int i = 0; i < n_size; i++)
   		{
   			drand48_r(&buffer, &num_aleatorio); //Gera um número entre 0 e 1
			sol_corrente[i] = (2.0 * num_aleatorio) - 1.0;
		}

		//Avalia o custo da solução corrente - X
		custo_sol_corrente = CSA_EvalCost(sol_corrente, n_size, func_num);
		custo_sol_melhor   = custo_sol_corrente;

		//Coloca o valor do custo da solução corrente no vetor com todos os custos de cada otimizador
		custos_correntes[my_rank] = custo_sol_corrente;

		#pragma omp barrier
		//Avaliação de gamma
		#pragma omp single
		{	
			//Calcula o custo máximo entre os otimizadores
			for (int i = 0; i < m_threads; ++i)
			{             
				if (i == 0 )
				{
					max_custo_corrente = custos_correntes[i];
				}
				else
				{
					if (custos_correntes[i] > max_custo_corrente)
					{
						max_custo_corrente = custos_correntes[i];
					}
				}
			}

			//Avalia GAMMA
			for (int i = 0; i < m_threads; ++i)
			{
				gamma += exp((custos_correntes[i] - max_custo_corrente)/tmp_ac);
			}
		}
		
		//Calcula a probabilidade de aceitação
		double e; 
		e = exp((custo_sol_corrente - max_custo_corrente)/tmp_ac); 
		A = e/gamma;
		probabilidades_A[my_rank] = A;

		//Iniciar o loop principal - critério de parada
		for (k = m_threads; k < parada;)
		{		
			//Gera as novas soluções - y
   			for (int i = 0; i < n_size; i++)
   			{
   				//srand48_r(my_rank*time(NULL),(&buffer)); //Gera semente
   				drand48_r(&buffer, &num_aleatorio); //Gera um número entre 0 e 1
				sol_nova[i] = fmod(sol_corrente[i] + (tmp_gen*(tan(PI * (num_aleatorio - 0.5)))), 1.0);
			}

			//Avalia o custo da nova solução - y
			custo_sol_nova = CSA_EvalCost(sol_nova, n_size, func_num);
			k++;

			//Joga a moeda
			//srand48_r(my_rank*time(NULL),(&buffer)); //Gera semente
			drand48_r(&buffer, &num_aleatorio); //Gera um número entre 0 e 1
			if (custo_sol_nova < custo_sol_corrente or A > num_aleatorio)
			{
				for (int i = 0; i < n_size; ++i)
				{
					sol_corrente[i] = sol_nova[i]; //Atualiza a solução corrente com o valor da solução nova
				}

				custo_sol_corrente        = custo_sol_nova;       //Atualiza o custo da solução corrente com o custo da solução nova
				custos_correntes[my_rank] = custo_sol_corrente;  //Atualiza o vetor com os custos correntes

				//Verifica se a nova solucao corrente eh a melhor
				if (custo_sol_corrente < custo_sol_melhor)
				{
					custo_sol_melhor = custo_sol_corrente;
				}
			}

			//Atualização
			#pragma omp barrier
			#pragma omp single 
			{
				//Calcula o custo máximo entre os otimizadores
				for (int i = 0; i < m_threads; ++i)
				{
					if (i == 0 )
					{
						max_custo_corrente = custos_correntes[i];
					}
					else
					{
						if (custos_correntes[i] > max_custo_corrente)
						{
							max_custo_corrente = custos_correntes[i];
						}
					}
				}

				//Reseta o GAMMA
				gamma = 0.0;
				
				//Avalia GAMMA
				for (int i = 0; i < m_threads; ++i)
				{
					gamma += exp((custos_correntes[i] - max_custo_corrente)/tmp_ac);
				}
			}

			//Calcula a probabilidade de aceitação
			e = exp((custo_sol_corrente - max_custo_corrente)/tmp_ac); 
			A = e/gamma;
			probabilidades_A[my_rank] = A;
			#pragma omp barrier
			
			#pragma omp single
			{
				//Somatorio das probabilidades de aceitação
				somatorio_A = 0;
				for (int i = 0; i < m_threads; ++i)
				{
					somatorio_A += pow(probabilidades_A[i], 2.0);
				}
				
				//Avalia a variância
				variancia = 0.0;
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
		}

		#pragma omp barrier
		//Stop the clock
		#pragma omp single
		{
			gettimeofday(&stop, NULL);
		}

		//Coloca as melhores soluções em um vetor de melhores soluções
		best_sol[my_rank] = custo_sol_melhor;
	}
	//Fim da região paralela

	//Calcula a média e a melhor solução
	for (int i = 0; i < m_threads; ++i)
	{
		media = media + best_sol[i];

		if (i == 0)
		{
			best = best_sol[i];
		}
		else
			if (best_sol[i] < best)
			{
				best = best_sol[i];
			}		
	}
	media = media/m_threads;

	//Abre o arquivo txt
   	char arquivo[] = "runtime_omp_csa.txt";
    FILE *fp;     
    fp = fopen(arquivo, "a");                 
       
    //Salva o tempo de execução na tabela
    if (fp == NULL) exit(1);
    else 
    {
		fprintf(fp, "%.3f\t", (stop.tv_sec + stop.tv_usec*1e-6)-(start.tv_sec + start.tv_usec*1e-6));
        fclose(fp);     
    }

    
	//Print
	printf("\n%d OTIMIZADORES, DIMENSÃO %d e UTILIZANDO A FUNÇÃO #%d\n", m_threads, n_size, func_num);
	printf("\n MELHOR CUSTO --> %1.2e\n", best);
	printf("\n MÉDIA DOS CUSTOS --> %1.2e\n", media);
	printf("\n");

	return 0; 
}  
/* main */
