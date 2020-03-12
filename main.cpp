#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include "Timer.h"
//#include<time.h>

#include "BRKGA.h"
#include "MTRand.h"

string retornar_sem_aspas(string palavra){ 
	string chars = "\"";
	for (char c : chars)
	{
		palavra.erase(remove(palavra.begin(), palavra.end(), c), palavra.end()); //remove todas as ocorrencias dos caracters selecionados
	}
	return palavra;
}

using namespace std;

int main(int argc, char*argv[])
{
    string aux = ".txt";
    int seed;

    vector<string> names(argc);
    cout << "Tamanho comando" << ' ' << argc << endl;

    for(int i = 0; i<argc ; i++){
      cout << argv[i] << ' ';
    }
    cout << endl;

    string cusin;
    for(int i = 0; i<argc ; i++){
       
        if(argv[i] == string("--sce")){
            names[0] = retornar_sem_aspas(argv[i+1]);
            
        }
          
        if(argv[i] == string("--his")){
            names[1] = retornar_sem_aspas(argv[i+1]);
            
        }
            
        if(argv[i] == string("--week")){
            names[2] = retornar_sem_aspas(argv[i+1]);
            
        }
            
        if(argv[i] == string("--cusIn")){
            cusin = string(argv[i]);
            names[5] = retornar_sem_aspas(argv[i+1]);
            names[5] = names[5] + aux;
            cout << "entrou string = cunIn" << endl;
        }
            
        if(argv[i] == string("--cusOut")){
            names[3] = retornar_sem_aspas(argv[i+1]);
            names[3] = names[3] + aux;  
           
        }
            
        if(argv[i] == string("--sol")){
            names[4] = retornar_sem_aspas(argv[i+1]); 
           
        }

        if(argv[i] == string("--rand")){
            names[6] = retornar_sem_aspas(argv[i+1]); 
        } 
    }
    seed = stoi(names[6]);
    
    MTRand mt(seed);
    
    BRKGA bkrga(35,10000,0.2,0.2,0.5,names,cusin,1,mt); // tamanho da populacao, numero de geracoes, porcentagem da populacao elitista, porcentagem de mutantes, probabilidade de herdar alelo do elitista
    
    return 0;
}


