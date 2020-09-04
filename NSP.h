#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <ctype.h>
#include <stdio.h>

using namespace std;

typedef struct svetor_qualificacao
{
  string nurse; 
  string contract; 
  int num_funcoes;
  vector<string> funcoes;   //bonificaçao, quando a bonificação é zero significa que o ponto é um hotel
} Vetor_qualificacao;

typedef struct svetor_sucessoes
{
  string turno; //coordenada X
  pair<int,int> x; //coordenada Y
} Vetor_sucessoes;

typedef struct svetor_contratos
{
  string contrato; 
  vector<pair<double,double>> min_max_1;
  vector<pair<int,int>> min_max; 
  int fim_semana, fim_semana_restricao;
} Vetor_contratos;

vector<int> retornar_pares(string); // Identifica os dois interios da string (1,52) e retorna um vetor com os dois números
int retornar_enfermeiro(string);


class Cenario {
public:
	//objeto para leitura de arquivo
	ifstream f;
	
	//nome do arquivo
	string cenario;
	string nomeArq;
    string s; 
	string aux;
    int n; // variavel para ler os numeros inteiros do arquivo de cenário
	int numero_enfermeiros;
	int num_semanas;
	

	//construtor e destrutor
	Cenario(string nomeArq);
	virtual ~Cenario();

    // vetor de sucessões de ss proibidas
	//vector<string> ss(5); // vetor com os 5 ss possíveis
	vector<string> vetor_turnos;
	vector<string> vetor_none;                          
    vector<string> vetor_early;
	vector<string> vetor_day;
	vector<string> vetor_late;
	vector<string> vetor_night;

    
	//vetor de qualificações
	vector<Vetor_qualificacao> vetor_qualificacao;
	vector<Vetor_sucessoes> vetor_sucessoes;
	vector<Vetor_contratos> vetor_contratos;

	// vetor de suc
};

class Week { 
public:
	//objeto para leitura de arquivo
	ifstream w;
	
	//nome do arquivo
	string nomeArq;
	string s; // passar s
	int tam_matriz;


	//construtor e destrutor
	Week(string nomeArq);
	virtual ~Week();
  
	//Matriz de requimentos minimo e ótimo
	vector<vector<pair<int,int>>> matriz_requerimentos;
	vector<vector<int>> matriz_requisicoes;
};

class History { 
public:
	//objeto para leitura de arquivo
	ifstream h;
	
	//nome do arquivo
	string nomeArq;
    string s; // ultimo s que foi alocado
	int n;
	 // variavel para ler os numeros inteiros da instancia
	//construtor e destrutor
	History(string nomeArq);
	virtual ~History();
  
	//Matriz de requimentos minimo e ótimo
	vector<vector<string>> matriz_borda_s;
	vector<vector<int>> matriz_borda_fracas;
};


class Custom { 
public:
	//objeto para leitura de arquivo
	ifstream cust;
	
	//nome do arquivo
	string nomeArq;
    string s;
	int n; // ultimo s que foi alocado
	 // variavel para ler os numeros inteiros da instancia
	//construtor e destrutor
	Custom(string nomeArq);
	virtual ~Custom();
	vector<pair<int,int>> dias_trabalhados;
  
};