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
	vector<string> funcoes;
} Vetor_qualificacao;

typedef struct svetor_sucessoes
{
	string turno;
	pair<int, int> x;
} Vetor_sucessoes;

typedef struct svetor_contratos
{
	string contrato;
	vector<pair<double, double>> min_max_1;
	vector<pair<int, int>> min_max;
	int fim_semana, fim_semana_restricao;
} Vetor_contratos;

vector<int> retornar_pares(string);
int retornar_enfermeiro(string);

class Cenario
{
public:
	ifstream f;


	string cenario;
	string nomeArq;
	string s;
	string aux;
	int n;
	int numero_enfermeiros;
	int num_semanas;

	Cenario(string nomeArq);
	virtual ~Cenario();

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
};

class Week
{
public:
	ifstream w;

	string nomeArq;
	string s;
	int tam_matriz;

	Week(string nomeArq);
	virtual ~Week();

	vector<vector<pair<int, int>>> matriz_requerimentos;
	vector<vector<int>> matriz_requisicoes;
};

class History
{
public:
	ifstream h;

	string nomeArq;
	string s;
	int n;

	History(string nomeArq);
	virtual ~History();

	vector<vector<string>> matriz_borda_s;
	vector<vector<int>> matriz_borda_fracas;
};

class Custom
{
public:
	ifstream cust;

	string nomeArq;
	string s;
	int n;

	Custom(string nomeArq);
	virtual ~Custom();
	vector<pair<int, int>> dias_trabalhados;
};