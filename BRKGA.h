#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <vector>
#include <ctime>
#include <stdlib.h>
#include "Timer.h"

#include "NSP.h"
#include "MTRand.h"

using namespace std;

typedef struct scompatibilidade
{
   vector<int> compatibilidade;
} Vetor_compatibilidade;


class BRKGA
{
public:
vector<vector<pair<double, double>>> VND_print(vector<vector<pair<double, double>>> & individuo_codificado);
void troca_populacao();
int multi_start(int &iteracoes_max, double &tempo);
vector<vector<pair<string,string>>> individuo_vnd;
Timer temp;
//vector<vector<pair<double,double>>> pai_elite;
//vector<vector<pair<double,double>>> pai_nao_elite;

vector<vector<pair<double, double>>> mutacao_2_cruzamento(vector<vector<pair<double, double>>> &pai_elite, vector<vector<pair<double, double>>> &pai_nao_elite);
vector<vector<pair<double,double>>> cruzamento_ponto_corte(vector<vector<pair<double, double>>> &pai_elite, vector<vector<pair<double, double>>> &pai_nao_elite);
bool primeira_melhora_2_aleatorio_vnd(vector<vector<pair<string, string>>> &individuo_decodificado);
vector<vector<pair<string,string>>> troca_qualificao_turno_aleatorio(vector<vector<pair<string,string>>> &individuo_decodificado);
vector<vector<pair<double,double>>> VND(vector<vector<pair<double,double>>> &individuo_codificado);
vector<vector<pair<string, string>>> decodificar_individuo_ponderado_2_sol(vector<vector<pair<double, double>>> &individuo_codificado, int &sol_alocacoes);
vector<vector<pair<double, double>>> codificador_ponderado_2(vector<vector<pair<string, string>>> &individuo_decodificado);
vector<vector<pair<string, string>>> decodificar_individuo_ponderado_2(vector<vector<pair<double, double>>> &individuo_codificado);
bool primeira_melhora_2_qualificacao_aleatorio (vector<vector<pair<string, string>>> &individuo_decodificado);
bool primeira_melhora_2_turno_aleatorio (vector<vector<pair<string, string>>> &individuo_decodificado);
vector<vector<pair<string, string>>> primeira_melhora_2_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado);
vector<vector<pair<string, string>>> primeira_melhora_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado);
bool primeira_melhora_2_folga_aleatorio (vector<vector<pair<string, string>>> &individuo_decodificado);
bool primeira_melhora_2_turno_blocos_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado);
bool primeira_melhora_2_troca_pares(vector<vector<pair<string, string>>> &individuo_decodificado);
//movimentos VND
vector<vector<pair<string, string>>> folga_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado);
vector<vector<pair<string, string>>> turno_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado);
vector<vector<pair<string, string>>> qualificacao_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado);
void solucao_gulosa();
void brincadeira();
vector<vector<pair<string, string>>> primeira_melhora_2(vector<vector<pair<string, string>>> &individuo_decodificado);
vector<vector<pair<string, string>>> decodificar_individuo_2(vector<vector<pair<double, double>>> &individuo_codificado);
vector<vector<pair<string,string>>> trocas(int i, int j, vector<vector<pair<string,string>>> &individuo_decodificado);
vector<vector<pair<string,string>>> trocas_2(int i, int j, vector<vector<pair<string,string>>> &individuo_decodificado);
int funcao_avaliacao_2(vector<vector<pair<string, string>>> &individuo_decodificado);
void solucao_gulosa_2();
  string nomeArq;
  string custom_out;
  int semana_1;
  int semana;
  int tam_populacao;
  int num_geracoes;
  int tam_elitista;
  int tam_mutante;
  int tam_nao_elitista;
  int tam_cruzamento;
  int iteracoes_max;
  int count_melhora = 0;
  double tempo;
  double prob_cruzamento;
  Cenario *cn;
  History *ht;
  Week *w0;
  Custom *cus;
  vector<int> dias_trabalhados;
  vector<vector<pair<string, string>>> decodificar_individuo_ponderado(vector<vector<pair<double, double>>> &individuo_codificado);
  vector<vector<pair<string, string>>> decodificar_individuo_2_sol(vector<vector<pair<double, double>>> &individuo_codificado, int &sol_alocacoes);
  void dividir_populacao_2(vector<vector<vector<pair<double, double>>>> &populacao);
  void porcentagem_contratos();
  vector<double> porcentagens;
  vector<string> sem_aspas;
  BRKGA(int,int,double,double,double,vector<string>,string,int,MTRand,double);
  string retornar_sem_aspas(string);
  vector<vector<pair<double,double>>> individuo_codificado;
  vector<vector<pair<string,string>>> individuo_decodificado;
  vector<vector<vector<pair<double,double>>>> populacao;
  vector<vector<vector<pair<double,double>>>> proxima_populacao;
  vector<vector<pair<double,double>>> criar_individuo();
  vector<vector<pair<string,string>>> decodificar_individuo_sol(vector<vector<pair<double,double>>> &individuo_codificado, int &sol_alocacoes);
  vector<vector<pair<string,string>>> decodificar_individuo(vector<vector<pair<double,double>>> &individuo_codificado);
  vector<vector<pair<string, string>>> decodificar_individuo_ponderado_sol(vector<vector<pair<double, double>>> &individuo_codificado, int &sol_alocacoes);
  int funcao_avaliacao(vector<vector<pair<string,string>>> &individuo_decodificado);
  int funcao_avaliacao_teste(vector<vector<pair<string,string>>> &individuo_decodificado);
  int funcao_avaliacao_minimo(vector<vector<pair<string, string>>> &individuo_decodificado, int);
  int funcao_avaliacao_minimo_completa(vector<vector<pair<string, string>>> &individuo_decodificado);
  int funcao_avaliacao_linha_coluna(vector<vector<pair<string, string>>> &individuo_decodificado, int linha, int coluna);
  vector<vector<pair<string,string>>> trocas_3(int i, int j, vector<vector<pair<string,string>>> &individuo_decodificado,int linha, int coluna);
  int fitness;
vector<vector<vector<pair<double,double>>>>  populacao_nao_elitista;
vector<vector<vector<pair<double,double>>>> populacao_elitista;
  void gerar_populacao();
  void print_individuo(vector<vector<pair<double,double>>> &individuo_codificado);
  void print_decodificado(vector<vector<pair<string,string>>> &individuo_decodificado);
  void print_populacao(vector<vector<vector<pair<double,double>>>> &populacao);
  void print_populacao_elitista(vector<vector<vector<pair<double,double>>>> &populacao_elitista);
  void print_populacao_nao_elitista(vector<vector<vector<pair<double,double>>>> &populacao_nao_elitista);
  vector<vector<pair<string, string>>> trocar_turno(vector<vector<pair<string, string>>> &individuo_decodificado);
  

  void fitness_populacao(vector<vector<vector<pair<double,double>>>> &populacao);
  vector<pair<int,int>> ordena_populacao_custo(vector<vector<vector<pair<double,double>>>> &populacao);
  void dividir_populacao(vector<vector<vector<pair<double,double>>>> &populacao);
  vector<vector<vector<pair<double,double>>>> vetor_mutante;
  vector<vector<pair<double,double>>> mutacao();
 vector<vector<pair<double,double>>> mutacao_2();
 void mutacao_3_vetor(int sem_melhora, int i);
  void mutacao_vetor(int sem_melhora);
 void  mutacao_2_vetor(int sem_melhora, int iteracao);
  vector<vector<pair<double,double>>> cruzamento(vector<vector<pair<double,double>>> &pai_elitista, vector<vector<pair<double,double>>> &pai_nao_elitista);
  vector<vector<pair<double,double>>> ciclo_evolutivo();
  vector<vector<pair<double,double>>> ciclo_evolutivo_ponderado(double &tempo);
  vector<vector<vector<pair<double, double>>>> mutacao_combinada_ponderada();
  vector<vector<vector<pair<double, double>>>> mutacao_combinada();
  vector<vector<pair<double, double>>> codificador_2(vector<vector<pair<string, string>>> &individuo_decodificado);
  
  // movimentos
  vector<vector<pair<string, string>>> primeira_melhora(vector<vector<pair<string, string>>> &individuo_decodificado);
  vector<vector<pair<string,string>>> procura_melhor_vizinho(vector<vector<pair<string,string>>> &individuo_decodificado);
  vector<vector<pair<string,string>>> melhor_melhora(vector<vector<pair<string,string>>> &individuo_decodificado);
  vector<vector<pair<string,string>>> trocar_qualificacao(vector<vector<pair<string,string>>> &individuo_decodificado);
  vector<vector<pair<string,string>>> trocar_qualificacao_turno(vector<vector<pair<string,string>>> &individuo_decodificado);
  vector<vector<pair<string,string>>> trocar_qualificacao_enfermeiro(vector<vector<pair<string, string>>> &individuo_decodificado, vector<int> sorteio, int posicao, int fitness, int i);
  vector<vector<pair<double,double>>> codificador(vector<vector<pair<string,string>>>&individuo_decodificado); 
  vector<vector<pair<double, double>>> codificador_ponderado(vector<vector<pair<string, string>>> &individuo_decodificado);
  void melhora_populacao_inicial_ponderado();
  void melhora_populacao_inicial();
  void melhora_populacao();
  vector<vector<pair<double, double>>> cruzamento_teste(vector<vector<pair<double, double>>> &pai_elite, vector<vector<pair<double, double>>> &pai_nao_elite);
  
  void matriz_compatibilidades();
  vector<Vetor_compatibilidade> vetor_compatibilidade;
  vector<vector<pair<double, double>>> VND_inicial(vector<vector<pair<double, double>>> &individuo_codificado);
  vector<vector<pair<double, double>>> VND_inicial_normal(vector<vector<pair<double, double>>> &individuo_codificado);
  vector<vector<pair<double, double>>> VND_normal(vector<vector<pair<double, double>>> &individuo_codificado);
  vector<vector<pair<string,string>>> h2(vector<vector<pair<string,string>>> &individuo_decodificado);
  bool troca_decoficador(vector<vector<pair<string,string>>> &individuo_decodificado);
  bool primeira_melhora_2_aleatorio_vnd_thalles(vector<vector<pair<string, string>>> &individuo_decodificado);
  vector<vector<vector<pair<double, double>>>> cruzamento_ponto_corte_2(vector<vector<pair<double, double>>> &pai_elite, vector<vector<pair<double, double>>> &pai_nao_elite);
  vector<vector<vector<pair<double, double>>>> resultado_cruzamento;
  //vector<vector<pair<double,double>>> troca_decoficador(vector<vector<pair<double,double>>> &individuo_codificado);
  //MT
  bool primeira_melhora_2_turno_blocos_aleatorio_2(vector<vector<pair<string, string>>> &individuo_decodificado);
  MTRand mt;
  int constante;
  vector<int> custo_colunas;
  double randzero();
void funcao_avaliacao_colunas(vector<vector<pair<string, string>>> &individuo_decodificado);
int funcao_avaliacao_linha(vector<vector<pair<string, string>>> &individuo_decodificado, int linha);
int funcao_avaliacao_linha_coluna_2(vector<vector<pair<string, string>>> &individuo_decodificado, int linha,vector<int> colunas);
bool primeira_melhora_2_turno_blocos_aleatorio_2_total(vector<vector<pair<string, string>>> &individuo_decodificado);
void matriz_compatibilidades_2();

  //destrutor
  virtual ~BRKGA();
  
}; 