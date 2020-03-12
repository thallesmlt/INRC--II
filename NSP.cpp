// Arquivo Responsável pela leitura dos arquivos do problema
#include "NSP.h"

int num_enfermeiros; // variavel global, número de enfermeiros da instancia




vector<double> retornar_pares_double(string palavra) // elimina () e , de uma palavra
{
	vector<double> conversao;
	conversao = vector<double>(2);
	string chars = "()"; // string que contém os caracteres que serão eliminados

	for (char c : chars)
	{
		palavra.erase(remove(palavra.begin(), palavra.end(), c), palavra.end()); //remove todas as ocorrencias dos caracters selecionados
	}

	size_t pos = 0;											//posicao inicial zero
	string delimiter = ",";									// a virgula delimita o tamanho do inteiro da esquerda e da direita
	string auxiliar;										// string para divisao
	while ((pos = palavra.find(delimiter)) != string::npos) //enquanto a posicao for diferente da posicao do delimitador
	{
		auxiliar = palavra.substr(0, pos);			//token recebe o conteudo da posicao atual da string
		palavra.erase(0, pos + delimiter.length()); //apaga a posicao da string inicial
	}
	double first = stod(auxiliar);
	double second = stod(palavra);
	conversao[0] = first;
	conversao[1] = second;
	return conversao;
}



Cenario::Cenario(string nomeArq) : nomeArq(nomeArq)
{
	//abrir o arquivo
	f = ifstream(nomeArq);

	//verificar se o arquivo foi aberto com sucesso
	if (!f.is_open())
		cout << "Erro ao abrir o arquivo!!!!  cenario" << endl;

	
	while (!f.eof())
	{
		f >> s;
		if( s == "SCENARIO"){
			f >> s;
			f >> cenario; 
		}

		if (s == "WEEKS")
		{
			f >> s;
			f >> num_semanas;
		}
		
		if (s == "SHIFT_TYPES")
		{
			for (int i = 0; i < 2; i++)
			{
				f >> s;
			}
			vetor_sucessoes = vector<Vetor_sucessoes>(4);
			for (int i = 0; i < 4; i++)
			{
				f >> vetor_sucessoes[i].turno;
				f >> s;
				vector<int> aux;
				aux = retornar_pares(s);
				vetor_sucessoes[i].x.first = aux[0];
				vetor_sucessoes[i].x.second = aux[1];
			}	
		}
												  
		if (s == "FORBIDDEN_SHIFT_TYPES_SUCCESSIONS") // inicia a leitura das sucessoes de sucessões proibidas
		{
			vetor_turnos = vector<string>(5);
			vetor_turnos[0] = "None";
			vetor_turnos[1] = "Early";
			vetor_turnos[2] = "Day";
			vetor_turnos[3] = "Late";
			vetor_turnos[4] = "Night";
			vetor_none = vetor_turnos;
			for(int i = 0 ; i < 4 ; i++){
				f >> s;
				f >> n;
				if(n>0){
				for(int j = 0; j < n; j++){
					f >> s;
					for(int h = 0; h < vetor_turnos.size(); h++){
						if(vetor_turnos[h] == s){
							vetor_turnos.erase(vetor_turnos.begin() + h); //apaga o s proibido do vetor de sucessoes
						}
					}
				}
			}
				if(i == 0){
					
					vetor_early = vector<string>(5-n);
					vetor_early = vetor_turnos;
				}
				if(i == 1){
					vetor_day = vector<string>(5-n);
					vetor_day = vetor_turnos;
				}
				if(i == 2){
					vetor_late = vector<string>(5-n);
					vetor_late = vetor_turnos;
				}
				if(i == 3){
					vetor_night = vector<string>(5-n);
					vetor_night = vetor_turnos;
				}
			
			vetor_turnos = vetor_none;
		}
		cout << "Vetor Late" << endl;
		for(int i = 0;i<vetor_late.size();i++){
			cout << vetor_late[i] << ' ';
		}
		cout << endl;
	}
	 // fim da sucessao de turnos proibidos

		if (s == "CONTRACTS")
		{
			f >> s;
			f >> n;
			vetor_contratos = vector<Vetor_contratos>(n);
			for (int i = 0; i < vetor_contratos.size(); i++)
			{
				f >> vetor_contratos[i].contrato;
				vetor_contratos[i].min_max = vector<pair<int, int>>(3);
				vetor_contratos[i].min_max_1 = vector<pair<double, double>>(1);
				for (int j = 0; j < 3; j++)
				{
					if (j == 0)
					{
						vector<double> aux;
						f >> s;
						aux = retornar_pares_double(s);
						vetor_contratos[i].min_max_1[j].first = aux[0];
						vetor_contratos[i].min_max_1[j].second = aux[1];
						vector<int> aux_1;
						aux_1 = retornar_pares(s);
						vetor_contratos[i].min_max[j].first = aux_1[0];
						vetor_contratos[i].min_max[j].second = aux_1[1];
					}
					else
					{
						vector<int> aux_1;
						f >> s;
						aux_1 = retornar_pares(s);
						vetor_contratos[i].min_max[j].first = aux_1[0];
						vetor_contratos[i].min_max[j].second = aux_1[1];
					}
				}
				f >> vetor_contratos[i].fim_semana;
				f >> vetor_contratos[i].fim_semana_restricao;
			}
		}

		if (s == "NURSES")
		{
			f >> s;				  // pula até chegar no número de enfermeiros da instancia
			f >> num_enfermeiros; // número de enfermeiros da instancia
			numero_enfermeiros = num_enfermeiros;
			vetor_qualificacao = vector<Vetor_qualificacao>(num_enfermeiros);

			for (int i = 0; i < num_enfermeiros; i++)
			{
				f >> vetor_qualificacao[i].nurse;
				f >> vetor_qualificacao[i].contract;
				f >> vetor_qualificacao[i].num_funcoes;
				vetor_qualificacao[i].funcoes = vector<string>(vetor_qualificacao[i].num_funcoes);

				if (vetor_qualificacao[i].num_funcoes == 1 || vetor_qualificacao[i].num_funcoes == 2 || vetor_qualificacao[i].num_funcoes == 3) // comparacao é com vetor_qualificacao[i].num_funcoes
				{
					f >> vetor_qualificacao[i].funcoes[0];
				}

				if (vetor_qualificacao[i].num_funcoes == 2 || vetor_qualificacao[i].num_funcoes == 3)
				{
					f >> vetor_qualificacao[i].funcoes[1];
				}
				if (vetor_qualificacao[i].num_funcoes == 3)
				{
					f >> vetor_qualificacao[i].funcoes[2];
				}
			}
		}	
	}
}

Cenario::~Cenario()
{
	f.close();
}

Week::Week(string nomeArq) : nomeArq(nomeArq)
{
	w = ifstream(nomeArq);

	//verificar se o arquivo foi aberto com sucesso
	if (!w.is_open())
		cout << "Erro ao abrir o arquivo! week" << endl;

	
	while (!w.eof())
	{			
		w >> s; 
		if (s == "REQUIREMENTS")
		{
			matriz_requerimentos = vector<vector<pair<int,int>>>(16,vector<pair<int,int>>(7));
			for (int i = 0; i < 16; i++) 
			{
				w >> s;
				w >> s;
				for (int j = 0; j < 7; j++)
				{
					w >> s;
					vector<int> aux;
					aux = retornar_pares(s);
					matriz_requerimentos[i][j].first = aux[0];
					matriz_requerimentos[i][j].second = aux[1];
				}
			}
		}

		if (s == "SHIFT_OFF_REQUESTS")
		{
			w >> s;
			w >> tam_matriz;
			matriz_requisicoes = vector<vector<int>>(tam_matriz, vector<int>(3));
			for (int i = 0; i < tam_matriz; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					w >> s;
					if(j == 0){
						matriz_requisicoes[i][j] = retornar_enfermeiro(s);
					}
					if(j == 1){
						if(s == "Any"){
							matriz_requisicoes[i][j] = 0;
						}
						if(s == "Early"){
							matriz_requisicoes[i][j] = 1;
						}
						if(s == "Day"){
							matriz_requisicoes[i][j] = 2;
						}
						if(s == "Late"){
							matriz_requisicoes[i][j] = 3;
						}
						if(s == "Night"){
							matriz_requisicoes[i][j] = 4;
						}
					}
					if(j == 2){
						if(s == "Mon"){
							matriz_requisicoes[i][j] = 0;
						}
						if(s == "Tue"){
							matriz_requisicoes[i][j] = 1;
						}
						if(s == "Wed"){
							matriz_requisicoes[i][j] = 2;
						}
						if(s == "Thu"){
							matriz_requisicoes[i][j] = 3;
						}
						if(s == "Fri"){
							matriz_requisicoes[i][j] = 4;
						}
						if(s == "Sat"){
							matriz_requisicoes[i][j] = 5;
						}
						if(s == "Sun"){
							matriz_requisicoes[i][j] = 6;
						}
					}
				}
			}
		}
	}
}

Week::~Week()
{
	w.close();
}

History::History(string nomeArq) : nomeArq(nomeArq)
{
	h = ifstream(nomeArq);

	
	if (!h.is_open())
		cout << "Erro ao abrir o arquivo! history" << endl;

	matriz_borda_s = vector<vector<string>>(num_enfermeiros, vector<string>(2));
	matriz_borda_fracas = vector<vector<int>>(num_enfermeiros, vector<int>(5));

	while (!h.eof())
	{			
		h >> s; 
		if (s == "NURSE_HISTORY")
		{
			for (int i = 0; i < num_enfermeiros; i++)
			{
				int count = 0; 
				for (int j = 0; j < 7; j++)
				{
					if (j == 0)
					{
						h >> matriz_borda_s[i][j];
						count++;
					}
					else if (j == 3)
					{
						h >> matriz_borda_s[i][1];
						count++;
					}
					else
					{
						h >> matriz_borda_fracas[i][j - count];
					}
				}
			}
		}
	}
}

History::~History()
{
	h.close();
}

vector<int> retornar_pares(string palavra)
{
	vector<int> conversao;
	conversao = vector<int>(2);
	string chars = "()"; // string que contém os caracteres que serão eliminados

	for (char c : chars)
	{
		palavra.erase(remove(palavra.begin(), palavra.end(), c), palavra.end()); //remove todas as ocorrencias dos caracters selecionados
	}

	size_t pos = 0;											//posicao inicial zero
	string delimiter = ",";									// a virgula delimita o tamanho do inteiro da esquerda e da direita
	string auxiliar;										// string para divisao
	while ((pos = palavra.find(delimiter)) != string::npos) //enquanto a posicao for diferente da posicao do delimitador
	{
		auxiliar = palavra.substr(0, pos);			//token recebe o conteudo da posicao atual da string
		palavra.erase(0, pos + delimiter.length()); //apaga a posicao da string inicial
	}
	int first = stoi(auxiliar);
	int second = stoi(palavra);
	conversao[0] = first;
	conversao[1] = second;
	return conversao;
}

int retornar_enfermeiro(string palavra)
{
	int conversao;
	string chars = "HNUCTR_"; // string que contém os caracteres que serão eliminados

	for (char c : chars)
	{
		palavra.erase(remove(palavra.begin(), palavra.end(), c), palavra.end()); //remove todas as ocorrencias dos caracters selecionados
	}
	conversao = stoi(palavra);	
	return conversao;
}


Custom::Custom(string nomeArq) : nomeArq(nomeArq)
{
	
	cust = ifstream(nomeArq);
	cout << nomeArq << endl;

	
	if (!cust.is_open())
		cout << "Erro ao abrir o Custom!!!! Custom" << endl;

	cust >> n;
	dias_trabalhados = vector<pair<int, int>>(n);
	cout << "dias trabalhados" << endl;
	for (int i = 0; i < n; i++)
	{
		cust >> dias_trabalhados[i].first;
		cout << dias_trabalhados[i].first << ' ';
		cust >> dias_trabalhados[i].second;
		cout << dias_trabalhados[i].second << endl;
	}
}

Custom::~Custom()
{
	cust.close();
}



