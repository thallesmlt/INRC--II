    #include "BRKGA.h"
#include <algorithm>

BRKGA::BRKGA(int TAM_POPULACAO, int NUM_GERACOES, double TAM_ELITISTA, double TAM_MUTANTE, double PROB_CRUZAMENTO, vector<string> names, string cusin, int iteracoes_max, MTRand mtrand)
{
    cout << endl;
    cout << "cenario:" << names[0] << endl;
    cout << "historico>" << names[1] << endl;
    cout << "semana:" << names[2] << endl;
    cout << "customout:" << names[3] << endl;
    cout << "solucao:" << names[4] << endl;
    if (cusin == string("--cusIn"))
    {
        cus = new Custom(names[5]);
        cout << "custonin:" << names[5] << endl;
        this->semana_1 = 0;
    }
    else
    {
        this->semana_1 = 1;
    }
    cn = new Cenario(names[0]);
    ht = new History(names[1]);
    w0 = new Week(names[2]);
    this->nomeArq = string(names[4]);
    this->custom_out = string(names[3]);

    this->tam_populacao = TAM_POPULACAO;
    this->num_geracoes = NUM_GERACOES;
    this->tam_elitista = (TAM_ELITISTA * TAM_POPULACAO);
    this->tam_nao_elitista = (TAM_POPULACAO - tam_elitista);
    this->tam_mutante = (TAM_MUTANTE * TAM_POPULACAO);
    this->prob_cruzamento = PROB_CRUZAMENTO;
    this->tam_cruzamento = TAM_POPULACAO - tam_elitista - tam_mutante;
    this->populacao_elitista = vector<vector<vector<pair<double, double>>>>(this->tam_elitista, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    this->populacao_nao_elitista = vector<vector<vector<pair<double, double>>>>(this->tam_nao_elitista, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    this->mt = mtrand;
    cout << "Tamanho populacoes" << endl;
    cout << tam_populacao << endl;
    cout << tam_elitista << endl;
    cout << tam_mutante << endl;
    cout << tam_cruzamento << endl;
    cout << tam_nao_elitista << endl;
    this->temp = temp;
    multi_start(iteracoes_max); // execução Do Ag e busca local

}
BRKGA::~BRKGA()
{
    delete cn;
    delete ht;
    delete w0;
    delete cus;
}

double BRKGA::randzero()
{
    return mt.rand();
}

void BRKGA::matriz_compatibilidades() // cria conjuntos de enfermeiros compativeis 
{
    int tamanho = 0;
    int count = 0;
    bool igual;
    int aux = 0;
    int num_vetor_compatibilidade = 0;
    vector<Vetor_compatibilidade> vetor_compatibilidade_aux;
    vetor_compatibilidade_aux = vector<Vetor_compatibilidade>(8);

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        cout << "Entrou" << ' ' << i << endl;

        if (aux == 0)
        {
            tamanho++;
            vetor_compatibilidade_aux[num_vetor_compatibilidade].compatibilidade.push_back(i);
            
            aux++;
        }
        else
        {
            
            if (cn->vetor_qualificacao[vetor_compatibilidade_aux[num_vetor_compatibilidade].compatibilidade[aux - 1]].num_funcoes == cn->vetor_qualificacao[i].num_funcoes)
            {

                for (int j = 0; j < cn->vetor_qualificacao[i].num_funcoes; j++)
                {

                    if (cn->vetor_qualificacao[vetor_compatibilidade_aux[num_vetor_compatibilidade].compatibilidade[aux - 1]].funcoes[j] != cn->vetor_qualificacao[i].funcoes[j])
                    {
                        igual = false;

                        break;
                    }
                    else
                    {
                        igual = true;
                    }
                }

                if (igual == false)
                {

                    aux = 0;
                    num_vetor_compatibilidade++;
                    i--;
                }
                else
                {
                    vetor_compatibilidade_aux[num_vetor_compatibilidade].compatibilidade.push_back(i);

                    aux++;
                }
            }
            else
            {
                aux = 0;
                num_vetor_compatibilidade++;
                i--;
            }
        }
    }

    
    vetor_compatibilidade = vector<Vetor_compatibilidade>(tamanho);
    

    for (int i = 0; i < tamanho; i++)
    {
        vetor_compatibilidade[i] = vetor_compatibilidade_aux[i];
    }
}

vector<vector<pair<double, double>>> BRKGA::criar_individuo() // cria uma solução representada por chaves aleatórias
{
    individuo_codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));


    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            individuo_codificado[i][j].first = randzero();
            individuo_codificado[i][j].second = randzero();
        }
    }
    return individuo_codificado;
}

vector<vector<pair<string, string>>> BRKGA::decodificar_individuo(vector<vector<pair<double, double>>> &individuo_codificado)
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int posicao = 0;

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_none.size());
                    individuo_decodificado[i][j].first = cn->vetor_none[posicao];
                }
                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_early.size());
                    individuo_decodificado[i][j].first = cn->vetor_early[posicao];
                }
                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_day.size());
                    individuo_decodificado[i][j].first = cn->vetor_day[posicao];
                }
                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_late.size());
                    individuo_decodificado[i][j].first = cn->vetor_late[posicao];
                }
                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_night.size());
                    individuo_decodificado[i][j].first = cn->vetor_night[posicao];
                }
            }
            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_none.size());
                    individuo_decodificado[i][j].first = cn->vetor_none[posicao];
                }
                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_early.size());
                    individuo_decodificado[i][j].first = cn->vetor_early[posicao];
                }
                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_day.size());
                    individuo_decodificado[i][j].first = cn->vetor_day[posicao];
                }
                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_late.size());
                    individuo_decodificado[i][j].first = cn->vetor_late[posicao];
                }
                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_night.size());
                    individuo_decodificado[i][j].first = cn->vetor_night[posicao];
                }
            }

            posicao = (cn->vetor_qualificacao[i].num_funcoes * individuo_codificado[i][j].second);

            if (posicao == 0)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[0];
            }
            if (posicao == 1)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[1];
            }
            if (posicao == 2)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[2];
            }
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::decodificar_individuo_2(vector<vector<pair<double, double>>> &individuo_codificado)
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int posicao = 0;

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            posicao = individuo_codificado[i][j].first * cn->vetor_none.size();
            individuo_decodificado[i][j].first = cn->vetor_none[posicao];

            posicao = (cn->vetor_qualificacao[i].num_funcoes * individuo_codificado[i][j].second);

            if (posicao == 0)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[0];
            }
            if (posicao == 1)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[1];
            }
            if (posicao == 2)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[2];
            }
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::decodificar_individuo_2_sol(vector<vector<pair<double, double>>> &individuo_codificado, int &sol_alocacoes)
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int posicao = 0;

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            posicao = individuo_codificado[i][j].first * cn->vetor_none.size();
            individuo_decodificado[i][j].first = cn->vetor_none[posicao];
            posicao = (cn->vetor_qualificacao[i].num_funcoes * individuo_codificado[i][j].second);

            if (individuo_decodificado[i][j].first != "None")
            {
                sol_alocacoes++;
            }

            if (posicao == 0)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[0];
            }
            if (posicao == 1)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[1];
            }
            if (posicao == 2)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[2];
            }
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::decodificar_individuo_sol(vector<vector<pair<double, double>>> &individuo_codificado, int &sol_alocacoes)
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int posicao = 0;

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_none.size());
                    individuo_decodificado[i][j].first = cn->vetor_none[posicao];
                }
                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_early.size());
                    individuo_decodificado[i][j].first = cn->vetor_early[posicao];
                }
                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_day.size());
                    individuo_decodificado[i][j].first = cn->vetor_day[posicao];
                }
                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_late.size());
                    individuo_decodificado[i][j].first = cn->vetor_late[posicao];
                }
                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_night.size());
                    individuo_decodificado[i][j].first = cn->vetor_night[posicao];
                }
            }
            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_none.size());
                    individuo_decodificado[i][j].first = cn->vetor_none[posicao];
                }
                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_early.size());
                    individuo_decodificado[i][j].first = cn->vetor_early[posicao];
                }
                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_day.size());
                    individuo_decodificado[i][j].first = cn->vetor_day[posicao];
                }
                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_late.size());
                    individuo_decodificado[i][j].first = cn->vetor_late[posicao];
                }
                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    posicao = (individuo_codificado[i][j].first * cn->vetor_night.size());
                    individuo_decodificado[i][j].first = cn->vetor_night[posicao];
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                sol_alocacoes++;
            }

            posicao = (cn->vetor_qualificacao[i].num_funcoes * individuo_codificado[i][j].second);

            if (posicao == 0)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[0];
            }
            if (posicao == 1)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[1];
            }
            if (posicao == 2)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[2];
            }
        }
    }
    return individuo_decodificado;
}

int BRKGA::funcao_avaliacao(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int s7 = 0;
    int sucessoes = 0;
    bool contem;
    int aux = 0;
    int count_dias_trabalhados = 0;
    int dias_global = 0;
    int i = 0;
    int count = 0;
    int turnos_consecutivos = 0;
    int count_folgas = 0;
    int dias_folgas = 0;
    string turno_anterior;
    int count_dias_consecutivos = 0;
    int dias_consecutivos = 0;
    int count_fim_semana = 0;
    int count_turnos = 0;
    int custo_total;
    int requisicoes = 0;
    int fim_semana_completo = 0;
    int requerimento_otimo = 0;
    int requerimento_minimo = 0;
    vector<pair<int, int>> turno_funcao;
    turno_funcao = vector<pair<int, int>>(16); // coloca todas as posições do vetor com valor 0
    for (int i = 0; i < 16; i++)
    {
        turno_funcao[i].first = 0;
        turno_funcao[i].second = 0;
    }

    for (i = 0; i < 7; i++)
    {
        for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
        {
            if (individuo_decodificado[j][i].first == "Early")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[0].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[1].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[2].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[3].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Day")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[4].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[5].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[6].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[7].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Late")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[8].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[9].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[10].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[11].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Night")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[12].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[13].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[14].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[15].first++;
                }
            }
            if (i == 5)
            {
                if (cn->vetor_qualificacao[j].contract == "FullTime" && cn->vetor_contratos[0].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[j].contract == "PartTime" && cn->vetor_contratos[1].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[j].contract == "HalfTime" && cn->vetor_contratos[2].fim_semana_restricao == 1)
                {

                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
                if (cn->vetor_qualificacao[j].contract == "20Percent" && cn->vetor_contratos.size() == 4 && cn->vetor_contratos[3].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
            }
        }

        for (int h = 0; h < 16; h++)
        {
            //fazer aqui o somatorio
            if (turno_funcao[h].first < w0->matriz_requerimentos[h][i].first)
            {
                requerimento_minimo = (requerimento_minimo + (w0->matriz_requerimentos[h][i].first - turno_funcao[h].first));
            }
            if (turno_funcao[h].first < w0->matriz_requerimentos[h][i].second)
            {
                requerimento_otimo = (requerimento_otimo + (w0->matriz_requerimentos[h][i].second - turno_funcao[h].first));
            }
            turno_funcao[h].first = 0;
        }
    }

    for (int i = 0; i < w0->matriz_requisicoes.size(); i++)
    {
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first != "None" && w0->matriz_requisicoes[i][1] == 0)
        {
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Early" && w0->matriz_requisicoes[i][1] == 1)
        {
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Day" && w0->matriz_requisicoes[i][1] == 2)
        {
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Late" && w0->matriz_requisicoes[i][1] == 3)
        {
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Night" && w0->matriz_requisicoes[i][1] == 4)
        {
            requisicoes++;
        }
    }

    for (int h = 0; h < custom_out.size(); h++)
    {
        if (isdigit(custom_out[h]))
        {
            aux = custom_out[h] - '0';
            aux = cn->num_semanas - aux;
            
        }
    } 

    // restricoes s2 e s3 ---- s2-peso 30 resolvida ---- s3 peso 30 a ser resolvida
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (individuo_decodificado[i][j].first == "None")
            {
                count_folgas++;
                if (j == 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][4] < cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + ht->matriz_borda_fracas[i][4];
                            }
                            if (ht->matriz_borda_fracas[i][4] >= cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + cn->vetor_contratos[h].min_max[2].second;
                            }
                            break;
                        }
                    }

                    // faz a subtração dos turnos que restarem da semana anterior
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].first)
                            {
                                count_dias_consecutivos = ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = cn->vetor_contratos[h].min_max[1].second;
                            }
                            break;
                        }
                    }
                }

                if (count_dias_consecutivos > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos < cn->vetor_contratos[h].min_max[1].first)
                            {
                                dias_consecutivos = dias_consecutivos + (cn->vetor_contratos[h].min_max[1].first - count_dias_consecutivos);
                            }
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                            count_dias_consecutivos = 0;
                            break;
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                       
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                    
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                                break;
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_consecutivos++;
                if (j == 0)
                {
                    count_folgas = ht->matriz_borda_fracas[i][4];
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + cn->vetor_contratos[h].min_max[1].second;
                            }
                            break;
                        }
                    }
                }

                if (count_folgas > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_folgas < cn->vetor_contratos[h].min_max[2].first)
                            {
                                dias_folgas = dias_folgas + (cn->vetor_contratos[h].min_max[2].first - count_folgas);
                            }
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                            }
                            count_folgas = 0;
                            break;
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                            break;
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_trabalhados++;
                if (j == 0)
                {
                    count_turnos++;
                    if (ht->matriz_borda_s[i][1] == individuo_decodificado[i][j].first)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + ht->matriz_borda_fracas[i][2];
                                }
                                if (ht->matriz_borda_fracas[i][2] >= cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + cn->vetor_sucessoes[h].x.second;
                                }
                                break;
                            }
                        }
                    }

                    if (ht->matriz_borda_s[i][1] != individuo_decodificado[i][j].first)
                    {
                        if (ht->matriz_borda_fracas[i][2] > 0)
                        {
                            for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                            {
                                if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                                {
                                    if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                    {
                                        turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first == turno_anterior)
                {
                    count_turnos++;
                    if (j == 6)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == turno_anterior)
                            {
                                if (count_turnos > cn->vetor_sucessoes[h].x.second)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                                }
                                break;
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first != turno_anterior && j != 0)
                {
                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {
                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                            break;
                        }
                    }
                    count_turnos = 1;
                }
                turno_anterior = individuo_decodificado[i][j].first;
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                if (j == 0)
                {
                    if (ht->matriz_borda_fracas[i][2] > 0)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                }
                                break;
                            }
                        }
                    }
                }
                if (count_turnos > 0)
                {

                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {

                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                            break;
                        }
                    }
                }

                count_turnos = 0;
                turno_anterior = individuo_decodificado[i][j].first;
            }
        }
        count_dias_consecutivos = 0;
        count_folgas = 0;
        count_turnos = 0;
        turno_anterior.erase();
        

        //restriçao s6
        if (semana_1 == 0)
        {
            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        //s7++;
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }

                    
                    int divisao_first = floor(double(cus->dias_trabalhados[i].first / aux));
                    int divisao_second = ceil(double(cus->dias_trabalhados[i].second / aux));

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                        //dias_global++;
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                        //dias_global++;
                    }
                    break;
                }
            }
        }

        if (semana_1 == 1)
        {

            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        //s7++;
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }
                    int divisao_first = floor(double(cn->vetor_contratos[t].min_max[0].first / cn->num_semanas));
                    int divisao_second = ceil(double(cn->vetor_contratos[t].min_max[0].second / cn->num_semanas));

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                        //dias_global++;
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                        //dias_global++;
                    }
                    break;
                }
            }
        }
        
        count_dias_trabalhados = 0;
    }

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            contem = false;
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }

            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }
        }
    }

    requerimento_otimo = (requerimento_otimo * 30);
    requisicoes = (requisicoes * 10);
    fim_semana_completo = (fim_semana_completo * 30);
    dias_consecutivos = (dias_consecutivos * 30);
    dias_folgas = (dias_folgas * 30);
    turnos_consecutivos = (turnos_consecutivos * 15);
    if (aux <= cn->num_semanas/2) // antes da metade do horizonte de planejamenti
    {
        sucessoes = (sucessoes * 60);
        requerimento_minimo = (requerimento_minimo * 60);
        s7 = (s7 *45);
        dias_global = (dias_global * 20);
    }
    else{ // depois da metade do horizonte de planejamento
        sucessoes = (sucessoes * 60);
    requerimento_minimo = (requerimento_minimo * 60);
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }
         

    
    custo_total = s7 + sucessoes + requerimento_minimo + requerimento_otimo + requisicoes + fim_semana_completo + dias_consecutivos + dias_folgas + turnos_consecutivos + dias_global;
    return custo_total;
}

void BRKGA::brincadeira()
{
    int aux = 0;
    porcentagem_contratos();
    melhora_populacao_inicial();
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < tam_populacao; j++)
        {
            individuo_decodificado = decodificar_individuo_2(populacao[j]);
            individuo_decodificado = primeira_melhora(individuo_decodificado);
            individuo_decodificado = procura_melhor_vizinho(individuo_decodificado);
            populacao[j] = codificador_2(individuo_decodificado);
        }
        if (i == 4)
        {
            for (int j = 0; j < tam_populacao; j++)
            {
                individuo_decodificado = decodificar_individuo_2(populacao[j]);
                print_decodificado(individuo_decodificado);
                aux = funcao_avaliacao(individuo_decodificado);
                cout << "Fitness = " << aux << endl;
            }
        }
    }
}

void BRKGA::print_individuo(vector<vector<pair<double, double>>> &individuo_codificado)
{
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            cout << "(" << individuo_codificado[i][j].first << ' ';
            cout << individuo_codificado[i][j].second << ")" << ' ';
        }
        cout << endl;
    }
}

void BRKGA::print_decodificado(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        cout << i << ' ';
        for (int j = 0; j < 7; j++)
        {
            cout << "(" << individuo_decodificado[i][j].first << ' ';
            cout << individuo_decodificado[i][j].second << ")" << ' ';
        }
        cout << endl;
    }
}

void BRKGA::gerar_populacao()
{
    populacao = vector<vector<vector<pair<double, double>>>>(this->tam_populacao, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    for (int i = 0; i < this->tam_populacao; i++)
    {
        populacao[i] = criar_individuo();
    }
}

void BRKGA::print_populacao(vector<vector<vector<pair<double, double>>>> &populacao)
{
    for (int i = 0; i < tam_populacao; i++)
    {
        print_individuo(populacao[i]);
        cout << endl;
    }
}

void BRKGA::print_populacao_elitista(vector<vector<vector<pair<double, double>>>> &populacao_elitista)
{
    for (int i = 0; i < populacao_elitista.size(); i++)
    {
        print_individuo(populacao_elitista[i]);
        cout << endl;
    }
}

void BRKGA::print_populacao_nao_elitista(vector<vector<vector<pair<double, double>>>> &populacao_nao_elitista)
{
    for (int i = 0; i < populacao_nao_elitista.size(); i++)
    {
        print_individuo(populacao_nao_elitista[i]);
        cout << endl;
    }
}

void BRKGA::fitness_populacao(vector<vector<vector<pair<double, double>>>> &populacao)
{
    vector<vector<pair<string, string>>> sol;
    sol = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));

    for (int i = 0; i < this->tam_populacao; i++)
    {
        sol = decodificar_individuo(populacao[i]);
        funcao_avaliacao(sol);
        sol.clear();
    }
}

vector<pair<int, int>> BRKGA::ordena_populacao_custo(vector<vector<vector<pair<double, double>>>> &populacao)
{
    vector<pair<int, int>> ordenacao(tam_populacao);
    vector<vector<pair<string, string>>> sol;
    sol = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    for (int i = 0; i < ordenacao.size(); i++)
    {
        sol = decodificar_individuo(populacao[i]);
        ordenacao[i].first = funcao_avaliacao(sol);
        ordenacao[i].second = i;
        sol.clear();
    }
    sort(ordenacao.begin(), ordenacao.end());
    return ordenacao;
}

void BRKGA::dividir_populacao(vector<vector<vector<pair<double, double>>>> &populacao)
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int aux = 0;
    vector<pair<int, int>> ordenacao(tam_populacao);
    vector<vector<pair<string, string>>> sol;
    sol = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));

    for (int i = 0; i < ordenacao.size(); i++)
    {
        sol = decodificar_individuo_2(populacao[i]);
        ordenacao[i].first = funcao_avaliacao(sol);
        ordenacao[i].second = i;
        //cout << ordenacao[i].first << ' ';
        //cout << ordenacao[i].second << endl;
        //sol.clear();
    }
    sort(ordenacao.begin(), ordenacao.end());

    for (int j = 0; j < tam_elitista; j++)
    {
        
        populacao_elitista[j] = populacao[ordenacao[j].second];
      
    }
    for (int j = 0; j < tam_nao_elitista; j++)
    {
        populacao_nao_elitista[j] = populacao[ordenacao[j + tam_elitista].second];
    }
    
}

void BRKGA::dividir_populacao_2(vector<vector<vector<pair<double, double>>>> &populacao)
{
    int custo;
    vector<pair<int, int>> ordenacao(tam_populacao);
    vector<vector<pair<string, string>>> sol;
    sol = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    for (int i = 0; i < ordenacao.size(); i++)
    {
        sol = decodificar_individuo_ponderado_2(populacao[i]);
        ordenacao[i].first = funcao_avaliacao(sol);
        ordenacao[i].second = i;
        
    }
    sort(ordenacao.begin(), ordenacao.end());
    
 

    for (int j = 0; j < tam_elitista; j++)
    {
        populacao_elitista[j] = populacao[ordenacao[j].second];

    }
    

    for (int j = 0; j < populacao_nao_elitista.size(); j++)
    {
        populacao_nao_elitista[j] = populacao[ordenacao[j+tam_elitista].second];

    }
}

vector<vector<pair<double, double>>> BRKGA::mutacao()
{
    
    
    individuo_codificado = criar_individuo();
    return individuo_codificado;

}

void BRKGA::mutacao_vetor(int sem_melhora)
{
    vetor_mutante = vector<vector<vector<pair<double, double>>>>(this->tam_mutante, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));

    for (int i = 0; i < vetor_mutante.size(); i++)
    {
        vetor_mutante[i] = criar_individuo();
    }
}

void BRKGA::mutacao_2_vetor(int sem_melhora, int iteracao)
{

    int count;
    int sorteio;
    vector<int> vetor_sorteio = vector<int>(tam_populacao);
    vector<int> vetor_sorteio_2 = vector<int>(cn->vetor_qualificacao.size());
    vector<int> vetor_sorteio_2_aux = vector<int>(cn->vetor_qualificacao.size());
    for (int i = 1; i < vetor_sorteio.size(); i++)
    {
        vetor_sorteio[i] = i;
    }

   
    for (int i = 0; i < vetor_sorteio_2.size(); i++)
    {
        vetor_sorteio_2[i] = i;
    }
    

    random_shuffle(vetor_sorteio.begin(), vetor_sorteio.end());

    if(sem_melhora < 40 ){
        constante = 1;
    }
    if(sem_melhora >= 40 && sem_melhora < 120){
        constante = 2;
    }
    else{
        constante = 2;  
    }
    if(sem_melhora % 60 == 0){
        constante = cn->vetor_qualificacao.size()/2;
    }
    
    int sorteio_tamanho;
    vetor_mutante = vector<vector<vector<pair<double, double>>>>(this->tam_mutante, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
   
    
    for (int i = 0; i < vetor_mutante.size(); i++)
    {
        random_shuffle(vetor_sorteio_2.begin(), vetor_sorteio_2.end());
        vetor_mutante[i] = populacao[vetor_sorteio[i]];
       
        for (int j = 0; j < constante; j++)
        {
            sorteio = (mt.randInt() % 7);
            
            vetor_mutante[i][vetor_sorteio_2[j]][sorteio].first = randzero();
            vetor_mutante[i][vetor_sorteio_2[j]][sorteio].second = randzero();
        }
    }
    for(int j = 0; j < vetor_mutante.size();j++){
        if (sem_melhora % 10 == 0 && sem_melhora < 40 /*&& iteracao > 1000*/)
        {
            vetor_mutante[j] = VND(vetor_mutante[j]);
        }
        else{
            if(sem_melhora % 10 == 0 && sem_melhora >= 40 /*&& iteracao > 1000*/){
                vetor_mutante[j] = VND(vetor_mutante[j]);
                
            }
            
        }
    }
      
}

void BRKGA::mutacao_3_vetor(int sem_melhora, int t)
{
    int sorteio_otimo;
    int limite;
    int sorteio, sorteio_individuo, sorteio_linha;
    vector<int> vetor_sorteio = vector<int>(tam_populacao);

    for (int i = 0; i < vetor_sorteio.size() -1; i++)
    {
        vetor_sorteio[i] = 1 + i;
        
    }

    random_shuffle(vetor_sorteio.begin(), vetor_sorteio.end());

    vetor_mutante = vector<vector<vector<pair<double, double>>>>(this->tam_mutante, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));


    
    for (int i = 0; i < tam_mutante; i++)
    {
        if(sem_melhora >= 50){
            limite = 2;
        }
        /*if(sem_melhora >= 150){
            limite = 8;
        }*/

        if(sem_melhora < 50){
            limite = 5;
        }
        
        vetor_mutante[i] = populacao[vetor_sorteio[i]];
        
        sorteio = (mt.randInt() % 7);
        
        for (int j = 0; j < limite; j++)
        {
            sorteio_linha = mt.randInt() % cn->vetor_qualificacao.size();
            vetor_mutante[i][sorteio_linha][sorteio].first = randzero();
            vetor_mutante[i][sorteio_linha][sorteio].second = randzero();
        }
        
        
        
     if(sem_melhora % 8 == 0 && sem_melhora <= 50)
        {
            vetor_mutante[i] = VND(vetor_mutante[i]);
        }

        if (sem_melhora % 8 == 0 && sem_melhora >= 50)
        {
            
            vetor_mutante[i] = VND(vetor_mutante[i]);
        }

    }
  
}

bool BRKGA::troca_decoficador(vector<vector<pair<string,string>>> &individuo_decodificado)
{
    int requerimentos = 0;
    bool melhora = false;
    int limite,tam_corte;
    individuo_codificado = codificador_ponderado_2(individuo_decodificado);
    int sorteio, sorteio_individuo, sorteio_linha,fitness,fitness_linha;
    vector<int> vetor_sorteio = vector<int>(cn->vetor_qualificacao.size());

    for (int i = 0; i < vetor_sorteio.size(); i++)
    {
        vetor_sorteio[i] = i;
    }
    
    
    vector<pair<double,double>> aux = vector<pair<double,double>>(7);
    vector<pair<double,double>> aux_2 = vector<pair<double,double>>(7);
    random_shuffle(vetor_sorteio.begin(), vetor_sorteio.end());
    

    sorteio = (mt.randInt() % 7);
    limite = 7 - sorteio;
    tam_corte = (mt.randInt() % limite);
    tam_corte++;
    tam_corte = 1;
   
    

    for (int i = 0; i < vetor_sorteio.size() - 1; i++)
    {
        
        
            individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);
            funcao_avaliacao_colunas(individuo_decodificado);
            for(int t = sorteio; t < sorteio + tam_corte; t++){
                requerimentos = requerimentos + custo_colunas[t];
                aux[t].first = individuo_codificado[vetor_sorteio[i]][sorteio].first;
                aux[t].second = individuo_codificado[vetor_sorteio[i]][sorteio].second;
                aux_2[t].first = individuo_codificado[vetor_sorteio[i + 1]][sorteio].first;
                aux_2[t].second = individuo_codificado[vetor_sorteio[i + 1]][sorteio].second;
            }
           
            fitness = requerimentos + funcao_avaliacao_linha(individuo_decodificado,i) + funcao_avaliacao_linha(individuo_decodificado, i + 1);
           
            requerimentos = 0;

            for(int t = sorteio; t < sorteio + tam_corte; t++){//fazer o swap
                individuo_codificado[vetor_sorteio[i]][sorteio].first = aux_2[t].first;
                individuo_codificado[vetor_sorteio[i]][sorteio].second = aux_2[t].second;
                individuo_codificado[vetor_sorteio[i + 1]][sorteio].first = aux[t].first;
                individuo_codificado[vetor_sorteio[i + 1]][sorteio].second = aux[t].second;
            }

            individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);
            funcao_avaliacao_colunas(individuo_decodificado);
            for(int t = sorteio; t < sorteio + tam_corte; t++){
                requerimentos = requerimentos = requerimentos + custo_colunas[t];  
            }
            
            
            fitness_linha = requerimentos + funcao_avaliacao_linha(individuo_decodificado,i) + funcao_avaliacao_linha(individuo_decodificado, i + 1);
            
            requerimentos = 0;
            if(fitness_linha < fitness){
                melhora = true;
                random_shuffle(vetor_sorteio.begin(), vetor_sorteio.end());
                    i = -1;
                    sorteio = (mt.randInt() % 7);
                    limite = 7 - sorteio;
                    tam_corte = (mt.randInt() % limite);
                    tam_corte++;
                    tam_corte = 1;
            }
            else{
                for(int t = sorteio; t < sorteio + tam_corte; t++){
                individuo_codificado[vetor_sorteio[i]][sorteio].first = aux[t].first;
                individuo_codificado[vetor_sorteio[i]][sorteio].second = aux[t].second;
                individuo_codificado[vetor_sorteio[i + 1]][sorteio].first = aux_2[t].first;
                individuo_codificado[vetor_sorteio[i + 1]][sorteio].second = aux_2[t].second;
                }
            } 
        
    }
    
    individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);
    return melhora;
}

vector<vector<pair<double, double>>> BRKGA::mutacao_2()
{

    int sorteio;
    vector<int> vetor_sorteio = vector<int>(this->tam_nao_elitista);
    vector<int> vetor_sorteio_2 = vector<int>(cn->vetor_qualificacao.size());
    vector<int> vetor_sorteio_2_aux = vector<int>(cn->vetor_qualificacao.size());
    for (int i = 0; i < vetor_sorteio.size(); i++)
    {
        vetor_sorteio[i] = tam_elitista + i;
        
    }

    
    for (int i = 0; i < vetor_sorteio_2.size(); i++)
    {
        vetor_sorteio_2[i] = i;
    }
 
    random_shuffle(vetor_sorteio_2.begin(), vetor_sorteio_2.end());
    
    for (int j = 0; j < (cn->vetor_qualificacao.size()); j++)
    {
        sorteio = (mt.randInt() % 7);


        individuo_codificado[vetor_sorteio_2[j]][sorteio].first = randzero();
        individuo_codificado[vetor_sorteio_2[j]][sorteio].second = randzero();

    }
    return individuo_codificado;
}

vector<vector<pair<double, double>>> BRKGA::mutacao_2_cruzamento(vector<vector<pair<double, double>>> &pai_elite, vector<vector<pair<double, double>>> &pai_nao_elite)
{
    int sorteio;
    int moeda_mutacao;
    vector<int> vetor_sorteio = vector<int>(populacao_nao_elitista.size());
    vector<int> vetor_sorteio_2 = vector<int>(cn->vetor_qualificacao.size());
    

    for (int i = 0; i < vetor_sorteio_2.size(); i++)
    {
        vetor_sorteio_2[i] = i;
    }

    int aux = 0;
    double moeda;

    for (int i = 0; i < pai_elite.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            moeda = randzero();
            if (moeda <= prob_cruzamento)
            {

                individuo_codificado[i][j].first = pai_elite[i][j].first;
                individuo_codificado[i][j].second = pai_elite[i][j].second;
            }
            else
            {
                individuo_codificado[i][j].first = pai_nao_elite[i][j].first;
                individuo_codificado[i][j].second = pai_nao_elite[i][j].second;
            }
        }
    }

    
    vetor_mutante = vector<vector<vector<pair<double, double>>>>(this->tam_mutante, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    moeda_mutacao = randzero();
    if (moeda <= 0.8)
    {
        for (int i = tam_elitista; i < tam_populacao; i++)
        {
            random_shuffle(vetor_sorteio_2.begin(), vetor_sorteio_2.end());
            for (int j = 0; j < (cn->vetor_qualificacao.size()); j++)
            {
                sorteio = (mt.randInt() % 7);

                individuo_codificado[vetor_sorteio_2[j]][sorteio].first = randzero();
                individuo_codificado[vetor_sorteio_2[j]][sorteio].second = randzero();
            }
        }
    }
    else
    {
        return individuo_codificado;
    }
}

vector<vector<vector<pair<double, double>>>> BRKGA::mutacao_combinada_ponderada()
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vetor_mutante = vector<vector<vector<pair<double, double>>>>(this->tam_mutante, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    for (int i = 0; i < tam_mutante; i++)
    {
       
        vetor_mutante[i] = criar_individuo();
    
    }
    return vetor_mutante;
}

vector<vector<vector<pair<double, double>>>> BRKGA::mutacao_combinada()
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vetor_mutante = vector<vector<vector<pair<double, double>>>>(this->tam_mutante, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    for (int i = 0; i < tam_mutante; i++)
    {
        if (i <= (tam_mutante / 2))
            vetor_mutante[i] = criar_individuo();
        else
        {
            vetor_mutante[i] = criar_individuo();
            individuo_decodificado = decodificar_individuo(vetor_mutante[i]);
            individuo_decodificado = trocar_qualificacao_turno(individuo_decodificado);
            vetor_mutante[i] = codificador(individuo_decodificado);
        }
    }
    return vetor_mutante;
}

vector<vector<pair<double, double>>> BRKGA::cruzamento(vector<vector<pair<double, double>>> &pai_elite, vector<vector<pair<double, double>>> &pai_nao_elite)
{
   
    vector<vector<pair<double, double>>> resultado;
    resultado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    double moeda;
    

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            moeda = randzero();
            if (moeda <= prob_cruzamento)
            {
                resultado[i][j].first = pai_elite[i][j].first;
                resultado[i][j].second = pai_elite[i][j].second;
            }
            else
            {
                resultado[i][j].first = pai_nao_elite[i][j].first;
                resultado[i][j].second = pai_nao_elite[i][j].second;
            } 
        }
    }
  
    return resultado;
}

vector<vector<pair<double, double>>> BRKGA::cruzamento_ponto_corte(vector<vector<pair<double, double>>> &pai_elite, vector<vector<pair<double, double>>> &pai_nao_elite)
{
    
  
    int corte;
 
    vector<vector<pair<double, double>>> resultado;

  
    
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        corte = mt.randInt() % 7;
        if (corte >= 3)
        {
              
            
            for (int j = corte + 1 ; j < 7; j++)
            {
                pai_elite[i][j].first = pai_nao_elite[i][j].first;
                pai_elite[i][j].second = pai_nao_elite[i][j].second;
               
            }
        }
        else
        {
            for (int j = 0; j <= corte; j++)
            {
                pai_elite[i][j].first = pai_nao_elite[i][j].first;
                pai_elite[i][j].second = pai_nao_elite[i][j].second;
               
            }
        }
      
    }
    
    return pai_elite;
    
}

vector<vector<vector<pair<double, double>>>> BRKGA::cruzamento_ponto_corte_2(vector<vector<pair<double, double>>> &pai_elite, vector<vector<pair<double, double>>> &pai_nao_elite)
{
 
    resultado_cruzamento = vector<vector<vector<pair<double, double>>>>(2, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    
    int corte;
  
    vector<vector<pair<double, double>>> resultado;
    resultado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    resultado_cruzamento[0] = pai_elite;
    resultado_cruzamento[1] = pai_nao_elite;    //resultado_2 = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    
    corte = mt.randInt() % 7;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
            
            for (int j = corte; j < 7; j++)
            {
                resultado_cruzamento[0][i][j].first = pai_nao_elite[i][j].first;
                resultado_cruzamento[0][i][j].second = pai_nao_elite[i][j].second;
                resultado_cruzamento[1][i][j].first =  pai_elite[i][j].first;
                resultado_cruzamento[1][i][j].second = pai_elite[i][j].second;
            }   
    }
    return resultado_cruzamento;
}


vector<vector<pair<double, double>>> BRKGA::ciclo_evolutivo()
{
    int count;
    int k;
    int aux = 0;
    int divisao = 50;
    int razao = 1;
    int sol_alocacoes = 0;
    int custo;
    int custo_aux;
    int melhor_custo;
    int sem_melhora = 1;
    int melhor_individuo;

    temp.start();
    clock_t t;
    t = clock();
    clock_t tempo;
    double resultado_tempo;
    vector<vector<pair<string, string>>> sol_decodificada;
    populacao_elitista = vector<vector<vector<pair<double, double>>>>(this->tam_elitista, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    populacao_nao_elitista = vector<vector<vector<pair<double, double>>>>(this->tam_nao_elitista, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    sol_decodificada = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<vector<pair<string, string>>> sol_decodificada_ponderada;
    sol_decodificada_ponderada = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<vector<pair<string, string>>> sol_decodificada_normal;
    sol_decodificada_normal = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    matriz_compatibilidades();
    int i;
    int elite, nao_elite;
    
    porcentagem_contratos();
    melhora_populacao_inicial();

    individuo_decodificado = decodificar_individuo_2(populacao[0]);
    melhor_individuo = funcao_avaliacao(individuo_decodificado);
    cout << "Melhor individuo Inicial" << melhor_individuo << endl;
    dividir_populacao(populacao);

    for ( i = 0; i < num_geracoes; i++)
    {
             mutacao_2_vetor(sem_melhora,i); 

            for (int j = 0; j < populacao_elitista.size(); j++)
            {
            populacao[j] = populacao_elitista[j];
            }

        
            for (int h = tam_elitista; h < (tam_elitista + tam_mutante); h++)
            {
                populacao[h] = vetor_mutante[h - tam_elitista];
            }

        
            for (int k = (tam_elitista + tam_mutante); k < tam_populacao; k++)
            {
                elite = (rand() % (tam_elitista));
                nao_elite = (rand() % (tam_nao_elitista));
                populacao[k] = cruzamento_ponto_corte(populacao_elitista[elite], populacao_nao_elitista[nao_elite]);
           }

        dividir_populacao(populacao);
        individuo_decodificado = decodificar_individuo_2(populacao_elitista[0]);
        melhor_custo = funcao_avaliacao(individuo_decodificado);
        cout << "Melhor solucao = " << melhor_custo << ' ' << "Sem melhora:" << sem_melhora << ' ' << i << endl;
        if (melhor_custo < melhor_individuo)
        {
            melhor_individuo = melhor_custo;
            sem_melhora = 1;
        }
        else
        {
            sem_melhora++;
        }
        

        tempo = clock();
        resultado_tempo = ((tempo - t)/(double)CLOCKS_PER_SEC);
         
        
    if(resultado_tempo >= 47.00){
        cout << "Acabou o tempo" << endl;
        cout << "Tempo" << ' ' << resultado_tempo << endl;
        temp.stop();
        return populacao_elitista[0];
        }
    }
}

void BRKGA::porcentagem_contratos()
{
    porcentagens = vector<double>(cn->vetor_contratos.size());
    for (int i = 0; i < cn->vetor_contratos.size(); i++)
    {
        porcentagens[i] = 1 - (((cn->vetor_contratos[i].min_max_1[0].first + cn->vetor_contratos[i].min_max_1[0].second) / 2) / (cn->num_semanas * 7));
    }
}

int BRKGA::multi_start(int &iteracoes_max)
{
    int num_iteracoes = 0;
    int melhor_individuo = 999999;
    int melhor_aux = 0;
    int sol_alocacoes = 0;
    int custo = 0;
    vector<vector<pair<string, string>>> sol_decodificada_ponderada;
    vector<vector<pair<double, double>>> solucao_final;
    vector<vector<pair<string, string>>> solucao_parcial_decodificada;
    while (num_iteracoes < iteracoes_max)
    {
        cout << "Entrou while" << endl;
        vector<vector<pair<double, double>>> solucao_parcial = ciclo_evolutivo_ponderado();
        //vector<vector<pair<double, double>>> solucao_parcial = ciclo_evolutivo();

        num_iteracoes++;
        cout << "Terminou ciclo evolutivo" << endl;
        solucao_parcial_decodificada = decodificar_individuo_ponderado_2(solucao_parcial);
        //solucao_parcial_decodificada = decodificar_individuo_2(solucao_parcial);

        melhor_aux = funcao_avaliacao(solucao_parcial_decodificada);
        if (melhor_aux < melhor_individuo)
        {
            melhor_individuo = melhor_aux;
            solucao_final = solucao_parcial;
        }
        cout << "Melhor individuo: " << melhor_individuo << endl;
    }
    cout << "Fim while" << endl;

    sol_decodificada_ponderada = decodificar_individuo_ponderado_2_sol(solucao_final, sol_alocacoes);
    //sol_decodificada_ponderada = decodificar_individuo_2_sol(solucao_final, sol_alocacoes);
    cout << "Chegou aqui!!!!" << endl;
    print_decodificado(sol_decodificada_ponderada);
    custo = funcao_avaliacao_teste(sol_decodificada_ponderada);
    cout << "Custo" << custo << endl;
    custo = funcao_avaliacao(sol_decodificada_ponderada);
    cout << "Custo_avaliacao_sem_print" << custo << endl;
    ofstream arquivo;
    arquivo.open(this->nomeArq);
    string s;

    string semana;

    for (int i = 0; i < nomeArq.length(); i++)
    {
        if (isdigit(nomeArq[i]))
        {
            semana = nomeArq[i];
        }
    }
    cout << "semana:" << ' ' << semana << endl;

    while (s != "fim")
    {
        arquivo << "SOLUTION" << endl;
        arquivo << semana << ' ';
        arquivo << cn->cenario << endl;
        arquivo << endl;
        arquivo << "ASSIGNMENTS = " << sol_alocacoes << endl;
        for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
        {
            for (int j = 0; j < 7; j++)
            {
                if (sol_decodificada_ponderada[i][j].first != "None")
                {
                    arquivo << cn->vetor_qualificacao[i].nurse << ' ';
                    if (j == 0)
                    {
                        arquivo << "Mon" << ' ';
                    }
                    if (j == 1)
                    {
                        arquivo << "Tue" << ' ';
                    }
                    if (j == 2)
                    {
                        arquivo << "Wed" << ' ';
                    }
                    if (j == 3)
                    {
                        arquivo << "Thu" << ' ';
                    }
                    if (j == 4)
                    {
                        arquivo << "Fri" << ' ';
                    }
                    if (j == 5)
                    {
                        arquivo << "Sat" << ' ';
                    }
                    if (j == 6)
                    {
                        arquivo << "Sun" << ' ';
                    }
                    arquivo << sol_decodificada_ponderada[i][j].first << ' ';
                    arquivo << sol_decodificada_ponderada[i][j].second << endl;
                }
            }
        }
        arquivo.close();
        s = "fim";
    }


    ofstream arquivo_custom;
    arquivo_custom.open(this->custom_out);
    dias_trabalhados = vector<int>(cn->vetor_qualificacao.size());
    string p;
    int aux, aux1;

    while (p != "fim")
    {
        arquivo_custom << cn->vetor_qualificacao.size() << endl;
        for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
        {
            dias_trabalhados[i] = 0;
            for (int j = 0; j < 7; j++)
            {
                if (sol_decodificada_ponderada[i][j].first != string("None"))
                {
                    dias_trabalhados[i]++;
                }
            }

            for (int h = 0; h < cn->vetor_contratos.size(); h++)
            {
                if (cn->vetor_contratos[h].contrato == cn->vetor_qualificacao[i].contract)
                {
                    if (this->semana_1 == 0)
                    {
                        aux = (cus->dias_trabalhados[i].first - dias_trabalhados[i]);
                        if (aux < 0)
                            aux = 0;

                        aux1 = (cus->dias_trabalhados[i].second - dias_trabalhados[i]);
                        if (aux1 < 0)
                            aux1 = 0;

                        arquivo_custom << aux << ' ';
                        arquivo_custom << aux1 << endl;
                    }
                    if (this->semana_1 == 1)
                    {
                        aux = (cn->vetor_contratos[h].min_max[0].first - dias_trabalhados[i]);
                        aux1 = (cn->vetor_contratos[h].min_max[0].second - dias_trabalhados[i]);
                        arquivo_custom << aux << ' ';
                        arquivo_custom << aux1 << endl;
                    }
                }
            }
        }
        arquivo_custom.close();
        p = "fim";
    }

    return melhor_individuo;
}

vector<vector<pair<double, double>>> BRKGA::ciclo_evolutivo_ponderado()
{
    Timer temp;
    temp.start();
    clock_t t;
    t = clock();
    clock_t tempo;
    int count;
    int k;
    int aux = 0;
    int divisao = 50;
    int razao = 1;
    int sol_alocacoes = 0;
    int custo;
    int custo_aux;
    int melhor_custo;
    int sem_melhora = 1;
    int melhor_individuo;
    vector<vector<pair<string, string>>> sol_decodificada;
    populacao_elitista = vector<vector<vector<pair<double, double>>>>(this->tam_elitista, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    populacao_nao_elitista = vector<vector<vector<pair<double, double>>>>(this->tam_nao_elitista, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    sol_decodificada = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<vector<pair<string, string>>> sol_decodificada_ponderada;
    sol_decodificada_ponderada = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<vector<pair<string, string>>> sol_decodificada_normal;
    sol_decodificada_normal = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    matriz_compatibilidades();
    int i;

    int elite, nao_elite;
    gerar_populacao();
    
    porcentagem_contratos();
    cout << "aqui" << endl;
    double resultado_tempo;
    dividir_populacao_2(populacao);
    individuo_decodificado = decodificar_individuo_ponderado_2(populacao[0]);
    melhor_individuo = funcao_avaliacao(individuo_decodificado);
    cout << "Melhor individuo Inicial" << melhor_individuo << endl;
    

    for (i = 0; i < num_geracoes; i++)
    {

            
             mutacao_2_vetor(sem_melhora,i);
            for (int j = 0; j < populacao_elitista.size(); j++)
            {
            populacao[j] = populacao_elitista[j];
            }

        
            for (int h = tam_elitista; h < (tam_elitista + tam_mutante); h++)
            {
                populacao[h] = vetor_mutante[h - tam_elitista];
            }

            for (int k = (tam_elitista + tam_mutante); k < tam_populacao; k++)
            {
                elite = (rand() % (tam_elitista));
                nao_elite = (rand() % (tam_nao_elitista));
                populacao[k] = cruzamento_ponto_corte(populacao_elitista[elite], populacao_nao_elitista[nao_elite]);
           }
  
        dividir_populacao_2(populacao);
        individuo_decodificado = decodificar_individuo_ponderado_2(populacao_elitista[0]);
        melhor_custo = funcao_avaliacao(individuo_decodificado);
        cout << "Melhor solucao = " << melhor_custo << ' ' << "Sem melhora:" << sem_melhora << ' ' << i << endl;
        if (melhor_custo < melhor_individuo)
        {
            melhor_individuo = melhor_custo;
            sem_melhora = 1;
        }
        else
        {
            sem_melhora++;
        }
       

        tempo = clock();
        resultado_tempo = ((tempo - t)/(double)CLOCKS_PER_SEC);
    
    if(resultado_tempo >= 370.00){ // alterar o tempo conforme o número de enfermeiros das instâncias
        cout << "Acabou o tempo" << endl;
        cout << "Tempo" << ' ' << resultado_tempo << endl;
        temp.stop();
        for(int i = 0; i < tam_populacao; i++){
            individuo_decodificado = decodificar_individuo_ponderado_2(populacao[i]);
            custo = funcao_avaliacao_teste(individuo_decodificado);
            cout << "Custo" << ' ' << custo << endl;
        }
        return populacao_elitista[0];
        }
    }   
}

vector<vector<pair<string, string>>> BRKGA::trocar_turno(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    string aux;
    int fit, fit_linha;
    fit = funcao_avaliacao(individuo_decodificado);
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {

            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "None")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_none[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_none.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }

                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "Earçy")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_early[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_early.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }

                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "Day")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_day[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_day.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }

                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "Late")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_late[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_late.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }

                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "Night")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_night[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_night.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }
            }

            if (j != 0)
            {

                if (individuo_decodificado[i][j - 1].first == "None")
                {

                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "None")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_none[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_none.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j - 1].first == "Early")
                {

                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "Earçy")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_early[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_early.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j - 1].first == "Day")
                {

                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "Day")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_day[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_day.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j - 1].first == "Late")
                {

                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "Late")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_late[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_late.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j - 1].first == "Night")
                {

                    aux = individuo_decodificado[i][j].first;
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (individuo_decodificado[i][j].first == "Night")
                            h++;
                        else
                        {
                            individuo_decodificado[i][j].first = cn->vetor_night[h];
                            fit_linha = funcao_avaliacao(individuo_decodificado);
                            if (fit_linha < fit)
                            {
                                fit = fit_linha;
                                h = cn->vetor_night.size();
                            }
                            else
                            {
                                individuo_decodificado[i][j].first = aux;
                            }
                        }
                    }
                }
            }
        }
    }

    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::trocar_qualificacao(vector<vector<pair<string, string>>> &individuo_decodificado)
{

    vector<vector<pair<string, string>>> sol_linha;
    vector<int> sorteio;
    int posicao;

    string aux;
    int fitness, fitness_linha;
    bool melhora;

    for (int i = 0; i < 7; i++)
    {
        sorteio = vector<int>(cn->vetor_qualificacao.size());
        for (int t = 0; t < cn->vetor_qualificacao.size(); t++) // preenche o vetor de sorteio
        {
            sorteio[t] = t;
        }
        for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
        {
            if (j == 0)
            {
                fitness = funcao_avaliacao_minimo(individuo_decodificado, i); 
            }

            
            if (sorteio.size() != 0)
            {
                posicao = (mt.randInt() % sorteio.size());


                for (int h = 0; h < cn->vetor_qualificacao[sorteio[posicao]].num_funcoes; h++)
                {
                    cout << cn->vetor_qualificacao[sorteio[posicao]].funcoes[h] << ' ';
                    
                    if (cn->vetor_qualificacao[sorteio[posicao]].funcoes[h] != individuo_decodificado[sorteio[posicao]][i].second)
                    {

                        aux = individuo_decodificado[sorteio[posicao]][i].second;
                        individuo_decodificado[sorteio[posicao]][i].second = cn->vetor_qualificacao[sorteio[posicao]].funcoes[h]; // faz a troca de qualificação
                        fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                        if (fitness_linha >= fitness)
                        { // troca deve ser desfeita
                            individuo_decodificado[sorteio[posicao]][i].second = aux;
                        }

                        if (fitness_linha < fitness)
                        {
                            fitness = fitness_linha;
                        }
                    }
                }
                sorteio.erase(sorteio.begin() + posicao);
            }
        }
        cout << endl;
    }

    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::trocar_qualificacao_turno(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    bool contem;
    vector<vector<pair<string, string>>> sol_linha;
    pair<string, string> aux1; // guardar o turno e qualificacao que uma dada posicao da matriz continha
    vector<int> sorteio;
    int posicao;

    string aux;
    int fitness, fitness_linha;
    bool melhora;

    for (int i = 0; i < 7; i++)
    {
        sorteio = vector<int>(cn->vetor_qualificacao.size());
        for (int t = 0; t < cn->vetor_qualificacao.size(); t++) // preenche o vetor de sorteio
        {
            sorteio[t] = t;
        }
        for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
        {
            contem = false;
            if (j == 0)
            {
                fitness = funcao_avaliacao_minimo(individuo_decodificado, i); // avalia o fitness da coluna
            }

            if (sorteio.size() != 0)
            {
                posicao = (mt.randInt() % sorteio.size());
                aux1.first = individuo_decodificado[sorteio[posicao]][i].first;

                if (i == 0)
                {
                    if (ht->matriz_borda_s[sorteio[posicao]][1] == "None")
                    {
                        for (int p = 0; p < cn->vetor_none.size(); p++) // verifica aquela posicao contem um turno que respeta a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_none[p])
                                contem = true;
                        }
                        if (contem == false)
                            aux1.first = cn->vetor_none[0];

                        for (int t = 1; t < cn->vetor_none.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_none[0];
                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_none[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            if (fitness_linha > fitness) 
                            {                            
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                                
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }

                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_none.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                    }

                    if (ht->matriz_borda_s[sorteio[posicao]][1] == "Early")
                    {
                        for (int p = 0; p < cn->vetor_early.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_early[p])
                                contem = true;
                        }
                        if (contem == false)
                            aux1.first = "Early";
                        for (int t = 1; t < cn->vetor_early.size(); t++) // t=0 sempre possui none como valor
                        {
                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_early[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                                
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }

                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_early.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                    }
                    if (ht->matriz_borda_s[sorteio[posicao]][1] == "Day")
                    {
                        for (int p = 0; p < cn->vetor_day.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_day[p])
                                contem = true;
                        }
                        if (contem == false)
                            aux1.first = "Day";
                        for (int t = 1; t < cn->vetor_day.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_day[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                                
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }

                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_day.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                    }
                    if (ht->matriz_borda_s[sorteio[posicao]][1] == "Late")
                    {
                        for (int p = 0; p < cn->vetor_late.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_late[p])
                                contem = true;
                        }
                        if (contem == false)
                            aux1.first = "Late";
                        for (int t = 1; t < cn->vetor_late.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_late[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                                
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }

                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_late.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                    }
                    if (ht->matriz_borda_s[sorteio[posicao]][1] == "Night")
                    {
                        for (int p = 0; p < cn->vetor_night.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_night[p])
                                contem = true;
                        }
                        if (contem == false)
                            aux1.first = "Night";
                        for (int t = 1; t < cn->vetor_night.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_night[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                                
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }

                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_night.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                    }
                }
                if (i != 0)
                {
                    //cout << "turno no individuo anteriro: " << individuo_decodificado[sorteio[posicao]][i - 1].first << ' ' << sorteio[posicao] << endl;
                    if (individuo_decodificado[sorteio[posicao]][i - 1].first == "None")
                    {
                        for (int p = 0; p < cn->vetor_none.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_none[p])
                                contem = true;
                        }
                        if (contem == false)
                            aux1.first = cn->vetor_none[0];
                        for (int t = 1; t < cn->vetor_none.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_none[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                               
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }
                             if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_early.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                    }
                    if (individuo_decodificado[sorteio[posicao]][i - 1].first == "Early")
                    {
                        for (int p = 0; p < cn->vetor_early.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_early[p])
                                contem = true;
                        }
                        if (contem == false)
                            
                            aux1.first = "Early";
                        for (int t = 1; t < cn->vetor_early.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_early[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                                
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }
                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_early.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                    }
                    if (individuo_decodificado[sorteio[posicao]][i - 1].first == "Day")
                    {
                        for (int p = 0; p < cn->vetor_day.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_day[p])
                                contem = true;
                        }
                        if (contem == false)
                            aux1.first = "Day";
                        for (int t = 1; t < cn->vetor_day.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_day[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                                
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }
                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_day.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                    
                    }
                    if (individuo_decodificado[sorteio[posicao]][i - 1].first == "Late")
                    {
                        for (int p = 0; p < cn->vetor_late.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_late[p])
                                contem = true;
                        }
                        if (contem == false)
                            
                            aux1.first = "Late";
                        for (int t = 1; t < cn->vetor_late.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_late[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                                
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }
                         
                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_late.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                        
                    }
                    if (individuo_decodificado[sorteio[posicao]][i - 1].first == "Night")
                    {
                        for (int p = 0; p < cn->vetor_night.size(); p++) // verifica aquela posicao contem um turno que respetai a sucessao de turnos
                        {
                            if (aux1.first == cn->vetor_night[p])
                                contem = true;
                        }
                        if (contem == false)
                            
                            aux1.first = "Night";
                        for (int t = 1; t < cn->vetor_night.size(); t++) // t=0 sempre possui none como valor
                        {

                            individuo_decodificado[sorteio[posicao]][i].first = cn->vetor_night[t]; // faz troca de turnos
                            fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                            
                            if (fitness_linha > fitness) // um requerimento minimo ja estava sendo cumprido e mudar o turno piorou a soluçao
                            {                            // troca deve ser desfeita
                                individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                            }

                            if (fitness_linha == fitness) // mudar o turno nao ocasiounou melhora nem piorar, entao as qualificações devem ser trocadas
                            {
                                individuo_decodificado = trocar_qualificacao_enfermeiro(individuo_decodificado, sorteio, posicao, fitness_linha, i);
                                fitness_linha = funcao_avaliacao_minimo(individuo_decodificado, i);
                               
                                if (fitness_linha == fitness)
                                { // se mesmo depois da troca de qualificao a solucao nao melhorou
                                    individuo_decodificado[sorteio[posicao]][i].first = aux1.first;
                                }
                            }
                          
                            if (fitness_linha < fitness) // trocar o turno melhorou a solucao, entao a qualificação ja estava correta
                            {
                                fitness = fitness_linha;
                                t = cn->vetor_night.size(); // se melhorou com a troca de turnos ou qualificacao, nao é mais necessario repetir o procedimento.
                            }
                        }
                        
                    }
                }
                sorteio.erase(sorteio.begin() + posicao);
            }
        }
    }

    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::troca_qualificao_turno_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    bool melhora = false;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<string> vetor_qualificacao;
    vector<string> vetor_turnos = vector<string>(5);
    vetor_turnos = cn->vetor_none;

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }
    random_shuffle(sorteio.begin(), sorteio.end());
    sorteio_aux = sorteio;

   
    while (sorteio.size() != 0)
    {
        melhora = false;
      

        i = sorteio[0].first;
        j = sorteio[0].second;
        
        fitness = funcao_avaliacao_minimo_completa(individuo_decodificado);
        atribuicoes.first = individuo_decodificado[i][j].first;
        atribuicoes.second = individuo_decodificado[i][j].second;

        for (int t = 1; t < cn->vetor_none.size(); t++)
        {
            random_shuffle(vetor_turnos.begin(), vetor_turnos.end());
            individuo_decodificado[i][j].first = vetor_turnos[t];
            for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
            {
                vetor_qualificacao = cn->vetor_qualificacao[i].funcoes;
                random_shuffle(vetor_qualificacao.begin(), vetor_qualificacao.end());
                individuo_decodificado[i][j].second = vetor_qualificacao[h];
                fitness_linha = funcao_avaliacao_minimo_completa(individuo_decodificado);
                if (fitness_linha < fitness)
                {
                    fitness = fitness_linha;
                    melhora = true;
                    break;
                }
            }
            if (melhora == true)
            {
                break;
            }
        }

        if (melhora == false)
        {
            individuo_decodificado[i][j].first = atribuicoes.first;
            individuo_decodificado[i][j].second = atribuicoes.second;
        }

        sorteio.erase(sorteio.begin());
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::trocar_qualificacao_enfermeiro(vector<vector<pair<string, string>>> &individuo_decodificado, vector<int> sorteio, int posicao, int fitness, int i)
{
    int fitness_linha;
    string aux;
    for (int h = 0; h < cn->vetor_qualificacao[sorteio[posicao]].num_funcoes; h++)
    {

        if (cn->vetor_qualificacao[sorteio[posicao]].funcoes[h] != individuo_decodificado[sorteio[posicao]][i].second)
        {

            aux = individuo_decodificado[sorteio[posicao]][i].second;
            individuo_decodificado[sorteio[posicao]][i].second = cn->vetor_qualificacao[sorteio[posicao]].funcoes[h]; // faz a troca de qualificação
            
            fitness_linha = funcao_avaliacao_minimo_completa(individuo_decodificado);
            
            if (fitness_linha >= fitness)
            { // troca deve ser desfeita
                individuo_decodificado[sorteio[posicao]][i].second = aux;
            }

            if (fitness_linha < fitness)
            {
                fitness = fitness_linha;
                break;
            }
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<double, double>>> BRKGA::codificador(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    vector<vector<pair<double, double>>> codificado;
    codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    double aux, aux1, aux2;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_none.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_none.size();
                        }
                    }
                }
                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_early.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_early.size();
                        }
                    }
                }
                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_day.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_day.size();
                        }
                    }
                }
                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_late.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_late.size();
                        }
                    }
                }
                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_night.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_night.size();
                        }
                    }
                }
            }

            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_none.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_none.size();
                        }
                    }
                }
                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_early.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_early.size();
                        }
                    }
                }
                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_day.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_day.size();
                        }
                    }
                }
                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_late.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_late.size();
                        }
                    }
                }
                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            aux1 = h * 1.0;
                            aux2 = cn->vetor_night.size() * 1.0;
                            aux = aux1 / aux2;
                            codificado[i][j].first = aux;
                            h = cn->vetor_night.size();
                        }
                    }
                }
            }

            for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++) // codificar a qualificacao
            {
                if (individuo_decodificado[i][j].second == cn->vetor_qualificacao[i].funcoes[h])
                {
                    aux1 = h * 1.0;
                    aux2 = cn->vetor_qualificacao[i].num_funcoes * 1.0;
                    aux = (aux1) / (aux2);
                    codificado[i][j].second = aux;
                    h = cn->vetor_qualificacao[i].num_funcoes;
                }
            }
        }
    }
    return codificado;
}

vector<vector<pair<double, double>>> BRKGA::codificador_2(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    vector<vector<pair<double, double>>> codificado;
    codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    double aux, aux1, aux2;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            for (int h = 0; h < cn->vetor_none.size(); h++)
            {
                if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                {
                    aux1 = h * 1.0;
                    aux2 = cn->vetor_none.size() * 1.0;
                    codificado[i][j].first = (aux1 / aux2);
                    h = cn->vetor_none.size();
                }
            }

            for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++) // codificar a qualificacao
            {
                if (individuo_decodificado[i][j].second == cn->vetor_qualificacao[i].funcoes[h])
                {
                    aux1 = h * 1.0;
                    aux2 = cn->vetor_qualificacao[i].num_funcoes * 1.0;
                    aux = (aux1) / (aux2);
                    codificado[i][j].second = aux;
                    h = cn->vetor_qualificacao[i].num_funcoes;
                }
            }
        }
    }
    return codificado;
}

vector<vector<pair<double, double>>> BRKGA::codificador_ponderado(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    double razao;
    vector<vector<pair<double, double>>> codificado;
    codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    double aux, aux1, aux2;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (j == 0)
            {
                for (int t = 0; t < cn->vetor_contratos.size(); t++)
                {
                    if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                    {
                        if (ht->matriz_borda_s[i][1] == "None")
                        {
                            for (int p = 0; p < cn->vetor_none.size(); p++)
                            {
                                if (cn->vetor_none[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao));
                                        
                                    }
                                    p = cn->vetor_none.size();
                                }
                            }
                        }

                        if (ht->matriz_borda_s[i][1] == "Early")
                        {
                            for (int p = 0; p < cn->vetor_early.size(); p++)
                            {
                                if (cn->vetor_early[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_early.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao));
                                    }
                                    p = cn->vetor_early.size();
                                }
                            }
                        }
                        if (ht->matriz_borda_s[i][1] == "Day")
                        {
                            for (int p = 0; p < cn->vetor_day.size(); p++)
                            {
                                if (cn->vetor_day[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_day.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao)); 
                                    }
                                    p = cn->vetor_day.size();
                                }
                            }
                        }
                        if (ht->matriz_borda_s[i][1] == "Late")
                        {
                            for (int p = 0; p < cn->vetor_late.size(); p++)
                            {
                                if (cn->vetor_late[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_late.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao));
                                    }
                                    p = cn->vetor_late.size();
                                }
                            }
                        }
                        if (ht->matriz_borda_s[i][1] == "Night")
                        {
                            for (int p = 0; p < cn->vetor_night.size(); p++)
                            {
                                if (cn->vetor_night[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_night.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao));
                                    }
                                    p = cn->vetor_night.size();
                                }
                            }
                        }
                    }
                }
            }

            if (j != 0)
            {
                for (int t = 0; t < cn->vetor_contratos.size(); t++)
                {
                    if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                    {
                        if (individuo_decodificado[i][j - 1].first == "None")
                        {
                            for (int p = 0; p < cn->vetor_none.size(); p++)
                            {
                                if (cn->vetor_none[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao));
                                    }
                                    p = cn->vetor_none.size();
                                }
                            }
                        }

                        if (individuo_decodificado[i][j - 1].first == "Early")
                        {
                            for (int p = 0; p < cn->vetor_early.size(); p++)
                            {
                                if (cn->vetor_early[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_early.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao));
                                    }
                                    p = cn->vetor_early.size();
                                }
                            }
                        }
                        if (individuo_decodificado[i][j - 1].first == "Day")
                        {
                            for (int p = 0; p < cn->vetor_day.size(); p++)
                            {
                                if (cn->vetor_day[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_day.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao)); 
                                    }
                                    p = cn->vetor_day.size();
                                }
                            }
                        }
                        if (individuo_decodificado[i][j - 1].first == "Late")
                        {
                            for (int p = 0; p < cn->vetor_late.size(); p++)
                            {
                                if (cn->vetor_late[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_late.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao));  
                                    }
                                    p = cn->vetor_late.size();
                                }
                            }
                        }
                        if (individuo_decodificado[i][j - 1].first == "Night")
                        {
                            for (int p = 0; p < cn->vetor_night.size(); p++)
                            {
                                if (cn->vetor_night[p] == individuo_decodificado[i][j].first)
                                {
                                    if (individuo_decodificado[i][j].first == "None")
                                    {
                                        codificado[i][j].first = porcentagens[t];
                                    }
                                    else
                                    {
                                        razao = ((1 - porcentagens[t]) / (cn->vetor_night.size() - 1));
                                        codificado[i][j].first = (porcentagens[t] + (p * razao));
                                    }
                                    p = cn->vetor_night.size();
                                }
                            }
                        }
                    }
                }
            }

            for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++) // codificar a qualificacao
            {
                if (individuo_decodificado[i][j].second == cn->vetor_qualificacao[i].funcoes[h])
                {
                    aux1 = h * 1.0;
                    aux2 = cn->vetor_qualificacao[i].num_funcoes * 1.0;
                    aux = (aux1) / (aux2);
                    codificado[i][j].second = aux;
                    h = cn->vetor_qualificacao[i].num_funcoes;
                }
            }
        }
    }
    return codificado;
}

vector<vector<pair<double, double>>> BRKGA::codificador_ponderado_2(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    double razao;
    vector<vector<pair<double, double>>> codificado;
    codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    double aux, aux1, aux2;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                {
                    for (int p = 0; p < cn->vetor_none.size(); p++)
                    {
                        if (cn->vetor_none[p] == individuo_decodificado[i][j].first)
                        {
                            if (individuo_decodificado[i][j].first == "None")
                            {
                                codificado[i][j].first = porcentagens[t];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                                codificado[i][j].first = (porcentagens[t] + (p * razao));
                                
                            }
                            break;
                        }
                    }
                    //break;
                }
            }
            for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++) // codificar a qualificacao
            {
                if (individuo_decodificado[i][j].second == cn->vetor_qualificacao[i].funcoes[h])
                {
                    aux1 = h * 1.0;
                    aux2 = cn->vetor_qualificacao[i].num_funcoes * 1.0;
                    aux = (aux1) / (aux2);
                    codificado[i][j].second = aux;
                    h = cn->vetor_qualificacao[i].num_funcoes;
                }
            }
        }
    }
    return codificado;
}

void BRKGA::melhora_populacao_inicial()
{
    int aux = 0;
    individuo_codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<vector<vector<pair<double, double>>>> solucoes;
    populacao = vector<vector<vector<pair<double, double>>>>(this->tam_populacao, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    solucoes = vector<vector<vector<pair<double, double>>>>(tam_elitista, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    for (int i = 0; i < tam_elitista; i++)
    {
        solucoes[i] = criar_individuo();
        individuo_decodificado = decodificar_individuo_2(solucoes[i]);

        for (int j = 0; j < 3; j++)
        {

            if (j == 0)
            {
                individuo_decodificado = troca_qualificao_turno_aleatorio(individuo_decodificado);
            }

        }
        solucoes[i] = codificador_2(individuo_decodificado);
       
    }

    vector<pair<int, int>> ordenacao(solucoes.size());
    int count = 0;
    for (int i = 0; i < ordenacao.size(); i++)
    {
        individuo_decodificado = decodificar_individuo_2(solucoes[i]);
        ordenacao[i].first = funcao_avaliacao(individuo_decodificado);
        if (ordenacao[i].first < 1000)
        {
            count++;
        }
        ordenacao[i].second = i;
    }
    cout << "Fitness populacao elitista inicial" << endl;
    sort(ordenacao.begin(), ordenacao.end());
    for (int i = 0; i < ordenacao.size(); i++)
    {
        cout << ordenacao[i].first << endl;
    }

    for (int i = 0; i < tam_elitista; i++)
    {
        populacao[i] = solucoes[ordenacao[i].second];
    }
}

void BRKGA::melhora_populacao_inicial_ponderado()
{
    int aux = 0;
    individuo_codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<vector<vector<pair<double, double>>>> solucoes;
    populacao = vector<vector<vector<pair<double, double>>>>(this->tam_populacao, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    solucoes = vector<vector<vector<pair<double, double>>>>(tam_populacao, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    for (int i = 0; i < tam_populacao; i++)
    {

        solucoes[i] = criar_individuo();
        individuo_decodificado = decodificar_individuo_ponderado_2(solucoes[i]);

        for (int j = 0; j < 3; j++)
        {

            if (j == 0)
            {
                individuo_decodificado = troca_qualificao_turno_aleatorio(individuo_decodificado);
                individuo_decodificado = trocar_qualificacao_turno(individuo_decodificado);
            }
        solucoes[i] = codificador_ponderado_2(individuo_decodificado);
    }

    vector<pair<int, int>> ordenacao(solucoes.size());
    int count = 0;
    for (int i = 0; i < ordenacao.size(); i++)
    {
        individuo_decodificado = decodificar_individuo_ponderado_2(solucoes[i]);
        ordenacao[i].first = funcao_avaliacao(individuo_decodificado);
        if (ordenacao[i].first < 1000)
        {
            count++;
        }
        ordenacao[i].second = i;
        
    }
    cout << "Fitness populacao elitista inicial" << endl;
    sort(ordenacao.begin(), ordenacao.end());
    for (int i = 0; i < ordenacao.size(); i++)
    {
        cout << ordenacao[i].first << endl;
    }

    
    for (int i = 0; i < tam_populacao; i++)
    {
        populacao[i] = solucoes[ordenacao[i].second];
    }
}

void BRKGA::troca_populacao()
{
    int aux = 0;
    individuo_codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<vector<vector<pair<double, double>>>> solucoes;
    populacao = vector<vector<vector<pair<double, double>>>>(this->tam_populacao, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    solucoes = vector<vector<vector<pair<double, double>>>>(tam_elitista - tam_populacao, vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7)));
    for (int i = 0; i < solucoes.size(); i++)
    {

        individuo_codificado = criar_individuo();

       

        individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);

        for (int j = 0; j < 3; j++)
        {

            if (j == 0)
            {

                individuo_decodificado = primeira_melhora_2_aleatorio(individuo_decodificado);
                individuo_decodificado = troca_qualificao_turno_aleatorio(individuo_decodificado);
                
            }


            individuo_decodificado = trocar_qualificacao_turno(individuo_decodificado);

            aux = funcao_avaliacao_minimo_completa(individuo_decodificado);
            
            if (aux == 0)
            {
                j = 3;
            }
        }

        individuo_decodificado = primeira_melhora_2_aleatorio(individuo_decodificado);
        solucoes[i] = codificador_ponderado_2(individuo_decodificado);
    }
    cout << "chegou aqui" << endl;
    vector<pair<int, int>> ordenacao(solucoes.size());
    int count = 0;
    for (int i = 0; i < ordenacao.size(); i++)
    {
        individuo_decodificado = decodificar_individuo_ponderado_2(solucoes[i]);
        ordenacao[i].first = funcao_avaliacao(individuo_decodificado);
        if (ordenacao[i].first < 1000)
        {
            count++;
        }
        ordenacao[i].second = i;
    }
    cout << "Fitness populacao elitista inicial" << endl;
    sort(ordenacao.begin(), ordenacao.end());
    for (int i = 0; i < ordenacao.size(); i++)
    {
        cout << ordenacao[i].first << endl;
    }

    
    for (int i = tam_elitista; i < tam_populacao; i++)
    {
        populacao[i] = solucoes[ordenacao[i].second];
    }
}

void BRKGA::melhora_populacao()
{
    int aux = 0;
    individuo_codificado = vector<vector<pair<double, double>>>(cn->vetor_qualificacao.size(), vector<pair<double, double>>(7));
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<vector<vector<pair<double, double>>>> solucoes;
    for (int i = 0; i < populacao.size(); i++)
    {
        
        individuo_decodificado = decodificar_individuo_ponderado_2(populacao[i]);
        for (int j = 0; j < 3; j++)
        {
            if (j == 0)
            {
                individuo_decodificado = primeira_melhora_2_aleatorio(individuo_decodificado);
            }
            individuo_decodificado = trocar_qualificacao_turno(individuo_decodificado);
            aux = funcao_avaliacao(individuo_decodificado);
            if (aux == 0)
            {
                break;
            } 
        }
        individuo_decodificado = primeira_melhora_2_aleatorio(individuo_decodificado);
        populacao[i] = codificador_ponderado_2(individuo_decodificado);
    }
}

// funcao de teste da funcao de avaliacao
int BRKGA::funcao_avaliacao_teste(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int s7 = 0;
    int requerimento_minimo = 0;
    int count_dias_trabalhados = 0;
    int dias_trabalhados = 0;
    int fim_semanas_trabalhados = 0;
    int count_dias_consecutivos = 0;
    int dias_consecutivos = 0;
    int turnos_consecutivos = 0;
    int count_turnos = 0;
    int count_folgas = 0;
    int dias_folgas = 0;
    int aux = 0;
    int i = 0;
    string turno_anterior;
    int custo_total;
    int requisicoes = 0;
    int fim_semana_completo = 0;
    int requerimento_otimo = 0;
    vector<pair<int, int>> turno_funcao;
    turno_funcao = vector<pair<int, int>>(16); // coloca todas as posições do vetor com valor 0
    for (int i = 0; i < 16; i++)
    {
        turno_funcao[i].first = 0;
        turno_funcao[i].second = 0;
    }

    for (i = 0; i < 7; i++)
    {
        for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
        {
            if (individuo_decodificado[j][i].first == "Early")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[0].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[1].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[2].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[3].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Day")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[4].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[5].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[6].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[7].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Late")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[8].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[9].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[10].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[11].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Night")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[12].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[13].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[14].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[15].first++;
                }
            }
            //s5 restricao
            if (i == 5)
            {
                if (cn->vetor_qualificacao[j].contract == "FullTime" && cn->vetor_contratos[0].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        cout << "fim_de_semana_completo" << ' ';
                        cout << j << ' ' << i << endl;
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        cout << "fim_de_semana_completo" << ' ';
                        cout << j << ' ' << i << endl;
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[j].contract == "PartTime" && cn->vetor_contratos[1].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        cout << "fim_de_semana_completo" << ' ';
                        cout << j << ' ' << i << endl;
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        cout << "fim_de_semana_completo" << ' ';
                        cout << j << ' ' << i << endl;
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[j].contract == "HalfTime" && cn->vetor_contratos[2].fim_semana_restricao == 1)
                {

                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        cout << "fim_de_semana_completo" << ' ';
                        cout << j << ' ' << i << endl;
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        cout << "fim_de_semana_completo" << ' ';
                        cout << j << ' ' << i << endl;
                        fim_semana_completo++;
                    }
                }
                if (cn->vetor_qualificacao[j].contract == "20Percent" && cn->vetor_contratos.size() == 4 && cn->vetor_contratos[3].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        cout << "fim_de_semana_completo" << ' ';
                        cout << j << ' ' << i << endl;
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        cout << "fim_de_semana_completo" << ' ';
                        cout << j << ' ' << i << endl;
                        fim_semana_completo++;
                    }
                }
            }
        }
        for (int h = 0; h < 16; h++)
        {
            //fazer aqui o somatorio
            if (turno_funcao[h].first < w0->matriz_requerimentos[h][i].first)
            {
                cout << "requerimento_minimo" << endl;
                requerimento_minimo = (requerimento_minimo + (w0->matriz_requerimentos[h][i].first - turno_funcao[h].first));
            }
            if (turno_funcao[h].first < w0->matriz_requerimentos[h][i].second)
            {
                requerimento_otimo = (requerimento_otimo + (w0->matriz_requerimentos[h][i].second - turno_funcao[h].first));
                cout << "requerimento_otimo" << ' ';
                cout << h << ' ';
                cout << i << ' ';
                cout << w0->matriz_requerimentos[h][i].second << ' ';
                cout << turno_funcao[h].first;
                cout << endl;
            }
            turno_funcao[h].first = 0;
        }
    }

    for (int i = 0; i < w0->matriz_requisicoes.size(); i++)
    {
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first != "None" && w0->matriz_requisicoes[i][1] == 0)
        {
            cout << "requisicoes" << ' ';
            cout << w0->matriz_requisicoes[i][0] << ' ';
            cout << w0->matriz_requisicoes[i][2] << ' ';
            cout << endl;
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Early" && w0->matriz_requisicoes[i][1] == 1)
        {
            cout << "requisicoes" << ' ';
            cout << w0->matriz_requisicoes[i][0] << ' ';
            cout << w0->matriz_requisicoes[i][2] << ' ';
            cout << endl;
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Day" && w0->matriz_requisicoes[i][1] == 2)
        {
            cout << "requisicoes" << ' ';
            cout << w0->matriz_requisicoes[i][0] << ' ';
            cout << w0->matriz_requisicoes[i][2] << ' ';
            cout << endl;
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Late" && w0->matriz_requisicoes[i][1] == 3)
        {
            cout << "requisicoes" << ' ';
            cout << w0->matriz_requisicoes[i][0] << ' ';
            cout << w0->matriz_requisicoes[i][2] << ' ';
            cout << endl;
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Night" && w0->matriz_requisicoes[i][1] == 4)
        {
            cout << "requisicoes" << ' ';
            cout << w0->matriz_requisicoes[i][0] << ' ';
            cout << w0->matriz_requisicoes[i][2] << ' ';
            cout << endl;
            requisicoes++;
        }
    }

    for (int h = 0; h < custom_out.size(); h++)
    {
        if (isdigit(custom_out[h]))
        {
            aux = custom_out[h] - '0';
            aux = cn->num_semanas - aux;
        }
    }

    // restricoes s2 e s3 ---- s2-peso 30 resolvida ---- s3 peso 30 a ser resolvida
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (individuo_decodificado[i][j].first == "None")
            {
                count_folgas++;
                if (j == 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][4] < cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + ht->matriz_borda_fracas[i][4];
                            }
                            if (ht->matriz_borda_fracas[i][4] >= cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + cn->vetor_contratos[h].min_max[2].second;
                            }
                        }
                    }

                    // faz a subtração dos turnos que restarem da semana anterior
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].first)
                            {
                                count_dias_consecutivos = ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_dias_consecutivos > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos < cn->vetor_contratos[h].min_max[1].first)
                            {
                                dias_consecutivos = dias_consecutivos + (cn->vetor_contratos[h].min_max[1].first - count_dias_consecutivos);
                            }
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                            count_dias_consecutivos = 0;
                            break;
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                                break;
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_consecutivos++;
                if (j == 0)
                {
                    count_folgas = ht->matriz_borda_fracas[i][4];
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_folgas > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_folgas < cn->vetor_contratos[h].min_max[2].first)
                            {
                                dias_folgas = dias_folgas + (cn->vetor_contratos[h].min_max[2].first - count_folgas);
                            }
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                            }
                            count_folgas = 0;
                            
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_trabalhados++;
                if (j == 0)
                {
                    count_turnos++;
                    if (ht->matriz_borda_s[i][1] == individuo_decodificado[i][j].first)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + ht->matriz_borda_fracas[i][2];
                                }
                                if (ht->matriz_borda_fracas[i][2] >= cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + cn->vetor_sucessoes[h].x.second;
                                }
                            }
                        }
                    }

                    if (ht->matriz_borda_s[i][1] != individuo_decodificado[i][j].first)
                    {
                        if (ht->matriz_borda_fracas[i][2] > 0)
                        {
                            for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                            {
                                if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                                {
                                    if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                    {
                                        turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                    }
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first == turno_anterior)
                {
                    count_turnos++;
                    if (j == 6)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == turno_anterior)
                            {
                                if (count_turnos > cn->vetor_sucessoes[h].x.second)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first != turno_anterior && j != 0)
                {
                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {
                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                    count_turnos = 1;
                }
                turno_anterior = individuo_decodificado[i][j].first;
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                if (j == 0)
                {
                    if (ht->matriz_borda_fracas[i][2] > 0)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                }
                            }
                        }
                    }
                }
                if (count_turnos > 0)
                {

                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {

                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                }

                count_turnos = 0;
                turno_anterior = individuo_decodificado[i][j].first;
            }
        }
        count_dias_consecutivos = 0;
        count_folgas = 0;
        count_turnos = 0;
        turno_anterior.erase();

        //s6
        if (semana_1 == 0)
        {
            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }

                    // int divisao_first = floor(double(cus->dias_trabalhados[i].first / aux));
                    int divisao_first = floor(double(cus->dias_trabalhados[i].first / aux));
                    int divisao_second = ceil(double(cus->dias_trabalhados[i].second / aux));

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_trabalhados = dias_trabalhados + (divisao_first - count_dias_trabalhados);
                        
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_trabalhados = dias_trabalhados + (count_dias_trabalhados - divisao_second);
                        
                    }
                    break;
                }
            }
        }

        if (semana_1 == 1)
        {

            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }
                    int divisao_first = floor(double(cn->vetor_contratos[t].min_max[0].first / cn->num_semanas));
                    
                    int divisao_second = ceil(double(cn->vetor_contratos[t].min_max[0].second / cn->num_semanas));

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_trabalhados = dias_trabalhados + (divisao_first - count_dias_trabalhados);
                        
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_trabalhados = dias_trabalhados + (count_dias_trabalhados - divisao_second);
                        
                    }
                    break;
                }
            }
        }
        count_dias_trabalhados = 0;
    }
    int consecutivos;
    int sucessoes = 0;
    bool contem;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            contem = false;
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }

            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }
        }
    }
    sucessoes = (sucessoes * 60);
    requerimento_minimo = (requerimento_minimo * 60);
    requerimento_otimo = (requerimento_otimo * 30);
    requisicoes = (requisicoes * 10);
    fim_semana_completo = (fim_semana_completo * 30);
    dias_consecutivos = (dias_consecutivos * 30);
    dias_folgas = (dias_folgas * 30);
    turnos_consecutivos = (turnos_consecutivos * 15);
    consecutivos = turnos_consecutivos + dias_consecutivos;
    if (aux <= cn->num_semanas/2)
    {
        s7 = (s7 * 30);
        dias_trabalhados = (dias_trabalhados * 20);
    }
    else
    {
        s7 = (s7 * 30);
        dias_trabalhados = (dias_trabalhados * 20);
    }

    cout << endl;
    cout << "S7: " << s7 << endl;
    cout << "sucessoes proibidas: " << sucessoes << endl;
    cout << "requerimento_minimo:" << requerimento_minimo << endl;
    cout << "requerimento_otimo:" << requerimento_otimo << endl;
    cout << "requisicoes:" << requisicoes << endl;
    cout << "fim_semana_completo:" << fim_semana_completo << endl;
    cout << "dias consecutivos:" << dias_consecutivos << endl;
    cout << "dias folgas:" << dias_folgas << endl;
    cout << "turnos consecutivos:" << turnos_consecutivos << endl;
    cout << "consecutivos: " << consecutivos << endl;
    cout << "dias_trabalhados: " << dias_trabalhados << endl;

    custo_total = requerimento_minimo + requerimento_otimo + requisicoes + fim_semana_completo + dias_consecutivos + dias_folgas + turnos_consecutivos + dias_trabalhados + sucessoes;
    return custo_total;
}

vector<vector<pair<string, string>>> BRKGA::decodificar_individuo_ponderado_sol(vector<vector<pair<double, double>>> &individuo_codificado, int &sol_alocacoes)
{
    
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int posicao = 0;
    double razao = 0;
    int aux, aux1;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (j == 0)
            {
                for (int t = 0; t < cn->vetor_contratos.size(); t++)
                {
                    if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                    {

                        if (ht->matriz_borda_s[i][1] == "None")
                        {

                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_none[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                                for (int h = 1; h < cn->vetor_none.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_none[h];
                                        break;
                                    }
                                }
                            }
                           
                        }
                        if (ht->matriz_borda_s[i][1] == "Early")
                        {

                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_early[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_early.size() - 1));
                                for (int h = 1; h < cn->vetor_early.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_early[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (ht->matriz_borda_s[i][1] == "Day")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_day[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_day.size() - 1));
                                for (int h = 1; h < cn->vetor_day.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_day[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (ht->matriz_borda_s[i][1] == "Late")
                        {

                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_late[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_late.size() - 1));
                                for (int h = 1; h < cn->vetor_late.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_late[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (ht->matriz_borda_s[i][1] == "Night")
                        {
                           
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_night[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_night.size() - 1));
                                
                                for (int h = 1; h < cn->vetor_night.size(); h++)
                                {

                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_night[h];
                                        break;
                                    }
                                }
                            }
                           
                        }
                    }
                }
            }
            if (j != 0)
            {
                for (int t = 0; t < cn->vetor_contratos.size(); t++)
                {
                    if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                    {
                        if (individuo_decodificado[i][j - 1].first == "None")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_none[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                                for (int h = 1; h < cn->vetor_none.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_none[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (individuo_decodificado[i][j - 1].first == "Early")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_early[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_early.size() - 1));
                                for (int h = 1; h < cn->vetor_early.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_early[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (individuo_decodificado[i][j - 1].first == "Day")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_day[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_day.size() - 1));
                                for (int h = 1; h < cn->vetor_day.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_day[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (individuo_decodificado[i][j - 1].first == "Late")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_late[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_late.size() - 1));
                                for (int h = 1; h < cn->vetor_late.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_late[h];
                                        break;
                                    }
                                }
                            }
                           
                        }
                        if (individuo_decodificado[i][j - 1].first == "Night")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_night[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_night.size() - 1));
                                for (int h = 1; h < cn->vetor_night.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_night[h];

                                        break;
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                sol_alocacoes++;
            }

            posicao = (cn->vetor_qualificacao[i].num_funcoes * individuo_codificado[i][j].second);

            if (posicao == 0)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[0];
            }
            if (posicao == 1)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[1];
            }
            if (posicao == 2)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[2];
            }
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::decodificar_individuo_ponderado(vector<vector<pair<double, double>>> &individuo_codificado)
{

    
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int posicao = 0;
    double razao = 0;
    int aux, aux1;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (j == 0)
            {
                for (int t = 0; t < cn->vetor_contratos.size(); t++)
                {
                    if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                    {

                        if (ht->matriz_borda_s[i][1] == "None")
                        {

                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_none[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                                for (int h = 1; h < cn->vetor_none.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_none[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (ht->matriz_borda_s[i][1] == "Early")
                        {

                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_early[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_early.size() - 1));
                                for (int h = 1; h < cn->vetor_early.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_early[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (ht->matriz_borda_s[i][1] == "Day")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_day[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_day.size() - 1));
                                for (int h = 1; h < cn->vetor_day.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_day[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (ht->matriz_borda_s[i][1] == "Late")
                        {

                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_late[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_late.size() - 1));
                                for (int h = 1; h < cn->vetor_late.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_late[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (ht->matriz_borda_s[i][1] == "Night")
                        {
                          
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_night[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_night.size() - 1));
                                // cout << razao << endl;
                                for (int h = 1; h < cn->vetor_night.size(); h++)
                                {

                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_night[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }

            if (j != 0)
            {
                for (int t = 0; t < cn->vetor_contratos.size(); t++)
                {
                    if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                    {
                        if (individuo_decodificado[i][j - 1].first == "None")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_none[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                                for (int h = 1; h < cn->vetor_none.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_none[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (individuo_decodificado[i][j - 1].first == "Early")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_early[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_early.size() - 1));
                                for (int h = 1; h < cn->vetor_early.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_early[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (individuo_decodificado[i][j - 1].first == "Day")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_day[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_day.size() - 1));
                                for (int h = 1; h < cn->vetor_day.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_day[h];
                                        break;
                                    }
                                }
                            }
                           
                        }
                        if (individuo_decodificado[i][j - 1].first == "Late")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_late[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_late.size() - 1));
                                for (int h = 1; h < cn->vetor_late.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_late[h];
                                        break;
                                    }
                                }
                            }
                            
                        }
                        if (individuo_decodificado[i][j - 1].first == "Night")
                        {
                            if (individuo_codificado[i][j].first <= porcentagens[t])
                            {
                                individuo_decodificado[i][j].first = cn->vetor_night[0];
                            }
                            else
                            {
                                razao = ((1 - porcentagens[t]) / (cn->vetor_night.size() - 1));
                                for (int h = 1; h < cn->vetor_night.size(); h++)
                                {
                                    if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                                    {
                                        individuo_decodificado[i][j].first = cn->vetor_night[h];

                                        break;
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }

            posicao = (cn->vetor_qualificacao[i].num_funcoes * individuo_codificado[i][j].second);

            if (posicao == 0)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[0];
            }
            if (posicao == 1)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[1];
            }
            if (posicao == 2)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[2];
            }
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::decodificar_individuo_ponderado_2(vector<vector<pair<double, double>>> &individuo_codificado)
{

    
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int posicao = 0;
    double razao = 0;
    int aux, aux1;

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                {

                    if (individuo_codificado[i][j].first <= porcentagens[t])
                    {
                        individuo_decodificado[i][j].first = cn->vetor_none[0];
                    }
                    else
                    {
                        razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                        for (int h = 1; h < cn->vetor_none.size(); h++)
                        {
                            if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                            {
                                individuo_decodificado[i][j].first = cn->vetor_none[h];
                                break;
                            }
                        }
                    }
                    break;
                }
            }
            posicao = (cn->vetor_qualificacao[i].num_funcoes * individuo_codificado[i][j].second);

            if (posicao == 0)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[0];
            }
            if (posicao == 1)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[1];
            }
            if (posicao == 2)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[2];
            }
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::decodificar_individuo_ponderado_2_sol(vector<vector<pair<double, double>>> &individuo_codificado, int &sol_alocacoes)
{

   
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int posicao = 0;
    double razao = 0;
    int aux, aux1;

    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_contratos[t].contrato == cn->vetor_qualificacao[i].contract)
                {

                    if (individuo_codificado[i][j].first <= porcentagens[t])
                    {
                        individuo_decodificado[i][j].first = cn->vetor_none[0];
                    }
                    else
                    {
                        razao = ((1 - porcentagens[t]) / (cn->vetor_none.size() - 1));
                        for (int h = 1; h < cn->vetor_none.size(); h++)
                        {
                            if (individuo_codificado[i][j].first <= ((h * razao) + porcentagens[t]))
                            {
                                individuo_decodificado[i][j].first = cn->vetor_none[h];
                                break;
                            }
                        }
                    }
                    break;
                }
            }
            if (individuo_decodificado[i][j].first != "None")
            {
                sol_alocacoes++;
            }
            posicao = (cn->vetor_qualificacao[i].num_funcoes * individuo_codificado[i][j].second);

            if (posicao == 0)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[0];
            }
            if (posicao == 1)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[1];
            }
            if (posicao == 2)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[2];
            }
        }
    }
    return individuo_decodificado;
}

int BRKGA::funcao_avaliacao_minimo(vector<vector<pair<string, string>>> &individuo_decodificado, int coluna)
{
    int custo = 0;
    int requerimento_minimo = 0;
    int requerimento_otimo = 0;
    vector<int> turno_funcao;
    turno_funcao = vector<int>(16);
    for (int i = 0; i < 16; i++)
    {
        turno_funcao[i] = 0;
    }

    for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
    {
        if (individuo_decodificado[j][coluna].first == "Early")
        {
            if (individuo_decodificado[j][coluna].second == "HeadNurse")
            {
                turno_funcao[0]++;
            }
            if (individuo_decodificado[j][coluna].second == "Nurse")
            {
                turno_funcao[1]++;
            }
            if (individuo_decodificado[j][coluna].second == "Caretaker")
            {
                turno_funcao[2]++;
            }
            if (individuo_decodificado[j][coluna].second == "Trainee")
            {
                turno_funcao[3]++;
            }
        }

        if (individuo_decodificado[j][coluna].first == "Day")
        {
            if (individuo_decodificado[j][coluna].second == "HeadNurse")
            {
                turno_funcao[4]++;
            }
            if (individuo_decodificado[j][coluna].second == "Nurse")
            {
                turno_funcao[5]++;
            }
            if (individuo_decodificado[j][coluna].second == "Caretaker")
            {
                turno_funcao[6]++;
            }
            if (individuo_decodificado[j][coluna].second == "Trainee")
            {
                turno_funcao[7]++;
            }
        }

        if (individuo_decodificado[j][coluna].first == "Late")
        {
            if (individuo_decodificado[j][coluna].second == "HeadNurse")
            {
                turno_funcao[8]++;
            }
            if (individuo_decodificado[j][coluna].second == "Nurse")
            {
                turno_funcao[9]++;
            }
            if (individuo_decodificado[j][coluna].second == "Caretaker")
            {
                turno_funcao[10]++;
            }
            if (individuo_decodificado[j][coluna].second == "Trainee")
            {
                turno_funcao[11]++;
            }
        }
        if (individuo_decodificado[j][coluna].first == "Night")
        {
            if (individuo_decodificado[j][coluna].second == "HeadNurse")
            {
                turno_funcao[12]++;
            }
            if (individuo_decodificado[j][coluna].second == "Nurse")
            {
                turno_funcao[13]++;
            }
            if (individuo_decodificado[j][coluna].second == "Caretaker")
            {
                turno_funcao[14]++;
            }
            if (individuo_decodificado[j][coluna].second == "Trainee")
            {
                turno_funcao[15]++;
            }
        }
    }

    for (int h = 0; h < 16; h++)
    {
        //fazer aqui o somatorio
        if (turno_funcao[h] < w0->matriz_requerimentos[h][coluna].first)
        {
            requerimento_minimo = (requerimento_minimo + (w0->matriz_requerimentos[h][coluna].first - turno_funcao[h]));
        }
        if (turno_funcao[h] < w0->matriz_requerimentos[h][coluna].second)
        {
            requerimento_otimo = (requerimento_otimo + (w0->matriz_requerimentos[h][coluna].second - turno_funcao[h]));
        }
        turno_funcao[h] = 0;
    }

    custo = (requerimento_minimo * 50) + (requerimento_otimo * 25);
    return custo;
}

int BRKGA::funcao_avaliacao_minimo_completa(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    custo_colunas = vector<int>(7);
    int custo = 0;
    int requerimento_minimo = 0;
    int requerimento_otimo = 0;
    vector<int> turno_funcao;
    turno_funcao = vector<int>(16);
    for (int i = 0; i < 16; i++)
    {
        turno_funcao[i] = 0;
    }

    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
        {
            if (individuo_decodificado[j][i].first == "Early")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[0]++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[1]++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[2]++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[3]++;
                }
            }

            {
                if (individuo_decodificado[j][i].first == "Day")
                {
                    if (individuo_decodificado[j][i].second == "HeadNurse")
                    {
                        turno_funcao[4]++;
                    }
                    if (individuo_decodificado[j][i].second == "Nurse")
                    {
                        turno_funcao[5]++;
                    }
                    if (individuo_decodificado[j][i].second == "Caretaker")
                    {
                        turno_funcao[6]++;
                    }
                    if (individuo_decodificado[j][i].second == "Trainee")
                    {
                        turno_funcao[7]++;
                    }
                }

                if (individuo_decodificado[j][i].first == "Late")
                {
                    if (individuo_decodificado[j][i].second == "HeadNurse")
                    {
                        turno_funcao[8]++;
                    }
                    if (individuo_decodificado[j][i].second == "Nurse")
                    {
                        turno_funcao[9]++;
                    }
                    if (individuo_decodificado[j][i].second == "Caretaker")
                    {
                        turno_funcao[10]++;
                    }
                    if (individuo_decodificado[j][i].second == "Trainee")
                    {
                        turno_funcao[11]++;
                    }
                }
                if (individuo_decodificado[j][i].first == "Night")
                {
                    if (individuo_decodificado[j][i].second == "HeadNurse")
                    {
                        turno_funcao[12]++;
                    }
                    if (individuo_decodificado[j][i].second == "Nurse")
                    {
                        turno_funcao[13]++;
                    }
                    if (individuo_decodificado[j][i].second == "Caretaker")
                    {
                        turno_funcao[14]++;
                    }
                    if (individuo_decodificado[j][i].second == "Trainee")
                    {
                        turno_funcao[15]++;
                    }
                }
            }

            
        }
        for (int h = 0; h < 16; h++)
        {
            //fazer aqui o somatorio
            if (turno_funcao[h] < w0->matriz_requerimentos[h][i].first)
            {
                requerimento_minimo = (requerimento_minimo + (w0->matriz_requerimentos[h][i].first - turno_funcao[h]));
            }
            if (turno_funcao[h] < w0->matriz_requerimentos[h][i].second)
            {
                requerimento_otimo = (requerimento_otimo + (w0->matriz_requerimentos[h][i].second - turno_funcao[h]));
            }
            turno_funcao[h] = 0;
        }
        
    }

    custo = (requerimento_minimo * 60) + (requerimento_otimo *30);
    return custo;
}

void BRKGA::funcao_avaliacao_colunas(vector<vector<pair<string, string>>> &individuo_decodificado)
{

    custo_colunas = vector<int>(7);
    int custo = 0;
    int requerimento_minimo = 0;
    int requerimento_otimo = 0;
    vector<int> turno_funcao;
    turno_funcao = vector<int>(16);
    for (int i = 0; i < 16; i++)
    {
        turno_funcao[i] = 0;
    }

    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
        {
            if (individuo_decodificado[j][i].first == "Early")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[0]++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[1]++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[2]++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[3]++;
                }
            }

            {
                if (individuo_decodificado[j][i].first == "Day")
                {
                    if (individuo_decodificado[j][i].second == "HeadNurse")
                    {
                        turno_funcao[4]++;
                    }
                    if (individuo_decodificado[j][i].second == "Nurse")
                    {
                        turno_funcao[5]++;
                    }
                    if (individuo_decodificado[j][i].second == "Caretaker")
                    {
                        turno_funcao[6]++;
                    }
                    if (individuo_decodificado[j][i].second == "Trainee")
                    {
                        turno_funcao[7]++;
                    }
                }

                if (individuo_decodificado[j][i].first == "Late")
                {
                    if (individuo_decodificado[j][i].second == "HeadNurse")
                    {
                        turno_funcao[8]++;
                    }
                    if (individuo_decodificado[j][i].second == "Nurse")
                    {
                        turno_funcao[9]++;
                    }
                    if (individuo_decodificado[j][i].second == "Caretaker")
                    {
                        turno_funcao[10]++;
                    }
                    if (individuo_decodificado[j][i].second == "Trainee")
                    {
                        turno_funcao[11]++;
                    }
                }
                if (individuo_decodificado[j][i].first == "Night")
                {
                    if (individuo_decodificado[j][i].second == "HeadNurse")
                    {
                        turno_funcao[12]++;
                    }
                    if (individuo_decodificado[j][i].second == "Nurse")
                    {
                        turno_funcao[13]++;
                    }
                    if (individuo_decodificado[j][i].second == "Caretaker")
                    {
                        turno_funcao[14]++;
                    }
                    if (individuo_decodificado[j][i].second == "Trainee")
                    {
                        turno_funcao[15]++;
                    }
                }
            }

           
        }
        for (int h = 0; h < 16; h++)
        {
            //fazer aqui o somatorio
            if (turno_funcao[h] < w0->matriz_requerimentos[h][i].first)
            {
                requerimento_minimo = (requerimento_minimo + (w0->matriz_requerimentos[h][i].first - turno_funcao[h]));
            }
            if (turno_funcao[h] < w0->matriz_requerimentos[h][i].second)
            {
                requerimento_otimo = (requerimento_otimo + (w0->matriz_requerimentos[h][i].second - turno_funcao[h]));
            }
            turno_funcao[h] = 0;
        }
         custo_colunas[i] = ((requerimento_minimo * 60) + (requerimento_otimo * 30));
         requerimento_minimo = 0;
         requerimento_otimo = 0;
    }

}

void BRKGA::solucao_gulosa()
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    int aux = 0;
    
    for (int i = cn->vetor_qualificacao.size() - 1; i >= 0; i--)
    {
        for (int j = 0; j < 7; j++)
        {
            
            individuo_decodificado = trocas_2(i, j, individuo_decodificado);
        }
    }

   
    for (int j = 0; j < 3; j++)
    {
       
        if (j == 0)
        {
            individuo_decodificado = primeira_melhora_2(individuo_decodificado);
            
        }
       
        individuo_decodificado = trocar_qualificacao_turno(individuo_decodificado);

        aux = funcao_avaliacao_minimo_completa(individuo_decodificado);
       
        if (aux == 0)
        {
           
            j = 3;
        }
    }
    individuo_decodificado = primeira_melhora_2(individuo_decodificado);
    individuo_decodificado = procura_melhor_vizinho(individuo_decodificado);
    aux = funcao_avaliacao(individuo_decodificado);
    cout << "Custo da solucao = " << aux << endl;

}

void BRKGA::solucao_gulosa_2()
{
    individuo_decodificado = vector<vector<pair<string, string>>>(cn->vetor_qualificacao.size(), vector<pair<string, string>>(7));
    vector<Vetor_qualificacao> auxiliar;
    int i = 0;
    int aux = 0;
    vector<int> contratos;
    contratos = vector<int>(cn->vetor_qualificacao.size());

    for (int h = 0; h < cn->vetor_qualificacao.size(); h++)
    {
        if (cn->vetor_qualificacao[h].contract == "FullTime")
        {
            
            contratos[i] = h;
            i++;
        }
    }

    for (int h = 0; h < cn->vetor_qualificacao.size(); h++)
    {
        if (cn->vetor_qualificacao[h].contract == "PartTime")
        {
            
            contratos[i] = h;
            i++;
        }
    }

    for (int h = 0; h < cn->vetor_qualificacao.size(); h++)
    {
        if (cn->vetor_qualificacao[h].contract == "HalfTime")
        {
            
            contratos[i] = h;
            i++;
        }
    }

    
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            
            individuo_decodificado = trocas_2(contratos[i], j, individuo_decodificado);
        }
    }
    print_decodificado(individuo_decodificado);
    aux = funcao_avaliacao_2(individuo_decodificado);
    cout << "Custo da solucao = " << aux << endl;
}

vector<vector<pair<string, string>>> BRKGA::trocas(int i, int j, vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int aux, fitness;
    fitness = funcao_avaliacao_minimo_completa(individuo_decodificado);
    
    for (int p = 1; p < cn->vetor_none.size(); p++)
    {
        individuo_decodificado[i][j].first = cn->vetor_none[p];
        for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
        {
            individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
            aux = funcao_avaliacao_minimo_completa(individuo_decodificado);

            if (aux < fitness)
            {
                return individuo_decodificado;
            }
        }
    }
    if (aux == fitness)
    {
        individuo_decodificado[i][j].first = "None";
        individuo_decodificado[i][j].second = "None";
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::trocas_2(int i, int j, vector<vector<pair<string, string>>> &individuo_decodificado)
{
    
    int aux, fitness;
    pair<string, string> atribuicoes;
    individuo_decodificado[i][j].first = "None";
    individuo_decodificado[i][j].second = "None";

    fitness = funcao_avaliacao(individuo_decodificado);
    

    for (int p = 0; p < cn->vetor_none.size(); p++)
    {
        individuo_decodificado[i][j].first = cn->vetor_none[p];
        for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
        {
            individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
            aux = funcao_avaliacao(individuo_decodificado);
   
            if (aux <= fitness)
            {
                fitness = aux;
                atribuicoes.first = individuo_decodificado[i][j].first;
                atribuicoes.second = individuo_decodificado[i][j].second;
            }
        }
    }

    individuo_decodificado[i][j].first = atribuicoes.first;
    individuo_decodificado[i][j].second = atribuicoes.second;
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::melhor_melhora(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    for (int i = 0; i < cn->vetor_contratos.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            procura_melhor_vizinho(individuo_decodificado);
        }
    }
}

vector<vector<pair<string, string>>> BRKGA::procura_melhor_vizinho(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    pair<string, string> salvar;
    pair<int, int> posicoes;
    bool melhora = true;
    while (melhora == true)
    {
        count++;
        melhora = false;

  
        for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
        {
            for (int j = 0; j < 7; j++)
            {

                fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
       
                salvar.first = individuo_decodificado[i][j].first;
                salvar.second = individuo_decodificado[i][j].second;
                
                for (int p = 0; p < cn->vetor_none.size(); p++)
                {
                    individuo_decodificado[i][j].first = cn->vetor_none[p];
                    for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
                    {
                        individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
                        fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
      
                        if (fitness_linha < fitness)
                        {
                            fitness = fitness_linha;
                            atribuicoes.first = individuo_decodificado[i][j].first;
                            atribuicoes.second = individuo_decodificado[i][j].second;
                            posicoes.first = i;
                            posicoes.second = j;
                            melhora = true;
                        }
                    }
                }

                individuo_decodificado[i][j].first = salvar.first;
                individuo_decodificado[i][j].second = salvar.second;   
            }
        }

        if (melhora == true)
        {
            
            individuo_decodificado[posicoes.first][posicoes.second].first = atribuicoes.first;
            individuo_decodificado[posicoes.first][posicoes.second].second = atribuicoes.second;
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::primeira_melhora(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    pair<string, string> atribuicoes;


    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {

            fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
        
            atribuicoes.first = individuo_decodificado[i][j].first;
            atribuicoes.second = individuo_decodificado[i][j].second;
            
            for (int p = 0; p < cn->vetor_none.size(); p++)
            {
                individuo_decodificado[i][j].first = cn->vetor_none[p];
                for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
                {
                    individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
                    fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
            
                    if (fitness_linha < fitness)
                    {
               
                        fitness = fitness_linha;
                      
                        atribuicoes.first = individuo_decodificado[i][j].first;
                        atribuicoes.second = individuo_decodificado[i][j].second;
                    }
                }
            }

            individuo_decodificado[i][j].first = atribuicoes.first;
            individuo_decodificado[i][j].second = atribuicoes.second;
        }
    }

    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::primeira_melhora_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int aux_i = 0;
    int aux_j = 0;
    int i, j;
    int fitness, fitness_linha, fitness_aux, aux;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    pair<string, string> atribuicoes;
    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }
    random_shuffle(sorteio.begin(), sorteio.end());
   

    while (sorteio.size() != 0)
    {
        i = sorteio[0].first;
        j = sorteio[0].second;
        fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, sorteio[0].first, sorteio[0].second);
     
        atribuicoes.first = individuo_decodificado[i][j].first;
        atribuicoes.second = individuo_decodificado[i][j].second;
        
        for (int p = 0; p < cn->vetor_none.size(); p++)
        {
            individuo_decodificado[i][j].first = cn->vetor_none[p];
            for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
                fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, sorteio[0].first, sorteio[0].second);
    
                if (fitness_linha < fitness)
                {
                   
                    fitness = fitness_linha;
                   
                    atribuicoes.first = individuo_decodificado[i][j].first;
                    atribuicoes.second = individuo_decodificado[i][j].second;
                }
            }
        }

        individuo_decodificado[i][j].first = atribuicoes.first;
        individuo_decodificado[i][j].second = atribuicoes.second;
        sorteio.erase(sorteio.begin());
    }

    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::primeira_melhora_2_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    bool melhora = true;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }
    sorteio_aux = sorteio;

    while (melhora == true)
    {
        sorteio = sorteio_aux;
        random_shuffle(sorteio.begin(), sorteio.end());
        melhora = false;
        fitness = funcao_avaliacao(individuo_decodificado);
        while (sorteio.size() != 0)
        {
            i = sorteio[0].first;
            j = sorteio[0].second;

            fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
   
            atribuicoes.first = individuo_decodificado[i][j].first;
            atribuicoes.second = individuo_decodificado[i][j].second;
            
            for (int p = 0; p < cn->vetor_none.size(); p++)
            {
                individuo_decodificado[i][j].first = cn->vetor_none[p];
                for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
                {
                    individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
                    fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
                
                    if (fitness_linha < fitness)
                    {
                        
                        fitness = fitness_linha;
                        
                        atribuicoes.first = individuo_decodificado[i][j].first;
                        atribuicoes.second = individuo_decodificado[i][j].second;
                        melhora = true;
                        h = cn->vetor_qualificacao[i].num_funcoes;
                        p = cn->vetor_none.size();
                    }
                }
            }
            individuo_decodificado[i][j].first = atribuicoes.first;
            individuo_decodificado[i][j].second = atribuicoes.second;
            if (melhora == true)
            {
                break;
            }
            sorteio.erase(sorteio.begin());
        }
    }

    return individuo_decodificado;
}

bool BRKGA::primeira_melhora_2_aleatorio_vnd(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux, fitnessAnterior, fitnessFinal;
    vector<string> vetor_none = vector<string>(5);
    vector<string> vetor_qualificacoes;
    vetor_none = cn->vetor_none;
    pair<string, string> atribuicoes;
    bool melhora = true;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }

    random_shuffle(sorteio.begin(), sorteio.end());
    melhora = false;
    for (int z = 0; z < sorteio.size(); z++)
    {
        i = sorteio[z].first;
        j = sorteio[z].second;

        fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);

        atribuicoes.first = individuo_decodificado[i][j].first;
        atribuicoes.second = individuo_decodificado[i][j].second;
        
        random_shuffle(vetor_none.begin(), vetor_none.end());
        for (int p = 0; p < cn->vetor_none.size(); p++)
        {
            individuo_decodificado[i][j].first = vetor_none[p];
            vetor_qualificacoes = vector<string>(cn->vetor_qualificacao[i].num_funcoes);
            for(int t = 0; t < cn->vetor_qualificacao[i].num_funcoes; t++){
                vetor_qualificacoes[t] = cn->vetor_qualificacao[i].funcoes[t];
            }
            random_shuffle(vetor_qualificacoes.begin(), vetor_qualificacoes.end());
            for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
            {
                individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
                individuo_decodificado[i][j].second = vetor_qualificacoes[h];
                fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);

                if (fitness_linha < fitness)
                {
                 

                    atribuicoes.first = individuo_decodificado[i][j].first;
                    atribuicoes.second = individuo_decodificado[i][j].second;
                    melhora = true;
                
                    random_shuffle(sorteio.begin(), sorteio.end());
                    z = 0;
               
                }
                else
                {
                    individuo_decodificado[i][j].first = atribuicoes.first;
                    individuo_decodificado[i][j].second = atribuicoes.second;
                }
            }
        }
    }
    return melhora;
}

bool BRKGA::primeira_melhora_2_aleatorio_vnd_thalles(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux, fitnessAnterior, fitnessFinal;
    vector<string> vetor_none = vector<string>(5);
    vector<string> vetor_qualificacoes;
    vetor_none = cn->vetor_none;
    pair<string, string> atribuicoes;
    bool melhora = true;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }

   
    fitnessAnterior = funcao_avaliacao(individuo_decodificado);

    random_shuffle(sorteio.begin(), sorteio.end());
    melhora = false;
    for (int z = 0; z < sorteio.size(); z++)
    {
        i = sorteio[0].first;
        j = sorteio[0].second;

        fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
        atribuicoes.first = individuo_decodificado[i][j].first;
        atribuicoes.second = individuo_decodificado[i][j].second;
        random_shuffle(vetor_none.begin(), vetor_none.end());
        for (int p = 0; p < cn->vetor_none.size(); p++)
        {
            individuo_decodificado[i][j].first = vetor_none[p];
            vetor_qualificacoes = cn->vetor_qualificacao[i].funcoes;
            random_shuffle(vetor_qualificacoes.begin(), vetor_qualificacoes.end());
            for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
            {
                
                individuo_decodificado[i][j].second = vetor_qualificacoes[h];
                fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);

                if (fitness_linha < fitness)
                {
                    atribuicoes.first = individuo_decodificado[i][j].first;
                    atribuicoes.second = individuo_decodificado[i][j].second;
                    melhora = true;
                    h = cn->vetor_qualificacao[i].num_funcoes;
                    p = cn->vetor_none.size();

                    return melhora;
                }
            }  
        }
        if(melhora == false){
            individuo_decodificado[i][j].first = atribuicoes.first;
            individuo_decodificado[i][j].second = atribuicoes.second;
        }
      
   sorteio.erase(sorteio.begin());
    }

}

bool BRKGA::primeira_melhora_2_turno_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    bool melhora = true;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }
    sorteio_aux = sorteio;



    random_shuffle(sorteio.begin(), sorteio.end());
    melhora = false;
    while (sorteio.size() != 0)
    {
        i = sorteio[0].first;
        j = sorteio[0].second;

     
        atribuicoes.first = individuo_decodificado[i][j].first;
        fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
        
        for (int p = 0; p < cn->vetor_none.size(); p++)
        {
            individuo_decodificado[i][j].first = cn->vetor_none[p];
            fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
            if (fitness_linha < fitness)
            {
                
                fitness = fitness_linha;
                
                atribuicoes.first = individuo_decodificado[i][j].first;

                melhora = true;

                p = cn->vetor_none.size();
            }
        }
        individuo_decodificado[i][j].first = atribuicoes.first;

        if (melhora == true)
        {
            break;
        }
        sorteio.erase(sorteio.begin());
    }

    return melhora;
}

bool BRKGA::primeira_melhora_2_qualificacao_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    bool melhora = true;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }
    sorteio_aux = sorteio;



    random_shuffle(sorteio.begin(), sorteio.end());
    melhora = false;
    while (sorteio.size() != 0)
    {
        i = sorteio[0].first;
        j = sorteio[0].second;
        fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);

        atribuicoes.second = individuo_decodificado[i][j].second;

        //cout << "fitness_aux = " << fitness_aux << endl;
        for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
        {
            individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
            fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);

            if (fitness_linha < fitness)
            {
                
                fitness = fitness_linha;
                
                atribuicoes.second = individuo_decodificado[i][j].second;
                melhora = true;
                h = cn->vetor_qualificacao[i].num_funcoes;
            }
        }
        individuo_decodificado[i][j].second = atribuicoes.second;

        if (melhora == true)
        {
            break;
        }
        sorteio.erase(sorteio.begin());
    }

    return melhora;
}

bool BRKGA::primeira_melhora_2_turno_blocos_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    bool melhora = false;
    int fitness, fitness_linha, inicio_corte, limite, tam_corte, fitnessAnterior, fitnessFinal;
    pair<string, string> atribuicoes;
    vector<vector<vector<pair<string, string>>>> individuos_antes_troca;
   
    individuos_antes_troca = vector<vector<vector<pair<string, string>>>>(2, vector<vector<pair<string, string>>>(1, vector<pair<string, string>>(7)));
   
    vector<Vetor_compatibilidade> compatibilade_aux = vector<Vetor_compatibilidade>(vetor_compatibilidade.size());
    compatibilade_aux = vetor_compatibilidade;

    random_shuffle(compatibilade_aux.begin(), compatibilade_aux.end());

    fitnessAnterior = funcao_avaliacao(individuo_decodificado);

    inicio_corte = (mt.randInt() % 7);
    limite = 7 - inicio_corte;
    tam_corte = (mt.randInt() % limite);
    tam_corte++;
    
    for (int i = 0; i < vetor_compatibilidade.size(); i++)
    {
        
       fitness = funcao_avaliacao(individuo_decodificado);
       bool melhora = false;
       
        random_shuffle(compatibilade_aux[i].compatibilidade.begin(), compatibilade_aux[i].compatibilidade.end());

  

        for (int j = 0; j < compatibilade_aux[i].compatibilidade.size() - 1; j++)
        {

            for (int k = j + 1; k < compatibilade_aux[i].compatibilidade.size(); k++)
            {
               

                for (int t = inicio_corte; t < inicio_corte + tam_corte; t++)
                { // salva as informacoes antes do movimento
                    
                    individuos_antes_troca[0][0][t].first = individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][t].first;
                    individuos_antes_troca[0][0][t].second = individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][t].second;
                    individuos_antes_troca[1][0][t].first = individuo_decodificado[compatibilade_aux[i].compatibilidade[k]][t].first;
                    individuos_antes_troca[1][0][t].second = individuo_decodificado[compatibilade_aux[i].compatibilidade[k]][t].second;

                    
                }

                for (int h = inicio_corte; h < inicio_corte + tam_corte; h++) // faz o swap
                {
                    individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][h].first = individuos_antes_troca[1][0][h].first;
                    individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][h].second = individuos_antes_troca[1][0][h].second;
                    individuo_decodificado[compatibilade_aux[i].compatibilidade[k]][h].first = individuos_antes_troca[0][0][h].first;
                    individuo_decodificado[compatibilade_aux[i].compatibilidade[k]][h].second = individuos_antes_troca[0][0][h].second;
                }

                fitness_linha = funcao_avaliacao(individuo_decodificado);
                if (fitness_linha < fitness)
                {
                    melhora = true;
                    fitness = fitness_linha;
                  
                    return melhora;
                }
                else
                { // desfaz o movimento
                    for (int t = inicio_corte; t < inicio_corte + tam_corte; t++)
                    {
                        individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][t].first = individuos_antes_troca[0][0][t].first;
                        individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][t].second = individuos_antes_troca[0][0][t].second;
                        individuo_decodificado[compatibilade_aux[i].compatibilidade[k]][t].first = individuos_antes_troca[1][0][t].first;
                        individuo_decodificado[compatibilade_aux[i].compatibilidade[k]][t].second = individuos_antes_troca[1][0][t].second;

                    }
                }
            }
        }
    }

    fitnessFinal = funcao_avaliacao(individuo_decodificado);

    if (fitnessFinal < fitnessAnterior)
        melhora = true;
    else
        melhora = false;

    return melhora;
}

bool BRKGA::primeira_melhora_2_turno_blocos_aleatorio_2(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int requerimentos;
    bool melhora = false;
    int fitness, fitness_linha, inicio_corte, limite, tam_corte, fitnessAnterior, fitnessFinal;
    pair<string, string> atribuicoes;
    vector<vector<vector<pair<string, string>>>> individuos_antes_troca;
    
    individuos_antes_troca = vector<vector<vector<pair<string, string>>>>(2, vector<vector<pair<string, string>>>(1, vector<pair<string, string>>(7)));
    vector<Vetor_compatibilidade> compatibilade_aux = vector<Vetor_compatibilidade>(vetor_compatibilidade.size());
    compatibilade_aux = vetor_compatibilidade;

    random_shuffle(compatibilade_aux.begin(), compatibilade_aux.end());

    inicio_corte = (mt.randInt() % 7);
    limite = 7 - inicio_corte;
    tam_corte = (mt.randInt() % limite);
    tam_corte++;

    
    for (int i = 0; i < vetor_compatibilidade.size(); i++)
    {
  
        
        random_shuffle(compatibilade_aux[i].compatibilidade.begin(), compatibilade_aux[i].compatibilidade.end());

        
         
        for (int j = 0; j < compatibilade_aux[i].compatibilidade.size() - 1; j++)
        {
            fitness = funcao_avaliacao(individuo_decodificado);
            
           for (int t = inicio_corte; t < inicio_corte + tam_corte; t++)
                { // salva as informacoes antes do mobiemtno
                    individuos_antes_troca[0][0][t].first = individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][t].first;
                    individuos_antes_troca[0][0][t].second = individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][t].second;
                    individuos_antes_troca[1][0][t].first = individuo_decodificado[compatibilade_aux[i].compatibilidade[j + 1]][t].first;
                    individuos_antes_troca[1][0][t].second = individuo_decodificado[compatibilade_aux[i].compatibilidade[j +1]][t].second; 
                }
            
             for (int h = inicio_corte; h < inicio_corte + tam_corte; h++) // faz o swap
                {
                    individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][h].first = individuos_antes_troca[1][0][h].first;
                    individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][h].second = individuos_antes_troca[1][0][h].second;
                    individuo_decodificado[compatibilade_aux[i].compatibilidade[j + 1]][h].first = individuos_antes_troca[0][0][h].first;
                    individuo_decodificado[compatibilade_aux[i].compatibilidade[j + 1]][h].second = individuos_antes_troca[0][0][h].second;
                }

                fitness_linha = funcao_avaliacao(individuo_decodificado);
                if (fitness_linha < fitness)
                {
                    melhora = true;
                    random_shuffle(compatibilade_aux.begin(), compatibilade_aux.end());
                    inicio_corte = (mt.randInt() % 7);
                    limite = 7 - inicio_corte;
                    tam_corte = (mt.randInt() % limite);
                    tam_corte++;
                    
                    i = -1;
                    j = -1;
                    break;
                    

                    
                }
                else
                { // desfaz o movimento
                    for (int t = inicio_corte; t < inicio_corte + tam_corte; t++)
                    {
                        
                        individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][t].first = individuos_antes_troca[0][0][t].first;
                        individuo_decodificado[compatibilade_aux[i].compatibilidade[j]][t].second = individuos_antes_troca[0][0][t].second;
                        individuo_decodificado[compatibilade_aux[i].compatibilidade[j + 1]][t].first = individuos_antes_troca[1][0][t].first;
                        individuo_decodificado[compatibilade_aux[i].compatibilidade[j + 1]][t].second = individuos_antes_troca[1][0][t].second;

                    }
                }
                   if(melhora){
           
                }
        }
     
    }


    return melhora;
}



bool BRKGA::primeira_melhora_2_troca_pares(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int fitness, fitness_linha, fitness_aux, aux;
    int posicao, posicao_2;
    int aux_i = 0;
    int aux_j = 0;
    int i, j;
    pair<string, string> atribuicoes;
    pair<string, string> atribuicoes_inicial;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 6);
    bool melhora;

    for (int i = 0; i < cn->vetor_qualificacao.size() * 6; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 5)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }

    random_shuffle(sorteio.begin(), sorteio.end());
    melhora = false;

    while (sorteio.size() != 0)
    {
        i = sorteio[0].first;
        j = sorteio[0].second;
        fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
        atribuicoes_inicial.first = individuo_decodificado[i][j].first;
        atribuicoes_inicial.second = individuo_decodificado[i][j].second;

        individuo_decodificado[i][j].first = individuo_decodificado[i][j + 1].first;
        individuo_decodificado[i][j].second = individuo_decodificado[i][j + 1].second;

        individuo_decodificado[i][j + 1].first = atribuicoes_inicial.first;
        individuo_decodificado[i][j + 1].second = atribuicoes_inicial.second;

        fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
        if (fitness_linha < fitness)
        {
            melhora = true;
            fitness = fitness_linha;
           
            return melhora;
        }
        else
        {
            individuo_decodificado[i][j + 1].first = individuo_decodificado[i][j].first;
            individuo_decodificado[i][j + 1].second = individuo_decodificado[i][j].second;
            individuo_decodificado[i][j].first = atribuicoes_inicial.first;
            individuo_decodificado[i][j].second = atribuicoes_inicial.second;
        }
        sorteio.erase(sorteio.begin());
    }
    return melhora;
}

vector<vector<pair<string, string>>> BRKGA::folga_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    bool melhora = true;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }
    sorteio_aux = sorteio;

    while (melhora == false)
    {
        sorteio = sorteio_aux;
        random_shuffle(sorteio.begin(), sorteio.end());
        melhora = false;

        i = sorteio[0].first;
        j = sorteio[0].second;
        fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
        atribuicoes.first = individuo_decodificado[i][j].first;

        

        individuo_decodificado[i][j].first = "None";
        fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
        if (fitness_linha < fitness)
        {
            fitness = fitness_linha;
            
            atribuicoes.first = individuo_decodificado[i][j].first;
            melhora = true;
        }

        individuo_decodificado[i][j].first = atribuicoes.first;

        if (melhora == true)
        {
            break;
        }
        sorteio.erase(sorteio.begin());
    }

    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::turno_aleatorio(vector<vector<pair<string, string>>> &individuo_vnd)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    bool melhora = true;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }
    sorteio_aux = sorteio;

    sorteio = sorteio_aux;
    random_shuffle(sorteio.begin(), sorteio.end());
    melhora = false;
    while (sorteio.size() != 0)
    {
        i = sorteio[0].first;
        j = sorteio[0].second;
        fitness = funcao_avaliacao_linha_coluna(individuo_vnd, i, j);
        atribuicoes.first = individuo_decodificado[i][j].first;

        for (int p = 1; p < cn->vetor_none.size(); p++)
        {
            individuo_decodificado[i][j].first = cn->vetor_none[p];
            fitness_linha = funcao_avaliacao_linha_coluna(individuo_vnd, i, j);
            if (fitness_linha < fitness)
            {
                fitness = fitness_linha;
                
                atribuicoes.first = individuo_decodificado[i][j].first;

                melhora = true;

                p = cn->vetor_none.size();
            }
        }
        individuo_decodificado[i][j].first = atribuicoes.first;

        if (melhora == true)
        {
            return individuo_decodificado;
        }
        sorteio.erase(sorteio.begin());
    }
    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::qualificacao_aleatorio(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int i, j;
    int aux_i = 0;
    int aux_j = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    bool melhora = true;
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    vector<pair<int, int>> sorteio_aux;
    sorteio_aux = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);

    for (int i = 0; i < cn->vetor_qualificacao.size() * 7; i++)
    {
        sorteio[i].first = aux_i;
        sorteio[i].second = aux_j;
        if (aux_j == 6)
        {
            aux_j = 0;
            aux_i++;
        }
        else
        {
            aux_j++;
        }
    }
    sorteio_aux = sorteio;


    sorteio = sorteio_aux;
    random_shuffle(sorteio.begin(), sorteio.end());
    melhora = false;
    while (sorteio.size() != 0)
    {
        i = sorteio[0].first;
        j = sorteio[0].second;
        fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
        atribuicoes.second = individuo_decodificado[i][j].second;

        for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
        {
            individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
            fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
            if (fitness_linha < fitness)
            {
                
                fitness = fitness_linha;
                atribuicoes.second = individuo_decodificado[i][j].second;
                melhora = true;
                h = cn->vetor_qualificacao[i].num_funcoes;
            }
        }
        individuo_decodificado[i][j].second = atribuicoes.second;

        if (melhora == true)
        {
            return individuo_decodificado;
        }
        sorteio.erase(sorteio.begin());
    }

    return individuo_decodificado;
}

vector<vector<pair<double, double>>> BRKGA::VND(vector<vector<pair<double, double>>> &individuo_codificado)
{
    
    int custo;
    bool melhora = false;
    int count = 0;
    int k = 1, r = 2;

    int moeda = mt.randInt() % 2;
    individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);
    
    while (k <= r)
    {
        if (moeda == 0)
        {
            if (k == 1)
            {
                
                melhora = primeira_melhora_2_aleatorio_vnd(individuo_decodificado);
                if (melhora == true)
                {
                    count++;
                }
            }

            if (k == 2)
            {
                melhora = primeira_melhora_2_turno_blocos_aleatorio_2(individuo_decodificado);
            }
        }
        else if (moeda == 1)
        {
            if (k == 2)
            {
                melhora = primeira_melhora_2_aleatorio_vnd(individuo_decodificado);
                if (melhora == true)
                {
                    count++;
                }
            }

            if (k == 1)
            {
                melhora = primeira_melhora_2_turno_blocos_aleatorio_2(individuo_decodificado);
            }
        }

        if (melhora == true)
        {
            count++;
            melhora = false;
            k = 1;
        }
        else
        {
            k++;
        }
    }
 
    individuo_codificado = codificador_ponderado_2(individuo_decodificado);
    return individuo_codificado;
}

vector<vector<pair<double, double>>> BRKGA::VND_inicial(vector<vector<pair<double, double>>> &individuo_codificado)
{
    
    int custo;
    bool melhora;
    int count = 0;
    int k = 1, r = 1;
    individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);
    while (k <= r)
    {
        if (k == 1)
        {
            melhora = primeira_melhora_2_aleatorio_vnd_thalles(individuo_decodificado);
            if (melhora == true)
            {
                count++; 
                k = 1;
            }
            else
            {
                k++;
            }
        }  
    }

    individuo_codificado = codificador_ponderado_2(individuo_decodificado);

    return individuo_codificado;
}

vector<vector<pair<double, double>>> BRKGA::VND_inicial_normal(vector<vector<pair<double, double>>> &individuo_codificado)
{
    int custo;
    bool melhora;
    int count = 0;
    int k = 1, r = 1;
    individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);
    while (k <= r)
    {
        if (k == 1)
        {
            melhora = primeira_melhora_2_aleatorio_vnd(individuo_decodificado);
            if (melhora == true)
            {
                count++;
                k = 1;
            }
            else
            {
                k++;
            }
        }
    }
    individuo_codificado = codificador_ponderado_2(individuo_decodificado);
    return individuo_codificado;
}

vector<vector<pair<string, string>>> BRKGA::h2(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    vector<int> enfermeiros = vector<int>(cn->vetor_qualificacao.size());
    vector<int> enfermeiros_aux = vector<int>(cn->vetor_qualificacao.size());
    int tamanho_requerimento = 0;
    bool contem = false;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        enfermeiros[i] = i;
    }
    enfermeiros_aux = enfermeiros;

    for (int j = 0; j < 7; j++)
    {
        tamanho_requerimento = 0;

        enfermeiros = enfermeiros_aux;
        random_shuffle(enfermeiros.begin(), enfermeiros.end());
        if (tamanho_requerimento < 4)
        {
            if (w0->matriz_requerimentos[0][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 0)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "HeadNurse")
                            {
                                cout << "DEBUG" << endl;
                                individuo_decodificado[enfermeiros[t]][j].first = "Early";
                                individuo_decodificado[enfermeiros[t]][j].second = "HeadNurse";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[1][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 1)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Nurse")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Early";
                                individuo_decodificado[enfermeiros[t]][j].second = "Nurse";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[2][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 2)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Caretaker")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Early";
                                individuo_decodificado[enfermeiros[t]][j].second = "Caretaker";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[3][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 3)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Trainee")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Early";
                                individuo_decodificado[enfermeiros[t]][j].second = "Trainee";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }
        }

        if (tamanho_requerimento >= 4 && tamanho_requerimento < 8)
        {
            if (w0->matriz_requerimentos[4][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 4)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "HeadNurse")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Day";
                                individuo_decodificado[enfermeiros[t]][j].second = "HeadNurse";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[5][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 5)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Nurse")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Day";
                                individuo_decodificado[enfermeiros[t]][j].second = "Nurse";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[6][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 6)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Caretaker")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Day";
                                individuo_decodificado[enfermeiros[t]][j].second = "Caretaker";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[7][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 7)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Trainee")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Day";
                                individuo_decodificado[enfermeiros[t]][j].second = "Trainee";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }
        }

        if (tamanho_requerimento >= 8 && tamanho_requerimento < 12)
        {
            if (w0->matriz_requerimentos[8][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 8)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "HeadNurse")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Late";
                                individuo_decodificado[enfermeiros[t]][j].second = "HeadNurse";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[9][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 9)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Nurse")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Late";
                                individuo_decodificado[enfermeiros[t]][j].second = "Nurse";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[10][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 10)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Caretaker")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Late";
                                individuo_decodificado[enfermeiros[t]][j].second = "Caretaker";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[11][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 11)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Trainee")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Late";
                                individuo_decodificado[enfermeiros[t]][j].second = "Trainee";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }
        }
        if (tamanho_requerimento >= 12 && tamanho_requerimento < 16)
        {
            contem = false;
            if (w0->matriz_requerimentos[12][j].first != 0)
            {
                if (tamanho_requerimento == 12)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "HeadNurse")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Night";
                                individuo_decodificado[enfermeiros[t]][j].second = "HeadNurse";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[13][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 13)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Nurse")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Night";
                                individuo_decodificado[enfermeiros[t]][j].second = "Nurse";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[14][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 14)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Caretaker")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Night";
                                individuo_decodificado[enfermeiros[t]][j].second = "Caretaker";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            tamanho_requerimento++;
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }

            if (w0->matriz_requerimentos[15][j].first != 0)
            {
                contem = false;
                if (tamanho_requerimento == 15)
                {
                    for (int t = 0; t < enfermeiros.size(); t++)
                    {
                        for (int j = 0; j < cn->vetor_qualificacao[enfermeiros[t]].num_funcoes; j++)
                        {
                            if (cn->vetor_qualificacao[enfermeiros[t]].funcoes[j] == "Trainee")
                            {
                                individuo_decodificado[enfermeiros[t]][j].first = "Night";
                                individuo_decodificado[enfermeiros[t]][j].second = "Trainee";
                                enfermeiros.erase(enfermeiros.begin() + enfermeiros[t]);
                                contem = true;
                                break;
                            }
                        }
                        if (contem == true)
                        {
                            break;
                        }
                    }
                }
            }
            else
            {
                tamanho_requerimento++;
            }
        }
    }
    return individuo_decodificado;
}

vector<vector<pair<double, double>>> BRKGA::VND_normal(vector<vector<pair<double, double>>> &individuo_codificado)
{
    //cout << endl;
    int custo;
    bool melhora;
    int count = 0;
    int k = 1, r = 1;
    individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);
    while (k <= r)
    {
        if (k == 1)
        {
            melhora = primeira_melhora_2_turno_blocos_aleatorio(individuo_decodificado);
            if (melhora == true)
            {
                count++;
                k = 1;
            }
            else
            {
                k++;
            }
        }
    }
    individuo_codificado = codificador_ponderado_2(individuo_decodificado);
    return individuo_codificado;
}

vector<vector<pair<double, double>>> BRKGA::VND_print(vector<vector<pair<double, double>>> &individuo_codificado)
{
    int custo;
    bool melhora;
    int count = 0;
    int k = 1, r = 2;
    individuo_decodificado = decodificar_individuo_ponderado_2(individuo_codificado);
    while (k <= r)
    {
        if (k == 1)
        {
            melhora = primeira_melhora_2_turno_blocos_aleatorio(individuo_decodificado);
            k++;
        }

        if (k == 2)
        {
            melhora = primeira_melhora_2_aleatorio_vnd(individuo_decodificado);
            count++;
            custo = funcao_avaliacao(individuo_decodificado);
            cout << "Custo k == 2: " << custo << endl;

            if (melhora == true)
            {
                k = 1;
            }

            else
            {
                k++;
            }
        }
    }
    individuo_codificado = codificador_ponderado_2(individuo_decodificado);
    return individuo_codificado;
}

vector<vector<pair<string, string>>> BRKGA::primeira_melhora_2(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int count = 0;
    int fitness, fitness_linha, fitness_aux, aux;
    pair<string, string> atribuicoes;
    bool melhora = true;

    while (melhora == true)
    {
        melhora = false;
        for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
        {
            for (int j = 0; j < 7; j++)
            {

                fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
                atribuicoes.first = individuo_decodificado[i][j].first;
                atribuicoes.second = individuo_decodificado[i][j].second;
                for (int p = 0; p < cn->vetor_none.size(); p++)
                {
                    individuo_decodificado[i][j].first = cn->vetor_none[p];
                    for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
                    {
                        individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
                        fitness_linha = funcao_avaliacao_linha_coluna(individuo_decodificado, i, j);
                        if (fitness_linha < fitness)
                        {
                            fitness = fitness_linha;
                            
                            atribuicoes.first = individuo_decodificado[i][j].first;
                            atribuicoes.second = individuo_decodificado[i][j].second;
                            melhora = true;
                            h = cn->vetor_qualificacao[i].num_funcoes;
                            p = cn->vetor_none.size();
                        }
                    }
                }
                individuo_decodificado[i][j].first = atribuicoes.first;
                individuo_decodificado[i][j].second = atribuicoes.second;
                if (melhora == true)
                {
                    i = cn->vetor_qualificacao.size();
                    j = 7;
                }
            }
        }
    }

    return individuo_decodificado;
}

vector<vector<pair<string, string>>> BRKGA::trocas_3(int i, int j, vector<vector<pair<string, string>>> &individuo_decodificado, int linha, int coluna)
{
    vector<pair<int, int>> sorteio;
    sorteio = vector<pair<int, int>>(cn->vetor_qualificacao.size() * 7);
    int aux, fitness;
    pair<string, string> atribuicoes;
    individuo_decodificado[i][j].first = "None";
    individuo_decodificado[i][j].second = "None";

    fitness = funcao_avaliacao_linha_coluna(individuo_decodificado, linha, coluna);

    for (int p = 0; p < cn->vetor_none.size(); p++)
    {
        individuo_decodificado[i][j].first = cn->vetor_none[p];
        for (int h = 0; h < cn->vetor_qualificacao[i].num_funcoes; h++)
        {
            individuo_decodificado[i][j].second = cn->vetor_qualificacao[i].funcoes[h];
            aux = funcao_avaliacao_linha_coluna(individuo_decodificado, linha, coluna);
            if (aux <= fitness)
            {
                fitness = aux;
                atribuicoes.first = individuo_decodificado[i][j].first;
                atribuicoes.second = individuo_decodificado[i][j].second;
            }
        }
    }

    individuo_decodificado[i][j].first = atribuicoes.first;
    individuo_decodificado[i][j].second = atribuicoes.second;
    return individuo_decodificado;
}

int BRKGA::funcao_avaliacao_2(vector<vector<pair<string, string>>> &individuo_decodificado)
{
    int s7 = 0;
    int aux = 0;
    int count_dias_trabalhados = 0;
    int dias_global = 0;
    int i = 0;
    int count = 0;
    int turnos_consecutivos = 0;
    int count_folgas = 0;
    int dias_folgas = 0;
    string turno_anterior;
    int count_dias_consecutivos = 0;
    int dias_consecutivos = 0;
    int count_fim_semana = 0;
    int count_turnos = 0;
    int custo_total;
    int requisicoes = 0;
    int fim_semana_completo = 0;
    int requerimento_otimo = 0;
    int requerimento_minimo = 0;
    vector<pair<int, int>> turno_funcao;
    turno_funcao = vector<pair<int, int>>(16); // coloca todas as posições do vetor com valor 0
    for (int i = 0; i < 16; i++)
    {
        turno_funcao[i].first = 0;
        turno_funcao[i].second = 0;
    }

    for (i = 0; i < 7; i++)
    {
        for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
        {
            if (individuo_decodificado[j][i].first == "Early")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[0].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[1].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[2].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[3].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Day")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[4].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[5].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[6].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[7].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Late")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[8].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[9].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[10].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[11].first++;
                }
            }
            if (individuo_decodificado[j][i].first == "Night")
            {
                if (individuo_decodificado[j][i].second == "HeadNurse")
                {
                    turno_funcao[12].first++;
                }
                if (individuo_decodificado[j][i].second == "Nurse")
                {
                    turno_funcao[13].first++;
                }
                if (individuo_decodificado[j][i].second == "Caretaker")
                {
                    turno_funcao[14].first++;
                }
                if (individuo_decodificado[j][i].second == "Trainee")
                {
                    turno_funcao[15].first++;
                }
            }
            if (i == 5)
            {
                if (cn->vetor_qualificacao[j].contract == "FullTime" && cn->vetor_contratos[0].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[j].contract == "PartTime" && cn->vetor_contratos[1].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[j].contract == "HalfTime" && cn->vetor_contratos[2].fim_semana_restricao == 1)
                {

                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
                if (cn->vetor_qualificacao[j].contract == "20Percent" && cn->vetor_contratos.size() == 4 && cn->vetor_contratos[3].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[j][i].first == "None" && individuo_decodificado[j][i + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[j][i].first != "None" && individuo_decodificado[j][i + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
            }
        }

        for (int h = 0; h < 16; h++)
        {
            //fazer aqui o somatorio
            if (turno_funcao[h].first < w0->matriz_requerimentos[h][i].first)
            {
                requerimento_minimo = (requerimento_minimo + (w0->matriz_requerimentos[h][i].first - turno_funcao[h].first));
            }
            if (turno_funcao[h].first < w0->matriz_requerimentos[h][i].second)
            {
                requerimento_otimo = (requerimento_otimo + (w0->matriz_requerimentos[h][i].second - turno_funcao[h].first));
            }
            turno_funcao[h].first = 0;
        }
    }

    for (int i = 0; i < w0->matriz_requisicoes.size(); i++)
    {
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first != "None" && w0->matriz_requisicoes[i][1] == 0)
        {
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Early" && w0->matriz_requisicoes[i][1] == 1)
        {
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Day" && w0->matriz_requisicoes[i][1] == 2)
        {
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Late" && w0->matriz_requisicoes[i][1] == 3)
        {
            requisicoes++;
        }
        if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Night" && w0->matriz_requisicoes[i][1] == 4)
        {
            requisicoes++;
        }
    }

    for (int h = 0; h < custom_out.size(); h++)
    {
        if (isdigit(custom_out[h]))
        {
            aux = custom_out[h] - '0';
            aux = cn->num_semanas - aux;
        }
    }

    // restricoes s2 e s3 ---- s2-peso 30 resolvida ---- s3 peso 30 a ser resolvida
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            if (individuo_decodificado[i][j].first == "None")
            {
                count_folgas++;
                if (j == 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][4] < cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + ht->matriz_borda_fracas[i][4];
                            }
                            if (ht->matriz_borda_fracas[i][4] >= cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + cn->vetor_contratos[h].min_max[2].second;
                            }
                        }
                    }

                    // faz a subtração dos turnos que restarem da semana anterior
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].first)
                            {
                                count_dias_consecutivos = ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_dias_consecutivos > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos < cn->vetor_contratos[h].min_max[1].first)
                            {
                                dias_consecutivos = dias_consecutivos + (cn->vetor_contratos[h].min_max[1].first - count_dias_consecutivos);
                            }
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                            count_dias_consecutivos = 0;
                            break;
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                                break;
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_consecutivos++;
                if (j == 0)
                {
                    count_folgas = ht->matriz_borda_fracas[i][4];
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_folgas > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_folgas < cn->vetor_contratos[h].min_max[2].first)
                            {
                                dias_folgas = dias_folgas + (cn->vetor_contratos[h].min_max[2].first - count_folgas);
                            }
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                            }
                            count_folgas = 0;
                            // break;
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_trabalhados++;
                if (j == 0)
                {
                    count_turnos++;
                    if (ht->matriz_borda_s[i][1] == individuo_decodificado[i][j].first)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + ht->matriz_borda_fracas[i][2];
                                }
                                if (ht->matriz_borda_fracas[i][2] >= cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + cn->vetor_sucessoes[h].x.second;
                                }
                            }
                        }
                    }

                    if (ht->matriz_borda_s[i][1] != individuo_decodificado[i][j].first)
                    {
                        if (ht->matriz_borda_fracas[i][2] > 0)
                        {
                            for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                            {
                                if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                                {
                                    if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                    {
                                        turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                    }
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first == turno_anterior)
                {
                    count_turnos++;
                    if (j == 6)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == turno_anterior)
                            {
                                if (count_turnos > cn->vetor_sucessoes[h].x.second)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first != turno_anterior && j != 0)
                {
                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {
                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                    count_turnos = 1;
                }
                turno_anterior = individuo_decodificado[i][j].first;
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                if (j == 0)
                {
                    if (ht->matriz_borda_fracas[i][2] > 0)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                }
                            }
                        }
                    }
                }
                if (count_turnos > 0)
                {

                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {

                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                }

                count_turnos = 0;
                turno_anterior = individuo_decodificado[i][j].first;
            }
        }
        count_dias_consecutivos = 0;
        count_folgas = 0;
        count_turnos = 0;
        turno_anterior.erase();
        //restrição turnos consecutivos

        //restriçao s6
        if (semana_1 == 0)
        {
            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                     if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }
                    int divisao_first = (cus->dias_trabalhados[i].first / aux);
                    int divisao_second = (cus->dias_trabalhados[i].second / aux);

                    if (divisao_first * aux != cus->dias_trabalhados[i].first)
                        divisao_first++;

                    if (divisao_second * aux != cus->dias_trabalhados[i].second)
                    {
                        divisao_second++;
                    }
                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                    }
                }
            }
        }

        if (semana_1 == 1)
        {
            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                     if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }
                    int divisao_first = (cn->vetor_contratos[t].min_max[0].first / cn->num_semanas);
                    if (divisao_first * cn->num_semanas != cn->vetor_contratos[t].min_max[0].first)
                        divisao_first++;

                    int divisao_second = (cn->vetor_contratos[t].min_max[0].second / cn->num_semanas);
                    if (divisao_second * cn->num_semanas != cn->vetor_contratos[t].min_max[0].second)
                        divisao_second++;

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                    }
                }
            }
        }
        count_dias_trabalhados = 0;
    }

    int sucessoes = 0;
    bool contem;
    for (int i = 0; i < cn->vetor_qualificacao.size(); i++)
    {
        for (int j = 0; j < 7; j++)
        {
            contem = false;
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }

            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }
        }
    }
    sucessoes = (sucessoes * 60);
    requerimento_minimo = (requerimento_minimo * 60);
    requerimento_otimo = (requerimento_otimo * 30);
    requisicoes = (requisicoes * 10);
    fim_semana_completo = (fim_semana_completo * 30);
    dias_consecutivos = (dias_consecutivos * 30);
    dias_folgas = (dias_folgas * 30);
    turnos_consecutivos = (turnos_consecutivos * 15);
    if (aux <= cn->num_semanas/2)
    {
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }
    else{
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }
        


    custo_total = s7 + sucessoes + requerimento_minimo + requerimento_otimo + requisicoes + fim_semana_completo + dias_consecutivos + dias_folgas + turnos_consecutivos + dias_global;
    return custo_total;
}

int BRKGA::funcao_avaliacao_linha(vector<vector<pair<string, string>>> &individuo_decodificado, int linha)
{

    // int j = coluna;
    int s7 = 0;
    int aux = 0;
    int count_dias_trabalhados = 0;
    int dias_global = 0;

    int count = 0;
    int turnos_consecutivos = 0;
    int count_folgas = 0;
    int dias_folgas = 0;
    string turno_anterior;
    int count_dias_consecutivos = 0;
    int dias_consecutivos = 0;
    int count_fim_semana = 0;
    int count_turnos = 0;
    int custo_total;
    int requisicoes = 0;
    int fim_semana_completo = 0;
    int requerimento_otimo = 0;
    int requerimento_minimo = 0;
    vector<int> turno_funcao;
    turno_funcao = vector<int>(16); // coloca todas as posições do vetor com valor 0
   

    for (int h = 0; h < custom_out.size(); h++)
    {
        if (isdigit(custom_out[h]))
        {
            aux = custom_out[h] - '0';
            aux = cn->num_semanas - aux;
        }
    }

    for (int i = 0; i < w0->matriz_requisicoes.size(); i++)
    {
        if (w0->matriz_requisicoes[i][0] == linha)
        {
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first != "None" && w0->matriz_requisicoes[i][1] == 0)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Early" && w0->matriz_requisicoes[i][1] == 1)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Day" && w0->matriz_requisicoes[i][1] == 2)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Late" && w0->matriz_requisicoes[i][1] == 3)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Night" && w0->matriz_requisicoes[i][1] == 4)
            {
                requisicoes++;
            }
        }
    }

    // restricoes s2 e s3 ---- s2-peso 30 resolvida ---- s3 peso 30 a ser resolvida
    for (int i = linha; i <= linha; i++)
    {
        for (int j = 0; j < 7; j++)
        {

            if (j == 5)
            {
                if (cn->vetor_qualificacao[i].contract == "FullTime" && cn->vetor_contratos[0].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[i].contract == "PartTime" && cn->vetor_contratos[1].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[i].contract == "HalfTime" && cn->vetor_contratos[2].fim_semana_restricao == 1)
                {

                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
                if (cn->vetor_qualificacao[i].contract == "20Percent" && cn->vetor_contratos.size() == 4 && cn->vetor_contratos[3].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                count_folgas++;
                if (j == 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][4] < cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + ht->matriz_borda_fracas[i][4];
                            }
                            if (ht->matriz_borda_fracas[i][4] >= cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + cn->vetor_contratos[h].min_max[2].second;
                            }
                        }
                    }

                    // faz a subtração dos turnos que restarem da semana anterior
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].first)
                            {
                                count_dias_consecutivos = ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_dias_consecutivos > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos < cn->vetor_contratos[h].min_max[1].first)
                            {
                                dias_consecutivos = dias_consecutivos + (cn->vetor_contratos[h].min_max[1].first - count_dias_consecutivos);
                            }
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                            count_dias_consecutivos = 0;
                            break;
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                                break;
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_consecutivos++;
                if (j == 0)
                {
                    count_folgas = ht->matriz_borda_fracas[i][4];
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_folgas > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_folgas < cn->vetor_contratos[h].min_max[2].first)
                            {
                                dias_folgas = dias_folgas + (cn->vetor_contratos[h].min_max[2].first - count_folgas);
                            }
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                            }
                            count_folgas = 0;
                            
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_trabalhados++;
                if (j == 0)
                {
                    count_turnos++;
                    if (ht->matriz_borda_s[i][1] == individuo_decodificado[i][j].first)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + ht->matriz_borda_fracas[i][2];
                                }
                                if (ht->matriz_borda_fracas[i][2] >= cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + cn->vetor_sucessoes[h].x.second;
                                }
                            }
                        }
                    }

                    if (ht->matriz_borda_s[i][1] != individuo_decodificado[i][j].first)
                    {
                        if (ht->matriz_borda_fracas[i][2] > 0)
                        {
                            for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                            {
                                if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                                {
                                    if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                    {
                                        turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                    }
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first == turno_anterior)
                {
                    count_turnos++;
                    if (j == 6)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == turno_anterior)
                            {
                                if (count_turnos > cn->vetor_sucessoes[h].x.second)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first != turno_anterior && j != 0)
                {
                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {
                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                    count_turnos = 1;
                }
                turno_anterior = individuo_decodificado[i][j].first;
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                if (j == 0)
                {
                    if (ht->matriz_borda_fracas[i][2] > 0)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                }
                            }
                        }
                    }
                }
                if (count_turnos > 0)
                {

                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {

                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                }

                count_turnos = 0;
                turno_anterior = individuo_decodificado[i][j].first;
            }
        }
        count_dias_consecutivos = 0;
        count_folgas = 0;
        count_turnos = 0;
        turno_anterior.erase();
        //restrição turnos consecutivos

        //restriçao s6

        if (semana_1 == 0)
        {

            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        //s7++;
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }

                    
                    int divisao_first = floor(double(cus->dias_trabalhados[i].first / aux));
                    int divisao_second = ceil(double(cus->dias_trabalhados[i].second / aux) );

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                        
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                        
                    }
                    break;
                }
            }
        }

        if (semana_1 == 1)
        {

            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {

                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        //s7++;
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }
                   
                    int divisao_first = floor(double(cn->vetor_contratos[t].min_max[0].first / cn->num_semanas));
                    int divisao_second = ceil(double(cn->vetor_contratos[t].min_max[0].second / cn->num_semanas) );

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                        
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                        
                    }
                    break;
                }
            }
        }
        count_dias_trabalhados = 0;
    }

    int sucessoes = 0;
    bool contem;
    for (int i = linha; i <= linha; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            contem = false;
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }

            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }
        }
    }


    sucessoes = (sucessoes * 60);
    requerimento_minimo = (requerimento_minimo * 60);
    requerimento_otimo = (requerimento_otimo * 30);
    requisicoes = (requisicoes * 10);
    fim_semana_completo = (fim_semana_completo * 30);
    dias_consecutivos = (dias_consecutivos * 30);
    dias_folgas = (dias_folgas * 30);
    turnos_consecutivos = (turnos_consecutivos * 15);
    if (aux <= cn->num_semanas/2)
    {
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }
    else{
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }

        

    custo_total = s7 + sucessoes + requisicoes + fim_semana_completo + dias_consecutivos + dias_folgas + turnos_consecutivos + dias_global;
    
    return custo_total;
}
//cout << endl;
int BRKGA::funcao_avaliacao_linha_coluna(vector<vector<pair<string, string>>> &individuo_decodificado, int linha, int coluna)
{

    
    int s7 = 0;
    int aux = 0;
    int count_dias_trabalhados = 0;
    int dias_global = 0;

    int count = 0;
    int turnos_consecutivos = 0;
    int count_folgas = 0;
    int dias_folgas = 0;
    string turno_anterior;
    int count_dias_consecutivos = 0;
    int dias_consecutivos = 0;
    int count_fim_semana = 0;
    int count_turnos = 0;
    int custo_total;
    int requisicoes = 0;
    int fim_semana_completo = 0;
    int requerimento_otimo = 0;
    int requerimento_minimo = 0;
    vector<int> turno_funcao;
    turno_funcao = vector<int>(16); // coloca todas as posições do vetor com valor 0
    for (int i = 0; i < 16; i++)
    {
        turno_funcao[i] = 0;
    }

    for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
    {
        if (individuo_decodificado[j][coluna].first == "Early")
        {
            if (individuo_decodificado[j][coluna].second == "HeadNurse")
            {
                turno_funcao[0]++;
            }
            if (individuo_decodificado[j][coluna].second == "Nurse")
            {
                turno_funcao[1]++;
            }
            if (individuo_decodificado[j][coluna].second == "Caretaker")
            {
                turno_funcao[2]++;
            }
            if (individuo_decodificado[j][coluna].second == "Trainee")
            {
                turno_funcao[3]++;
            }
        }

        if (individuo_decodificado[j][coluna].first == "Day")
        {
            if (individuo_decodificado[j][coluna].second == "HeadNurse")
            {
                turno_funcao[4]++;
            }
            if (individuo_decodificado[j][coluna].second == "Nurse")
            {
                turno_funcao[5]++;
            }
            if (individuo_decodificado[j][coluna].second == "Caretaker")
            {
                turno_funcao[6]++;
            }
            if (individuo_decodificado[j][coluna].second == "Trainee")
            {
                turno_funcao[7]++;
            }
        }

        if (individuo_decodificado[j][coluna].first == "Late")
        {
            if (individuo_decodificado[j][coluna].second == "HeadNurse")
            {
                turno_funcao[8]++;
            }
            if (individuo_decodificado[j][coluna].second == "Nurse")
            {
                turno_funcao[9]++;
            }
            if (individuo_decodificado[j][coluna].second == "Caretaker")
            {
                turno_funcao[10]++;
            }
            if (individuo_decodificado[j][coluna].second == "Trainee")
            {
                turno_funcao[11]++;
            }
        }
        if (individuo_decodificado[j][coluna].first == "Night")
        {
            if (individuo_decodificado[j][coluna].second == "HeadNurse")
            {
                turno_funcao[12]++;
            }
            if (individuo_decodificado[j][coluna].second == "Nurse")
            {
                turno_funcao[13]++;
            }
            if (individuo_decodificado[j][coluna].second == "Caretaker")
            {
                turno_funcao[14]++;
            }
            if (individuo_decodificado[j][coluna].second == "Trainee")
            {
                turno_funcao[15]++;
            }
        }
    }

    for (int h = 0; h < 16; h++)
    {
        //fazer aqui o somatorio
        if (turno_funcao[h] < w0->matriz_requerimentos[h][coluna].first)
        {
            requerimento_minimo = (requerimento_minimo + (w0->matriz_requerimentos[h][coluna].first - turno_funcao[h]));
        }

        if (turno_funcao[h] < w0->matriz_requerimentos[h][coluna].second)
        {
            requerimento_otimo = (requerimento_otimo + (w0->matriz_requerimentos[h][coluna].second - turno_funcao[h]));
        }
        turno_funcao[h] = 0;
    }

    for (int h = 0; h < custom_out.size(); h++)
    {
        if (isdigit(custom_out[h]))
        {
            aux = custom_out[h] - '0';
            aux = cn->num_semanas - aux;
        }
    }

    for (int i = 0; i < w0->matriz_requisicoes.size(); i++)
    {
        if (w0->matriz_requisicoes[i][0] == linha)
        {
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first != "None" && w0->matriz_requisicoes[i][1] == 0)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Early" && w0->matriz_requisicoes[i][1] == 1)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Day" && w0->matriz_requisicoes[i][1] == 2)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Late" && w0->matriz_requisicoes[i][1] == 3)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Night" && w0->matriz_requisicoes[i][1] == 4)
            {
                requisicoes++;
            }
        }
    }

    // restricoes s2 e s3 ---- s2-peso 30 resolvida ---- s3 peso 30 a ser resolvida
    for (int i = linha; i <= linha; i++)
    {
        for (int j = 0; j < 7; j++)
        {

            if (j == 5)
            {
                if (cn->vetor_qualificacao[i].contract == "FullTime" && cn->vetor_contratos[0].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[i].contract == "PartTime" && cn->vetor_contratos[1].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[i].contract == "HalfTime" && cn->vetor_contratos[2].fim_semana_restricao == 1)
                {

                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
                if (cn->vetor_qualificacao[i].contract == "20Percent" && cn->vetor_contratos.size() == 4 && cn->vetor_contratos[3].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                count_folgas++;
                if (j == 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][4] < cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + ht->matriz_borda_fracas[i][4];
                            }
                            if (ht->matriz_borda_fracas[i][4] >= cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + cn->vetor_contratos[h].min_max[2].second;
                            }
                        }
                    }

                    // faz a subtração dos turnos que restarem da semana anterior
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].first)
                            {
                                count_dias_consecutivos = ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_dias_consecutivos > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos < cn->vetor_contratos[h].min_max[1].first)
                            {
                                dias_consecutivos = dias_consecutivos + (cn->vetor_contratos[h].min_max[1].first - count_dias_consecutivos);
                            }
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                            count_dias_consecutivos = 0;
                            break;
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                                break;
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_consecutivos++;
                if (j == 0)
                {
                    count_folgas = ht->matriz_borda_fracas[i][4];
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_folgas > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_folgas < cn->vetor_contratos[h].min_max[2].first)
                            {
                                dias_folgas = dias_folgas + (cn->vetor_contratos[h].min_max[2].first - count_folgas);
                            }
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                            }
                            count_folgas = 0;
                            
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_trabalhados++;
                if (j == 0)
                {
                    count_turnos++;
                    if (ht->matriz_borda_s[i][1] == individuo_decodificado[i][j].first)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + ht->matriz_borda_fracas[i][2];
                                }
                                if (ht->matriz_borda_fracas[i][2] >= cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + cn->vetor_sucessoes[h].x.second;
                                }
                            }
                        }
                    }

                    if (ht->matriz_borda_s[i][1] != individuo_decodificado[i][j].first)
                    {
                        if (ht->matriz_borda_fracas[i][2] > 0)
                        {
                            for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                            {
                                if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                                {
                                    if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                    {
                                        turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                    }
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first == turno_anterior)
                {
                    count_turnos++;
                    if (j == 6)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == turno_anterior)
                            {
                                if (count_turnos > cn->vetor_sucessoes[h].x.second)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first != turno_anterior && j != 0)
                {
                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {
                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                    count_turnos = 1;
                }
                turno_anterior = individuo_decodificado[i][j].first;
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                if (j == 0)
                {
                    if (ht->matriz_borda_fracas[i][2] > 0)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                }
                            }
                        }
                    }
                }
                if (count_turnos > 0)
                {

                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {

                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                }

                count_turnos = 0;
                turno_anterior = individuo_decodificado[i][j].first;
            }
        }
        count_dias_consecutivos = 0;
        count_folgas = 0;
        count_turnos = 0;
        turno_anterior.erase();
        //restrição turnos consecutivos

        //restriçao s6

        if (semana_1 == 0)
        {

            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        //s7++;
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }

                    
                    int divisao_first = floor(double(cus->dias_trabalhados[i].first / aux));
                    int divisao_second = ceil(double(cus->dias_trabalhados[i].second / aux) );

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                        
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                        
                    }
                    break;
                }
            }
        }

        if (semana_1 == 1)
        {

            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {

                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        //s7++;
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }
                    int divisao_first = floor(double(cn->vetor_contratos[t].min_max[0].first / cn->num_semanas));
                    int divisao_second = ceil(double(cn->vetor_contratos[t].min_max[0].second / cn->num_semanas) );

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                        
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                       
                    }
                    break;
                }
            }
        }
        count_dias_trabalhados = 0;
    }

    int sucessoes = 0;
    bool contem;
    for (int i = linha; i <= linha; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            contem = false;
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }

            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }
        }
    }

    
    sucessoes = (sucessoes * 60);
    requerimento_minimo = (requerimento_minimo * 60);
    requerimento_otimo = (requerimento_otimo * 30);
    requisicoes = (requisicoes * 10);
    fim_semana_completo = (fim_semana_completo * 30);
    dias_consecutivos = (dias_consecutivos * 30);
    dias_folgas = (dias_folgas * 30);
    turnos_consecutivos = (turnos_consecutivos * 15);
    if (aux <= cn->num_semanas/2)
    {
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }
    else{
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }
        

    custo_total = s7 + sucessoes + requerimento_minimo + requerimento_otimo + requisicoes + fim_semana_completo + dias_consecutivos + dias_folgas + turnos_consecutivos + dias_global;
    //cout << endl;
    return custo_total;
}

int BRKGA::funcao_avaliacao_linha_coluna_2(vector<vector<pair<string, string>>> &individuo_decodificado, int linha, vector<int> colunas)
{

    // int j = coluna;
    int s7 = 0;
    int aux = 0;
    int count_dias_trabalhados = 0;
    int dias_global = 0;

    int count = 0;
    int turnos_consecutivos = 0;
    int count_folgas = 0;
    int dias_folgas = 0;
    string turno_anterior;
    int count_dias_consecutivos = 0;
    int dias_consecutivos = 0;
    int count_fim_semana = 0;
    int count_turnos = 0;
    int custo_total;
    int requisicoes = 0;
    int fim_semana_completo = 0;
    int requerimento_otimo = 0;
    int requerimento_minimo = 0;
    vector<int> turno_funcao;
    turno_funcao = vector<int>(16); // coloca todas as posições do vetor com valor 0
    for (int i = 0; i < 16; i++)
    {
        turno_funcao[i] = 0;
    }

for(int i = 0; i < colunas.size(); i++){
    for (int j = 0; j < cn->vetor_qualificacao.size(); j++)
    {
        if (individuo_decodificado[j][colunas[i]].first == "Early")
        {
            if (individuo_decodificado[j][colunas[i]].second == "HeadNurse")
            {
                turno_funcao[0]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Nurse")
            {
                turno_funcao[1]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Caretaker")
            {
                turno_funcao[2]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Trainee")
            {
                turno_funcao[3]++;
            }
        }

        if (individuo_decodificado[j][colunas[i]].first == "Day")
        {
            if (individuo_decodificado[j][colunas[i]].second == "HeadNurse")
            {
                turno_funcao[4]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Nurse")
            {
                turno_funcao[5]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Caretaker")
            {
                turno_funcao[6]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Trainee")
            {
                turno_funcao[7]++;
            }
        }

        if (individuo_decodificado[j][colunas[i]].first == "Late")
        {
            if (individuo_decodificado[j][colunas[i]].second == "HeadNurse")
            {
                turno_funcao[8]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Nurse")
            {
                turno_funcao[9]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Caretaker")
            {
                turno_funcao[10]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Trainee")
            {
                turno_funcao[11]++;
            }
        }
        if (individuo_decodificado[j][colunas[i]].first == "Night")
        {
            if (individuo_decodificado[j][colunas[i]].second == "HeadNurse")
            {
                turno_funcao[12]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Nurse")
            {
                turno_funcao[13]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Caretaker")
            {
                turno_funcao[14]++;
            }
            if (individuo_decodificado[j][colunas[i]].second == "Trainee")
            {
                turno_funcao[15]++;
            }
        }
    }

    for (int h = 0; h < 16; h++)
    {
        //fazer aqui o somatorio
        if (turno_funcao[h] < w0->matriz_requerimentos[h][colunas[i]].first)
        {
            requerimento_minimo = (requerimento_minimo + (w0->matriz_requerimentos[h][colunas[i]].first - turno_funcao[h]));
        }

        if (turno_funcao[h] < w0->matriz_requerimentos[h][colunas[i]].second)
        {
            requerimento_otimo = (requerimento_otimo + (w0->matriz_requerimentos[h][colunas[i]].second - turno_funcao[h]));
        }
        turno_funcao[h] = 0;
    }

}
    
    for (int h = 0; h < custom_out.size(); h++)
    {
        if (isdigit(custom_out[h]))
        {
            aux = custom_out[h] - '0';
            aux = cn->num_semanas - aux;
        }
    }

    for (int i = 0; i < w0->matriz_requisicoes.size(); i++)
    {
        if (w0->matriz_requisicoes[i][0] == linha)
        {
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first != "None" && w0->matriz_requisicoes[i][1] == 0)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Early" && w0->matriz_requisicoes[i][1] == 1)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Day" && w0->matriz_requisicoes[i][1] == 2)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Late" && w0->matriz_requisicoes[i][1] == 3)
            {
                requisicoes++;
            }
            if (individuo_decodificado[w0->matriz_requisicoes[i][0]][w0->matriz_requisicoes[i][2]].first == "Night" && w0->matriz_requisicoes[i][1] == 4)
            {
                requisicoes++;
            }
        }
    }

    // restricoes s2 e s3 ---- s2-peso 30 resolvida ---- s3 peso 30 a ser resolvida
    for (int i = linha; i <= linha; i++)
    {
        for (int j = 0; j < 7; j++)
        {

            if (j == 5)
            {
                if (cn->vetor_qualificacao[i].contract == "FullTime" && cn->vetor_contratos[0].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[i].contract == "PartTime" && cn->vetor_contratos[1].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }

                if (cn->vetor_qualificacao[i].contract == "HalfTime" && cn->vetor_contratos[2].fim_semana_restricao == 1)
                {

                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
                if (cn->vetor_qualificacao[i].contract == "20Percent" && cn->vetor_contratos.size() == 4 && cn->vetor_contratos[3].fim_semana_restricao == 1)
                {
                    if (individuo_decodificado[i][j].first == "None" && individuo_decodificado[i][j + 1].first != "None")
                    {
                        fim_semana_completo++;
                    }
                    if (individuo_decodificado[i][j].first != "None" && individuo_decodificado[i][j + 1].first == "None")
                    {
                        fim_semana_completo++;
                    }
                }
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                count_folgas++;
                if (j == 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][4] < cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + ht->matriz_borda_fracas[i][4];
                            }
                            if (ht->matriz_borda_fracas[i][4] >= cn->vetor_contratos[h].min_max[2].second)
                            {
                                count_folgas = count_folgas + cn->vetor_contratos[h].min_max[2].second;
                            }
                        }
                    }

                    // faz a subtração dos turnos que restarem da semana anterior
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].first)
                            {
                                count_dias_consecutivos = ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_dias_consecutivos > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos < cn->vetor_contratos[h].min_max[1].first)
                            {
                                dias_consecutivos = dias_consecutivos + (cn->vetor_contratos[h].min_max[1].first - count_dias_consecutivos);
                            }
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                            count_dias_consecutivos = 0;
                            break;
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                                break;
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_consecutivos++;
                if (j == 0)
                {
                    count_folgas = ht->matriz_borda_fracas[i][4];
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (ht->matriz_borda_fracas[i][3] < cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + ht->matriz_borda_fracas[i][3];
                            }
                            if (ht->matriz_borda_fracas[i][3] >= cn->vetor_contratos[h].min_max[1].second)
                            {
                                count_dias_consecutivos = count_dias_consecutivos + cn->vetor_contratos[h].min_max[1].second;
                            }
                        }
                    }
                }

                if (count_folgas > 0)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_folgas < cn->vetor_contratos[h].min_max[2].first)
                            {
                                dias_folgas = dias_folgas + (cn->vetor_contratos[h].min_max[2].first - count_folgas);
                            }
                            if (count_folgas > cn->vetor_contratos[h].min_max[2].second)
                            {
                                dias_folgas = dias_folgas + (count_folgas - cn->vetor_contratos[h].min_max[2].second);
                            }
                            count_folgas = 0;
                            
                        }
                    }
                }

                if (j == 6)
                {
                    for (int h = 0; h < cn->vetor_contratos.size(); h++)
                    {
                        if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[h].contrato)
                        {
                            if (count_dias_consecutivos > cn->vetor_contratos[h].min_max[1].second)
                            {
                                dias_consecutivos = dias_consecutivos + (count_dias_consecutivos - cn->vetor_contratos[h].min_max[1].second);
                            }
                        }
                    }
                }
            }

            if (individuo_decodificado[i][j].first != "None")
            {
                count_dias_trabalhados++;
                if (j == 0)
                {
                    count_turnos++;
                    if (ht->matriz_borda_s[i][1] == individuo_decodificado[i][j].first)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + ht->matriz_borda_fracas[i][2];
                                }
                                if (ht->matriz_borda_fracas[i][2] >= cn->vetor_sucessoes[h].x.second)
                                {
                                    count_turnos = count_turnos + cn->vetor_sucessoes[h].x.second;
                                }
                            }
                        }
                    }

                    if (ht->matriz_borda_s[i][1] != individuo_decodificado[i][j].first)
                    {
                        if (ht->matriz_borda_fracas[i][2] > 0)
                        {
                            for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                            {
                                if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                                {
                                    if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                    {
                                        turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                    }
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first == turno_anterior)
                {
                    count_turnos++;
                    if (j == 6)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == turno_anterior)
                            {
                                if (count_turnos > cn->vetor_sucessoes[h].x.second)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                                }
                            }
                        }
                    }
                }

                if (individuo_decodificado[i][j].first != turno_anterior && j != 0)
                {
                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {
                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                    count_turnos = 1;
                }
                turno_anterior = individuo_decodificado[i][j].first;
            }

            if (individuo_decodificado[i][j].first == "None")
            {
                if (j == 0)
                {
                    if (ht->matriz_borda_fracas[i][2] > 0)
                    {
                        for (int h = 0; h < cn->vetor_sucessoes.size(); h++)
                        {
                            if (cn->vetor_sucessoes[h].turno == ht->matriz_borda_s[i][1])
                            {
                                if (ht->matriz_borda_fracas[i][2] < cn->vetor_sucessoes[h].x.first)
                                {
                                    turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - ht->matriz_borda_fracas[i][2]);
                                }
                            }
                        }
                    }
                }
                if (count_turnos > 0)
                {

                    for (int h = 0; h < cn->vetor_sucessoes.size(); h++) // se o enfermeiro trabalhar e o historico for menos do que o minimo
                    {
                        if (cn->vetor_sucessoes[h].turno == turno_anterior)
                        {

                            if (count_turnos < cn->vetor_sucessoes[h].x.first)
                            {
                                turnos_consecutivos = turnos_consecutivos + (cn->vetor_sucessoes[h].x.first - count_turnos);
                            }
                            if (count_turnos > cn->vetor_sucessoes[h].x.second)
                            {
                                turnos_consecutivos = turnos_consecutivos + (count_turnos - cn->vetor_sucessoes[h].x.second);
                            }
                        }
                    }
                }

                count_turnos = 0;
                turno_anterior = individuo_decodificado[i][j].first;
            }
        }
        count_dias_consecutivos = 0;
        count_folgas = 0;
        count_turnos = 0;
        turno_anterior.erase();
        //restrição turnos consecutivos

        //restriçao s6

        if (semana_1 == 0)
        {

            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {
                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }

                    
                    int divisao_first = floor(double(cus->dias_trabalhados[i].first / aux));
                    int divisao_second = ceil(double(cus->dias_trabalhados[i].second / aux) );

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                        
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                        
                    }
                    break;
                }
            }
        }

        if (semana_1 == 1)
        {

            for (int t = 0; t < cn->vetor_contratos.size(); t++)
            {
                if (cn->vetor_qualificacao[i].contract == cn->vetor_contratos[t].contrato)
                {

                    if (individuo_decodificado[i][5].first != "None" || individuo_decodificado[i][6].first != "None")
                    {
                        //s7++;
                        if (ht->matriz_borda_fracas[i][1] + 1 > cn->vetor_contratos[t].fim_semana)
                        {
                            s7++;
                        }
                    }
                   
                    int divisao_first = floor(double(cn->vetor_contratos[t].min_max[0].first / cn->num_semanas));
                    int divisao_second = ceil(double(cn->vetor_contratos[t].min_max[0].second / cn->num_semanas) );

                    if (count_dias_trabalhados < divisao_first)
                    {
                        dias_global = dias_global + (divisao_first - count_dias_trabalhados);
                        
                    }
                    if (count_dias_trabalhados > divisao_second)
                    {
                        dias_global = dias_global + (count_dias_trabalhados - divisao_second);
                        
                    }
                    break;
                }
            }
        }
        count_dias_trabalhados = 0;
    }

    int sucessoes = 0;
    bool contem;
    for (int i = linha; i <= linha; i++)
    {
        for (int j = 0; j < 7; j++)
        {
            contem = false;
            if (j == 0)
            {
                if (ht->matriz_borda_s[i][1] == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (ht->matriz_borda_s[i][1] == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }

            if (j != 0)
            {
                if (individuo_decodificado[i][j - 1].first == "None")
                {
                    for (int h = 0; h < cn->vetor_none.size(); h++)
                    {
                        if (cn->vetor_none[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_none.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Early")
                {
                    for (int h = 0; h < cn->vetor_early.size(); h++)
                    {
                        if (cn->vetor_early[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_early.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Day")
                {
                    for (int h = 0; h < cn->vetor_day.size(); h++)
                    {
                        if (cn->vetor_day[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_day.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Late")
                {
                    for (int h = 0; h < cn->vetor_late.size(); h++)
                    {
                        if (cn->vetor_late[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_late.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }

                if (individuo_decodificado[i][j - 1].first == "Night")
                {
                    for (int h = 0; h < cn->vetor_night.size(); h++)
                    {
                        if (cn->vetor_night[h] == individuo_decodificado[i][j].first)
                        {
                            contem = true;
                            h = cn->vetor_night.size();
                        }
                    }
                    if (contem == false)
                        sucessoes++;
                }
            }
        }
    }

    
    sucessoes = (sucessoes * 60);
    requerimento_minimo = (requerimento_minimo * 60);
    requerimento_otimo = (requerimento_otimo * 30);
    requisicoes = (requisicoes * 10);
    fim_semana_completo = (fim_semana_completo * 30);
    dias_consecutivos = (dias_consecutivos * 30);
    dias_folgas = (dias_folgas * 30);
    turnos_consecutivos = (turnos_consecutivos * 15);
    if (aux <= cn->num_semanas/2)
    {
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }
    else{
        s7 = (s7 * 30);
        dias_global = (dias_global * 20);
    }
        
 

    custo_total = s7 + sucessoes + requerimento_minimo + requerimento_otimo + requisicoes + fim_semana_completo + dias_consecutivos + dias_folgas + turnos_consecutivos + dias_global;
    
    return custo_total;
}