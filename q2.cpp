
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <chrono>

using namespace std;

struct Conexao {
    int origem;
    int destino;
    int custo;
};

struct DetalhesGrafo {
    int numVertices;
    vector<vector<int>> custos;
    unordered_map<int, int> demandas;
};

struct Solucao {
    vector<int> melhorRota;
    int menorCusto;
};

const int CAPACIDADE_VEICULO = 15;

DetalhesGrafo lerArquivoGrafo(const string& nomeArquivo) {
    DetalhesGrafo grafo;
    ifstream arquivo(nomeArquivo);

    if (!arquivo.is_open()) {
        cerr << "Erro ao abrir o arquivo!" << endl;
        exit(1);
    }

    arquivo >> grafo.numVertices;
    grafo.custos.resize(grafo.numVertices, vector<int>(grafo.numVertices, numeric_limits<int>::max()));
    grafo.demandas[0] = 0;

    for (int i = 1; i < grafo.numVertices; ++i) {
        int vertice, demanda;
        arquivo >> vertice >> demanda;
        grafo.demandas[vertice] = demanda;
    }

    int numConexoes;
    arquivo >> numConexoes;

    for (int i = 0; i < numConexoes; ++i) {
        int origem, destino, custo;
        arquivo >> origem >> destino >> custo;
        grafo.custos[origem][destino] = custo;
        grafo.custos[destino][origem] = custo; // Considera-se que o grafo é não direcionado
    }

    arquivo.close();
    return grafo;
}

int encontrarMaisProximo(const DetalhesGrafo& grafo, const vector<bool>& visitado, int ultimo, int cargaRestante) {
    int distanciaMinima = numeric_limits<int>::max();
    int maisProximo = -1;

    for (int i = 1; i < grafo.numVertices; ++i) {
        if (!visitado[i] && grafo.custos[ultimo][i] < distanciaMinima && grafo.demandas.at(i) <= cargaRestante) {
            distanciaMinima = grafo.custos[ultimo][i];
            maisProximo = i;
        }
    }

    return maisProximo;
}

Solucao heuristicaInsercaoMaisProxima(const DetalhesGrafo& grafo) {
    int n = grafo.numVertices;
    vector<bool> visitado(n, false);
    vector<int> rota;
    int custoTotal = 0;
    int cargaAtual = 0;

    rota.push_back(0);
    visitado[0] = true;
    int ultimo = 0;

    for (int i = 1; i < n; ++i) {
        int maisProximo = encontrarMaisProximo(grafo, visitado, ultimo, CAPACIDADE_VEICULO - cargaAtual);
        if (maisProximo != -1) {
            rota.push_back(maisProximo);
            visitado[maisProximo] = true;
            cargaAtual += grafo.demandas.at(maisProximo);
            custoTotal += grafo.custos[ultimo][maisProximo];
            ultimo = maisProximo;
        } else {
            custoTotal += grafo.custos[ultimo][0];
            rota.push_back(0);
            cargaAtual = 0;
            ultimo = 0;
            i--;  // Reavaliar a inserção de nós restantes
        }
    }

    if (ultimo != 0) {
        custoTotal += grafo.custos[ultimo][0];
        rota.push_back(0);
    }

    Solucao solucao;
    solucao.melhorRota = rota;
    solucao.menorCusto = custoTotal;

    return solucao;
}

int main() {
    auto inicio = std::chrono::high_resolution_clock::now();
    string nomeArquivo = "./grafo.txt";
    int numVertices;
    unordered_map<int, int> demandasVertices;
    vector<Conexao> listaConexoes;

    DetalhesGrafo grafo = lerArquivoGrafo(nomeArquivo);
    numVertices = grafo.numVertices;
    demandasVertices = grafo.demandas;

    // Convertendo a matriz de custos para a lista de conexões
    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            if (grafo.custos[i][j] != numeric_limits<int>::max() && i != j) {
                Conexao conexao;
                conexao.origem = i;
                conexao.destino = j;
                conexao.custo = grafo.custos[i][j];
                listaConexoes.push_back(conexao);
            }
        }
    }

    vector<int> vertices;
    for (int i = 1; i < numVertices; ++i) {
        vertices.push_back(i);
    }

    int capacidade = 15;

    Solucao melhorSolucao = heuristicaInsercaoMaisProxima(grafo);

    cout << "Melhores rotas:" << endl;
    int numeroRota = 1;
    vector<int> rotaAtual;
    for (int i = 0; i < melhorSolucao.melhorRota.size(); ++i) {
        if (melhorSolucao.melhorRota[i] == 0 && !rotaAtual.empty()) {
            cout << numeroRota << ": [0, ";
            for (size_t j = 0; j < rotaAtual.size(); ++j) {
                cout << rotaAtual[j];
                if (j < rotaAtual.size() - 1) cout << ", ";
            }
            cout << ", 0]" << endl;
            rotaAtual.clear();
            numeroRota++;
        } else if (melhorSolucao.melhorRota[i] != 0) {
            rotaAtual.push_back(melhorSolucao.melhorRota[i]);
        }
    }

    cout << "Menor custo: " << melhorSolucao.menorCusto << endl;

    auto fim = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duracao = fim - inicio;
    std::cout << "Tempo de execução: " << duracao.count() << " segundos" << std::endl;
    std::cout << "Número de nós: " << numVertices << "\n";

    return 0;
}
