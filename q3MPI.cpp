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
#include <omp.h>
#include <mpi.h>

using namespace std;

struct Aresta {
    int origem;
    int destino;
    int custo;
};

void lerArquivoGrafo(const string& nomeArquivo, int& numVertices, unordered_map<int, int>& demandasVertices, vector<Aresta>& listaConexoes) {
    ifstream arquivo(nomeArquivo);
    if (!arquivo.is_open()) {
        cerr << "Erro ao abrir o arquivo!" << endl;
        exit(1);
    }

    arquivo >> numVertices;
    demandasVertices[0] = 0;

    for (int i = 1; i < numVertices; ++i) {
        int vertice, demanda;
        arquivo >> vertice >> demanda;
        demandasVertices[vertice] = demanda;
    }

    int numConexoes;
    arquivo >> numConexoes;

    for (int i = 0; i < numConexoes; ++i) {
        Aresta aresta;
        arquivo >> aresta.origem >> aresta.destino >> aresta.custo;
        listaConexoes.push_back(aresta);
    }

    arquivo.close();
}

vector<vector<int>> gerarPermutacoes(vector<int>& vertices) {
    vector<vector<int>> permutacoes;
    sort(vertices.begin(), vertices.end());
    do {
        permutacoes.push_back(vertices);
    } while (next_permutation(vertices.begin(), vertices.end()));
    return permutacoes;
}

vector<vector<int>> verificarCapacidade(vector<vector<int>>& rotas, const unordered_map<int, int>& demandasVertices, int capacidade) {
    vector<vector<int>> rotasValidas;

    for (auto& rota : rotas) {
        int cargaAtual = 0;
        bool rotaValida = true;

        for (int vertice : rota) {
            if (vertice == 0) {
                cargaAtual = 0;
            } else {
                cargaAtual += demandasVertices.at(vertice);
                if (cargaAtual > capacidade) {
                    rotaValida = false;
                    break;
                }
            }
        }
        if (rotaValida) {
            rotasValidas.push_back(rota);
        }
    }

    return rotasValidas;
}

int calcularCusto(const vector<int>& rota, const unordered_map<int, unordered_map<int, int>>& mapaCustoConexoes) {
    int custoTotal = 0;
    for (size_t i = 0; i < rota.size() - 1; ++i) {
        int origem = rota[i];
        int destino = rota[i + 1];
        auto it = mapaCustoConexoes.find(origem);
        if (it != mapaCustoConexoes.end()) {
            auto it2 = it->second.find(destino);
            if (it2 != it->second.end()) {
                custoTotal += it2->second;
            }
        }
    }
    return custoTotal;
}

vector<vector<int>> validarRotas(vector<vector<int>>& permutacoes, const vector<Aresta>& listaConexoes) {
    vector<vector<int>> rotasValidas;

    for (auto& permutacao : permutacoes) {
        if (permutacao[0] != 0) {
            permutacao.insert(permutacao.begin(), 0);
        }
        if (permutacao.back() != 0) {
            permutacao.push_back(0);
        }

        bool rotaValida = true;
        for (size_t i = 0; i < permutacao.size() - 1; ++i) {
            int origem = permutacao[i];
            int destino = permutacao[i + 1];
            bool arestaExiste = false;
            for (const auto& aresta : listaConexoes) {
                if (aresta.origem == origem && aresta.destino == destino) {
                    arestaExiste = true;
                    break;
                }
            }
            if (!arestaExiste) {
                bool rotaValida1 = false;
                bool rotaValida2 = false;
                for (const auto& aresta : listaConexoes) {
                    if (aresta.origem == origem && aresta.destino == 0) {
                        rotaValida1 = true;
                    }
                    if (aresta.origem == 0 && aresta.destino == destino) {
                        rotaValida2 = true;
                    }
                }
                if (!rotaValida1 || !rotaValida2) {
                    rotaValida = false;
                    break;
                } else {
                    permutacao.insert(permutacao.begin() + i + 1, 0);
                }
            }
        }
        if (rotaValida) {
            rotasValidas.push_back(permutacao);
        }
    }

    return rotasValidas;
}

unordered_map<int, unordered_map<int, int>> construirMapaCustos(const vector<Aresta>& listaConexoes) {
    unordered_map<int, unordered_map<int, int>> mapaCustoConexoes;
    for (const auto& aresta : listaConexoes) {
        mapaCustoConexoes[aresta.origem][aresta.destino] = aresta.custo;
    }
    return mapaCustoConexoes;
}

pair<vector<int>, int> resolverVRP(vector<int>& vertices, const unordered_map<int, int>& demandasVertices, const vector<Aresta>& listaConexoes, int capacidade) {
    vector<int> melhorRota;
    int menorCusto = numeric_limits<int>::max();

    auto mapaCustoConexoes = construirMapaCustos(listaConexoes);

    vector<vector<int>> permutacoes = gerarPermutacoes(vertices);

    permutacoes = validarRotas(permutacoes, listaConexoes);

    permutacoes = verificarCapacidade(permutacoes, demandasVertices, capacidade);

    #pragma omp parallel for
    for (int i = 0; i < permutacoes.size(); ++i) {
        int custo = calcularCusto(permutacoes[i], mapaCustoConexoes);
        #pragma omp critical
        if (custo < menorCusto) {
            menorCusto = custo;
            melhorRota = permutacoes[i];
        }
    }

    return {melhorRota, menorCusto};
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto inicio = std::chrono::high_resolution_clock::now();

    int numVertices;
    unordered_map<int, int> demandasVertices;
    vector<Aresta> listaConexoes;

    if (rank == 0) {
        lerArquivoGrafo("grafo.txt", numVertices, demandasVertices, listaConexoes);
    }

    MPI_Bcast(&numVertices, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        demandasVertices[0] = 0;
        for (int i = 1; i < numVertices; ++i) {
            demandasVertices[i] = 0;
        }
    }

    MPI_Bcast(&demandasVertices[0], numVertices, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&listaConexoes[0], listaConexoes.size() * sizeof(Aresta), MPI_BYTE, 0, MPI_COMM_WORLD);

    vector<int> vertices;
    for (int i = 1; i < numVertices; ++i) {
        vertices.push_back(i);
    }

    int capacidade = 15;

    auto [melhorRota, menorCusto] = resolverVRP(vertices, demandasVertices, listaConexoes, capacidade);

    int globalMenorCusto;
    MPI_Reduce(&menorCusto, &globalMenorCusto, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Melhores rotas:" << endl;
        int numeroRota = 1;
        vector<int> rotaAtual;
        for (int i = 0; i < melhorRota.size(); ++i) {
            if (melhorRota[i] == 0 && !rotaAtual.empty()) {
                cout << numeroRota << ": [0, ";
                for (size_t j = 0; j < rotaAtual.size(); ++j) {
                    cout << rotaAtual[j];
                    if (j < rotaAtual.size() - 1) cout << ", ";
                }
                cout << ", 0]" << endl;
                rotaAtual.clear();
                numeroRota++;
            } else if (melhorRota[i] != 0) {
                rotaAtual.push_back(melhorRota[i]);
            }
        }

        cout << "Menor custo: " << globalMenorCusto << endl;

        auto fim = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duracao = fim - inicio;
        std::cout << "Tempo de execução: " << duracao.count() << " segundos" << std::endl;
        std::cout << "Número de nós: " << numVertices << "\n";
    }

    MPI_Finalize();
    return 0;
}
