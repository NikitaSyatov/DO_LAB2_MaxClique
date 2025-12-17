#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <ctime>
#include <random>
#include <unordered_set>
#include <algorithm>
#include <stack>
#include <set>
using namespace std;

#define PACKAGE_FOLDER "../test_packages/"

class MaxCliqueProblem
{
private:
    vector<unordered_set<int>> neighbour_sets;
    vector<int> best_clique;
    vector<vector<bool>> adj_matrix;
    int n; // количество вершин
    
    // Строим матрицу смежности для быстрого доступа
    void buildAdjMatrix() {
        if (!adj_matrix.empty()) return;
        
        adj_matrix.assign(n, vector<bool>(n, false));
        for (int i = 0; i < n; i++) {
            for (int j : neighbour_sets[i]) {
                adj_matrix[i][j] = true;
            }
        }
    }
    
    // Простой жадный алгоритм (используется как базовая эвристика)
    vector<int> simpleGreedyClique() {
        vector<int> clique;
        vector<bool> used(n, false);
        
        // Начинаем с вершины с максимальной степенью
        int max_degree = -1;
        int start_vertex = 0;
        
        for (int i = 0; i < n; i++) {
            if ((int)neighbour_sets[i].size() > max_degree) {
                max_degree = neighbour_sets[i].size();
                start_vertex = i;
            }
        }
        
        clique.push_back(start_vertex);
        used[start_vertex] = true;
        
        // Пытаемся добавить другие вершины
        for (int i = 0; i < n; i++) {
            if (used[i]) continue;
            
            bool can_add = true;
            for (int v : clique) {
                if (!adj_matrix[i][v]) {
                    can_add = false;
                    break;
                }
            }
            
            if (can_add) {
                clique.push_back(i);
                used[i] = true;
            }
        }
        
        return clique;
    }
    
    // Основной алгоритм: рандомизированный жадный с использованием iterations и randomization
    void randomizedGreedySearch(int iterations, int randomization) {
        if (n == 0) return;
        
        buildAdjMatrix(); // строим матрицу один раз
        
        mt19937 rng(static_cast<unsigned int>(time(0)));
        
        for (int iter = 0; iter < iterations; ++iter) {
            vector<int> clique;
            vector<int> candidates(n);
            for (int i = 0; i < n; i++) candidates[i] = i;
            
            // Случайное перемешивание вершин
            shuffle(candidates.begin(), candidates.end(), rng);
            
            // Основной цикл построения клики
            while (!candidates.empty()) {
                int last = candidates.size() - 1;
                
                // Выбираем случайную вершину из первых randomization кандидатов
                int rnd = 0;
                if (randomization > 1) {
                    uniform_int_distribution<int> dist(0, min(randomization - 1, last));
                    rnd = dist(rng);
                }
                
                int vertex = candidates[rnd];
                clique.push_back(vertex);
                
                // Фильтруем кандидатов: оставляем только тех, кто соединен с выбранной вершиной
                vector<int> new_candidates;
                for (int candidate : candidates) {
                    if (candidate == vertex) continue;
                    if (adj_matrix[vertex][candidate]) {
                        new_candidates.push_back(candidate);
                    }
                }
                
                candidates = new_candidates;
                
                // Перемешиваем оставшихся кандидатов для следующей итерации
                if (!candidates.empty() && randomization > 1) {
                    shuffle(candidates.begin(), candidates.end(), rng);
                }
            }
            
            // Улучшаем найденную клику: пытаемся добавить пропущенные вершины
            improveClique(clique);
            
            // Обновляем лучший результат
            if (clique.size() > best_clique.size()) {
                best_clique = clique;
            }
        }
    }
    
    // Улучшение клики: пытаемся добавить все возможные вершины
    void improveClique(vector<int>& clique) {
        vector<bool> in_clique(n, false);
        for (int v : clique) in_clique[v] = true;
        
        bool improved = true;
        while (improved) {
            improved = false;
            
            for (int v = 0; v < n; v++) {
                if (in_clique[v]) continue;
                
                bool can_add = true;
                for (int u : clique) {
                    if (!adj_matrix[v][u]) {
                        can_add = false;
                        break;
                    }
                }
                
                if (can_add) {
                    clique.push_back(v);
                    in_clique[v] = true;
                    improved = true;
                }
            }
        }
    }
    
    // Метод ветвей и границ с использованием параметров
    void branchAndBoundSearch(int max_nodes, int randomization) {
        if (n == 0) return;
        
        buildAdjMatrix();
        
        // Начальное решение через жадный алгоритм
        best_clique = simpleGreedyClique();
        
        // Структура для узла поиска
        struct Node {
            vector<int> clique;
            vector<int> candidates;
            int bound; // верхняя оценка
            
            Node(vector<int> c, vector<int> cand) : clique(c), candidates(cand) {
                bound = clique.size() + candidates.size();
            }
        };
        
        // Создаем начальные кандидаты (все вершины)
        vector<int> all_vertices(n);
        for (int i = 0; i < n; i++) all_vertices[i] = i;
        
        // Сортируем по убыванию степени
        sort(all_vertices.begin(), all_vertices.end(), [&](int a, int b) {
            return neighbour_sets[a].size() > neighbour_sets[b].size();
        });
        
        stack<Node> stk;
        stk.push(Node({}, all_vertices));
        
        int nodes_processed = 0;
        mt19937 rng(static_cast<unsigned int>(time(0)));
        
        while (!stk.empty() && nodes_processed < max_nodes) {
            Node node = stk.top();
            stk.pop();
            nodes_processed++;
            
            // Отсечение по оценке
            if (node.bound <= (int)best_clique.size()) {
                continue;
            }
            
            if (node.candidates.empty()) {
                if (node.clique.size() > best_clique.size()) {
                    best_clique = node.clique;
                }
                continue;
            }
            
            // Выбор вершины для ветвления с использованием randomization
            int selected_idx = 0;
            if (randomization > 1 && node.candidates.size() > 1) {
                int k = min(randomization, (int)node.candidates.size());
                uniform_int_distribution<int> dist(0, k - 1);
                selected_idx = dist(rng);
            }
            
            int selected_vertex = node.candidates[selected_idx];
            
            // Ветвь 1: включаем вершину в клику
            vector<int> new_clique = node.clique;
            new_clique.push_back(selected_vertex);
            
            // Формируем новых кандидатов
            vector<int> new_candidates;
            for (int v : node.candidates) {
                if (v == selected_vertex) continue;
                
                bool compatible = true;
                for (int u : new_clique) {
                    if (!adj_matrix[v][u]) {
                        compatible = false;
                        break;
                    }
                }
                
                if (compatible) {
                    new_candidates.push_back(v);
                }
            }
            
            Node with_node(new_clique, new_candidates);
            if (with_node.bound > (int)best_clique.size()) {
                if (new_candidates.empty()) {
                    if (new_clique.size() > best_clique.size()) {
                        best_clique = new_clique;
                    }
                } else {
                    stk.push(with_node);
                }
            }
            
            // Ветвь 2: исключаем вершину из клики
            vector<int> candidates_without;
            for (int v : node.candidates) {
                if (v != selected_vertex) {
                    candidates_without.push_back(v);
                }
            }
            
            Node without_node(node.clique, candidates_without);
            if (without_node.bound > (int)best_clique.size() && !candidates_without.empty()) {
                stk.push(without_node);
            }
        }
    }
    
    // Адаптивный выбор алгоритма в зависимости от размера графа
    void adaptiveSearch(int iterations, int randomization, double time_limit = 9.5) {
        clock_t start = clock();
        
        // Для очень маленьких графов (< 50 вершин) используем полный поиск
        if (n < 50) {
            branchAndBoundSearch(iterations * 100, randomization);
        }
        // Для средних графов (50-200 вершин) комбинируем оба подхода
        else if (n < 200) {
            // Первая половина итераций - рандомизированный жадный
            randomizedGreedySearch(iterations / 2, randomization);
            
            // Вторая половина - ветви и границы, если осталось время
            if ((double)(clock() - start) / CLOCKS_PER_SEC < time_limit / 2) {
                branchAndBoundSearch(iterations * 50, randomization);
            }
        }
        // Для больших графов (> 200 вершин) только рандомизированный жадный
        else {
            randomizedGreedySearch(iterations, randomization);
        }
        
        // Финальное улучшение
        if (!best_clique.empty()) {
            improveClique(best_clique);
        }
    }

public:
    void ReadGraphFile(string filename)
    {
        ifstream fin(filename);
        if (!fin.is_open()) {
            cout << "Error: Cannot open file " << filename << endl;
            return;
        }
        
        string line;
        int vertices = 0, edges = 0;
        bool got_header = false;
        
        while (getline(fin, line))
        {
            if (line.empty()) continue;
            
            if (line[0] == 'c') {
                continue;
            }

            stringstream line_input(line);
            char command;
            if (line[0] == 'p') {
                string type;
                line_input >> command >> type >> vertices >> edges;
                neighbour_sets.resize(vertices);
                n = vertices;
                cout << "Graph: " << vertices << " vertices, " << edges << " edges" << endl;
                got_header = true;
            }
            else if (got_header) {
                int start, finish;
                line_input >> command >> start >> finish;
                if (start > 0 && start <= vertices && finish > 0 && finish <= vertices) {
                    neighbour_sets[start - 1].insert(finish - 1);
                    neighbour_sets[finish - 1].insert(start - 1);
                }
            }
        }
        fin.close();
        
        // Проверяем, что граф загружен
        int edge_count = 0;
        for (const auto& s : neighbour_sets) {
            edge_count += s.size();
        }
        cout << "Loaded graph with " << vertices << " vertices and " << edge_count/2 << " edges" << endl;
    }
    
    void FindClique(int iterations, int randomization=3)
    {
        best_clique.clear();
        adj_matrix.clear();
        
        if (neighbour_sets.empty()) {
            cout << "Error: Graph is empty!" << endl;
            return;
        }
        
        n = neighbour_sets.size();
        cout << "Searching for max clique in graph with " << n << " vertices" << endl;
        cout << "Using " << iterations << " iterations, randomization = " << randomization << endl;
        
        clock_t start = clock();
        
        // Адаптивный выбор алгоритма
        adaptiveSearch(iterations, randomization, 9.5);
        
        double elapsed = double(clock() - start) / CLOCKS_PER_SEC;
        cout << "Time: " << elapsed << " seconds, clique size: " << best_clique.size() << endl;
    }

    const vector<int>& GetClique()
    {
        return best_clique;
    }

    bool Check()
    {
        if (best_clique.empty()) {
            cout << "Empty clique" << endl;
            return false;
        }
        
        // Проверяем индексы
        for (int v : best_clique) {
            if (v < 0 || v >= n) {
                cout << "Vertex " << v << " is out of bounds!" << endl;
                return false;
            }
        }
        
        // Проверяем на дубликаты
        set<int> unique_vertices(best_clique.begin(), best_clique.end());
        if (unique_vertices.size() != best_clique.size()) {
            cout << "Duplicate vertices found!" << endl;
            return false;
        }
        
        // Проверяем ребра
        for (size_t i = 0; i < best_clique.size(); i++) {
            for (size_t j = i + 1; j < best_clique.size(); j++) {
                int v1 = best_clique[i];
                int v2 = best_clique[j];
                if (neighbour_sets[v1].find(v2) == neighbour_sets[v1].end()) {
                    cout << "Missing edge between " << v1 << " and " << v2 << endl;
                    return false;
                }
            }
        }
        
        return true;
    }
};

int main()
{
    // int iterations;
    // cout << "Number of iterations: ";
    // cin >> iterations;
    // int randomization = 3;
    // cout << "Randomization level (1=deterministic, >1=randomized): ";
    // cin >> randomization;
    
    vector<string> files = { "C125.9.clq", "johnson8-2-4.clq", "johnson16-2-4.clq", "MANN_a9.clq", "MANN_a27.clq",
        "p_hat1000-1.clq", "keller4.clq", "hamming8-4.clq", "brock200_1.clq", "brock200_2.clq", "brock200_3.clq", "brock200_4.clq",
        "gen200_p0.9_44.clq", "gen200_p0.9_55.clq", "brock400_1.clq", "brock400_2.clq", "brock400_3.clq", "brock400_4.clq",
        "MANN_a45.clq", "sanr400_0.7.clq", "p_hat1000-2.clq", "p_hat500-3.clq", "p_hat1500-1.clq", "p_hat300-3.clq", "san1000.clq",
        "sanr200_0.9.clq" };
    vector<int> iterations = {700000, 1, 1, 1000, 500000,
         900000, 20000, 20000, 100000, 1000000, 1000000, 1000000,
         500000, 200000, 5000000, 2000000, 3000000, 3000000,
         100000, 4000000, 13000000, 4000000, 13000000, 2000000, 13000000,
         1000000
         };
    
    ofstream fout("clique.csv");
    fout << "File; Clique; Time (sec)\n";
    
    int count = 25;
    string file = files[count];
    // for (string file : files)
    // {
        MaxCliqueProblem problem;
        string full_path = PACKAGE_FOLDER + file;
        cout << "\n========================================" << endl;
        cout << "Processing file: " << full_path << endl;
        
        problem.ReadGraphFile(full_path);
        
        clock_t start = clock();
        problem.FindClique(iterations[count]);
        
        bool is_valid = problem.Check();
        
        if (!is_valid) {
            cout << "*** WARNING: incorrect clique ***" << endl;
            fout << "*** WARNING: incorrect clique ***\n";
        }
        
        double elapsed = double(clock() - start) / CLOCKS_PER_SEC;
        fout << file << "; " << problem.GetClique().size() << "; " << elapsed << '\n';
        
        cout << "Result: clique size = " << problem.GetClique().size() 
             << ", time = " << elapsed << " seconds" << endl;
        cout << "========================================\n" << endl;

        count++;
    // }
    
    fout.close();
    return 0;
}
