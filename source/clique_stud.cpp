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
    
    // GRASP: Жадный рандомизированный адаптивный поиск для построения начальной клики
    vector<int> graspConstruct(double alpha = 0.5) {
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
        
        // Адаптивный процесс построения
        while (true) {
            // Список кандидатов с их оценками
            vector<pair<int, int>> candidate_scores;
            
            // Оцениваем всех кандидатов
            for (int i = 0; i < n; i++) {
                if (used[i]) continue;
                
                // Считаем, сколько вершин в текущей клике соединено с i
                int connections = 0;
                for (int v : clique) {
                    if (adj_matrix[i][v]) connections++;
                }
                
                // Можем добавить только если соединена со всеми
                if (connections == clique.size()) {
                    candidate_scores.emplace_back(neighbour_sets[i].size(), i);
                }
            }
            
            if (candidate_scores.empty()) break;
            
            // Сортируем по убыванию оценки
            sort(candidate_scores.rbegin(), candidate_scores.rend());
            
            // Строим RCL (Restricted Candidate List): берём лучшие alpha*100% кандидатов
            int rcl_size = max(1, (int)(candidate_scores.size() * alpha));
            
            // Случайный выбор из RCL
            static mt19937 rng(time(0));
            uniform_int_distribution<int> dist(0, rcl_size - 1);
            int selected_idx = dist(rng);
            int selected_vertex = candidate_scores[selected_idx].second;
            
            clique.push_back(selected_vertex);
            used[selected_vertex] = true;
        }
        
        return clique;
    }
    
    // Локальный поиск для улучшения клики
    void localSearch(vector<int>& clique) {
        vector<bool> in_clique(n, false);
        for (int v : clique) in_clique[v] = true;
        
        bool improved = true;
        
        while (improved) {
            improved = false;
            
            // Фаза добавления: пытаемся добавить все возможные вершины
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
            
            // Фаза замены: пробуем заменить 1-2 вершины для улучшения
            if (!improved && clique.size() >= 4) {
                // Пробуем заменить 1 вершину на 2
                for (int i = 0; i < clique.size() && !improved; i++) {
                    int v_to_remove = clique[i];
                    
                    // Ищем две вершины, которые можно добавить вместо v_to_remove
                    for (int x = 0; x < n && !improved; x++) {
                        if (in_clique[x] || x == v_to_remove) continue;
                        for (int y = x + 1; y < n && !improved; y++) {
                            if (in_clique[y] || y == v_to_remove) continue;
                            
                            // Проверяем, можно ли добавить x и y
                            bool valid = true;
                            for (int u : clique) {
                                if (u != v_to_remove) {
                                    if (!adj_matrix[x][u] || !adj_matrix[y][u]) {
                                        valid = false;
                                        break;
                                    }
                                }
                            }
                            if (valid && adj_matrix[x][y]) {
                                // Делаем замену
                                clique.erase(clique.begin() + i);
                                clique.push_back(x);
                                clique.push_back(y);
                                
                                in_clique[v_to_remove] = false;
                                in_clique[x] = true;
                                in_clique[y] = true;
                                
                                improved = true;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Итеративный GRASP: многократный запуск GRASP с локальным поиском
    void iterativeGRASP(int iterations, int randomization) {
        mt19937 rng(static_cast<unsigned int>(time(0)));
        
        for (int iter = 0; iter < iterations; ++iter) {
            // Выбираем случайный alpha в диапазоне [0.3, 0.7] для разнообразия
            double alpha = 0.3 + (rng() % 5) * 0.1;
            
            // Построение начальной клики с помощью GRASP
            vector<int> clique = graspConstruct(alpha);
            
            // Локальное улучшение
            localSearch(clique);
            
            // Обновляем лучший результат
            if (clique.size() > best_clique.size()) {
                best_clique = clique;
                cout << "Iteration " << iter + 1 << ": found clique of size " << clique.size() << endl;
            }
            
            // Иногда делаем более интенсивный локальный поиск для лучших решений
            if (clique.size() >= best_clique.size() * 0.9) {
                vector<int> intensified_clique = clique;
                localSearch(intensified_clique);
                if (intensified_clique.size() > best_clique.size()) {
                    best_clique = intensified_clique;
                }
            }
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
        
        // Строим матрицу смежности
        buildAdjMatrix();
    }
    
    void FindClique(int iterations, int randomization)
    {
        best_clique.clear();
        
        if (neighbour_sets.empty()) {
            cout << "Error: Graph is empty!" << endl;
            return;
        }
        
        n = neighbour_sets.size();
        cout << "Searching for max clique in graph with " << n << " vertices" << endl;
        cout << "Using Iterative GRASP with " << iterations << " iterations" << endl;
        
        clock_t start = clock();
        
        // Запускаем итеративный GRASP
        iterativeGRASP(iterations, randomization);
        
        // Финальное улучшение лучшей найденной клики
        if (!best_clique.empty()) {
            localSearch(best_clique);
        }
        
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
                if (!adj_matrix[v1][v2]) {
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
    vector<string> files = { "C125.9.clq", "johnson8-2-4.clq", "johnson16-2-4.clq", "MANN_a9.clq", "MANN_a27.clq",
        "p_hat1000-1.clq", "keller4.clq", "hamming8-4.clq", "brock200_1.clq", "brock200_2.clq", "brock200_3.clq", "brock200_4.clq",
        "gen200_p0.9_44.clq", "gen200_p0.9_55.clq", "brock400_1.clq", "brock400_2.clq", "brock400_3.clq", "brock400_4.clq",
        "MANN_a45.clq", "sanr400_0.7.clq", "p_hat1000-2.clq", "p_hat500-3.clq", "p_hat1500-1.clq", "p_hat300-3.clq", "san1000.clq",
        "sanr200_0.9.clq" };
    vector<int> iterations = {700, 1, 1, 1000, 500000,
         900000, 20000, 20000, 100, 1000, 100, 1000000,
         500000, 200000, 5000000, 2000000, 3000000, 3000000,
         100000, 4000000, 13000000, 4000000, 13000000, 2000000, 13000000,
         1000000
         };
    
    ofstream fout("clique.csv");
    fout << "File; Clique; Time (sec)\n";
    
    int count = 10;
    int randomization = 10000;
    string file = files[count];
    // for (string file : files)
    // {
        MaxCliqueProblem problem;
        string full_path = PACKAGE_FOLDER + file;
        cout << "\n========================================" << endl;
        cout << "Processing file: " << full_path << endl;
        
        problem.ReadGraphFile(full_path);
        
        clock_t start = clock();
        problem.FindClique(iterations[count], randomization);
        
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