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
#include <queue>
#include <list>
#include <cmath>
using namespace std;

#define PACKAGE_FOLDER "../test_packages/"

class MaxCliqueProblem
{
private:
    vector<unordered_set<int>> neighbour_sets;
    vector<int> best_clique;
    vector<vector<bool>> adj_matrix;
    int n;
    
    void buildAdjMatrix() {
        if (!adj_matrix.empty()) return;
        
        adj_matrix.assign(n, vector<bool>(n, false));
        for (int i = 0; i < n; i++) {
            for (int j : neighbour_sets[i]) {
                adj_matrix[i][j] = true;
            }
        }
    }
    
    // ================== GRASP (Greedy Randomized Adaptive Search Procedure) ==================
    vector<int> graspConstruct(int randomization, double alpha = 0.5) {
        vector<int> clique;
        vector<bool> used(n, false);
        
        // RCL (Restricted Candidate List)
        vector<pair<int, int>> candidate_scores;
        
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
            candidate_scores.clear();
            
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
            
            // Строим RCL: берём лучшие alpha*100% кандидатов
            int rcl_size = max(1, (int)(candidate_scores.size() * alpha));
            if (randomization > 1) {
                rcl_size = min(rcl_size, randomization);
            }
            
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
    
    // Локальный поиск для GRASP
    void graspLocalSearch(vector<int>& clique, int max_iter = 1000) {
        vector<bool> in_clique(n, false);
        for (int v : clique) in_clique[v] = true;
        
        bool improved = true;
        int iter = 0;
        
        while (improved && iter < max_iter) {
            improved = false;
            iter++;
            
            // Пытаемся добавить вершины
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
            
            // Если не удалось добавить, пробуем замену 1 на 2
            if (!improved && clique.size() >= 3) {
                for (int i = 0; i < clique.size() && !improved; i++) {
                    for (int j = i + 1; j < clique.size() && !improved; j++) {
                        int v1 = clique[i];
                        int v2 = clique[j];
                        
                        // Ищем две вершины, которые можно добавить вместо v1 и v2
                        for (int x = 0; x < n && !improved; x++) {
                            if (in_clique[x] || x == v1 || x == v2) continue;
                            for (int y = x + 1; y < n && !improved; y++) {
                                if (in_clique[y] || y == v1 || y == v2) continue;
                                
                                // Проверяем, можно ли добавить x и y
                                bool valid = true;
                                for (int u : clique) {
                                    if (u != v1 && u != v2) {
                                        if (!adj_matrix[x][u] || !adj_matrix[y][u]) {
                                            valid = false;
                                            break;
                                        }
                                    }
                                }
                                if (valid && adj_matrix[x][y]) {
                                    // Делаем замену
                                    clique.erase(find(clique.begin(), clique.end(), v1));
                                    clique.erase(find(clique.begin(), clique.end(), v2));
                                    clique.push_back(x);
                                    clique.push_back(y);
                                    
                                    in_clique[v1] = false;
                                    in_clique[v2] = false;
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
    }
    
    // ================== Tabu Search ==================
    vector<int> tabuSearch(const vector<int>& initial_clique, int tabu_tenure, int max_iter, 
                          int max_time_seconds, clock_t start_time) {
        vector<int> current_clique = initial_clique;
        vector<int> best_clique_local = initial_clique;
        
        vector<int> tabu_list(n, 0); // Когда истекает запрет для вершины
        int iter = 0;
        
        while (iter < max_iter && 
               double(clock() - start_time) / CLOCKS_PER_SEC < max_time_seconds) {
            iter++;
            
            // Генерируем соседние решения
            vector<pair<int, vector<int>>> neighbors; // (оценка, клика)
            
            // 1. Добавление вершины
            vector<bool> in_current(n, false);
            for (int v : current_clique) in_current[v] = true;
            
            for (int v = 0; v < n; v++) {
                if (in_current[v]) continue;
                if (tabu_list[v] > iter) continue; // Вершина в табу-листе
                
                bool can_add = true;
                for (int u : current_clique) {
                    if (!adj_matrix[v][u]) {
                        can_add = false;
                        break;
                    }
                }
                
                if (can_add) {
                    vector<int> new_clique = current_clique;
                    new_clique.push_back(v);
                    neighbors.emplace_back(new_clique.size(), new_clique);
                }
            }
            
            // 2. Удаление вершины (только если клика большая)
            if (current_clique.size() > 10) {
                for (int v : current_clique) {
                    if (tabu_list[v] > iter) continue;
                    
                    vector<int> new_clique;
                    for (int u : current_clique) {
                        if (u != v) new_clique.push_back(u);
                    }
                    neighbors.emplace_back(new_clique.size(), new_clique);
                }
            }
            
            // 3. Замена вершины
            if (current_clique.size() >= 3) {
                for (int i = 0; i < current_clique.size(); i++) {
                    int v = current_clique[i];
                    if (tabu_list[v] > iter) continue;
                    
                    for (int x = 0; x < n; x++) {
                        if (in_current[x] || tabu_list[x] > iter) continue;
                        
                        // Проверяем, можно ли заменить v на x
                        bool valid = true;
                        for (int u : current_clique) {
                            if (u != v && !adj_matrix[x][u]) {
                                valid = false;
                                break;
                            }
                        }
                        
                        if (valid) {
                            vector<int> new_clique = current_clique;
                            new_clique[i] = x;
                            neighbors.emplace_back(new_clique.size(), new_clique);
                        }
                    }
                }
            }
            
            if (neighbors.empty()) {
                // Нет допустимых соседей, сбрасываем табу-лист
                fill(tabu_list.begin(), tabu_list.end(), 0);
                continue;
            }
            
            // Выбираем лучшего соседа (допускаем ухудшение)
            sort(neighbors.rbegin(), neighbors.rend());
            vector<int> next_clique = neighbors[0].second;
            
            // Обновляем табу-лист для изменённых вершин
            vector<bool> in_next(n, false);
            for (int v : next_clique) in_next[v] = true;
            
            for (int v : current_clique) {
                if (!in_next[v]) {
                    // Вершина удалена - запрещаем добавлять обратно
                    tabu_list[v] = iter + tabu_tenure + rand() % 3; // Небольшая вариация
                }
            }
            for (int v : next_clique) {
                if (!in_current[v]) {
                    // Новая вершина - запрещаем удалять
                    tabu_list[v] = iter + tabu_tenure / 2 + rand() % 2;
                }
            }
            
            current_clique = next_clique;
            
            // Обновляем лучшее решение
            if (current_clique.size() > best_clique_local.size()) {
                best_clique_local = current_clique;
            }
        }
        
        return best_clique_local;
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
        srand(time(0));
        
        if (neighbour_sets.empty()) {
            cout << "Error: Graph is empty!" << endl;
            return;
        }
        
        cout << "Searching for max clique in graph with " << n << " vertices" << endl;
        cout << "Using GRASP + Tabu Search with " << iterations << " iterations, randomization = " << randomization << endl;
        
        clock_t start = clock();
        double time_limit = 9.5; // секунд
        
        // Этап 1: GRASP для построения начального решения
        vector<int> initial_clique;
        double grasp_time_limit = time_limit * 0.3;
        clock_t grasp_start = clock();
        
        while (double(clock() - grasp_start) / CLOCKS_PER_SEC < grasp_time_limit) {
            vector<int> clique = graspConstruct(randomization, 0.3 + (rand() % 7) * 0.1);
            graspLocalSearch(clique, 100);
            
            if (clique.size() > initial_clique.size()) {
                initial_clique = clique;
                if (initial_clique.size() > best_clique.size()) {
                    best_clique = initial_clique;
                }
            }
        }
        
        cout << "GRASP found clique of size: " << initial_clique.size() << endl;
        
        // Этап 2: Tabu Search для улучшения решения
        double remaining_time = time_limit - double(clock() - start) / CLOCKS_PER_SEC;
        if (remaining_time > 0.1 && initial_clique.size() > 0) {
            vector<int> tabu_clique = tabuSearch(initial_clique, 
                min(10, n/20), iterations, remaining_time, start);
            if (tabu_clique.size() > best_clique.size()) {
                best_clique = tabu_clique;
            }
            cout << "Tabu Search found: " << tabu_clique.size() << endl;
        }
        
        // Финальное локальное улучшение
        graspLocalSearch(best_clique, 200);
        
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
    int randomization = 10000;
    
    vector<string> files = { 
        "C125.9.clq", "johnson8-2-4.clq", "johnson16-2-4.clq", "MANN_a9.clq", "MANN_a27.clq",
        "p_hat1000-1.clq", "keller4.clq", "hamming8-4.clq", "brock200_1.clq", "brock200_2.clq", 
        "brock200_3.clq", "brock200_4.clq", "gen200_p0.9_44.clq", "gen200_p0.9_55.clq", 
        "brock400_1.clq", "brock400_2.clq", "brock400_3.clq", "brock400_4.clq",
        "MANN_a45.clq", "sanr400_0.7.clq", "p_hat1000-2.clq", "p_hat500-3.clq", 
        "p_hat1500-1.clq", "p_hat300-3.clq", "san1000.clq", "sanr200_0.9.clq" 
    };
    vector<int> iterations = {700000, 1, 1, 1000, 500000,
         900000, 20000, 20000, 100000, 1000000, 1000000, 1000000,
         500000, 200000, 5000000, 2000000, 3000000, 3000000,
         100000, 4000000, 13000000, 4000000, 13000000, 2000000, 13000000,
         1000000
         };
    
    ofstream fout("clique.csv");
    fout << "File; Clique; Time (sec)\n";
    
    int count = 0;
    for (string file : files)
    {
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
    }
    
    fout.close();
    return 0;
}