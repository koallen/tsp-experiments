/*
 * This can get 20.507495 points
 */
#include <iostream>
#include <cstring>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>

#define MAX 1000
#define TIME_LIMIT 1990000000
#define pf pair<float, float>
#define pi pair<int, int>
#define endl "\n"

using namespace std;

vector<int> neighbor[MAX]; // stores the 2 neighbors of every city
pf points[MAX];            // stores the coordinates of all the cities
int used[MAX] = {0};       // degree of every city
bool in_tour[MAX];         // whether a city is already in the tour, used when constructing the tour
int counter;               // keeping track at how many edges we have found
int parent[MAX];           // for checking whether a cycle is formed
int tour[MAX];             // stores the complete tour
int backup[MAX];           // a backup of the current best tour
int distMatrix[MAX][MAX];  // distance matrix

/*
 * Finds the parent of a city using the union find algorithm
 *
 * @param x the city which you want to get its parent
 */
int find_set(int x)
{
    if (x != parent[x])
        parent[x] = find_set(parent[x]);
    return parent[x];
}

/*
 * Checks whether an edge can be added to the tour fragment
 *
 * @param idxs a pair of indexes containing the two cities of an edge
 * @param N    total number of cities
 */
bool check(pi idxs, int N)
{
    if (used[idxs.first] < 2 && used[idxs.second] < 2)
    {
        int first_parent = find_set(idxs.first), second_parent = find_set(idxs.second);
        if (counter < N)
        {
            if (first_parent != second_parent)
            {
                parent[first_parent] = parent[second_parent];
                ++used[idxs.first];
                ++used[idxs.second];
                ++counter;
                return true;
            }
        } else if (counter == N) {
            if (first_parent == second_parent)
            {
                ++used[idxs.first];
                ++used[idxs.second];
                ++counter;
                return true;
            }
        }
        return false;
    }
    return false;
}

/*
 * Calculates the distance between 2 cities, rounded to the nearest integer
 *
 * @param idxA The index of the first city
 * @param idxB The index of the second city
 */
int dist(int idxA, int idxB)
{
    pf pointA = points[idxA];
    pf pointB = points[idxB];
    return sqrt(pow(pointA.first - pointB.first, 2) + pow(pointA.second - pointB.second, 2)) + 0.5;
}

/*
 * Generates the graph, which contains all the edges (undirected) in non-decreasing order
 *
 * @param N     total number of cities
 * @param graph the vector which contains all the edges (unsorted)
 */
void generateDistances(int N, vector<pair<int, pi> > &graph)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < i; ++j)
            if (i != j)
                graph.push_back(pair<int, pi>(dist(i, j), pi(i, j)));
    sort(graph.begin(), graph.end());
}

/*
 * Generates the whole distance matrix
 *
 * @param N     total number of cities
 */
void generateDistancesMatrix(int N)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < i; ++j)
            if (i != j)
            {
                int distance = dist(i, j);
                distMatrix[i][j] = distance;
                distMatrix[j][i] = distance;
            }
}

/*
 * Generates the final tour starting from city with index 0 and store it in an array
 *
 * @param N total number of cities
 */
void generateTour(int N)
{
    int i = 0, current = 0;
    in_tour[0] = true;
    tour[0] = 0;
    while (i < N)
    {
        tour[i] = current;
        for (int j = 0; j < neighbor[current].size(); ++j)
            if (!in_tour[neighbor[current][j]])
            {
                in_tour[neighbor[current][j]] = true;
                current = neighbor[current][j];
                break;
            }
        ++i;
    }
}

int calculateTour(int *tour, int N)
{
    int tour_distance = 0;
    for (int i = 0; i < N - 1; ++i)
        tour_distance += distMatrix[tour[i]][tour[i+1]];
    tour_distance += distMatrix[tour[N-1]][tour[0]];
    return tour_distance;
}

int deltaEval(int i, int k, int *tour, int current_best, int N)
{
    int new_distance;
    if (k == N - 1)
        new_distance = current_best + distMatrix[tour[i-1]][tour[k]] + distMatrix[tour[i]][tour[0]] - distMatrix[tour[i-1]][tour[i]] - distMatrix[tour[k]][tour[0]];
    else
        new_distance = current_best + distMatrix[tour[i-1]][tour[k]] + distMatrix[tour[i]][tour[k+1]] - distMatrix[tour[i-1]][tour[i]] - distMatrix[tour[k]][tour[k+1]];
    return new_distance;
}

void swap(int a, int b, int *tour)
{
    int temp;
    for (int i = 0; i < (b - a + 1) / 2; ++i)
    {
        temp = tour[a + i];
        tour[a + i] = tour[b - i];
        tour[b - i] = temp;
    }
}

int main()
{
    ios::sync_with_stdio(false); // optimize I/O
    auto t1 = chrono::high_resolution_clock::now(); // start timer

    vector<pair<int, pi> > graph;
    int N;
    float x, y;

    counter = 1;

    cin >> N;
    for (int i = 0; i < N; ++i)
        parent[i] = i;
    for (int i = 0; i < N; ++i)
    {
        cin >> x >> y;
        points[i] = make_pair(x, y);
    }
    generateDistances(N, graph);
    generateDistancesMatrix(N);

    // multi fragment algorithm
    for (int i = 0; i < graph.size(); ++i)
        if (check(graph[i].second, N))
        {
            neighbor[graph[i].second.first].push_back(graph[i].second.second);
            neighbor[graph[i].second.second].push_back(graph[i].second.first);
            if (counter > N) break;
        }

    // construct tour
    generateTour(N);

    // 2-opt local search
    int best_tour_dist = calculateTour(tour, N); // best tour for the 1st round
    int best_tour_dist_2;                        // best tour for the 2nd round
    int new_distance;

START_AGAIN:
    for (int i = 1; i < N - 1; ++i)
        for (int k = i + 1; k < N; ++k)
        {
            auto t2 = chrono::high_resolution_clock::now();
            chrono::duration<int64_t,nano> elapsed = t2 - t1;
            if (elapsed.count() > TIME_LIMIT)
                goto END1;
            new_distance = deltaEval(i, k, tour, best_tour_dist, N);
            if (new_distance < best_tour_dist)
            {
                swap(i, k, tour);
                best_tour_dist = new_distance;
                goto START_AGAIN;
            }
        }

    // backup the current best tour and randomly shuffle the tour to restart
    memcpy(backup, tour, sizeof(int)*N);
    random_shuffle(begin(tour), begin(tour)+N);
    best_tour_dist_2 = calculateTour(tour, N);

SECOND_ROUND:
    for (int i = 1; i < N - 1; ++i)
        for (int k = i + 1; k < N; ++k)
        {
            auto t3 = chrono::high_resolution_clock::now();
            chrono::duration<int64_t,nano> elapsed_2 = t3 - t1;
            if (elapsed_2.count() > TIME_LIMIT)
                goto END2;
            new_distance = deltaEval(i, k, tour, best_tour_dist_2, N);
            if (new_distance < best_tour_dist_2)
            {
                swap(i, k, tour);
                best_tour_dist_2 = new_distance;
                goto SECOND_ROUND;
            }
        }
    goto END2;

END1:
    for (int i = 0; i < N; ++i)
        cout << tour[i] << endl;
    goto END;

END2:
    if (best_tour_dist < best_tour_dist_2)
    {
        for (int i = 0; i < N; ++i)
            cout << backup[i] << endl;
    } else {
        for (int i = 0; i < N; ++i)
            cout << tour[i] << endl;
    }

END:
    return 0;
}
