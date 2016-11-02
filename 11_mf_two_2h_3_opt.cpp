/*
 * This can get 20.507495 points
 */
#include <iostream>
#include <cstring>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <climits>
#include <random>

#define MAX 1000
#define TIME_LIMIT 1990000000
#define pf pair<float, float>
#define pi pair<int, int>
#define endl "\n"

#define timestamp chrono::high_resolution_clock::time_point
#define duration chrono::duration<int64_t,nano>
#define TIMER(t) auto t = chrono::high_resolution_clock::now();
#define DURATION(ela, t1, t2) ela = t2 - t1;

using namespace std;

vector<int> neighbor[MAX]; // stores the 2 neighbors of every city
pf points[MAX];            // stores the coordinates of all the cities
int used[MAX] = {0};       // degree of every city
bool in_tour[MAX];         // whether a city is already in the tour, used when constructing the tour
int counter;               // keeping track at how many edges we have found
int parent[MAX];           // for checking whether a cycle is formed
int distMatrix[MAX][MAX];  // distance matrix
vector<pi> neighbors[MAX]; // sorted K nearest neighbors of each city
int tour[MAX];
int backup[MAX];
int position[MAX];            // keeps track of where a city is
random_device rd;
mt19937 eng(rd());

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
                graph.push_back(pair<int, pi>(distMatrix[i][j], pi(i, j)));
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
                neighbors[i].push_back(pi(distance, j));
                neighbors[j].push_back(pi(distance, i));
            }
    for (int i = 0; i < N; ++i)
        sort(neighbors[i].begin(), neighbors[i].end());
}

/*
 * Generates the final tour starting from city with index 0 and store it in an array
 *
 * @param N total number of cities
 */
void generateTour(int *tour, int N)
{
    int i = 0, current = 0;
    in_tour[0] = true;
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

void reverse(int start, int end, int N, int *tour)
{
    int num = ((start <= end ? end - start : (end + N) - start) + 1) / 2;
    int temp, currentHead = start, currentTail = end;
    for (int i = 0; i < num; ++i)
    {
        // swap 2 elements
        temp = tour[currentHead];
        tour[currentHead] = tour[currentTail];
        tour[currentTail] = temp;
        // update city position
        position[tour[currentHead]] = currentHead;
        position[tour[currentTail]] = currentTail;
        // update head and tail
        currentHead = (currentHead + 1) % N;
        currentTail = ((currentTail + N) - 1) % N;
    }
}

int deltaTwoH(int i, int k, int *tour, int current_best, int N)
{
    int new_distance = current_best - distMatrix[tour[i-1]][tour[i]] - distMatrix[tour[i]][tour[i+1]] - distMatrix[tour[k]][tour[k+1]];
    new_distance += distMatrix[tour[i-1]][tour[i+1]] + distMatrix[tour[k]][tour[i]] + distMatrix[tour[i]][tour[k+1]];
    return new_distance;
}

void swapTwoH(int a, int b, int *tour)
{
    int temp = tour[a];
    for (int i = a; i < b; ++i)
        tour[i] = tour[i+1];
    tour[b] = temp;
}

void updateCityPosition(int *tour, int N)
{
    for (int i = 0; i < N; ++i)
        position[tour[i]] = i;
}

void MF(vector<pair<int, pi> > &graph, int N)
{
    counter = 1;
    // multi fragment algorithm
    for (int i = 0; i < graph.size(); ++i)
        if (check(graph[i].second, N))
        {
            neighbor[graph[i].second.first].push_back(graph[i].second.second);
            neighbor[graph[i].second.second].push_back(graph[i].second.first);
            if (counter > N) break;
        }
}

void init(int &N, int &K, vector<pair<int, pi> > &graph)
{
    float x, y;

    cin >> N;
    K = N / 17;

    for (int i = 0; i < N; ++i)
        parent[i] = i;

    for (int i = 0; i < N; ++i)
    {
        cin >> x >> y;
        points[i] = make_pair(x, y);
    }

    generateDistancesMatrix(N);
    generateDistances(N, graph);
}

void initTour(int N, vector<pair<int, pi> > &graph, int *tour)
{
    MF(graph, N);
    generateTour(tour, N);
    updateCityPosition(tour, N);
}

/*
 * Optimize with 2-opt exchange until local optima
 */
void twoOpt(int *tour, int &bestTourDist, const int N, const int K)
{
    int u, v, x, y;
    int pos_u, pos_v, pos_x, pos_y;
    bool localOptima;

    // 2-opt local search
    do {
        localOptima = true;
        for (int pos_u = 0; pos_u < N; ++pos_u)
        {
            pos_v = (pos_u + 1) % N;
            u = tour[pos_u];
            v = tour[pos_v];
            for (int k = 0; k < K; ++k)
            {
                pos_x = position[neighbors[u][k].second];
                pos_y = (pos_x + 1) % N;
                x = tour[pos_x];
                y = tour[pos_y];

                if (v == x || y == u)
                    continue;

                // TODO: minimize runtime with min and max

                if (distMatrix[u][x] + distMatrix[v][y] < distMatrix[u][v] + distMatrix[x][y])
                {
                    bestTourDist = deltaEval(pos_u, pos_x, tour, bestTourDist, N);
                    reverse(pos_v, pos_x, N, tour);
                    localOptima = false;
                    break;
                }
            }
        }
    } while (!localOptima);
}

void twoHOpt(int *tour, int &best_tour_dist, int N)
{
    // TODO: look into faster implementation of 2.5 opt
    int new_distance;
    // 2.5 opt
    for (int i = 1; i < N - 3; ++i)
        for (int k = i + 2; k < N - 1; ++k)
        {
            new_distance = deltaTwoH(i, k, tour, best_tour_dist, N);
            if (new_distance < best_tour_dist)
            {
                swapTwoH(i, k, tour);
                best_tour_dist = new_distance;
                break;
            }
        }
}

void orderEdges(int &A, int &B, int &C, int &D, int &E, int &F,
                int U, int V, int W, int X, int Y, int Z)
{
    E = Y;
    F = Z;

    if ((W < U && U < Y) || (Y < W && W < U) || (U < Y && Y < W))
    {
        A = W;
        B = X;
        C = U;
        D = V;
    } else {
        A = U;
        B = V;
        C = W;
        D = X;
    }
}

/*
 * Optimize with 3-opt exchange until local optima
 */
void threeOpt(int *tour, int &bestTourDist, const int N, const int K, timestamp t1)
{
    int u, v, a, b, x, y;
    int pos_u, pos_v, pos_a, pos_b, pos_x, pos_y;
    int pos_1, pos_2, pos_3, pos_4, pos_5, pos_6;
    int one, two, three, four, five, six;
    bool localOptima;
    duration e1;

    // 3-opt local search
    do {
        localOptima = true;
        for (int pos_u = 0; pos_u < N; ++pos_u)
        {
            TIMER(t2)
            DURATION(e1, t1, t2)
            if (e1.count() > TIME_LIMIT)
                return;

            pos_v = (pos_u + 1) % N;
            u = tour[pos_u];
            v = tour[pos_v];

            for (int i = 0; i < K; ++i) // find b that's closer to u
            {
                pos_b = position[neighbors[u][i].second];
                pos_a = (pos_b + N - 1) % N;
                a = tour[pos_a];
                b = tour[pos_b];

                // TODO: use max and min to reduce runtime
                if (v == a || u == a)
                    continue;

                for (int j = 0; j < K; ++j)
                {
                    pos_y = position[neighbors[v][j].second];
                    pos_x = (pos_y + N - 1) % N;
                    x = tour[pos_x];
                    y = tour[pos_y];

                    // TODO: use max and min to reduce runtime
                    if (x == b || y == b || y == a || x == u || x == v)
                        continue;

                    int dist = distMatrix[u][v] + distMatrix[a][b] + distMatrix[x][y];
                    orderEdges(pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, pos_u, pos_v, pos_a, pos_b, pos_x, pos_y);
                    one = tour[pos_1]; two = tour[pos_2]; three = tour[pos_3];
                    four = tour[pos_4]; five = tour[pos_5]; six = tour[pos_6];
                    if (distMatrix[one][four] + distMatrix[two][six] + distMatrix[three][five] < dist)
                    {
                        reverse(pos_6, pos_1, N, tour);
                        reverse(pos_4, pos_5, N, tour);
                        localOptima = false;
                        goto NEXT_ROUND;
                    } else if (distMatrix[two][four] + distMatrix[five][one] + distMatrix[six][three] < dist) {
                        reverse(pos_6, pos_1, N, tour);
                        reverse(pos_2, pos_3, N, tour);
                        localOptima = false;
                        goto NEXT_ROUND;
                    } else if (distMatrix[one][three] + distMatrix[two][five] + distMatrix[four][six] < dist) {
                        reverse(pos_2, pos_3, N, tour);
                        reverse(pos_4, pos_5, N, tour);
                        localOptima = false;
                        goto NEXT_ROUND;
                    } else if (distMatrix[two][five] + distMatrix[four][one] + distMatrix[six][three] < dist) {
                        reverse(pos_1, pos_6, N, tour);
                        reverse(pos_2, pos_3, N, tour);
                        reverse(pos_4, pos_5, N, tour);
                        localOptima = false;
                        goto NEXT_ROUND;
                    }
                }
            }
NEXT_ROUND:
            continue;
        }
    } while (!localOptima);
}

/*
 * Uses local search to optimize TSP until time runs our
 */
void TSP(int *tour, const int N, const int K, timestamp &t1)
{
    int bestTourDist; // best tour for the 1st round
    int bestTourDist2 = INT_MAX;
    bool backedup = false;
    int start, finish;
    uniform_int_distribution<> rand_N(N / 4, N / 3);

    // do optimization
    twoOpt(tour, bestTourDist, N, K);
    twoHOpt(tour, bestTourDist, N);
    threeOpt(tour, bestTourDist, N, K, t1);
    bestTourDist = calculateTour(tour, N);

    // local retour
    memcpy(backup, tour, sizeof(int) * N);
    backedup = true;
    TIMER(t2)
    duration e1;
    DURATION(e1, t1, t2)
    while (e1.count() < TIME_LIMIT)
    {
        start = rand_N(eng);
        finish = rand_N(eng);
        random_shuffle(tour + start, tour + N - finish);
        //random_shuffle(begin(tour), end(tour));
        twoOpt(tour, bestTourDist2, N, K);
        twoHOpt(tour, bestTourDist2, N);
        threeOpt(tour, bestTourDist2, N, K, t1);
        bestTourDist2 = calculateTour(tour, N);
        if (bestTourDist2 < bestTourDist)
        {
            memcpy(backup, tour, sizeof(int) * N);
            bestTourDist = bestTourDist2;
        }
        TIMER(t2)
        DURATION(e1, t1, t2)
    }

    // print tour
    if (backedup)
        for (int i = 0; i < N; ++i)
            cout << backup[i] << endl;
    else
        for (int i = 0; i < N; ++i)
            cout << tour[i] << endl;
}

int main()
{
    // initialization
    TIMER(t1);
    int N, K;
    vector<pair<int, pi> > graph;
    init(N, K, graph);

    // run TSP
    initTour(N, graph, tour);
    TSP(tour, N, K, t1);

    return 0;
}
