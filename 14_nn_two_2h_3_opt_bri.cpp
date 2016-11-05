/*
 * This can get 20.507495 points
 */
#include <iostream>
#include <cstring>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>

#define MAX 1000
#define TIME_LIMIT 1985000000
#define pf pair<float, float>
#define pi pair<int, int>
#define endl "\n"

#define timestamp chrono::high_resolution_clock::time_point
#define duration chrono::duration<int64_t,nano>
#define TIMER(t) auto t = chrono::high_resolution_clock::now();
#define DURATION(ela, t1, t2) ela = t2 - t1;

using namespace std;

pf points[MAX];            // stores the coordinates of all the cities
bool used[MAX];       // degree of every city
int distMatrix[MAX][MAX];  // distance matrix
vector<pi> neighbors[MAX]; // sorted K nearest neighbors of each city
uint16_t tour[MAX];
uint16_t backup[MAX];
int position[MAX];            // keeps track of where a city is
random_device rd;
mt19937 eng(rd());

/*
 * Calculates the distance between 2 cities, rounded to the nearest integer
 *
 * @param idxA The index of the first city
 * @param idxB The index of the second city
 */
inline int dist(int idxA, int idxB)
{
    pf pointA = points[idxA];
    pf pointB = points[idxB];
    return sqrt(pow(pointA.first - pointB.first, 2) + pow(pointA.second - pointB.second, 2)) + 0.5;
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

int calculateTour(uint16_t *tour, int N)
{
    int tour_distance = 0;
    for (int i = 0; i < N - 1; ++i)
        tour_distance += distMatrix[tour[i]][tour[i+1]];
    tour_distance += distMatrix[tour[N-1]][tour[0]];
    return tour_distance;
}

inline void reverse(int start, int end, int N, uint16_t *tour)
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

inline void swapTwoH(int a, int b, uint16_t *tour)
{
    int temp = tour[a];
    for (int i = a; i < b; ++i)
        tour[i] = tour[i+1];
    tour[b] = temp;
}

inline void updateCityPosition(uint16_t *tour, int N, int &max)
{
    for (int i = 0; i < N; ++i)
    {
        max = ::max(max, distMatrix[tour[i]][tour[(i+1)%N]]);
        position[tour[i]] = i;
    }
}

inline void init(int &N, int &K)
{
    float x, y;

    cin >> N;
    K = N / 17;

    for (int i = 0; i < N; ++i)
    {
        cin >> x >> y;
        points[i] = make_pair(x, y);
    }

    generateDistancesMatrix(N);
}

void NN(uint16_t *tour, int N)
{
    int best;
    tour[0] = 0;
    used[0] = true;
    for (int i = 1; i < N; ++i)
    {
        best = -1;
        for (int j = 0; j < N; ++j)
        {
            if (!used[j] && (best == -1 || distMatrix[tour[i-1]][j] < distMatrix[tour[i-1]][best]))
                best = j;
        }
        tour[i] = best;
        used[best] = true;
    }
}

void findMin(int &min, int N)
{
    for (int i = 0; i < N; ++i)
        min = ::min(neighbors[i][0].first, min);
}

inline void initTour(int N, uint16_t *tour, int &max, int &min)
{
    findMin(min, N);
    NN(tour, N);
    updateCityPosition(tour, N, max);
}

/*
 * Optimize with 2-opt exchange until local optima
 */
void twoOpt(uint16_t *tour, const int N, const int K, timestamp &t1, int &max, int &min)
{
    int u, v, x, y;
    int pos_u, pos_v, pos_x, pos_y;
    bool localOptima;
    duration e1;

    // 2-opt local search
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
            for (int k = 0; k < K; ++k)
            {
                pos_x = position[neighbors[u][k].second];
                pos_y = (pos_x + 1) % N;
                x = tour[pos_x];
                y = tour[pos_y];

                if (v == x || y == u)
                    continue;

                if (distMatrix[u][x] + min > distMatrix[u][v] + max)
                    break;

                if (distMatrix[u][x] + distMatrix[v][y] < distMatrix[u][v] + distMatrix[x][y])
                {
                    reverse(pos_v, pos_x, N, tour);
                    max = ::max(max, ::max(distMatrix[u][x], distMatrix[v][y]));
                    localOptima = false;
                    break;
                }
            }
        }
    } while (!localOptima);
}

void twoHOpt(uint16_t *tour, int N, timestamp &t1)
{
    // TODO: look into faster implementation of 2.5 opt
    int u, v, x, y, z;
    duration e1;
    bool localOptima;
    // 2.5 opt
    do {
        localOptima = true;
        for (int i = 1; i < N - 3; ++i)
        {
            TIMER(t2)
            DURATION(e1, t1, t2)
            if (e1.count() > TIME_LIMIT)
                return;

            x = tour[i-1];
            y = tour[i];
            z = tour[i+1];

            for (int k = i + 2; k < N - 1; ++k)
            {
                u = tour[k];
                v = tour[k+1];
                if (distMatrix[x][y] + distMatrix[y][z] + distMatrix[u][v] > distMatrix[x][z] + distMatrix[u][y] + distMatrix[y][v])
                {
                    swapTwoH(i, k, tour);
                    localOptima = false;
                    break;
                }
            }
        }
    } while (!localOptima);
}

inline void orderEdges(int &A, int &B, int &C, int &D, int &E, int &F,
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
void threeOpt(uint16_t *tour, const int N, const int K, timestamp t1, int &max, int &min)
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

                if (v == a || u == a)
                    continue;

                if (distMatrix[u][b] + 2 * min > distMatrix[u][v] + 2 * max)
                    break;

                if (distMatrix[u][b] + 2 * min > distMatrix[u][v] + distMatrix[a][b] + max)
                    continue;

                for (int j = 0; j < K; ++j)
                {
                    pos_y = position[neighbors[v][j].second];
                    pos_x = (pos_y + N - 1) % N;
                    x = tour[pos_x];
                    y = tour[pos_y];

                    if (x == b || y == b || y == a || x == u || x == v)
                        continue;

                    if (distMatrix[u][b] + distMatrix[v][y] + min > distMatrix[u][v] + distMatrix[a][b] + max)
                        break;

                    int dist = distMatrix[u][v] + distMatrix[a][b] + distMatrix[x][y];
                    orderEdges(pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, pos_u, pos_v, pos_a, pos_b, pos_x, pos_y);
                    one = tour[pos_1]; two = tour[pos_2]; three = tour[pos_3];
                    four = tour[pos_4]; five = tour[pos_5]; six = tour[pos_6];
                    if (distMatrix[one][four] + distMatrix[two][six] + distMatrix[three][five] < dist)
                    {
                        reverse(pos_6, pos_1, N, tour);
                        reverse(pos_4, pos_5, N, tour);
                        localOptima = false;
                        max = ::max(max, ::max(distMatrix[one][four], ::max(distMatrix[two][six], distMatrix[three][five])));
                        goto NEXT_ROUND;
                    } else if (distMatrix[two][four] + distMatrix[five][one] + distMatrix[six][three] < dist) {
                        reverse(pos_6, pos_1, N, tour);
                        reverse(pos_2, pos_3, N, tour);
                        localOptima = false;
                        max = ::max(max, ::max(distMatrix[two][four], ::max(distMatrix[five][one], distMatrix[six][three])));
                        goto NEXT_ROUND;
                    } else if (distMatrix[one][three] + distMatrix[two][five] + distMatrix[four][six] < dist) {
                        reverse(pos_2, pos_3, N, tour);
                        reverse(pos_4, pos_5, N, tour);
                        localOptima = false;
                        max = ::max(max, ::max(distMatrix[one][three], ::max(distMatrix[two][five], distMatrix[four][six])));
                        goto NEXT_ROUND;
                    } else if (distMatrix[two][five] + distMatrix[four][one] + distMatrix[six][three] < dist) {
                        reverse(pos_1, pos_6, N, tour);
                        reverse(pos_2, pos_3, N, tour);
                        reverse(pos_4, pos_5, N, tour);
                        localOptima = false;
                        max = ::max(max, ::max(distMatrix[two][five], ::max(distMatrix[four][one], distMatrix[six][three])));
                        goto NEXT_ROUND;
                    }
                }
            }
NEXT_ROUND:
            continue;
        }
    } while (!localOptima);
}

inline void doubleBridge(uint16_t *tour, const int N)
{
    uniform_int_distribution<int> pickCity(1, N/4);
    int temp, swap1, swap2;
    int A = pickCity(eng);
    int AtoB = pickCity(eng);
    int B = A + AtoB;
    int C = N - AtoB;

    for (int i = 0; i < AtoB; ++i)
    {
        swap1 = A + i;
        swap2 = C + i;
        temp = tour[swap1];
        tour[swap1] = tour[swap2];
        tour[swap2] = temp;
    }
}

/*
 * Uses local search to optimize TSP until time runs our
 */
void TSP(uint16_t *tour, const int N, const int K, timestamp &t1, int &max, int &min)
{
    int bestTourDist; // best tour for the 1st round
    int bestTourDist2;
    bool backedup = false;

    // do optimization
    twoOpt(tour, N, K, t1, max, min);
    twoHOpt(tour, N, t1);
    threeOpt(tour, N, K, t1, max, min);
    bestTourDist = calculateTour(tour, N);

    // local retour
    memcpy(backup, tour, sizeof(uint16_t) * N);
    backedup = true;
    TIMER(t2)
    duration e1;
    DURATION(e1, t1, t2)
    while (e1.count() < TIME_LIMIT)
    {
        // purturbation
        if (N >= 8)
            doubleBridge(tour, N);
        else
            random_shuffle(tour, tour+N);
        // update max
        for (int i = 0; i < N; ++i) {
            max = ::max(max, distMatrix[tour[i]][tour[(i + 1) % N]]);
        }
        // rerun local search
        twoOpt(tour, N, K, t1, max, min);
        twoHOpt(tour, N, t1);
        threeOpt(tour, N, K, t1, max, min);
        bestTourDist2 = calculateTour(tour, N);
        if (bestTourDist2 < bestTourDist)
        {
            memcpy(backup, tour, sizeof(uint16_t) * N);
            bestTourDist = bestTourDist2;
        } else {
            memcpy(tour, backup, sizeof(uint16_t) * N);
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
    int N, K, max = 0, min;
    init(N, K);

    // run TSP
    initTour(N, tour, max, min);
    TSP(tour, N, K, t1, max, min);

    return 0;
}
