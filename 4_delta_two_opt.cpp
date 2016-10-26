/*
 * This can get 18.769668 points in 1.99s
 */
#include <iostream>
#include <cstring>
#include <cmath>
#include <climits>
#include <chrono>

#define MAX 1000
#define TIME_LIMIT 1990000000
#define pf pair<float, float>

using namespace std;

pf points[MAX];
bool used[MAX];
int tour[MAX];
int new_tour[MAX];
int graph[MAX][MAX];

int dist(int idxA, int idxB)
{
    pf pointA = points[idxA];
    pf pointB = points[idxB];
    return sqrt(pow(pointA.first - pointB.first, 2) + pow(pointA.second - pointB.second, 2)) + 0.5;
}

void generateDistances(int N)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < i; ++j)
            if (i != j)
            {
                int distance = dist(i, j);
                graph[i][j] = distance;
                graph[j][i] = distance;
            }
}

int calculateTour(int *tour, int N)
{
    int tour_distance = 0;
    for (int i = 0; i < N - 1; ++i)
        tour_distance += graph[tour[i]][tour[i+1]];
    tour_distance += graph[tour[N-1]][tour[0]];
    return tour_distance;
}

int deltaEval(int i, int k, int *tour, int current_best, int N)
{
    int new_distance;
    if (k == N - 1)
        new_distance = current_best + graph[tour[i-1]][tour[k]] + graph[tour[i]][tour[0]] - graph[tour[i-1]][tour[i]] - graph[tour[k]][tour[0]];
    else
        new_distance = current_best + graph[tour[i-1]][tour[k]] + graph[tour[i]][tour[k+1]] - graph[tour[i-1]][tour[i]] - graph[tour[k]][tour[k+1]];
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
    auto t1 = chrono::high_resolution_clock::now();

    int N, best;
    float x, y;

    cin >> N;
    for (int i = 0; i < N; ++i)
    {
        cin >> x >> y;
        points[i] = make_pair(x, y);
    }

    generateDistances(N);

    // NN algorithm
    tour[0] = 0;
    used[0] = true;
    for (int i = 1; i < N; ++i)
    {
        best = -1;
        for (int j = 0; j < N; ++j)
            if (!used[j] && (best == -1 || dist(tour[i-1], j) < dist(tour[i-1], best)))
                best = j;
        tour[i] = best;
        used[best] = true;
    }

    int best_tour_dist = calculateTour(tour, N);
    int new_distance;

START_AGAIN:
    for (int i = 1; i < N - 1; ++i)
        for (int k = i + 1; k < N; ++k)
        {
            auto t2 = chrono::high_resolution_clock::now();
            chrono::duration<int64_t,nano> elapsed = t2 - t1;
            if (elapsed.count() > TIME_LIMIT)
                goto END;
            //swap(i, k, N, tour, new_tour);
            //new_distance = calculateTour(new_tour, N);
            //cout << new_distance << endl;
            new_distance = deltaEval(i, k, tour, best_tour_dist, N);
            //cout << new_distance << endl;
            if (new_distance < best_tour_dist)
            {
                swap(i, k, tour);
                //memcpy(tour, new_tour, sizeof(int) * N);
                best_tour_dist = new_distance;
                goto START_AGAIN;
            }
        }

END:
    for (int i = 0; i < N; ++i)
        cout << tour[i] << endl;
}

