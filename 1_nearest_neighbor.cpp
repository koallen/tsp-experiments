/*
 * This can get 3.051487 points
 */
#include <iostream>
#include <cstring>
#include <cmath>

#define MAX 1000
#define pf pair<float, float>

using namespace std;

pf points[MAX];
bool used[MAX];
int tour[MAX];

int dist(int idxA, int idxB)
{
    pf pointA = points[idxA];
    pf pointB = points[idxB];
    return sqrt(pow(pointA.first - pointB.first, 2) + pow(pointA.second - pointB.second, 2));
}

int main()
{
    int N, best;
    float x, y;

    cin >> N;
    for (int i = 0; i < N; ++i)
    {
        cin >> x >> y;
        points[i] = make_pair(x, y);
    }

    // NN algorithm
    tour[0] = 0;
    used[0] = true;
    for (int i = 1; i < N; ++i)
    {
        best = -1;
        for (int j = 0; j < N; ++j)
        {
            if (!used[j] && (best == -1 || dist(tour[i-1], j) < dist(tour[i-1], best)))
                best = j;
        }
        tour[i] = best;
        used[best] = true;
    }

    for (int i = 0; i < N; ++i)
        cout << tour[i] << endl;
}
