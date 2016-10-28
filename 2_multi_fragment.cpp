/*
 * This can get 5.793946 points
 */
#include <iostream>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

#define MAX 1000
#define pf pair<float, float>
#define pi pair<int, int>

using namespace std;

vector<int> neighbor[MAX];
pf points[MAX];
int used[MAX] = {0};
bool tour[MAX];
int counter;
int parent[MAX];

int find_set(int x)
{
    if (x != parent[x])
        parent[x] = find_set(parent[x]);
    return parent[x];
}

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

int dist(int idxA, int idxB)
{
    pf pointA = points[idxA];
    pf pointB = points[idxB];
    return sqrt(pow(pointA.first - pointB.first, 2) + pow(pointA.second - pointB.second, 2)) + 0.5;
}

void generateDistances(int N, vector<pair<int, pi> > &graph)
{
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < i; ++j)
            if (i != j)
                graph.push_back(pair<int, pi>(dist(i, j), pi(i, j)));
    sort(graph.begin(), graph.end());
}

void generateTour(int N)
{
    int i = 0, current = 0;
    tour[0] = true;
    while (i < N)
    {
        cout << current << endl;
        for (int j = 0; j < neighbor[current].size(); ++j)
            if (!tour[neighbor[current][j]])
            {
                tour[neighbor[current][j]] = true;
                current = neighbor[current][j];
                break;
            }
        ++i;
    }
}

int main()
{
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
}
