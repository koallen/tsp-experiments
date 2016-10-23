/*
 * This can get 1.212971 points
 */
#include <iostream>

#define MAX 1000

using namespace std;

// the distance matrix
int dist[MAX][MAX];

int main()
{
    int N;
    float x, y;

    cin >> N;
    for (int i = 0; i < N; ++i)
        cin >> x >> y;

    for (int i = 0; i < N; ++i)
        cout << i << endl;
}
