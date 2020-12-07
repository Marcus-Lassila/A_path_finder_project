/* 
This program solves an ice maze of in the format:

map = "  #  x # \n
           #   #\n
       #  #  S #\n
             #  \n
        E  #    "

S: Starting point (slippery ground).
E: Exit point (non-slippery ground).
 : Slippery ground.
x: Non-slippery ground.
#: An obstacle which we cannot occupy or pass through.

The map is assumed to be surronded by walls, and we cannot go outside of the map.

We start at the starting point S and we should reach the exit point in as few steps as possible.
We can move in the four directions up, down, left and right, but if we move to slippery ground, 
we keep sliding in that direction until we hit a wall, an obstacle or non-slippery ground. 

The program uses Dijkstra's algorithm to find the shortest path from the start to exit in terms of
steps taken not accounting for the distance traveled by sliding. If there are multiple shortest 
paths in terms of steps taken, we tiebreak by selecting the path with the overall shortest distance.
If there are still are multiple shortest paths we just select one of them.

The program then prints how the character denoted by '&' travels from S to E in the maze, step by step
in a sort of animation in the terminal, having a poor frame rate. The program is made for a windows 
operating system, and the size of the terminal have to be adjusted manually.

The first map below is a replica of the original ice maze found in a Pokemon game, and the second is a randomly
generated map. 
*/

#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <windows.h>

using namespace std;

struct Vertex { 
    /* Structure used in the priority queue for Dijkstra's algorithm. */

    int steps; // Number of steps to get to vertex from a given initial vertex.
    int distance; // Distance to vertex from a given initial vertex.
    pair<int,int> position; // 2D coordinates of vertex.
    Vertex(int s , int d, pair<int,int> pos) : steps(s), distance(d), position(pos){}
    bool operator < (const struct Vertex& other_vertex) const { //Overwrite "less than" operator to compare steps.
        return (steps < other_vertex.steps);
    }
};

bool is_valid_position(int i, int j, const vector<vector<char>> &grid) {
    /* Return true if position (i, j) can be occupied, false otherwise. */
    int M = grid.size();
    int N = grid[0].size();
    return ((i >= 0) && (i < M) && (j >= 0) && (j < N) && (grid[i][j] != '#'));
}

bool is_slipery(int i, int j, const vector<vector<char>> &grid) {
    /* Return true if position (i, j) is slippery, false otherwise. */
    return (grid[i][j] == ' ') || (grid[i][j] == 'S');
}

vector<pair<int,int>> neighbours(const pair<int,int> &pos, const vector<vector<char>> &grid) {
    /* Return the vertices connected to pos by one step, assuming we may slide some distance. */
    vector<pair<int,int>> neighs; // Vector of neighbours to be returned.
    int i = pos.first;
    int j = pos.second;
    int k = 1;
    if (is_valid_position(i-k, j, grid)) { // If we may go up, there is a neighbouring vertex in that direction.
        pair<int,int> u;
        if (is_slipery(i-k, j, grid)) { // If we slide, count how long we slide.
            while (is_valid_position(i-k-1, j, grid) && is_slipery(i-k-1, j, grid)) {
                k++;
            }
            if (is_valid_position(i-k-1, j, grid))
                u = make_pair(i-k-1,j);
            else
                u = make_pair(i-k,j);
        }
        else { // If we do not slide, the neighbouring vertex is just one unit distance above.
            u = make_pair(i-1,j);
        }
        neighs.push_back(u); // Add the neighbour to neighs.
        k = 1;
    }
    if (is_valid_position(i+k, j, grid)) { // Find the neighbouring vertex in the downwards direction.
        pair<int,int> d;
        if (is_slipery(i+k, j, grid)) { // If we slide, count how long we slide.
            while (is_valid_position(i+k+1, j, grid) && is_slipery(i+k+1, j, grid)) {
                k++;
            }
            if (is_valid_position(i+k+1, j, grid))
                d = make_pair(i+k+1,j);
            else
                d = make_pair(i+k,j);
        }
        else {  // If we do not slide, the neighbouring vertex is just one unit distance below.
            d = make_pair(i+1,j);
        }
        neighs.push_back(d); // Add the neighbour to neighs.
        k = 1;
    }
    if (is_valid_position(i, j-k, grid)) { // Find the neighbouring vertex in the left direction.
        pair<int,int> l;
        if (is_slipery(i, j-k, grid)) { // If we slide, count how long we slide.
            while (is_valid_position(i, j-k-1, grid) && is_slipery(i, j-k-1, grid)) {
                k++;
            }
            if (is_valid_position(i, j-k-1, grid))
                l = make_pair(i,j-k-1);
            else
                l = make_pair(i,j-k);
        }
        else { // If we do not slide, the neighbouring vertex is just one unit distance to the left.
            l = make_pair(i,j-1);
        }
        neighs.push_back(l); // Add the neighbour to neighs.
        k = 1;
    }
    if (is_valid_position(i, j+k, grid)) { // Find the neighbouring vertex in the right direction.
        pair<int,int> r;
        if (is_slipery(i, j+k, grid)) { // If we slide, count how long we slide.
            while (is_valid_position(i, j+k+1, grid) && is_slipery(i, j+k+1, grid)) {
                k++;
            }
            if (is_valid_position(i, j+k+1, grid))
                r = make_pair(i,j+k+1);
            else
                r = make_pair(i,j+k);
        }
        else { // If we do not slide, the neighbouring vertex is just one unit distance to the right.
            r = make_pair(i,j+1);
        }
        neighs.push_back(r); // Add the neighbour to neighs.
        k = 1;
    }
    return neighs;
}

vector<pair<int,int>> maze_solver(const string &map) {
    /* Solve ice maze using Dijkstra's algorithm and represent each position in the maze 
    as a vertex defined in struct Vertex. Takes an ice maze map in the form of a string. */
    
    vector<vector<char>> grid; // To store map as a 2D vector of chars for convinience.
    vector<char> row;
    
    pair<int,int> start; // Start coordinates.
    pair<int,int> end; // End coordinates.

    int i = 0, j = -1;
    for (auto x : map) {
        j++;
        if (x != '\n') {
            row.push_back(x);
            if (x == 'S') start = make_pair(i,j); // This is the start position.
            if (x == 'E') end = make_pair(i,j); // This is the exit position.
        }
        else {
            grid.push_back(row); // Add row to grid.
            row.clear(); // Clear before filling up the next row.
            i++;
            j = -1;
        }
        
    }
    grid.push_back(row); // Add the last row.

    int M = grid.size(); // Number of rows.
    int N = grid[0].size(); // Number of columns.
    int inf = numeric_limits<int>::max(); // "Infinity" for Dijkstra's algorithm.

    vector<vector<int>> steps(M, vector<int>(N, inf)); // steps[i][j] holds least number of steps from start to (i,j).
    steps[start.first][start.second] = 0;

    vector<vector<int>> dist(M, vector<int>(N, inf)); // dist[i][j] holds distance from start to (i,j).
    dist[start.first][start.second] = 0;

    priority_queue<Vertex> pq; // Priority queue for Dijkstra's algorithm.
    pq.push(Vertex(0, 0, start));

    pair<int,int> successor[M][N]; // Map to store connection between points in the path.
    successor[start.first][start.second] = make_pair(-1,-1);

    vector<pair<int,int>> neighs; 
    
    while (!pq.empty()) {
        int s = pq.top().steps; // Number of steps (u,d,l,r) from start.
        int d = pq.top().distance; // Distance from start.
        pair<int,int> pos = pq.top().position; // Current position.
        pq.pop();

        neighs = neighbours(pos, grid); // Neighbouring vertices.
        for (auto v : neighs) {
            int w = abs(v.first - pos.first) + abs(v.second - pos.second);
            if (steps[v.first][v.second] > s + 1) { // There is a way to v using less steps.
                steps[v.first][v.second] = s + 1; // Relax least number of steps to v.
                dist[v.first][v.second] = d + w; // Relax distance to v.
                successor[v.first][v.second] = pos; // Update connection to previous point.
                pq.push(Vertex(s+1, d+w, v)); // Add vertex to priority queue.
            }
            else {
                if (steps[v.first][v.second] == s + 1 && dist[v.first][v.second] > d + w) {
                    // There is way to v using the same number of steps but having shorter distance.
                    dist[v.first][v.second] = d + w; // Relax distance to v.
                    successor[v.first][v.second] = pos; // Update connection to previous point.
                    pq.push(Vertex(s+1, d+w, v)); // Add vertex to priority queue.
                }
            }
        }
    }
    if (steps[end.first][end.second] == inf) return {}; // Exit could not be reached, return empty vector. 
  
    // Find path from successor array, going backwards from the end.
    vector<pair<int,int>> path;
    pair<int,int> pos = end, prev, tmp;
    int dx, dy;
    while ((pos != start)) { // Until we reach the start.
        prev = successor[pos.first][pos.second];
        dx = pos.second - prev.second;
        dy = pos.first - prev.first;
        if (dx == 0) { // We went up or down.
            tmp.second = pos.second;
            if (dy > 0) { // We went down.
                for (int i = 0; i < dy; i++) { // Go down and add visited coordinates to path.
                    tmp.first = pos.first - i;
                    path.push_back(tmp);
                }
            }
            else {  // We went up.
                for (int i = 0; i > dy; i--) { // Go up and add visited coordinates to path.
                    tmp.first = pos.first - i;
                    path.push_back(tmp);
                }
            }
        }
        else { // We went left or right.
            tmp.first = pos.first;
            if (dx > 0) { // We went right.
                for (int i = 0; i < dx; i++) { // Go right and add visited coordinates to path.
                    tmp.second = pos.second - i;
                    path.push_back(tmp);
                }
            }
            else { // We went left.
                for (int i = 0; i > dx; i--) { // Go left and add visited coordinates to path.
                    tmp.second = pos.second - i;
                    path.push_back(tmp);
                }
            }
        }
        pos = prev; // Repeat for the previous vertex.
    }
    path.push_back(start); // Add the start to the path.
    reverse(path.begin(), path.end()); // Since we went backwards from end to start, reverse path.
    return path;
}

int flatten(int M, int N, pair<int,int> &pos) {
    /* Return the "flattened" index in the string map corresponding to the coordinate pos.
    M and N are the dimensions of map. */
    
    // If pos is outside of map, return -1.
    if (0 > pos.first || M <= pos.first || 0 > pos.second || N <= pos.second) 
        return -1;
    return (N+1)*pos.first + pos.second;
}

void print_path(const vector<pair<int,int>> &path, string &map) {
    /* Prints the path taken in the ice maze as an animation in the terminal, having a poor frame rate... */

    int N = 0; 
    while (map[N] != '\n') { // Find number of columns in map.
        N++;
    }
    int M = 1;
    for (int i = 0; i < map.size(); i++) { // Find number of rows in map.
        if (map[i] == '\n') M++;
    }

    int i;
    int j = -1; 
    char hold; // To hold the character of the position we visit, in order to reset it in the next frame.
    cout << string(20, '\n'); // Print some empty lines to clear terminal.
    cout << map; // Print map.
    Sleep(2000); // Paus for two seconds.
    cout << string(20, '\n'); // Print some empty lines to clear terminal.
    for (auto x : path) {
        if (j != -1) // If there is something to hold.
            map[j] = hold; 
        i = flatten(M, N, x); // Convert coordinate in path to flattened string map index.
        hold = map[i]; // Hold character.
        map[i] = '&'; // Our position in the maze is represented by &.
        j = i;
        cout << map;
        Sleep(400); // Paus for 0.4 seconds.
        cout << string(20, '\n'); // Print some empty lines to clear terminal.
        if (hold == 'E') { // If we have exited the maze, print maze once again.
            map[i] = 'E';
            cout << map;
        }
    }
}

int main() {

    // Define some ice maze maps.
    string map1 = "        #     \n"
                  "   #          \n"
                  "         #    \n"
                  " #            \n"
                  "#       #     \n"
                  "             #\n"
                  "      #       \n"
                  "  #          E\n"
                  "             #\n"
                  "       #      \n"
                  "     #   #    \n"
                  "              \n"
                  "#############S";

    string map2 = "    ##  # #   #     # # x   #  ##\n"
                  "      ##      #   #####  #    #  \n"
                  "             x#   x ##  #     x  \n"
                  "    #      x     #   ### #  #  S#\n"
                  "         E# # x  # #  x  #       \n"
                  "#  # # #        #  #         x # \n"
                  "      #    ##      #        #    \n"
                  " #  #           #      #         \n"
                  "## x #     #  #     #       # #  \n"
                  "#     #     #     x#    #  #  #  \n"
                  "          #      ##    x    #    \n"
                  "            #   ##    #  x x     \n"
                  "  #x   #  # ##   ##   x  x##  #  \n"
                  " x#     x ##  # # #   # #  ##    \n"
                  "#   x         #          x##   ##\n"
                  "  ##  x  xx #         x  #       ";

    // Solve the ice maze and find the minimum step path.
    vector<pair<int,int>> path1 = maze_solver(map1);
    vector<pair<int,int>> path2 = maze_solver(map2);

    // Print the solution as a continous sequence of steps in the terminal.
    if (path1.size() == 0 && path2.size() == 0)
        cout << "Neither of the mazes could be solved!\n";
    if (path1.size() != 0) {
        print_path(path1, map1);
        Sleep(3000);
    }
    if (path2.size() != 0)
        print_path(path2, map2);
    return 0;
}