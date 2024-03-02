#include "DG.h"
#include <iostream>
#include <ctime>

int main()
{
    int N_int; // number of interior faces
    int N_bound; // number of boundary faces
    int N_elem; // number of elements 
    int N_curv; // number of curved elements
    int N_nodes; // number of nodes
    string files[6]; //name of the files containing mesh and matrices

    int p = 0;
    int mesh = 0;
    int max_iter= 130000;
    bool free_stream = false;
    double tolerance = 1.0e-7;
    double CFL=0.8;
    if(mesh==0){
        N_int = 109; 
        N_bound = 28; 
        N_elem = 82;
        N_curv = 21;
        N_nodes =  203;
        files[0] ="bump0_curv";
        files[1] = "../curve_mesh/I2E0.txt";
        files[2] = "../curve_mesh/B2E0.txt"; 
        files[3] = "../curve_mesh/Bn0.txt";
        files[4] = "../curve_mesh/In0.txt";
        files[5] = "../curve_mesh/Area0.txt";
    }
    clock_t begin = clock();
    solve_GD(max_iter, free_stream, tolerance, N_nodes, N_elem, N_curv, N_int, N_bound, files, p, CFL);
    clock_t end = clock();
    cout << "Total time: " << double(end-begin)/CLOCKS_PER_SEC <<endl;

    return 0;
}


