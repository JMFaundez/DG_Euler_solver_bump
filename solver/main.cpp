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

    int p = 2;
    int mesh = 2;
    int max_iter= 130000;
    bool free_stream = false;
    double tolerance = 1.0e-7;
    double CFL=1.3;
    if(mesh==0){
        N_int = 132; 
        N_bound = 36; 
        N_elem = 100;
        N_curv = 28;
        N_nodes =  265;
        files[0] ="bump0_curv";
        files[1] = "./Meshes/I2E0.txt";
        files[2] = "./Meshes/B2E0.txt"; 
        files[3] = "./Meshes/Bn0.txt";
        files[4] = "./Meshes/In0.txt";
        files[5] = "./Meshes/Area0.txt";
    }
    if(mesh==1){
        N_int = 564; 
        N_bound = 72; 
        N_elem = 400;
        N_curv = 56;
        N_nodes =  629;
        files[0] ="bump1_curv";
        files[1] = "./Meshes/I2E1.txt";
        files[2] = "./Meshes/B2E1.txt"; 
        files[3] = "./Meshes/Bn1.txt";
        files[4] = "./Meshes/In1.txt";
        files[5] = "./Meshes/Area1.txt";
    }
    if(mesh==2){
        N_int = 2328; 
        N_bound = 144; 
        N_elem = 1600;
        N_curv = 112;
        N_nodes =  1657;
        files[0] ="bump2_curv";
        files[1] = "./Meshes/I2E2.txt";
        files[2] = "./Meshes/B2E2.txt"; 
        files[3] = "./Meshes/Bn2.txt";
        files[4] = "./Meshes/In2.txt";
        files[5] = "./Meshes/Area2.txt";
    }
    clock_t begin = clock();
    solve_GD(max_iter, free_stream, tolerance, N_nodes, N_elem, N_curv, N_int, N_bound, files, p, CFL);
    clock_t end = clock();
    cout << "Total time: " << double(end-begin)/CLOCKS_PER_SEC <<endl;

    return 0;
}


