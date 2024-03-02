#include "DG.h"

void read_file(double V[], int El[], int Ec[],int I2E[], int B2E[],
                double In[],double Bn[],double P[],double Area[],
                double dl_int[], double dl_bound[],
                string namefile[6],int N_int, int N_bound,
                int N_elem, int N_curv, int N_nodes, double I_mid[], double B_mid[]){
    
    ifstream grifile(namefile[0]);
    ifstream I2Efile(namefile[1]);
    ifstream B2Efile(namefile[2]);
    ifstream Bnfile(namefile[3]);
    ifstream Infile(namefile[4]);
    ifstream Areafile(namefile[5]);
    int N[3]; // number of nodes, elements and dimension
    for(int i=0;i<3;i++){
        grifile >> N[i];
    }
    if(N[1]!=N_elem){
        cout <<"check number of elements" << endl;
        cout << N_elem << N[1] << endl;
    }
    // coordinates of the nodes
    for(int i=0;i<N[0];i++){
        grifile >> V[i*N[2]];
        grifile >> V[i*N[2]+1];
    }
    int NB;
    // number of boundaires
    grifile >> NB;
    //  names of the boundaries
    string* names = new string[NB];
    // number of faces per boundary
    int* N_bf = new int[NB];
    int extra;
    // boundaries
    int** B = new int*[NB];
    for(int i=0;i<NB;i++){
        grifile >> N_bf[i];
        grifile >> extra;
        grifile >> names[i];
        B[i] = new int[N_bf[i]*N[2]];
        for(int j=0;j<N_bf[i];j++){
            grifile >> B[i][j*N[2]];
            grifile >> B[i][j*N[2]+1];
        }   
    }
    string name_el;
    grifile >> extra;
    grifile >> extra;
    grifile >> name_el;
    for(int i=0;i<N_elem-N_curv;i++){
        grifile >> El[3*i];
        grifile >> El[3*i+1];
        grifile >> El[3*i+2];
        // 0-based
        El[i*3] = El[i*3]-1;
        El[i*3+1] = El[i*3+1]-1;
        El[i*3+2] = El[i*3+2]-1;
    }
    grifile >> extra;
    grifile >> extra;
    grifile >> name_el;
    for(int i=0;i<N_curv;i++){
        for(int j=0;j<10;j++){
            grifile >> Ec[10*i+j];
            Ec[10*i+j] = Ec[10*i+j]-1;
        }
    }
    int node1, node2;
    int elem_L;
    int index_clock[3] = {2,0,1};
    int index_anti[3] = {1,2,0};
    int index_clock_c[3] = {9,0,3};
    int index_anti_c[3] = {3,9,0};
    // Read I2E 
    for(int i=0;i<N_int;i++){
        I2Efile >> I2E[i*4];
        I2Efile >> I2E[i*4+1];
        I2Efile >> I2E[i*4+2];
        I2Efile >> I2E[i*4+3];
        // 0-based
        I2E[i*4] = I2E[i*4]-1;
        I2E[i*4+1] = I2E[i*4+1]-1;
        I2E[i*4+2] = I2E[i*4+2]-1;
        I2E[i*4+3] = I2E[i*4+3]-1;
        if(I2E[i*4]<N_elem-N_curv){
            node1 = El[3*I2E[i*4]+index_anti[I2E[i*4+1]]];
            node2 = El[3*I2E[i*4]+index_clock[I2E[i*4+1]]];
        }
        else{
            node1 = Ec[10*(I2E[i*4]-(N_elem-N_curv))+index_anti_c[I2E[i*4+1]]];
            node2 = Ec[10*(I2E[i*4]-(N_elem-N_curv))+index_clock_c[I2E[i*4+1]]];
        }
        dl_int[i] = sqrt(pow(V[2*node1]-V[2*node2],2)+pow(V[2*node1+1]-V[2*node2+1],2));
        I_mid[2*i] = (V[2*node1]+V[2*node2])/2;
        I_mid[2*i+1] = (V[2*node1+1]+V[2*node2+1])/2;
        P[I2E[i*4]]+= dl_int[i];
        P[I2E[i*4+2]]+= dl_int[i];
    }
    
    // Read B2E 
    for(int i=0;i<N_bound;i++){
        B2Efile >> B2E[i*3];
        B2Efile >> B2E[i*3+1];
        B2Efile >> B2E[i*3+2];
        // 0-based 
        B2E[i*3] = B2E[i*3]-1;
        B2E[i*3+1] = B2E[i*3+1]-1;
        B2E[i*3+2] = B2E[i*3+2]-1;
        if(B2E[i*3]<N_elem-N_curv){
            node1 = El[3*B2E[i*3]+index_anti[B2E[i*3+1]]];
            node2 = El[3*B2E[i*3]+index_clock[B2E[i*3+1]]];
        }
        else{
            node1 = Ec[10*(B2E[i*3]-(N_elem-N_curv))+index_anti_c[B2E[i*3+1]]];
            node2 = Ec[10*(B2E[i*3]-(N_elem-N_curv))+index_clock_c[B2E[i*3+1]]];

        }
        dl_bound[i] = sqrt(pow(V[2*node1]-V[2*node2],2)+pow(V[2*node1+1]-V[2*node2+1],2));
        B_mid[2*i] = (V[2*node1]+V[2*node2])/2;
        B_mid[2*i+1] = (V[2*node1+1]+V[2*node2+1])/2;
        P[B2E[i*3]]+= dl_bound[i];
    }
    
    // Read In
    for(int i=0;i<N_int;i++){
        Infile >> In[i*2];
        Infile >> In[i*2+1];\
    }
    // Read Bn 
    for(int i=0;i<N_bound;i++){
        Bnfile >> Bn[i*2];
        Bnfile >> Bn[i*2+1];\
    }
    
    // Read Area
    for(int i=0;i<N_elem;i++){
        Areafile >> Area[i];
    }
    
    // Delete 
    delete [] names;
    delete [] N_bf;
    for(int i=0;i<NB;i++){
        delete [] B[i];
    }
    grifile.close();
    I2Efile.close();
    B2Efile.close();
    Infile.close();
    Bnfile.close();
    Areafile.close();
}