#include "DG.h"

void print_matrix(double Ma[], int start, int n_row, int n_col){
    for(int i=0 ; i<n_row; i++ ){
        for(int j=0 ; j<n_col; j++){
            cout << Ma[start*n_col*n_row+ i*n_col +j] << '\t';
        }
        cout << endl;
    }
}

void print_matrix_int(int Ma[], int start, int n_row, int n_col){
    for(int i=0 ; i<n_row; i++ ){
        for(int j=0 ; j<n_col; j++){
            cout << Ma[start*n_col*n_row+ i*n_col +j] << '\t';
        }
        cout << endl;
    }
}

void set_zero(double M[], int elem){
    for(int i=0;i<elem;i++){
        M[i] = 0;
    }
}

void invert_matrix(double M[], double M_inv[], int n_row, int init){
    Mat<double> MR;
    MR.zeros(n_row,n_row);
    for(int i=0; i<n_row; i++){
        for(int j=0 ; j<n_row ; j++){
            MR(i,j) = M[init*n_row*n_row + i*n_row + j];
        }
    }
    Mat<double> MI;
    MI.zeros(n_row,n_row);
    MI = arma::inv(MR);
    for(int i=0; i<n_row; i++){
        for(int j=0 ; j<n_row ; j++){
            M_inv[init*n_row*n_row + i*n_row + j] = MI(i,j);
        }
    }
}
