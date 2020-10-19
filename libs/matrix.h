#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>

template <typename Ty>
class Matrix {
   public:
    int row;
    int colum;
    std::vector<std::vector<Ty>> value;

    Matrix(int _row, int _colum) {
        row = _row;
        colum = _colum;
        for (int i = 0; i < row; i++) {
            std::vector<Ty> tmpv;
            for (int j = 0; j < colum; j++) {
                tmpv.push_back(0);
            }
            value.push_back(tmpv);
        };
    };

    ~Matrix() {};

    std::vector<Ty>& operator[](int n){
        return value[n];
    }

    std::vector<int> size(){
        std::vector<int> size = {row, colum};
        return size;
    }

    void T(){
        int tmp = row;
        row = colum;
        colum = tmp;

        std::vector<std::vector<Ty>>TM;
        for(int i=0; i<row; i++){
            std::vector<Ty>TMtemp;
            for(int j=0; j<colum; j++){
                TMtemp.push_back(value[j][i]);
            }
            TM.push_back(TMtemp);
        }

        value = TM;
    }
};

template <typename Ty>
Matrix<Ty>multiply(Matrix<Ty> &F, Matrix<Ty> &xv){
    auto FM = F.value;
    auto xvm = xv.value;

    if(F.colum != xv.row){
        std::cout << "Matrix Shape Error" << std::endl;
        return Matrix<double>(0,0);
    }
    else{
        Matrix<Ty> Ans(F.row, xv.colum);
        for(int k = 0; k<F.row; k++){
            for(int i=0; i<xv.colum; i++){
                for(int j=0; j<xv.row; j++){
                    Ans.value[k][i] += FM[k][j]*xvm[j][i];
                }
            }
        }

        return Ans;
    }
}

template <typename Ty>
Matrix<Ty> Identity(int n){
    Matrix<Ty> E(n, n);
    for(int i=0; i<n; i++){
        E[i][i] = 1.0;
    }

    return E;
}

#endif