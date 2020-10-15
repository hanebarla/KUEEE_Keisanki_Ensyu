#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>

template <typename T>
class Matrix {
   public:
    int row;
    int colum;
    std::vector<std::vector<T>> value;

    Matrix(int _row, int _colum) {
        row = _row;
        colum = _colum;
        for (int i = 0; i < row; i++) {
            std::vector<T> tmpv;
            for (int j = 0; j < colum; j++) {
                tmpv.push_back(0);
            }
            value.push_back(tmpv);
        };
    };

    ~Matrix() {};

    std::vector<int> size(){
        std::vector<int> size = {row, colum};
        return size;
    }

    void T(){
        int tmp = row;
        row = colum;
        colum = tmp;
    }
};

template <typename T>
Matrix<T>multiply(Matrix<T> &F, Matrix<T> &xv){
    auto FM = F.value;
    auto xvm = xv.value;

    if(F.colum != xv.row){
        std::cout << "Matrix Shape Error" << std::endl;
        return Matrix<double>(0,0);
    }
    else{
        Matrix<T> Ans(F.row, xv.colum);
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

#endif