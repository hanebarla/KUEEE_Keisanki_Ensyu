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

    Matrix(int _row, int _colum, Ty n) {
        row = _row;
        colum = _colum;
        for (int i = 0; i < row; i++) {
            std::vector<Ty> tmpv;
            for (int j = 0; j < colum; j++) {
                tmpv.push_back(n);
            }
            value.push_back(tmpv);
        };
    };

    ~Matrix() {};

    Matrix<Ty>& operator*=(Ty n){
        for(int i=0; i<row; i++){
            for(int j=0; j<colum; j++){
                value[i][j] *= n;
            }
        }

        return *this;
    }

    Matrix<Ty>& operator/=(Ty n){
        for(int i=0; i<row; i++){
            for(int j=0; j<colum; j++){
                value[i][j] /= n;
            }
        }

        return *this;
    }

    Matrix<Ty>& operator+=(Ty n){
        for(int i=0; i<row; i++){
            for(int j=0; j<colum; j++){
                value[i][j] += n;
            }
        }

        return *this;
    }

    Matrix<Ty>& operator-=(Ty n){
        for(int i=0; i<row; i++){
            for(int j=0; j<colum; j++){
                value[i][j] += n;
            }
        }

        return *this;
    }

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

// 行列とスカラーの四則演算
template <typename Ty>
Matrix<Ty> operator+(const Matrix<Ty>& M, Ty n){
    Matrix<Ty> Ans(M.row, M.colum);
    for(int i=0; i<M.row; i++){
        for(int j=0; j<M.colum; j++){
            Ans[i][j] = M[i][j] + n;
        }
    }

    return Ans;
}

template <typename Ty>
Matrix<Ty> operator+(Ty n, const Matrix<Ty>& M){
    return M + n;
}

template <typename Ty>
Matrix<Ty> operator*(const Matrix<Ty>& M, Ty n){
    Matrix<Ty> Ans(M.row, M.colum);
    for(int i=0; i<M.row; i++){
        for(int j=0; j<M.colum; j++){
            Ans[i][j] = M[i][j] * n;
        }
    }

    return Ans;
}

template <typename Ty>
Matrix<Ty> operator*(Ty n, const Matrix<Ty>& M){
    return M * n;
}

template <typename Ty>
Matrix<Ty> operator-(const Matrix<Ty>& M, Ty n){
    Matrix<Ty> Ans(M.row, M.colum);
    for(int i=0; i<M.row; i++){
        for(int j=0; j<M.colum; j++){
            Ans[i][j] = M[i][j] - n;
        }
    }

    return Ans;
}

template <typename Ty>
Matrix<Ty> operator-(Ty n, const Matrix<Ty>& M){
    Matrix<Ty> Ans(M.row, M.colum);
    for(int i=0; i<M.row; i++){
        for(int j=0; j<M.colum; j++){
            Ans[i][j] = n - M[i][j];
        }
    }

    return Ans;
}

template <typename Ty>
Matrix<Ty> operator/(const Matrix<Ty>& M, Ty n){
    Matrix<Ty> Ans(M.row, M.colum);
    for(int i=0; i<M.row; i++){
        for(int j=0; j<M.colum; j++){
            Ans[i][j] = M[i][j] / n;
        }
    }

    return Ans;
}

// 行列の積 (後でSUUMAの実装)
template <typename Ty>
Matrix<Ty>dot(const Matrix<Ty> &F, const Matrix<Ty> &xv){
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

// 行列とベクターの積
template <typename Ty>
std::vector<Ty>dot(Matrix<Ty> &F, const std::vector<Ty> &xv){
    int xsize = xv.size();

    if(F.colum != xsize){
        std::cout << "Matrix Shape Error" << std::endl;
        std::vector<Ty> a;
        return a;
    }
    else{
        std::vector<Ty> Ans(F.row, 0);
        for(int k = 0; k<F.row; k++){
            for(int j=0; j<F.colum; j++){
                Ans[k] += F[k][j] * xv[j];
            }
        }

        return Ans;
    }
}

// 単位行列の代入
template <typename Ty>
Matrix<Ty> Identity(int n){
    Matrix<Ty> E(n, n);
    for(int i=0; i<n; i++){
        E[i][i] = 1.0;
    }

    return E;
}

#endif