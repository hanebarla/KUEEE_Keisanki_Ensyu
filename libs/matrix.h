#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <utility>
#include <vector>

template <typename Ty>
class Matrix {
   public:
    int row;
    int colum;
    std::vector<std::vector<Ty>> value;

    // コンストラクタ
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

    // デストラクタ
    ~Matrix(){};

    Matrix<Ty>& operator=(const std::vector<std::vector<Ty>>& vin) {
        row = int(vin.size());
        colum = int(vin[0].size());
        for (int i = 0; i < row; i++) {
            value[i] = vin[i];
        }

        return *this;
    }

    // スカラーとの四則演算
    Matrix<Ty>& operator*=(Ty n) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < colum; j++) {
                value[i][j] *= n;
            }
        }

        return *this;
    }

    // スカラーとのブロードキャスト計算
    Matrix<Ty>& operator/=(Ty n) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < colum; j++) {
                value[i][j] /= n;
            }
        }

        return *this;
    }

    Matrix<Ty>& operator+=(Ty n) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < colum; j++) {
                value[i][j] += n;
            }
        }

        return *this;
    }

    Matrix<Ty>& operator-=(Ty n) {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < colum; j++) {
                value[i][j] += n;
            }
        }

        return *this;
    }

    // 配列のアクセス
    std::vector<Ty>& operator[](int n) { return value[n]; }

    // シェイプの確認
    std::vector<int> size() {
        std::vector<int> size = {row, colum};
        return size;
    }

    // 転置
    void T() {
        int tmp = row;
        row = colum;
        colum = tmp;

        std::vector<std::vector<Ty>> TM;
        for (int i = 0; i < row; i++) {
            std::vector<Ty> TMtemp;
            for (int j = 0; j < colum; j++) {
                TMtemp.push_back(value[j][i]);
            }
            TM.push_back(TMtemp);
        }

        value = TM;
    }
};

// 行列の出力
template <typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T> v);

// フラット化
template <typename Ty>
Matrix<Ty> flatten(const Matrix<Ty>& M);

//リシェイプ
template <typename Ty>
Matrix<Ty> reshape(const Matrix<Ty>& M, int r, int c);

// 行列とスカラーの四則演算
template <typename Ty>
Matrix<Ty> operator+(const Matrix<Ty>& M, Ty n);

template <typename Ty>
Matrix<Ty> operator+(Ty n, const Matrix<Ty>& M);

template <typename Ty>
Matrix<Ty> operator*(const Matrix<Ty>& M, Ty n);

template <typename Ty>
Matrix<Ty> operator*(Ty n, const Matrix<Ty>& M);

template <typename Ty>
Matrix<Ty> operator-(const Matrix<Ty>& M, Ty n);

template <typename Ty>
Matrix<Ty> operator-(Ty n, const Matrix<Ty>& M);

template <typename Ty>
Matrix<Ty> operator/(const Matrix<Ty>& M, Ty n);

// 行列の積 (後でSUUMAの実装)
template <typename Ty>
Matrix<Ty> dot(const Matrix<Ty>& F, const Matrix<Ty>& xv);

// 行列とベクターの積
template <typename Ty>
std::vector<Ty> dot(Matrix<Ty>& F, const std::vector<Ty>& xv);

// 単位行列の代入
template <typename Ty>
Matrix<Ty> Identity(int n);

// 連続値の行列
template <typename Ty>
Matrix<Ty> Arrange(int n);

// LU分解
template <typename Ty>
Matrix<Ty> LUdecom(const Matrix<Ty>& M, std::vector<int>& idx);

// LとUを取り出す
template <typename Ty>
std::vector<Matrix<Ty>> LU(Matrix<Ty>& MC);

// callicurate
template <typename Ty>
std::vector<Ty> LUans(Matrix<Ty>& MC, std::vector<int>& idx,
                      std::vector<Ty>& bin);

// aggrigate LU
template <typename Ty>
std::vector<Ty> MatEqu(const Matrix<Ty>& A, std::vector<Ty>& b);

// 逆行列を求める
template <typename Ty>
Matrix<Ty> Inv(const Matrix<Ty>& M);

// 前身消去
template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Forward_easure(const Matrix<Ty>& M,
                                                      const std::vector<Ty>& b,
                                                      std::vector<int>& pbidx);
#endif