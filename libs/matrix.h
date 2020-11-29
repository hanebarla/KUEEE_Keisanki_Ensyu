#ifndef MATRIX_H
#define MATRIX_H
#include <omp.h>

#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "utils.h"

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

    void zeros() {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < colum; j++) {
                value[i][j] = 0;
            }
        }
    }
};

// 行列の出力
template <typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T> v) {
    stream << "[";
    for (int i = 0; i < v.row; i++) {
        if (i != 0) {
            stream << "\n ";
        }
        stream << "[ ";
        for (int j = 0; j < v.colum; j++) {
            stream << std::setw(10) << std::setprecision(2) << v.value[i][j]
                   << ", ";
        }
        stream << "]";
    }
    stream << "]\n";

    return stream;
}

// フラット化
template <typename Ty>
Matrix<Ty> flatten(const Matrix<Ty>& M) {
    int R = M.row;
    int C = M.colum;
    Matrix<Ty> flat(1, R * C);

    for (int r = 0; r < R; r++) {
        for (int c = 0; c < C; c++) {
            flat[0][r * C + c] = M.value[r][c];
        }
    }
    return flat;
}

//リシェイプ
template <typename Ty>
Matrix<Ty> reshape(const Matrix<Ty>& M, int r, int c) {
    if (r * c != M.row * M.colum) throw std::range_error("Shape size error.");
    Matrix<Ty> R(r, c);
    auto F = flatten(M);

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            R[i][j] = F[0][i * c + j];
        }
    }

    return R;
}

// 行列とスカラーの四則演算
// 行列 + スカラー (ブロードキャスト)
template <typename Ty>
Matrix<Ty> operator+(const Matrix<Ty>& Mo, Ty n) {
    auto M = Mo;
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = M[i][j] + n;
        }
    }

    return Ans;
}

// スカラー + 行列 ( ブロードキャスト)
template <typename Ty>
Matrix<Ty> operator+(Ty n, const Matrix<Ty>& M) {
    return M + n;
}

// 行列 * スカラー (ブロードキャスト)
template <typename Ty>
Matrix<Ty> operator+(const Matrix<Ty>& Ao, const Matrix<Ty>& Bo) {
    auto A = Ao;
    auto B = Bo;
    Matrix<Ty> Ans(A.row, A.colum);

    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.colum; j++) {
            Ans[i][j] = A[i][j] + B[i][j];
        }
    }

    return Ans;
}

// 行列 * スカラー (ブロードキャスト)
template <typename Ty>
Matrix<Ty> operator*(const Matrix<Ty>& Mo, Ty n) {
    auto M = Mo;
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = M[i][j] * n;
        }
    }

    return Ans;
}

// スカラー * 行列 (ブロードキャスト)
template <typename Ty>
Matrix<Ty> operator*(Ty n, const Matrix<Ty>& M) {
    return M * n;
}

// アダマール積
template <typename Ty>
Matrix<Ty> operator*(const Matrix<Ty>& Ao, const Matrix<Ty>& Bo) {
    auto A = Ao;
    auto B = Bo;
    Matrix<Ty> Ans(A.row, A.colum);

    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.colum; j++) {
            Ans[i][j] = A[i][j] * B[i][j];
        }
    }

    return Ans;
}

// 行列 - スカラー (ブロードキャスト)
template <typename Ty>
Matrix<Ty> operator-(const Matrix<Ty>& Mo, Ty n) {
    auto M = Mo;
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = M[i][j] - n;
        }
    }

    return Ans;
}

// スカラー - 行列 (ブロードキャスト)
template <typename Ty>
Matrix<Ty> operator-(Ty n, const Matrix<Ty>& Mo) {
    auto M = Mo;
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = n - M[i][j];
        }
    }

    return Ans;
}

// 行列 - 行列
template <typename Ty>
Matrix<Ty> operator-(const Matrix<Ty>& Ao, const Matrix<Ty>& Bo) {
    auto A = Ao;
    auto B = Bo;
    Matrix<Ty> Ans(A.row, A.colum);

    for (int i = 0; i < A.row; i++) {
        for (int j = 0; j < A.colum; j++) {
            Ans[i][j] = A[i][j] - B[i][j];
        }
    }

    return Ans;
}

// 行列 / スカラー (ブロードキャスト)
template <typename Ty>
Matrix<Ty> operator/(const Matrix<Ty>& Mo, Ty n) {
    auto M = Mo;
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = M[i][j] / n;
        }
    }

    return Ans;
}

// 行列の積 (後でSUUMAの実装)
template <typename Ty>
Matrix<Ty> dot(const Matrix<Ty>& F, const Matrix<Ty>& xv) {
    auto FM = F.value;
    auto xvm = xv.value;

    if (F.colum != xv.row) {
        throw std::range_error("Don't mache the shapes.");
    } else {
        Matrix<Ty> Ans(F.row, xv.colum);
        for (int k = 0; k < F.row; k++) {
            for (int i = 0; i < xv.colum; i++) {
                for (int j = 0; j < xv.row; j++) {
                    Ans.value[k][i] += FM[k][j] * xvm[j][i];
                }
            }
        }

        return Ans;
    }
}

// 行列とベクターの積
template <typename Ty>
std::vector<Ty> dot(Matrix<Ty>& F, const std::vector<Ty>& xv) {
    int xsize = xv.size();

    if (F.colum != xsize) throw std::range_error("Don't mache the shapes.");

    std::vector<Ty> Ans(F.row, 0);
    for (int k = 0; k < F.row; k++) {
        for (int j = 0; j < F.colum; j++) {
            Ans[k] += F[k][j] * xv[j];
        }
    }

    return Ans;
}

// 単位行列の代入
template <typename Ty>
Matrix<Ty> Identity(int n) {
    Matrix<Ty> E(n, n);
    for (int i = 0; i < n; i++) {
        E[i][i] = 1.0;
    }

    return E;
}

// 連続値の行列
template <typename Ty>
Matrix<Ty> Arrange(int n) {
    Matrix<Ty> Ar(1, n);
    for (int i = 0; i < n; i++) {
        Ar[0][i] = i;
    }

    return Ar;
}

// 残差ノルム
template <typename Ty>
Ty Resnorm(const Matrix<Ty>& Mo, const std::vector<Ty>& bo,
           const std::vector<Ty>& xko) {
    auto M = Mo;
    auto b = bo;
    auto xk = xko;
    Ty resleng = 0;
    Ty ansleng = 0;

    std::vector<Ty> res = dot(M, xk) - b;
    for (int i = 0; i < M.row; i++) {
        resleng += res[i] * res[i];
        ansleng += b[i] * b[i];
    }

    return Ty(resleng / ansleng);
}
// 残差ノルム(pair)
template <typename Ty>
Ty Resnorm(const std::pair<Matrix<Ty>, std::vector<Ty>>& pr,
           const std::vector<Ty>& xko) {
    auto M = pr.first;
    auto b = pr.second;
    auto xk = xko;
    Ty resleng = 0;
    Ty ansleng = 0;

    std::vector<Ty> res = dot(M, xk) - b;
    for (int i = 0; i < M.row; i++) {
        resleng += res[i] * res[i];
        ansleng += b[i] * b[i];
    }

    return Ty(resleng / ansleng);
}

// LU分解
template <typename Ty>
Matrix<Ty> LUdecom(const Matrix<Ty>& M, std::vector<int>& idx) {
    auto MC = M;
    int R = M.row;
    int C = M.colum;
    int imax = 0;
    Ty temp, dum, big, sum;
    std::vector<Ty> scale(R, 0);

    if (R != C) throw std::range_error("Not Square matirix.");

    // スケーリング
    for (int i = 0; i < R; i++) {
        big = 0;
        for (int j = 0; j < C; j++) {
            temp = fabsl(MC[i][j]);
            if (temp > big) big = temp;
        }

        if (big == 0.0) throw std::range_error("Max score is Zero.");

        scale[i] = 1.0 / big;
    }

    // 分解
    for (int j = 0; j < C; j++) {
        for (int i = 0; i < j; i++) {
            sum = MC[i][j];
            for (int k = 0; k < i; k++) {
                sum -= MC[i][k] * MC[k][j];
            }
            MC[i][j] = sum;
        }
        big = 0;

        for (int i = j; i < R; i++) {
            sum = MC[i][j];
            for (int k = 0; k < j; k++) {
                sum -= MC[i][k] * MC[k][j];
            }
            MC[i][j] = sum;

            dum = scale[i] * fabsl(sum);
            if (dum >= big) {
                big = dum;
                imax = i;
            }
        }

        if (j != imax) {
            for (int k = 0; k < C; k++) {
                dum = MC[imax][k];
                MC[imax][k] = MC[j][k];
                MC[j][k] = dum;
            }
            scale[imax] = scale[j];
        }

        idx[j] = imax;

        if (MC[j][j] == 0) MC[j][j] = 1e-20;

        if (j != C - 1) {
            dum = 1 / MC[j][j];
            for (int i = j + 1; i < R; i++) MC[i][j] *= dum;
        }
    }

    return MC;
}

// LとUを取り出す
template <typename Ty>
std::vector<Matrix<Ty>> LU(Matrix<Ty>& MC) {
    int R = MC.row;
    int C = MC.colum;

    if (R != C) throw std::range_error("Don't mache the shapes.");

    auto L = Identity<Ty>(R);
    Matrix<Ty> U(R, C);

    for (int i = 0; i < R; i++) {
        for (int j = 0; j < C; j++) {
            if (i <= j) {
                U[i][j] = MC[i][j];
            } else {
                L[i][j] = MC[i][j];
            }
        }
    }
    std::vector<Matrix<Ty>> LU = {L, U};

    return LU;
}

// callicurate
template <typename Ty>
std::vector<Ty> LUans(Matrix<Ty>& MC, std::vector<int>& idx,
                      std::vector<Ty>& bin) {
    auto b = bin;
    int R = MC.row;
    int C = MC.colum;
    int diognal = -1;
    int tmp_idx;
    Ty sum;

    for (int i = 0; i < R; i++) {
        tmp_idx = idx[i];
        sum = b[tmp_idx];
        b[tmp_idx] = b[i];
        if (diognal != -1) {
            for (int j = diognal - 1; j < i; j++) {
                sum -= MC[i][j] * b[j];
            }
        } else if (sum != 0) {
            diognal = i;
        }
        b[i] = sum;
    }

    for (int i = R - 1; i >= 0; i--) {
        sum = b[i];
        for (int j = i + 1; j < C; j++) {
            sum -= MC[i][j] * b[j];
        }
        b[i] = sum / MC[i][i];
    }

    return b;
}

// aggrigate LU
template <typename Ty>
std::vector<Ty> MatEqu(const Matrix<Ty>& A, std::vector<Ty>& b) {
    int R = A.row;
    int C = A.colum;

    if (R != C) throw std::range_error("Don't mache the shapes.");

    std::vector<int> idx(R, 0);
    auto MC = LUdecom(A, idx);
    auto x = LUans(MC, idx, b);

    return x;
}

// 逆行列を求める
template <typename Ty>
Matrix<Ty> Inv(const Matrix<Ty>& M) {
    int n = M.row;
    auto A = M;
    Matrix<Ty> inv(n, n);

    for (int i = 0; i < n; i++) {
        std::vector<Ty> b_inv(n, 0);
        b_inv[i] = 1;

        auto inv_c = MatEqu(M, b_inv);

        for (int j = 0; j < n; j++) {
            inv[j][i] = inv_c[j];
        }
    }

    return inv;
}

// 前身消去
template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Forward_easure(const Matrix<Ty>& M,
                                                      const std::vector<Ty>& b,
                                                      std::vector<int>& pbidx) {
    int R = M.row;
    int C = M.colum;
    int idx = pbidx.size();
    Ty tmax = 0;
    int tmaxidx = 0;
    auto FE = M;
    auto FB = b;

    // 初期化
    if (R != C) throw std::range_error("Not Square");
    if (R != idx) throw std::range_error("Dont Match the row size");
    for (int i = 0; i < idx; i++) {
        pbidx[i] = i;
    }

    // 一列ずつ
    for (int j = 0; j < (C - 1); j++) {
        // その列での一番大きい数字とそれがある行のindexを保存(ピボット選択)
        tmax = FE[j][j];
        tmaxidx = j;
        for (int i = j; i < R; i++) {
            if (tmax < FE[i][j]) {
                tmax = FE[i][j];
                tmaxidx = i;
            }
        }

        if (tmax == 0) throw std::range_error("Max is 0.");

        // j行目と最大行目を交換
        if (tmaxidx != j) {
            std::vector<Ty> maxrow = FE[tmaxidx];
            FE[tmaxidx] = FE[j];
            FE[j] = maxrow;
            Ty FBmax = FB[tmaxidx];
            FB[tmaxidx] = FB[j];
            FB[j] = FBmax;
            int temp_max_idx = pbidx[tmaxidx];
            pbidx[tmaxidx] = pbidx[j];
            pbidx[j] = temp_max_idx;
        }

        for (int i = j + 1; i < R; i++) {
            Ty m = Ty(FE[i][j] / FE[j][j]);
            for (int k = j; k < C; k++) {
                FE[i][k] -= m * FE[j][k];
            }
            FB[i] -= m * FB[j];
        }
    }

    std::pair<Matrix<Ty>, std::vector<Ty>> FEpair(FE, FB);

    return FEpair;
}
// 前進消去(pair)
template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Forward_easure_pair(
    const std::pair<Matrix<Ty>, std::vector<Ty>> pr, std::vector<int>& pbidx) {
    return Forward_easure<Ty>(pr.first, pr.second, pbidx);
}

// 後退代入
template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Backward_subsitution(
    const Matrix<Ty>& M, const std::vector<Ty>& b, std::vector<int>& pbidx) {
    int R = M.row;
    int C = M.colum;
    int idx = pbidx.size();
    Ty tmax = 0;
    int tmaxidx = 0;
    auto BS = M;
    auto Bb = b;

    if (R != C) throw std::range_error("Not Square");
    if (R != idx) throw std::range_error("Dont Match the row size");

    for (int i = (R - 1); i >= 0; i--) {
        for (int j = (C - 1); j > i; j--) {
            Bb[i] -= BS[i][j] * Bb[j];
        }
        if (BS[i][i] == 0) throw std::range_error("Diagonal is 0.");
        Bb[i] /= BS[i][i];
    }

    std::pair<Matrix<Ty>, std::vector<Ty>> BSpair(BS, Bb);

    return BSpair;
}
// 後退代入(pair)
template <typename Ty>
std::pair<Matrix<Ty>, std::vector<Ty>> Backward_subsitution_pair(
    const std::pair<Matrix<Ty>, std::vector<Ty>> pr, std::vector<int>& pbidx) {
    return Backward_subsitution(pr.first, pr.second, pbidx);
}

// DLU分解
template <typename Ty>
void DLUdecom(const Matrix<Ty>& Mo, Matrix<Ty>& D, Matrix<Ty>& L,
              Matrix<Ty>& U) {
    auto M = Mo;
    int R = M.row;
    int C = M.colum;
    D.zeros();
    L.zeros();
    U.zeros();

    if (R != C) throw std::range_error("Size Error");

    for (int i = 0; i < R; i++) {
        D[i][i] = M[i][i];

        for (int j = i + 1; j < C; j++) {
            L[j][i] = M[j][i];
            U[i][j] = M[i][j];
        }
    }
}

// Jacobi Step
template <typename Ty>
std::vector<Ty> Jacobi_Step(const Matrix<Ty>& M, const std::vector<Ty>& b,
                            std::vector<Ty>& xk) {
    std::vector<Ty> xkp = xk;
    auto D = M;
    auto L = M;
    auto U = M;
    DLUdecom(M, D, L, U);

    auto Dinv = Inv(D);
    auto B = -1.0 * dot(Dinv, (L + U));
    std::vector<Ty> c = dot(Dinv, b);
    std::vector<Ty> Bxk = dot(B, xk);

    return Bxk + c;
}
// ヤコビ法(pair)
template <typename Ty>
std::vector<Ty> Jacobi_Step_pair(
    const std::pair<Matrix<Ty>, std::vector<Ty>> pr, std::vector<Ty>& xk) {
    return Jacobi_Step(pr.first, pr.second, xk);
}

// ガウスザイデル法
template <typename Ty>
std::vector<Ty> Gauss_Seidel_Step(const Matrix<Ty>& M, const std::vector<Ty>& b,
                                  std::vector<Ty>& xk) {
    std::vector<Ty> xkp = xk;
    auto D = M;
    auto L = M;
    auto U = M;
    DLUdecom(M, D, L, U);

    auto DLinv = Inv(D + L);
    auto B = -1.0 * dot(DLinv, U);
    std::vector<Ty> c = dot(DLinv, b);
    std::vector<Ty> Bxk = dot(B, xk);

    return Bxk + c;
}
// ガウスザイデル法(pair)
template <typename Ty>
std::vector<Ty> Gauss_Seidel_Step_pair(
    const std::pair<Matrix<Ty>, std::vector<Ty>> pr, std::vector<Ty>& xk) {
    return Gauss_Seidel_Step(pr.first, pr.second, xk);
}

// SOR法
template <typename Ty>
std::vector<Ty> SOR_Step(const Matrix<Ty>& M, const std::vector<Ty>& b,
                                  std::vector<Ty>& xk, const Ty omega) {
    std::vector<Ty> xkp = xk;
    auto D = M;
    auto L = M;
    auto U = M;
    DLUdecom(M, D, L, U);

    auto Momega = (D + omega * L);
    auto Nomega = ((1.0 - omega) * D - omega * U);
    auto Momegainv = Inv(Momega);
    auto B = dot(Momegainv, Nomega);
    std::vector<Ty> c = omega * dot(Momegainv, b);
    std::vector<Ty> Bxk = dot(B, xk);

    return Bxk + c;
}
// SOR法(pair)
template <typename Ty>
std::vector<Ty> SOR_Step_pair(
    const std::pair<Matrix<Ty>, std::vector<Ty>> pr, std::vector<Ty>& xk,
    const Ty omega) {
    return SOR_Step(pr.first, pr.second, xk, omega);
}

// 共役勾配法
// 引数と返り値
template <typename Ty>
struct _CGReturn {
    std::vector<Ty> x;
    std::vector<Ty> r;
    std::vector<Ty> p;
    Ty rho;
};
// 共役勾配法
template <typename Ty>
_CGReturn<Ty> Conjugate_Gradient(const Matrix<Ty>& Mo, const _CGReturn<Ty>& cgin, int step){
    Matrix<Ty> M = Mo;
    std::vector<Ty> x = cgin.x;
    std::vector<Ty> r = cgin.r;
    std::vector<Ty> p_m = cgin.p;
    Ty rho_m = cgin.rho;
    std::vector<Ty> p = r;
    Ty beta = 0;

    _CGReturn<Ty> cgre = {x, r, p, 0.0};

    Ty rho = dot(r, r);
    if (step == 0){
        p = r;
    }else{
        beta = rho / rho_m;
        p = r + beta * p_m;
    }

    auto q = dot(M, p);
    auto alpha = rho / (dot(p, q));
    x = x + alpha * p;
    r = r - alpha * q;

    cgre.x = x;
    cgre.r = r;
    cgre.p = p;
    cgre.rho = rho;

    return cgre;
}

#endif