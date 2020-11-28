#include "matrix.h"

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
template <typename Ty>
Matrix<Ty> operator+(const Matrix<Ty>& M, Ty n) {
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = M[i][j] + n;
        }
    }

    return Ans;
}

template <typename Ty>
Matrix<Ty> operator+(Ty n, const Matrix<Ty>& M) {
    return M + n;
}

template <typename Ty>
Matrix<Ty> operator*(const Matrix<Ty>& M, Ty n) {
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = M[i][j] * n;
        }
    }

    return Ans;
}

template <typename Ty>
Matrix<Ty> operator*(Ty n, const Matrix<Ty>& M) {
    return M * n;
}

template <typename Ty>
Matrix<Ty> operator-(const Matrix<Ty>& M, Ty n) {
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = M[i][j] - n;
        }
    }

    return Ans;
}

template <typename Ty>
Matrix<Ty> operator-(Ty n, const Matrix<Ty>& M) {
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = n - M[i][j];
        }
    }

    return Ans;
}

template <typename Ty>
Matrix<Ty> operator/(const Matrix<Ty>& M, Ty n) {
    Matrix<Ty> Ans(M.row, M.colum);
    for (int i = 0; i < M.row; i++) {
        for (int j = 0; j < M.colum; j++) {
            Ans[i][j] = M[i][j] / n;
        }
    }

    return Ans;
}

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
            stream << std::setw(10) << std::setprecision(2) << v.value[i][j] << ", ";
        }
        stream << "]";
    }
    stream << "]\n";

    return stream;
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
            for (int j = diognal-1; j < i; j++) {
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
std::pair<Matrix<Ty>, std::vector<Ty>> Forward_easure(const Matrix<Ty>& M, const std::vector<Ty>& b, std::vector<int>& pbidx){
    int R = M.row;
    int C = M.colum;
    int idx = pbidx.size();
    Ty tmax = 0;
    int tmaxidx = 0;
    auto FE = M;
    auto FB = b;

    if(R!=C) throw std::range_error("Not Square");
    if(R!=idx) throw std::range_error("Dont Match the row size");
    for(int i=0; i<idx; i++){
        pbidx[i] = i;
    }

    for(int j=0; j<C-1; j++){
        tmax = FE[j][j];
        for(int i=j; i<R; i++){
            if (tmax < FE[i][j]){
                tmax = FE[i][j];
                tmaxidx = i;
            }
        }

        if (tmax==0) throw std::range_error("Max is 0.");

        if(tmaxidx!=j){
            std::vector<Ty> maxrow = FE[tmaxidx];
            FE[tmaxidx] = FE[j];
            FE[j] = maxrow;
            Ty FBmax = FB[tmaxidx];
            FB[tmaxidx] = FB[j];
            FB[j] = FB[tmaxidx];
        }

        for(int i=j+1; j<R; j++){
            double m = FE[i][j] / FE[j][j];
            m = Ty(m);
            for(int k=j+1; k<C; k++){
                FE[i][k] -= m*FE[j][k];
            }
            b[i] -= m*b[j];
        }
    }

    std::pair<Matrix<Ty>, std::vector<Ty>> FEpair(FE, FB);

    return FEpair;
}