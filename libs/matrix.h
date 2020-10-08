#ifndef MATRIX_H
#define MATRIX_H

template <typename T>
class Matrix {
   public:
    int row;
    int colum;
    T **Arr;

    Matrix(int _row, int _colum) {
        row = _row;
        colum = _colum;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < colum; j++) Arr[i][j] = 0;
        };
    };

    ~Matrix() {
        std::cout << "Delete" << std::endl;
        for (int i = 0; i < row; i++) {
            delete[] Arr[i];
        }
        delete[] Arr;
    };
};

#endif