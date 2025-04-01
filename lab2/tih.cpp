#include <iostream>  
#include <vector>    
#include <iomanip>   
#include <cmath>     

using namespace std;  

// Функция для выполнения матричного умножения (C = A * B)
vector<vector<double>> matrixMultiply(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int rowsA = A.size();     // Получаем количество строк в матрице A
    int colsA = A[0].size();    // Получаем количество столбцов в матрице A (предполагается, что все строки имеют одинаковую длину)
    int colsB = B[0].size();    // Получаем количество столбцов в матрице B

    // Создаем результирующую матрицу C размером rowsA x colsB, заполненную нулями
    vector<vector<double>> result(rowsA, vector<double>(colsB, 0.0));

    // Выполняем матричное умножение
    for (int i = 0; i < rowsA; ++i) {        // Перебираем строки матрицы A
        for (int j = 0; j < colsB; ++j) {    // Перебираем столбцы матрицы B
            for (int k = 0; k < colsA; ++k) { // Перебираем столбцы матрицы A (или строки матрицы B)
                result[i][j] += A[i][k] * B[k][j]; // Вычисляем элемент C[i][j] как сумму произведений элементов A[i][k] и B[k][j]
            }
        }
    }
    return result; // Возвращаем результирующую матрицу C
}

// Функция для транспонирования матрицы
vector<vector<double>> matrixTranspose(const vector<vector<double>>& A) {
    int rows = A.size();     // Получаем количество строк в матрице A
    int cols = A[0].size();    // Получаем количество столбцов в матрице A

    // Создаем результирующую матрицу размером cols x rows, заполненную нулями
    vector<vector<double>> result(cols, vector<double>(rows, 0.0));

    // Выполняем транспонирование
    for (int i = 0; i < rows; ++i) {      // Перебираем строки матрицы A
        for (int j = 0; j < cols; ++j) {  // Перебираем столбцы матрицы A
            result[j][i] = A[i][j]; // Меняем местами индексы: элемент A[i][j] становится элементом result[j][i]
        }
    }
    return result; // Возвращаем транспонированную матрицу
}

// Функция для сложения двух матриц (C = A + B)
vector<vector<double>> matrixAdd(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int rows = A.size();     // Получаем количество строк в матрице A
    int cols = A[0].size();    // Получаем количество столбцов в матрице A

    // Создаем результирующую матрицу C размером rows x cols, заполненную нулями
    vector<vector<double>> result(rows, vector<double>(cols, 0.0));

    // Выполняем сложение
    for (int i = 0; i < rows; ++i) {      // Перебираем строки матрицы A и B
        for (int j = 0; j < cols; ++j) {  // Перебираем столбцы матрицы A и B
            result[i][j] = A[i][j] + B[i][j]; // Складываем соответствующие элементы матриц A и B
        }
    }
    return result; // Возвращаем результирующую матрицу C
}

// Функция для умножения матрицы на скаляр
vector<vector<double>> scalarMultiply(const vector<vector<double>>& A, double scalar) {
    int rows = A.size();     // Получаем количество строк в матрице A
    int cols = A[0].size();    // Получаем количество столбцов в матрице A

    // Создаем результирующую матрицу размером rows x cols, заполненную нулями
    vector<vector<double>> result(rows, vector<double>(cols, 0.0));

    // Выполняем умножение
    for (int i = 0; i < rows; ++i) {      // Перебираем строки матрицы A
        for (int j = 0; j < cols; ++j) {  // Перебираем столбцы матрицы A
            result[i][j] = A[i][j] * scalar; // Умножаем каждый элемент матрицы A на скаляр
        }
    }
    return result; // Возвращаем результирующую матрицу
}

// Функция для создания единичной матрицы
vector<vector<double>> createIdentityMatrix(int size) {
    // Создаем единичную матрицу размером size x size, заполненную нулями
    vector<vector<double>> result(size, vector<double>(size, 0.0));

    // Заполняем главную диагональ единицами
    for (int i = 0; i < size; ++i) {
        result[i][i] = 1.0; // Элемент result[i][i] равен 1
    }
    return result; // Возвращаем единичную матрицу
}

// Функция для выполнения умножения матрицы на вектор (y = A * x)
vector<double> matrixVectorMultiply(const vector<vector<double>>& A, const vector<double>& x) {
    int rowsA = A.size();     // Получаем количество строк в матрице A
    int colsA = A[0].size();    // Получаем количество столбцов в матрице A
    int sizeX = x.size();     // Получаем размерность вектора x

    // Создаем результирующий вектор размером rowsA, заполненный нулями
    vector<double> result(rowsA, 0.0);

    // Выполняем умножение
    for (int i = 0; i < rowsA; ++i) {      // Перебираем строки матрицы A
        for (int j = 0; j < colsA; ++j) {  // Перебираем столбцы матрицы A (или элементы вектора x)
            result[i] += A[i][j] * x[j]; // Вычисляем элемент result[i] как сумму произведений элементов A[i][j] и x[j]
        }
    }
    return result; // Возвращаем результирующий вектор
}

// Функция для решения системы линейных уравнений методом Гаусса (Ax = b)
vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();  // Получаем размерность матрицы A

    // Прямой ход (forward elimination)
    for (int i = 0; i < n; ++i) {
        // Поиск опорного элемента (pivot element) в столбце i
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) { // Ищем строку с максимальным по модулю элементом в столбце i
                maxRow = k; // Запоминаем индекс этой строки
            }
        }

        // Перестановка строк 
        swap(A[i], A[maxRow]); // Меняем местами строки i и maxRow в матрице A
        swap(b[i], b[maxRow]); // Меняем местами элементы i и maxRow в векторе b

        // Обнуление элементов ниже опорного
        for (int k = i + 1; k < n; ++k) { // Перебираем строки ниже опорной
            double c = -A[k][i] / A[i][i]; // Вычисляем коэффициент для обнуления элемента A[k][i]
            for (int j = i; j < n; ++j) {   // Перебираем элементы в строке k, начиная с элемента A[k][i]
                if (i == j) {
                    A[k][j] = 0; // Обнуляем элемент A[k][i] (из-за округления может быть неточно 0)
                } else {
                    A[k][j] += c * A[i][j]; // Вычитаем из элемента A[k][j] элемент A[i][j], умноженный на коэффициент c
                }
            }
            b[k] += c * b[i]; // Модифицируем соответствующий элемент вектора b
        }
    }

    // Обратный ход (back substitution)
    vector<double> x(n, 0); // Создаем вектор решения x, заполненный нулями
    for (int i = n - 1; i >= 0; --i) {  // Перебираем строки в обратном порядке
        x[i] = b[i] / A[i][i]; // Вычисляем элемент x[i]
        for (int k = i - 1; k >= 0; --k) { // Перебираем строки выше текущей
            b[k] -= A[k][i] * x[i]; // Модифицируем элементы вектора b
        }
    }

    return x; // Возвращаем вектор решения x
}

// Функция для выполнения регуляризации Тихонова
vector<double> tikhonovRegularization(const vector<vector<double>>& A, const vector<double>& b, double lambda) {
    // 1. Вычисляем A_transpose * A (A^T * A)
    vector<vector<double>> At = matrixTranspose(A); // Транспонируем матрицу A
    vector<vector<double>> AtA = matrixMultiply(At, A); // Умножаем транспонированную матрицу на исходную

    // 2. Создаем lambda * I (λ * I), где I - единичная матрица
    vector<vector<double>> lambdaI = scalarMultiply(createIdentityMatrix(A[0].size()), lambda); // Создаем единичную матрицу и умножаем ее на скаляр lambda

    // 3. Вычисляем A_transpose * A + lambda * I (A^T * A + λ * I)
    vector<vector<double>> AtA_lambdaI = matrixAdd(AtA, lambdaI); // Складываем матрицы (A^T * A) и (λ * I)

    // 4. Вычисляем A_transpose * b (A^T * b)
    vector<double> Atb = matrixVectorMultiply(At, b); // Умножаем транспонированную матрицу на вектор b

    // 5. Решаем систему (A_transpose * A + lambda * I) * x = A_transpose * b методом Гаусса
    vector<double> x = gaussianElimination(AtA_lambdaI, Atb); // Решаем систему линейных уравнений методом Гаусса

    return x; // Возвращаем вектор решения x
}

int main() {
    // Задаем исходную матрицу A и вектор b
    vector<vector<double>> A = {
        {3, 1, -1},
        {2, 4, 1},
        {1, -1, 3}
    };

    vector<double> b = {6, 9, 4};

    // Задаем параметр регуляризации lambda
    double lambda =2.449* 0.00001;

    // Выполняем регуляризацию Тихонова
    vector<double> x = tikhonovRegularization(A, b, lambda); // Вызываем функцию регуляризации Тихонова

    // Выводим решение
    cout << "Solution x:\n";
    for (double val : x) { // Перебираем элементы вектора x
        cout << fixed << setprecision(6) << val << " "; // Выводим каждый элемент с фиксированной точностью 6 знаков после запятой
    }
    cout << endl; // Переходим на новую строку

    return 0; // Завершаем программу
}