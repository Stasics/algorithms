#include <iostream>  
#include <vector>    
#include <cmath>     
#include <iomanip>   

using namespace std;  

// Функция для вывода матрицы в консоль
void printMatrix(const vector<vector<double>>& matrix) {
    // const vector<vector<double>>& matrix: константная ссылка на матрицу (вектор векторов double), которую нужно вывести
    int rows = matrix.size();      // Получаем количество строк в матрице
    int cols = matrix[0].size();     // Получаем количество столбцов в матрице (предполагается, что все строки имеют одинаковую длину)
    for (int i = 0; i < rows; ++i) { // Перебираем строки матрицы
        for (int j = 0; j < cols; ++j) { // Перебираем столбцы матрицы
            cout << fixed << setprecision(2) << matrix[i][j] << " "; // Выводим элемент matrix[i][j] с фиксированной точкой и точностью 2 знака после запятой, добавляем пробел
        }
        cout << endl; // После вывода строки переходим на новую строку
    }
}

// Функция для решения системы линейных уравнений методом Гаусса-Жордана
vector<double> gaussJordan(vector<vector<double>> A, vector<double> b) {
    // vector<vector<double>> A: матрица коэффициентов системы уравнений (передается по значению, чтобы функция могла ее изменять)
    // vector<double> b: вектор свободных членов системы уравнений (передается по значению, чтобы функция могла его изменять)

    int n = A.size(); // Определяем размерность системы (количество уравнений или неизвестных)

    // Объединяем A и b в расширенную матрицу
    vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1)); // Создаем расширенную матрицу размером n x (n+1)
    for (int i = 0; i < n; ++i) { // Перебираем строки матрицы A
        for (int j = 0; j < n; ++j) { // Перебираем столбцы матрицы A
            augmentedMatrix[i][j] = A[i][j]; // Копируем элементы матрицы A в расширенную матрицу
        }
        augmentedMatrix[i][n] = b[i]; // Добавляем элементы вектора b в конец каждой строки расширенной матрицы
    }

    // Прямой и обратный ход метода Гаусса-Жордана
    for (int i = 0; i < n; ++i) {  // Перебираем столбцы (ведущие элементы)

        // Ищем строку с максимальным по модулю элементом в текущем столбце (ниже текущей строки)
        int maxRow = i; // Изначально предполагаем, что текущая строка содержит максимальный элемент
        for (int k = i + 1; k < n; ++k) { // Перебираем строки ниже текущей
            if (abs(augmentedMatrix[k][i]) > abs(augmentedMatrix[maxRow][i])) { // Если находим элемент, больше по модулю, чем текущий максимум
                maxRow = k; // Запоминаем индекс строки с максимальным элементом
            }
        }

        // Перестановка строк (если необходимо)
        if (maxRow != i) { // Если нашли строку с большим элементом, чем в текущей строке
            swap(augmentedMatrix[i], augmentedMatrix[maxRow]); // Меняем местами текущую строку и строку с максимальным элементом
            cout << "Row swap " << i << " and " << maxRow << endl; // Выводим информацию о перестановке строк
            printMatrix(augmentedMatrix); // Выводим текущее состояние расширенной матрицы
            cout << endl; // Выводим пустую строку для разделения этапов
        }

        // Нормализация ведущей строки (делаем ведущий элемент равным 1)
        double pivot = augmentedMatrix[i][i]; // Получаем значение ведущего элемента
        if (abs(pivot) < 1e-9) { // Если ведущий элемент близок к нулю (с учетом погрешности вычислений)
            cerr << "Matrix is singular. No unique solution." << endl; // Выводим сообщение об ошибке: матрица вырождена, нет единственного решения
            return vector<double>(n, NAN); // Возвращаем вектор, заполненный NaN (Not a Number), чтобы обозначить отсутствие решения
        }

        for (int j = i; j <= n; ++j) { // Перебираем элементы в текущей строке (включая элемент вектора b)
            augmentedMatrix[i][j] /= pivot; // Делим каждый элемент в строке на ведущий элемент, чтобы сделать его равным 1
        }
        cout << "Normalize row " << i << endl; // Выводим сообщение о нормализации строки
        printMatrix(augmentedMatrix); // Выводим текущее состояние расширенной матрицы
        cout << endl; // Выводим пустую строку для разделения этапов

        // Обнуление элементов выше и ниже ведущего элемента
        for (int k = 0; k < n; ++k) { // Перебираем все строки
            if (k != i) { // Если текущая строка не является ведущей
                double factor = augmentedMatrix[k][i]; // Вычисляем коэффициент, необходимый для обнуления элемента в текущем столбце
                for (int j = i; j <= n; ++j) { // Перебираем элементы в текущей строке (включая элемент вектора b)
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j]; // Вычитаем из текущей строки ведущую строку, умноженную на коэффициент
                }
                cout << "Eliminate row " << k << " using row " << i << endl; // Выводим сообщение об обнулении строки
                printMatrix(augmentedMatrix); // Выводим текущее состояние расширенной матрицы
                cout << endl; // Выводим пустую строку для разделения этапов
            }
        }
    }

    // Извлекаем решение из расширенной матрицы
    vector<double> solution(n); // Создаем вектор для хранения решения
    for (int i = 0; i < n; ++i) { // Перебираем строки
        solution[i] = augmentedMatrix[i][n]; // Копируем последний элемент каждой строки (элемент вектора b) в вектор решения
    }
    return solution; // Возвращаем вектор решения
}

int main() {
    // Задаем матрицу коэффициентов A
    vector<vector<double>> A = {
        {3, 1, -1},
        {2, 4, 1},
        {100, -1, 3}
    };

    // Задаем вектор свободных членов b
    vector<double> b = {6, 9, 4};

    // Выводим исходные данные
    cout << "Original Matrix A:" << endl;
    printMatrix(A); // Выводим матрицу A
    cout << "Vector b:" << endl;
    for (double val : b) { // Перебираем элементы вектора b
        cout << fixed << setprecision(2) << val << " "; // Выводим каждый элемент с фиксированной точкой и точностью 2 знака после запятой
    }
    cout << endl << endl; // Переходим на новую строку и выводим пустую строку для разделения этапов

    // Решаем систему уравнений методом Гаусса-Жордана
    vector<double> solution = gaussJordan(A, b); // Вызываем функцию gaussJordan для решения системы уравнений

    // Выводим решение
    cout << "Solution:" << endl;
    for (int i = 0; i < solution.size(); ++i) { // Перебираем элементы вектора решения
        cout << "x_" << i + 1 << " = " << fixed << setprecision(2) << solution[i] << endl; // Выводим каждый элемент решения с фиксированной точкой и точностью 2 знака после запятой
    }

    return 0; // Завершаем программу
}