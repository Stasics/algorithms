#include <iostream>  // Подключение библиотеки для ввода и вывода данных (например, std::cout, std::cin)
#include <vector>    // Подключение библиотеки для работы с динамическими массивами (векторами)
#include <cmath>     // Подключение библиотеки для математических функций (например, pow)
#include <algorithm> // Подключение библиотеки для работы с алгоритмами (например, reverse)

using namespace std; // Использование стандартного пространства имен, чтобы не писать std:: перед функциями

// Функция для вычисления коэффициентов базисного полинома l_j(x)
vector<double> computeBasisPolynomial(const vector<double>& x, int j) {
    int n = x.size(); // Количество узлов интерполяции
    vector<double> coefficients(n, 0.0); // Вектор для хранения коэффициентов полинома l_j(x)
    coefficients[0] = 1.0; // Начальное значение: l_j(x) = 1

    // Цикл для построения базисного полинома l_j(x)
    for (int m = 0; m < n; ++m) {
        if (m != j) { // Пропускаем узел x_j, так как l_j(x_j) = 1
            // Умножаем текущий полином на (x - x_m) / (x_j - x_m)
            vector<double> temp(n, 0.0); // Временный вектор для хранения текущих коэффициентов
            for (int k = 0; k < n; ++k) {
                temp[k] = coefficients[k]; // Копируем текущие коэффициенты во временный вектор
            }
            // Обновляем коэффициенты полинома
            for (int k = 1; k < n; ++k) {
                coefficients[k] = temp[k - 1] - x[m] * temp[k]; // Умножение на (x - x_m)
            }
            coefficients[0] = -x[m] * temp[0]; // Обновляем нулевой коэффициент
            double denominator = x[j] - x[m]; // Вычисляем знаменатель (x_j - x_m)
            // Делим все коэффициенты на знаменатель
            for (int k = 0; k < n; ++k) {
                coefficients[k] /= denominator;
            }
        }
    }

    return coefficients; // Возвращаем коэффициенты базисного полинома l_j(x)
}

// Функция для вывода многочлена (от старшей к младшей)
void printPolynomial(vector<double> coefficients) {
    int n = coefficients.size(); // Количество коэффициентов
    reverse(coefficients.begin(), coefficients.end()); // Разворачиваем вектор коэффициентов (чтобы старшая степень шла первой)
    bool firstTerm = true; // Флаг для отслеживания первого члена многочлена

    // Цикл для вывода многочлена
    for (int i = 0; i < n; ++i) {
        int power = n - 1 - i; // Степень текущего члена

        if (coefficients[i] != 0) { // Пропускаем нулевые коэффициенты
            if (!firstTerm && coefficients[i] > 0) { // Добавляем "+" перед положительными членами (кроме первого)
                cout << " + ";
            } else if (coefficients[i] < 0) { // Добавляем "-" перед отрицательными членами
                cout << " - ";
            }

            cout << abs(coefficients[i]); // Выводим абсолютное значение коэффициента

            if (power > 0) { // Если степень больше нуля, добавляем "x"
                cout << "x";
                if (power > 1) { // Если степень больше 1, добавляем степень
                    cout << "^" << power;
                }
            }
            firstTerm = false; // После первого члена сбрасываем флаг
        }
    }
    cout << endl; // Переход на новую строку после вывода многочлена
}

int main() {
    // Входные данные: узлы интерполяции и значения функции
    vector<double> x = {1, 2, 3, 4, 5};  // Узлы интерполяции
    vector<double> y = {7.2, 5.9, 4.9, 4.0, 3.2};  // Значения функции

    int n = x.size(); // Количество узлов интерполяции
    vector<double> finalCoefficients(n, 0.0); // Вектор для хранения коэффициентов итогового многочлена

    // Вычисляем итоговый многочлен
    for (int j = 0; j < n; ++j) {
        vector<double> basisPoly = computeBasisPolynomial(x, j); // Вычисляем базисный полином l_j(x)
        for (int k = 0; k < n; ++k) {
            finalCoefficients[k] += y[j] * basisPoly[k]; // Добавляем вклад базисного полинома в итоговый многочлен
        }
    }

    // Выводим итоговый многочлен
    cout << "Lagrange: ";
    printPolynomial(finalCoefficients);

    // Вычисляем значение многочлена в точке x = 2.5
    double point = 1.5; // Точка, в которой вычисляем значение многочлена
    double interpolatedValue = 0.0; // Переменная для хранения результата
    vector<double> reversedCoefficients = finalCoefficients; // Создаем копию коэффициентов для разворота
    reverse(reversedCoefficients.begin(), reversedCoefficients.end()); // Разворачиваем коэффициенты (чтобы старшая степень шла первой)

    // Вычисляем значение многочлена в точке point
    for (int i = 0; i < n; ++i) {
        interpolatedValue += reversedCoefficients[i] * pow(point, n - 1 - i); // Суммируем члены многочлена
    }
    cout << "answer in " << point << " is " << interpolatedValue << endl; // Выводим результат

    return 0; // Успешное завершение программы
}