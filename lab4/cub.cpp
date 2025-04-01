#include <iostream>  // Подключение библиотеки для ввода/вывода
#include <vector>    // Подключение библиотеки для работы с векторами
#include <cmath>     // Подключение библиотеки для математических функций
#include <iomanip>   // Подключение библиотеки для форматирования вывода
#include <functional>

using namespace std; // Использование стандартного пространства имен

// Функции для интерполяции
double func_1(double x) { return exp(x); }  // Экспонента
double func_2(double x) { return sinh(x); } // Гиперболический синус
double func_3(double x) { return cosh(x); } // Гиперболический косинус
double func_4(double x) { return sin(x); }  // Синус
double func_5(double x) { return cos(x); }  // Косинус
double func_6(double x) { return log(x); }  // Натуральный логарифм
double func_7(double x) { return exp(-x); } // Экспонента с отрицательным аргументом

// Метод прогонки для решения трехдиагональной системы
vector<double> tridiagonalSolve(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) {
    int n = a.size();  // Размер системы
    vector<double> m(n);  // Вектор для хранения решения
    vector<double> lambda(n), mu(n);  // Векторы для прогоночных коэффициентов

    // Прямая прогонка
    lambda[0] = -c[0] / b[0];  // Начальное значение lambda
    mu[0] = d[0] / b[0];  // Начальное значение mu
    for (int i = 1; i < n; ++i) {
        double denominator = b[i] + a[i] * lambda[i - 1];  // Вычисление знаменателя
        lambda[i] = -c[i] / denominator;  // Вычисление lambda[i]
        mu[i] = (d[i] - a[i] * mu[i - 1]) / denominator;  // Вычисление mu[i]
    }

    // Обратная прогонка
    m[n - 1] = mu[n - 1];  // Начальное значение решения
    for (int i = n - 2; i >= 0; --i) {
        m[i] = lambda[i] * m[i + 1] + mu[i];  // Вычисление решения
    }

    return m;  // Возврат решения
}

// Построение кубического сплайна
vector<vector<double>> cubicSpline(const vector<double>& x, const vector<double>& y) {
    int n = x.size() - 1;  // Количество интервалов
    vector<double> h(n);  // Вектор для хранения длин интервалов
    for (int i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i];  // Вычисление длины каждого интервала
    }

    // Векторы для системы уравнений
    vector<double> a(n - 1), b(n - 1), c(n - 1), d(n - 1);

    for (int i = 1; i < n; ++i) {
        a[i - 1] = h[i - 1] / 6.0;  // Коэффициент a
        b[i - 1] = (h[i - 1] + h[i]) / 3.0;  // Коэффициент b
        c[i - 1] = h[i] / 6.0;  // Коэффициент c
        d[i - 1] = (y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1];  // Коэффициент d
    }

    // Решение системы методом прогонки
    vector<double> m = tridiagonalSolve(a, b, c, d);
    m.insert(m.begin(), 0.0); // m0 = 0
    m.push_back(0.0); // mn = 0

    // Вычисление коэффициентов сплайна
    vector<vector<double>> coefficients(n, vector<double>(4));
    for (int i = 0; i < n; ++i) {
        coefficients[i][0] = y[i];  // Коэффициент a
        coefficients[i][1] = (y[i + 1] - y[i]) / h[i] - h[i] * (m[i + 1] + 2 * m[i]) / 6.0;  // Коэффициент b
        coefficients[i][2] = m[i] / 2.0;  // Коэффициент c
        coefficients[i][3] = (m[i + 1] - m[i]) / (6.0 * h[i]);  // Коэффициент d
    }

    return coefficients;  // Возврат коэффициентов сплайна
}

// Вычисление значения сплайна в точке
double evaluateSpline(double x, const vector<double>& x_values, const vector<vector<double>>& coefficients) {
    int n = x_values.size() - 1;  // Количество интервалов
    int i = 0;
    while (i < n && x > x_values[i + 1]) {  // Поиск интервала, в который попадает x
        ++i;
    }
    double dx = x - x_values[i];  // Разность между x и началом интервала
    return coefficients[i][0] + coefficients[i][1] * dx + coefficients[i][2] * dx * dx + coefficients[i][3] * dx * dx * dx;  // Вычисление значения сплайна
}

int main() {
    double a = 1.00, b = 1.20, h = 0.04;  // Определение интервала и шага
    vector<double> x_values;  // Вектор для хранения значений x
    for (double x = a; x <= b + 1e-9; x += h) {
        x_values.push_back(x);  // Заполнение вектора x_values
    }

    // Выбор функции
    int function_choice;
    cout << "Choose a function to interpolate (1-7):" << endl;
    cin >> function_choice;
    if (function_choice < 1 || function_choice > 7) {
        cerr << "Invalid function choice. Exiting." << endl;
        return 1;  // Завершение программы с ошибкой
    }

    vector<function<double(double)>> functions = { func_1, func_2, func_3, func_4, func_5, func_6, func_7 };  // Вектор функций
    function<double(double)> selected_function = functions[function_choice - 1];  // Выбор функции

    vector<double> y_values;  // Вектор для хранения значений y
    for (double x : x_values) {
        y_values.push_back(selected_function(x));  // Заполнение вектора y_values
    }

    // Построение сплайна
    vector<vector<double>> coefficients = cubicSpline(x_values, y_values);

    // Точки для вычисления
    vector<double> evaluation_points = { 1.05, 1.09, 1.13, 1.15, 1.17 };
    cout << fixed << setprecision(6);  // Установка точности вывода
    for (double x : evaluation_points) {
        double spline_value = evaluateSpline(x, x_values, coefficients);  // Вычисление значения сплайна
        cout << "Spline value at x = " << x << ": " << spline_value << endl;
        cout << "Actual value at x = " << x << ": " << selected_function(x) << endl;
        cout << "Error is: " << abs(spline_value - selected_function(x)) << endl;  // Вычисление ошибки
    }

    return 0;  // Завершение программы
}