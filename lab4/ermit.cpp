#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip> // Для форматирования вывода

using namespace std;

// Структура для хранения коэффициентов сплайна
struct HermiteCoeff {
    double a, b, c, d;
};

// Функция для построения сплайна Эрмита
vector<HermiteCoeff> buildHermiteSpline(const vector<double>& x, const vector<double>& y, const vector<double>& y_prime) {
    int n = x.size() - 1; // Количество интервалов
    vector<HermiteCoeff> coefficients(n);

    for (int i = 0; i < n; ++i) {
        double h = x[i + 1] - x[i]; // Длина интервала
        double dy = y[i + 1] - y[i]; // Разность значений функции

        coefficients[i].a = y[i];
        coefficients[i].b = y_prime[i];
        coefficients[i].c = (3 * dy - h * (y_prime[i + 1] + 2 * y_prime[i])) / (h * h);
        coefficients[i].d = (-2 * dy + h * (y_prime[i + 1] + y_prime[i])) / (h * h * h);
    }

    return coefficients;
}

// Функция для вычисления значения сплайна в точке x
double evaluateHermiteSpline(double x, const vector<double>& x_values, const vector<HermiteCoeff>& coefficients) {
    int n = x_values.size() - 1;
    int i = 0;

    // Находим интервал, в который попадает x
    while (i < n && x > x_values[i + 1]) {
        ++i;
    }

    double dx = x - x_values[i]; // Разность между x и началом интервала
    const HermiteCoeff& coeff = coefficients[i];

    // Вычисляем значение сплайна
    return coeff.a + coeff.b * dx + coeff.c * dx * dx + coeff.d * dx * dx * dx;
}

// Функция для вывода сплайна
void printHermiteSpline(const vector<double>& x, const vector<HermiteCoeff>& coefficients) {
    int n = x.size() - 1;

    cout << "Hermite Spline:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "Interval [" << x[i] << ", " << x[i + 1] << "]:" << endl;
        cout << "S(x) = "
             << coefficients[i].a << " + "
             << coefficients[i].b << " * (x - " << x[i] << ") + "
             << coefficients[i].c << " * (x - " << x[i] << ")^2 + "
             << coefficients[i].d << " * (x - " << x[i] << ")^3" << endl;

        // Вывод значений сплайна в нескольких точках на интервале
        cout << "Values on the interval:" << endl;
        for (double t = x[i]; t <= x[i + 1]; t += 0.1) {
            double value = evaluateHermiteSpline(t, x, coefficients);
            cout << "S(" << t << ") = " << value << endl;
        }
        cout << endl;
    }
}

int main() {
    // Пример данных
    vector<double> x = {0.0, 1.0, 2.0, 3.0}; // Узлы интерполяции
    vector<double> y = {0.0, 1.0, 0.0, -1.0}; // Значения функции в узлах
    vector<double> y_prime = {1.0, 0.0, -1.0, 0.0}; // Значения производных в узлах

    // Построение сплайна Эрмита
    vector<HermiteCoeff> coefficients = buildHermiteSpline(x, y, y_prime);

    // Вывод построенного сплайна
    printHermiteSpline(x, coefficients);

    // Вычисление значения сплайна в точке
    double point = 0.5;
    double value = evaluateHermiteSpline(point, x, coefficients);

    cout << "Hermite spline value at x = " << point << " is " << value << endl;

    return 0;
}