#include <iostream>  // Подключение библиотеки iostream для ввода и вывода данных.
#include <cmath>     // Подключение библиотеки cmath для математических функций, таких как exp, abs, и nan.

// Функция f(x) - определяет функцию, корень которой мы хотим найти.
double f(double x) {
    return exp(-x * x) - (x - 1) * (x - 1); // Возвращает значение функции e^(-x^2) - (x - 1)^2.
}

// Численная производная функции f'(x) - вычисляет производную с использованием приращения.
double numerical_derivative(double x, double h = 1e-5) {
    return (f(x + h) - f(x)) / h; // Возвращает приближенное значение производной.
}

// Метод Ньютона - реализует метод Ньютона (или метод касательных) для нахождения корня функции.
double newton_method(double x0, double epsilon) {
    double x = x0;          // Инициализация переменной x начальным приближением x0.
    double delta;           // Переменная для хранения изменения x на каждой итерации.
    int iterations = 0;     // Переменная для подсчета количества итераций.

    do {
        double fx = f(x);   // Вычисление значения функции f(x) в текущей точке x.
        double dfx = numerical_derivative(x);  // Вычисление значения производной f'(x) в текущей точке x.

        if (dfx == 0) {
            std::cerr << "Производная равна нулю. Метод не может быть применен." << std::endl; // Вывод сообщения об ошибке, если производная равна нулю.
            return NAN;      // Возвращает NAN (Not a Number), если производная равна нулю, т.к. деление на ноль невозможно.
        }

        delta = fx / dfx;  // Вычисление изменения x (delta) по формуле: delta = f(x) / f'(x).
        x = x - delta;     // Обновление значения x: x = x - delta.
        iterations++;      // Увеличение счетчика итераций.

    } while (std::abs(delta) >= epsilon); // Продолжает цикл до тех пор, пока абсолютное значение delta не станет меньше epsilon (заданной точности).

    std::cout << "iterations: " << iterations << std::endl; // Выводит количество итераций, потребовавшихся для нахождения корня.
    return x; // Возвращает найденный корень.
}

int main() {
    double x0 = -1.0;         // Начальное приближение для поиска корня.
    double epsilon = 0.0001;   // Желаемая точность нахождения корня.

    double root = newton_method(x0, epsilon); // Вызов функции newton_method для нахождения корня.

    if (!std::isnan(root)) {
        std::cout << "root: " << root << std::endl; // Вывод найденного корня, если он был найден (не является NAN).
    } else {
        std::cout << "error" << std::endl;  // Вывод сообщения, если корень не был найден (возвращено NAN).
    }

    return 0; // Успешное завершение программы.
}