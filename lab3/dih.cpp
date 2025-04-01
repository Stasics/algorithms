#include <iostream>
#include <cmath>
#include <iomanip>

// Функция, корень которой мы ищем
double f(double x) {
    return exp(-x * x) - (x - 1) * (x - 1);
}

// Метод дихотомии
double dichotomy(double a, double b, double epsilon) {
    double c;
    int iterations = 0; // Добавлено: Счетчик итераций
    while ((b - a) / 2 > epsilon) { // Пока интервал не станет меньше заданной точности
        c = (a + b) / 2; // Середина интервала
        if (f(c) == 0.0) { // Если нашли точный корень
            break;
        } else if (f(a) * f(c) < 0) { // Если корень в левой половине
            b = c;
        } else { // Если корень в правой половине
            a = c;
        }
        iterations++; // Добавлено: Увеличиваем счетчик на каждой итерации
    }
    std::cout << "iterations " << iterations << std::endl; // Добавлено: Выводим количество итераций
    return (a + b) / 2; // Возвращаем середину последнего интервала
}

int main() {
    double a = 1.0; // Начало интервала
    double b = 2.0; // Конец интервала
    double epsilon = 0.0001; // Точность
    double root = dichotomy(a, b, epsilon); // Находим корень

    // Выводим результат с точностью до 6 знаков после запятой
    std::cout << "answer: " << std::fixed << std::setprecision(4) << root << std::endl;
    return 0;
}