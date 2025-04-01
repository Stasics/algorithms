#include <iostream>
#include <cmath>
#include <iomanip>

// Функция, минимум которой мы ищем
double f(double x) {
    return (x + 1) * (x + 1); // Пример: f(x) = (x + 1)^2
}

// Метод обратной квадратичной интерполяции
double inverseQuadraticInterpolation(double x0, double x1, double x2, double epsilon) {
    double x3; // Новое значение x
    int iterations = 0; // Счетчик итераций

    while (true) {
        // Вычисляем новое значение x с использованием формулы обратной квадратичной интерполяции
        double numerator = x0 * x0 * (f(x1) - f(x2)) + x1 * x1 * (f(x2) - f(x0)) + x2 * x2 * (f(x0) - f(x1));
        double denominator = 2 * (x0 * (f(x1) - f(x2)) + x1 * (f(x2) - f(x0)) + x2 * (f(x0) - f(x1)));

        // Проверка на деление на ноль
        if (denominator == 0) {
            std::cerr << "0 in denom" << std::endl;
            return NAN;
        }

        x3 = numerator / denominator;

        // Проверка на достижение точности
        if (std::abs(x3 - x2) < epsilon) {
            std::cout << "minimum " << iterations << " iterations" << std::endl;
            return x3;
        }

        // Обновляем точки для следующей итерации
        x0 = x1;
        x1 = x2;
        x2 = x3;

        iterations++;

        // Защита от бесконечного цикла
        if (iterations > 1000) {
            std::cerr << "max iter." << std::endl;
            return x3;
        }
    }
}

int main() {
    double x0, x1, x2; // Три начальные точки
    double epsilon;    // Точность

    // Ввод трех значений x
    std::cout << " x0, x1, x2: ";
    std::cin >> x0 >> x1 >> x2;

    // Ввод точности
    std::cout << " (epsilon): ";
    std::cin >> epsilon;

    // Поиск минимума
    double result = inverseQuadraticInterpolation(x0, x1, x2, epsilon);

    // Вывод результата
    if (!std::isnan(result)) {
        std::cout << std::fixed << std::setprecision(10);
        std::cout << " x = " << result << std::endl;
        std::cout << " f(x) = " << f(result) << std::endl;
    }

    return 0;
}