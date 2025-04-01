#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

// Определяемая функция (пример)
double f(double x) {
    return  1/(1 - 0.49 * pow(sin(x), 2));
}

// Метод Монте-Карло для вычисления интеграла
double MonteCarlo(double a, double b, int n) {
    // a - нижний предел интегрирования
    // b - верхний предел интегрирования
    // n - количество случайных точек (проб)

    double max_f = 0.0;      // Максимальное значение функции на интервале [a, b]
    double min_f = 0.0;      // Минимальное значение функции на интервале [a, b]
    double height;           // Высота прямоугольника

    // Оценка максимума и минимума функции
    for (int i = 0; i < 1000; ++i) {
        double x = a + (double)rand() / RAND_MAX * (b - a);
        max_f = max(max_f, f(x));
        min_f = min(min_f, f(x));
    }

    height = max_f - min_f;

    int hits_positive = 0;   // Количество "попаданий" выше оси X
    int hits_negative = 0;   // Количество "попаданий" ниже оси X
    double x, y;

    // Генерируем случайные точки и считаем "попадания"
    for (int i = 0; i < n; ++i) {
        x = a + (double)rand() / RAND_MAX * (b - a);
        y = min_f + (double)rand() / RAND_MAX * height;

        if (y <= f(x)) {
            if (y >= 0) {
                hits_positive++;
            } else {
                hits_negative++;
            }
        }
    }

    // Вычисляем площадь прямоугольника и приближенное значение интеграла
    double area = (b - a) * height;
    double positive_area = area * (double)hits_positive / n;
    double negative_area = area * (double)hits_negative / n;

    double integral = positive_area - negative_area;
    return integral;
}

int main() {
    // Инициализация генератора случайных чисел
    srand(time(0));

    double a, b;          // Пределы интегрирования
    double n;          // Количество точек для каждой итерации
    int k;               // Количество итераций (вводится с клавиатуры)

    a = 0.0;
    b = 5.0;

    cout << "(k): ";
    cin >> k;  // Считываем количество итераций с клавиатуры

    cout << "(n): ";
    cin >> n;  // Считываем количество точек с клавиатуры

    double sum_s = 0.0;         // Сумма значений интеграла
    double sum_s_squared = 0.0;   // Сумма квадратов значений интеграла

    // Выполняем k итераций
    for (int i = 0; i < k; ++i) {
        // Вычисляем интеграл на текущей итерации
        double current_integral = MonteCarlo(a, b, n);

        // Обновляем суммы
        sum_s += current_integral;                  // Суммируем значения интеграла
        sum_s_squared += current_integral * current_integral;  // Суммируем квадраты значений интеграла

        cout << i + 1 << ": integral = " << current_integral << endl;
    }

    double average_integral = sum_s / k;          // Вычисляем среднее значение интеграла
    double average_squared = sum_s_squared / k;

    double standard_deviation = sqrt(average_squared - average_integral * average_integral); //Оцениваем стандартное отклонение

    cout << "integral: " << average_integral << endl;
    cout << "dev: " << standard_deviation << endl;

    return 0;
}