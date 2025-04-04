#include <iostream> // Подключение библиотеки iostream для ввода и вывода данных
#include <cmath>   // Подключение библиотеки cmath для математических функций (cos, sin, abs)
#include <vector>  // Подключение библиотеки vector для работы с динамическими массивами
#include <locale>

using namespace std; // Использование стандартного пространства имен std, чтобы не писать std:: перед cout, cin и т.д.

const double EPSILON = 1e-6; // Определение константы EPSILON, представляющей желаемую точность решения. 1e-6 эквивалентно 0.000001.
const int MAX_ITER = 100;    // Определение константы MAX_ITER, задающей максимальное количество итераций для метода Ньютона.


// Функции системы (определяем систему нелинейных уравнений)
double f1(double x, double y) {
    return cos(x + 0.5) + y - 0.8; // Первая функция системы уравнений
}

double f2(double x, double y) {
    return sin(y) + 2 * x - 1.6;   // Вторая функция системы уравнений
}

// Численные частные производные (метод конечных разностей)
double df1_dx(double x, double y, double h = 1e-5) {
    return (f1(x + h, y) - f1(x , y)) / h; // Производная f1 по x
}

double df1_dy(double x, double y, double h = 1e-5) {
    return (f1(x, y + h) - f1(x, y)) / h; // Производная f1 по y
}

double df2_dx(double x, double y, double h = 1e-5) {
    return (f2(x + h, y) - f2(x , y)) / h; // Производная f2 по x
}

double df2_dy(double x, double y, double h = 1e-5) {
    return (f2(x, y + h) - f2(x, y )) / h; // Производная f2 по y
}

// Решение системы линейных уравнений методом Гаусса
vector<double> solve_linear_system(vector<vector<double>> A, vector<double> B) {
    int n = A.size(); // Получаем размерность матрицы A (количество уравнений)
    for (int i = 0; i < n; i++) {
        // Поиск максимального элемента в столбце (для повышения устойчивости метода)
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[max_row][i])) {
                max_row = k; // Находим строку с максимальным элементом в текущем столбце
            }
        }
        // Перестановка строк (если необходимо)
        swap(A[i], A[max_row]); // Меняем местами текущую строку и строку с максимальным элементом
        swap(B[i], B[max_row]); // Аналогично меняем местами соответствующие элементы вектора B

        // Приведение к верхнетреугольному виду
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i]; // Вычисляем коэффициент для обнуления элементов под диагональю
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j]; // Вычитаем из текущей строки строку, умноженную на коэффициент
            }
            B[k] -= factor * B[i]; // Обновляем вектор B аналогично
        }
    }

    // Обратный ход (нахождение решения)
    vector<double> X(n); // Создаем вектор X для хранения решения
    for (int i = n - 1; i >= 0; i--) {
        X[i] = B[i]; // Присваиваем элементу X[i] соответствующее значение из B[i]
        for (int j = i + 1; j < n; j++) {
            X[i] -= A[i][j] * X[j]; // Вычитаем из X[i] сумму произведений элементов A[i][j] и X[j]
        }
        X[i] /= A[i][i]; // Делим X[i] на соответствующий диагональный элемент
    }

    return X; // Возвращаем вектор решения
}

// Метод Ньютона–Рафсона для решения системы нелинейных уравнений
vector<double> newton_raphson(double x0, double y0) {
    double x = x0, y = y0; // Инициализируем начальное приближение
    int iter = 0;           // Счетчик итераций

    while (iter < MAX_ITER) {
        // Вычисляем матрицу Якоби и вектор F
        vector<vector<double>> J = {
            {df1_dx(x, y), df1_dy(x, y)}, // Первая строка матрицы Якоби (численные частные производные f1)
            {df2_dx(x, y), df2_dy(x, y)}  // Вторая строка матрицы Якоби (численные частные производные f2)
        };
        vector<double> F = {-f1(x, y), -f2(x, y)}; // Вектор F, содержащий значения функций с противоположным знаком

        // Решаем систему линейных уравнений J * ΔX = F
        vector<double> delta = solve_linear_system(J, F); // Решаем систему линейных уравнений методом Гаусса

        // Обновляем решение
        x += delta[0]; // Обновляем значение x
        y += delta[1]; // Обновляем значение y

        // Проверяем условие завершения (достижение заданной точности)
        if (abs(delta[0]) < EPSILON && abs(delta[1]) < EPSILON) {
            break; // Выходим из цикла, если достигнута требуемая точность
        }

        iter++; // Увеличиваем счетчик итераций
    }

    return {x, y}; // Возвращаем найденное решение
}

int main() {
    // Устанавливаем локаль для поддержки русского языка
    setlocale(LC_ALL, "Russian");

    double x0 = 1.0, y0 = 1.0; // Начальное приближение для x и y
    vector<double> solution = newton_raphson(x0, y0); // Вызываем метод Ньютона-Рафсона для нахождения решения

    wcout << L"x = " << solution[0] << L", y = " << solution[1] << endl; // Выводим найденное решение

    // Проверка решения: подставляем найденные x и y в уравнения системы
    double f1_value = f1(solution[0], solution[1]); // Вычисляем значение первой функции
    double f2_value = f2(solution[0], solution[1]); // Вычисляем значение второй функции

    wcout << L"f1(x, y) = " << f1_value << endl; // Выводим значение f1(x, y)
    wcout << L"f2(x, y) = " << f2_value << endl; // Выводим значение f2(x, y)

    // Проверяем, насколько близки значения к нулю
    if (abs(f1_value) < EPSILON && abs(f2_value) < EPSILON) {
        wcout << L"Решение корректно: значения функций близки к нулю." << endl;
    } else {
        wcout << L"Решение некорректно: значения функций не близки к нулю." << endl;
    }

    return 0; // Успешное завершение программы
}