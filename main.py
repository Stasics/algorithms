# Импорт необходимых библиотек
import numpy as np
import matplotlib.pyplot as plt


# Определение функции косинуса 5x
def cosine_5x(t):
    """
    Генерирует косинусоидальный сигнал cos(5x)
    :param t: массив временных точек
    :return: массив значений сигнала
    """
    return np.cos(5 * np.asarray(t))


# Реализация метода Симпсона для численного интегрирования с контролем точности
def simpson_integral(f, a, b, n=1000, tol=1e-6, max_iter=10):
    """
    Вычисляет определенный интеграл методом Симпсона с контролем точности
    :param f: интегрируемая функция
    :param a: нижний предел интегрирования
    :param b: верхний предел интегрирования
    :param n: начальное количество интервалов разбиения
    :param tol: требуемая точность (по умолчанию 1e-6)
    :param max_iter: максимальное количество итераций уточнения
    :return: значение интеграла, достигнутая точность
    """
    if n % 2 != 0:
        n += 1

    # Первое вычисление
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = f(x)
    integral_prev = h / 3 * (y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-1:2]))

    # Итеративное уточнение
    for _ in range(max_iter):
        n *= 2  # Удваиваем количество интервалов
        h = (b - a) / n
        x_new = np.linspace(a, b, n + 1)
        # Используем уже вычисленные точки
        y_new = np.zeros(n + 1)
        y_new[::2] = y  # Четные индексы - старые точки
        y_new[1::2] = f(x_new[1::2])  # Нечетные - новые точки

        integral = h / 3 * (y_new[0] + y_new[-1] + 4 * np.sum(y_new[1:-1:2]) + 2 * np.sum(y_new[2:-1:2]))

        # Проверка точности
        error = abs(integral - integral_prev)
        if error < tol:
            return integral, error

        integral_prev = integral
        y = y_new

    return integral_prev, error


# Класс для работы с рядом Фурье с контролем точности
class FourierSeries:
    def __init__(self, func, period, num_harmonics, integral_tol=1e-6):
        """
        Инициализация ряда Фурье
        :param func: исходная периодическая функция
        :param period: период функции
        :param num_harmonics: количество учитываемых гармоник
        :param integral_tol: точность интегрирования (по умолчанию 1e-6)
        """
        self.func = func
        self.T = period
        self.N = num_harmonics
        self.omega = 2 * np.pi / period
        self.integral_tol = integral_tol

        self._compute_coefficients()
        self._compute_errors()

    def _compute_coefficients(self):
        """Вычисление коэффициентов ряда Фурье с контролем точности"""
        # Вычисляем a0 с оценкой точности
        self.a0, self.a0_error = simpson_integral(
            self.func, -self.T / 2, self.T / 2, tol=self.integral_tol)
        self.a0 /= self.T

        self.an = np.zeros(self.N)
        self.bn = np.zeros(self.N)
        self.an_errors = np.zeros(self.N)
        self.bn_errors = np.zeros(self.N)

        for n in range(1, self.N + 1):
            # Вычисляем an с оценкой точности
            self.an[n - 1], self.an_errors[n - 1] = simpson_integral(
                lambda t: self.func(t) * np.cos(n * self.omega * t),
                -self.T / 2, self.T / 2, tol=self.integral_tol)
            self.an[n - 1] *= 2 / self.T

            # Вычисляем bn с оценкой точности
            self.bn[n - 1], self.bn_errors[n - 1] = simpson_integral(
                lambda t: self.func(t) * np.sin(n * self.omega * t),
                -self.T / 2, self.T / 2, tol=self.integral_tol)
            self.bn[n - 1] *= 2 / self.T

    def _compute_errors(self):
        """Вычисление оценок ошибок для коэффициентов"""
        self.total_error = np.sqrt(
            (self.a0_error / self.T) ** 2 +
            np.sum((2 / self.T * self.an_errors) ** 2) +
            np.sum((2 / self.T * self.bn_errors) ** 2))

    def evaluate(self, t, return_error=False):
        """
        Вычисление значения ряда Фурье в заданных точках
        :param t: массив моментов времени
        :param return_error: если True, возвращает также оценку ошибки
        :return: массив значений ряда Фурье (и ошибки, если требуется)
        """
        t = np.asarray(t)
        result = self.a0 / 2 * np.ones_like(t)

        for n in range(1, self.N + 1):
            result += (self.an[n - 1] * np.cos(n * self.omega * t) +
                       self.bn[n - 1] * np.sin(n * self.omega * t))

        if return_error:
            # Оценка ошибки - сумма ошибок всех компонент
            error_estimate = (abs(self.a0_error / 2) +
                              np.sum(np.abs(self.an_errors)) +
                              np.sum(np.abs(self.bn_errors)))
            return result, error_estimate
        return result

    def get_coefficients_info(self):
        """Возвращает информацию о коэффициентах и их точности"""
        info = f"Fourier Series Coefficients (T={self.T}, N={self.N}):\n"
        info += f"a0 = {self.a0:.6f} ± {self.a0_error:.2e}\n"
        for n in range(1, self.N + 1):
            info += (f"a{n} = {self.an[n - 1]:.6f} ± {self.an_errors[n - 1]:.2e}, "
                     f"b{n} = {self.bn[n - 1]:.6f} ± {self.bn_errors[n - 1]:.2e}\n")
        info += f"Total estimated error: {self.total_error:.2e}"
        return info


# Основная функция программы
def main():
    # Параметры анализа
    T = 2 * np.pi  # Период для cos(5x) (2π/5 * 5 = 2π)
    N = 6
    t_start = -7
    t_end = 7
    integral_tol = 1e-8  # Точность интегрирования

    # Создаем объект ряда Фурье с контролем точности
    fourier = FourierSeries(cosine_5x, T, N, integral_tol)

    # Выводим информацию о коэффициентах
    print(fourier.get_coefficients_info())

    # Генерируем данные для построения графиков
    t = np.linspace(t_start, t_end, 1000)
    y_original = cosine_5x(t)
    y_fourier, y_error = fourier.evaluate(t, return_error=True)

    # Построение графиков
    plt.figure(figsize=(12, 6))

    # Основные графики
    plt.plot(t, y_original, label='cos(5x)', linewidth=2)
    plt.plot(t, y_fourier, label=f'Ряд Фурье ({N} гармоник)', linewidth=1.5)

    # Область ошибки
    plt.fill_between(t, y_fourier - y_error, y_fourier + y_error,
                     color='gray', alpha=0.2, label='Оценка ошибки')

    plt.title('Разложение cos(5x) в ряд Фурье', fontsize=14)
    plt.xlabel('Время', fontsize=12)
    plt.ylabel('Амплитуда', fontsize=12)
    plt.legend(fontsize=12)

    # Настройка сетки и осей
    plt.xticks(np.arange(t_start, t_end + 0.1, 0.5))
    plt.yticks(np.arange(-1.2, 1.3, 0.1))
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.ylim(-1.2, 1.2)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()