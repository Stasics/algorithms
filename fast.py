# Импорт необходимых библиотек
import numpy as np  # Для численных операций и работы с массивами
import matplotlib.pyplot as plt  # Для визуализации данных
from scipy.fft import fft, fftfreq  # Для быстрого преобразования Фурье (FFT)


# Функция генерации прямоугольного сигнала (меандра)
def square_wave(t, period=2.0, duty_cycle=0.5):
    """
    Генерирует прямоугольный сигнал с заданным периодом и скважностью
    t - массив временных точек
    period - период сигнала
    duty_cycle - коэффициент заполнения (0.5 для симметричного меандра)
    """
    # Вычисляем фазу сигнала (нормированную от 0 до 1)
    phase = (np.asarray(t) % period) / period
    # Возвращаем 1.0 когда фаза меньше duty_cycle, иначе -1.0
    return np.where(phase < duty_cycle, 1.0, -1.0)

"""def cosine_5x(t):
    Функция косинуса с частотой 5 Г
    return np.cos(5 * 2 * np.pi * t)  # 5 - частота в Гц, 2π - перевод в радиан"""
# Реализация численного интегрирования методом Симпсона
def simpson_integral(f, a, b, n=1000, tol=1e-6, max_iter=10):
    """
    Вычисляет интеграл функции f от a до b с контролем точности
    Алгоритм:
    1. Разбивает интервал на n подынтервалов
    2. Применяет формулу Симпсона
    3. Удваивает количество интервалов пока не достигнет заданной точности
    """
    # Метод Симпсона требует четного числа интервалов
    if n % 2 != 0:
        n += 1

    # Первое вычисление интеграла
    h = (b - a) / n  # шаг интегрирования
    x = np.linspace(a, b, n + 1)  # точки разбиения
    y = f(x)  # значения функции в точках
    # Формула Симпсона (1-4-2-4-...-2-4-1)
    integral_prev = h / 3 * (y[0] + y[-1] + 4 * np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-1:2]))

    # Итеративное уточнение результата
    for _ in range(max_iter):
        n *= 2  # удваиваем число интервалов
        h = (b - a) / n
        x_new = np.linspace(a, b, n + 1)

        # Используем уже вычисленные точки для оптимизации
        y_new = np.zeros(n + 1)
        y_new[::2] = y  # старые точки
        y_new[1::2] = f(x_new[1::2])  # новые точки

        # Пересчет интеграла
        integral = h / 3 * (y_new[0] + y_new[-1] + 4 * np.sum(y_new[1:-1:2]) + 2 * np.sum(y_new[2:-1:2]))

        # Проверка достижения требуемой точности
        error = abs(integral - integral_prev)
        if error < tol:
            return integral, error

        integral_prev = integral
        y = y_new

    return integral_prev, error


# Класс для работы с рядом Фурье
class FourierSeries:
    def __init__(self, func, period, num_harmonics, integral_tol=1e-6):
        """
        Инициализация ряда Фурье для заданной функции
        func - анализируемая функция
        period - период функции
        num_harmonics - количество учитываемых гармоник
        integral_tol - точность интегрирования
        """
        self.func = func
        self.T = period
        self.N = num_harmonics
        self.omega = 2 * np.pi / period  # круговая частота
        self.integral_tol = integral_tol

        self._compute_coefficients()  # Вычисление коэффициентов Фурье
        self._compute_errors()  # Оценка ошибок

    def _compute_coefficients(self):
        """Вычисление коэффициентов ряда Фурье"""
        # Постоянная составляющая (a0)
        self.a0, self.a0_error = simpson_integral(
            self.func, -self.T / 2, self.T / 2, tol=self.integral_tol)
        self.a0 /= self.T  # нормировка

        # Коэффициенты для гармоник
        self.an = np.zeros(self.N)
        self.bn = np.zeros(self.N)
        self.an_errors = np.zeros(self.N)
        self.bn_errors = np.zeros(self.N)

        for n in range(1, self.N + 1):
            # Косинусные коэффициенты (an)
            self.an[n - 1], self.an_errors[n - 1] = simpson_integral(
                lambda t: self.func(t) * np.cos(n * self.omega * t),
                -self.T / 2, self.T / 2, tol=self.integral_tol)
            self.an[n - 1] *= 2 / self.T  # нормировка

            # Синусные коэффициенты (bn)
            self.bn[n - 1], self.bn_errors[n - 1] = simpson_integral(
                lambda t: self.func(t) * np.sin(n * self.omega * t),
                -self.T / 2, self.T / 2, tol=self.integral_tol)
            self.bn[n - 1] *= 2 / self.T  # нормировка

    def _compute_errors(self):
        """Вычисление суммарной ошибки коэффициентов"""
        # Вычисляем общую ошибку как корень из суммы квадратов ошибок
        self.total_error = np.sqrt(
            (self.a0_error / self.T)**2 +
            np.sum((2 / self.T * self.an_errors)**2) +
            np.sum((2 / self.T * self.bn_errors)**2))

    def evaluate(self, t, return_error=False):
        """Вычисление значения ряда Фурье в заданных точках"""
        t = np.asarray(t)
        result = self.a0 / 2 * np.ones_like(t)  # Постоянная составляющая

        # Добавление вкладов гармоник
        for n in range(1, self.N + 1):
            result += (self.an[n - 1] * np.cos(n * self.omega * t) +
                      self.bn[n - 1] * np.sin(n * self.omega * t))

        if return_error:
            error_estimate = (abs(self.a0_error / 2) +
                            np.sum(np.abs(self.an_errors)) +
                            np.sum(np.abs(self.bn_errors)))
            return result, error_estimate
        return result

    def get_coefficients_info(self):
        """Возвращает строку с информацией о коэффициентах"""
        info = f"Fourier Series Coefficients (T={self.T}, N={self.N}):\n"
        info += f"a0 = {self.a0:.6f} ± {self.a0_error:.2e}\n"
        for n in range(1, self.N + 1):
            info += (f"a{n} = {self.an[n-1]:.6f} ± {self.an_errors[n-1]:.2e}, "
                    f"b{n} = {self.bn[n-1]:.6f} ± {self.bn_errors[n-1]:.2e}\n")
        info += f"Total estimated error: {self.total_error:.2e}"
        return info

# Класс для анализа FFT
class FFTAnalyzer:
    def __init__(self, func, period, num_samples=1024):
        """
        Инициализация FFT анализатора
        func - анализируемая функция
        period - период функции
        num_samples - количество отсчетов
        """
        self.func = func
        self.T = period
        self.N = num_samples

        # Дискретизация сигнала (2 периода для уменьшения краевых эффектов)
        self.t = np.linspace(0, 2 * period, 2 * num_samples, endpoint=False)
        self.signal = func(self.t)

        # Вычисление FFT
        self.fft_result = fft(self.signal)
        self.frequencies = fftfreq(2 * num_samples, d=period / num_samples)
        self.amplitudes = 2 / (2 * num_samples) * np.abs(self.fft_result)  # нормировка

    def reconstruct_signal(self, t, max_harmonics):
        """Реконструкция сигнала по ограниченному числу гармоник"""
        reconstructed = np.zeros_like(t, dtype=float)
        base_freq = 1 / self.T  # основная частота

        for k in range(len(self.fft_result)):
            freq = self.frequencies[k]
            harmonic_num = round(freq / base_freq)  # номер гармоники

            if 0 < harmonic_num <= max_harmonics:
                amplitude = self.amplitudes[k]  # амплитуда
                phase = np.angle(self.fft_result[k])  # фаза
                reconstructed += amplitude * np.cos(2 * np.pi * freq * t + phase)

        return reconstructed


# Основная функция сравнения методов
def compare_methods():
    # Параметры анализа
    T = 6.0  # период сигнала
    N_harmonics = 20  # количество гармоник
    t_start = -7  # начало временного интервала
    t_end = 7  # конец временного интервала

    # Создание временной оси
    t = np.linspace(t_start, t_end, 1000)
    y_original = square_wave(t, T)  # исходный сигнал

    # Инициализация анализаторов
    fourier = FourierSeries(lambda t: square_wave(t, T), T, N_harmonics)
    fft_analyzer = FFTAnalyzer(lambda t: square_wave(t, T), T)

    # Реконструкция сигналов
    y_fourier = fourier.evaluate(t)  # ряд Фурье
    y_fft = fft_analyzer.reconstruct_signal(t, N_harmonics)  # FFT

    # Построение графиков
    plt.figure(figsize=(14, 6))

    # График сигналов
    plt.subplot(1, 2, 1)
    plt.plot(t, y_original, label='Оригинальный сигнал', linewidth=2)
    plt.plot(t, y_fourier, '--', label=f'Ряд Фурье ({N_harmonics} гармоник)')
    plt.plot(t, y_fft, ':', label=f'FFT ({N_harmonics} гармоник)')
    plt.title('Сравнение методов реконструкции')
    plt.xlabel('Время')
    plt.ylabel('Амплитуда')
    plt.legend()
    plt.grid(True)
    plt.xticks(np.arange(t_start, t_end + 0.1, 0.5))  # шаг 0.5 по оси X
    plt.yticks(np.arange(-1.2, 1.3, 0.1))  # шаг 0.1 по оси Y

    # График спектров
    plt.subplot(1, 2, 2)
    # Спектр Фурье
    harmonics = range(1, N_harmonics + 1)
    fourier_amps = [np.sqrt(fourier.an[n - 1] ** 2 + fourier.bn[n - 1] ** 2) for n in harmonics]
    plt.stem(harmonics, fourier_amps, 'b', markerfmt='bo', label='Ряд Фурье')

    # Спектр FFT
    fft_harmonics = [round(freq * fft_analyzer.T) for freq in fft_analyzer.frequencies]
    fft_amps = []
    for h in harmonics:
        mask = [fh == h for fh in fft_harmonics]
        if any(mask):
            fft_amps.append(np.mean(fft_analyzer.amplitudes[mask]))

    plt.stem(harmonics, fft_amps, 'r', markerfmt='ro', label='FFT')
    plt.title('Сравнение гармоник')
    plt.xlabel('Номер гармоники')
    plt.ylabel('Амплитуда')
    plt.legend()
    plt.grid(True)
    plt.xticks(range(0, N_harmonics + 1, 2))  # шаг 2 по оси X
    plt.yticks(np.arange(0, 1.3, 0.1))  # шаг 0.1 по оси Y

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    compare_methods()