import numpy as np
import argparse


class RabinKarp:
    """
    Класс, реализующий поиск подстроки в строке с помощью алгоритма Рабина-Карпа"""

    def __init__(self, string, substring, x=2, q=5):
        """

        Args:
            string (str): строка, в которой ищем
            substring (str): подстрока, которую ищем
            x (int): множитель для кольцевой функции хэширования
            q (int): делитель (mod) для кольцевой функции хэширования
        """
        self.x = x
        self.q = q
        self.string = string.lower()
        self.substring = substring.lower()

    def get_letter(self, char):
        """
        Получаем ASCII код символа

        Args:
            char (str): символ

        Returns:
            int: ASCII код
        """
        return ord(char)

    def hash(self, sub, letter=None, i=None):
        """
        Функция хэширования

        Args:
            sub (int): значение хэшфункции, если letter != None, иначе - str, для которой получаем хэш
            letter (str): буква, которую хотим добавить в хэш
            i (int): позиция первого символа подстроки у sub

        Returns:
            int: хэш
        """
        if letter is None:
            return self._hash(sub)

        # вычисяем результат по формуле получения кольцевого хэша
        result = (sub - self.hash(self.string[i]) * self.x ** (len(self.substring) - 1)) * self.x + self._hash(letter)
        return np.mod(result, self.q)

    def _hash(self, sub):
        """
        Вычисление функции хэширования

        Args:
            sub (str): подстрока, для которой получаем хэш

        Returns:
            int: хэш
        """
        len_sub = len(sub)

        # список из кодов букв sub
        letters = list(map(self.get_letter, list(sub)))

        # степени, в которые будем возводить x
        powers = [_ for _ in range(1, len_sub + 1)]

        # получаем сумму для дальнейшей формулы
        hash_sum = np.sum(np.array([a * self.x ** (len_sub - i) for a, i in zip(letters, powers)]))
        return np.mod(hash_sum, self.q)

    def search(self):
        """
        Поиск подстроки в строке

        Returns:
            int: положение в строке
        """
        substring_len = len(self.substring)

        # получаем первый хэш для подстроки и части строки
        hpattern = self.hash(self.substring)
        hs = self.hash(self.string[:substring_len])

        # вычисляем, до какого числа итерируемся
        total = len(self.string) - substring_len
        i = 0

        while i <= total:

            # если хэши совпали, проверяем, совпали ли строки (мб коллизия)
            if hs == hpattern:
                if self.string[i:substring_len + i] == self.substring:
                    print(f'Substring found starting on {i}\n')
                    return i

            # если ещё не выходим за длину строки, то продолжаем искать
            if i + substring_len < len(self.string):
                hs = self.hash(hs, self.string[i + substring_len], i)
            i += 1
        print('Substring not found\n')
        return -1


def parse_args(args=None):
    """ Функция, которая создает парсер для аргументов командной строки, и потом сама парсит аргументы и результат
    возвращает.
    Args:
        args: Параметры, которые необходимо распарсить вместо аргументов командной строки.
    Returns:
        Namespace: с параметрами, которые смогли распарсить.
    """

    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=True, description="Alignment class")
    parser.add_argument("-string", "--string", required=True, default=None, type=str,
                        help="Str to search in")
    parser.add_argument("-substring", "--substring", required=True, default=None, type=str,
                        help="substring to search.")
    parser.add_argument("-x", "--x", required=False, type=int, default=2,
                        help="X for hashing")

    parser.add_argument("-q", "--q", required=False, type=int, default=5,
                        help="Q for hashing")

    args = parser.parse_args(args)
    return args


def run(args=None):
    """
    Функция парсит аргументы командной строки и передает их в класс RabinKarp.
    Args:
        args: Аргументы, которые надо парсить. Если None, то парсим аргументы командной строки.
    """
    # Парсим аргументы командной строки или переданную строку с аргументами
    args = parse_args(args)

    print(f"Alignment started with parameters: {vars(args)}")

    # Создаем экземпляр класса
    # Передаем аргументы командной строки как параметры в метод __init__
    rk = RabinKarp(args.string, args.substring, x=args.x, q=args.q)

    # получаем матрицу
    rk.search()


if __name__ == "__main__":
    run()

