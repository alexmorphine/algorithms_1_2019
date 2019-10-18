import argparse


class Table:
    """ Класс для таблицы длин максимальных префиксов-суффиксов """

    def __init__(self, word):
        """
        Args:
            word (str): подстрока, которую ищем
        """
        self.word = word

        # сама таблица
        self.pf = self.prefix_function()

    def prefix_function(self):
        """
        Расчёт таблицы

        Returns:
            list: заполненный список с максимальными длинами
        """
        lenght = len(self.word)

        # инициализируем список
        prefix = [0] * lenght

        # на первом месте оставим 0
        for i in range(1, lenght):

            # смотрим на предыдущий префикс
            previous = prefix[i - 1]

            # если нам есть что проверять, то смотрим на совпадение текущего символа с символом, следующим за
            # максимальным текущим суффиксом
            while (previous > 0) & (self.word[i] != self.word[previous]):

                # если неудачно, то делаем шаг назад (получаем более короткую длину)
                previous = prefix[previous - 1]

            # если символ после макс длины и текущий совпадают
            if self.word[i] == self.word[previous]:

                # то увеличиваем макс длину суффикса
                previous += 1

            # записываем получившуюся длину суффикса для текущей позиции
            prefix[i] = previous
        return prefix

    def __repr__(self):
        """
        Красивая печать таблицы

        Returns:
            str: текст для вывода
        """
        word = ' '.join(list(self.word))
        table = ' '.join(map(str, self.pf))
        repr_str = f"{word}\n{table}"
        return repr_str

    def __getitem__(self, key):
        """
        Делаем из объекта класса контейнер, чтобы сразу к таблице обращаться

        Args:
            key (int): позиция в массиве

        Returns:
            int: максимальная длина суффикса
        """
        return self.pf[key]


class KnuthMorrisPratt:
    """
        Класс, реализующий поиск подстроки в строке с помощью алгоритма Кнута-Морриса-Пратта"""

    def __init__(self, string, substring):
        """

        Args:
            string (str): строка, в которой ищем
            substring (str): подстрока, которую ищем
        """
        self.string = string.lower()
        self.substring = substring.lower()

        # таблица суффиксов
        self.table = Table(substring)

    def search(self):
        """
        Поиск подстроки в строке

        Returns:
            int: позиция первого вхождения
        """

        # счётчик по подстроке
        match = 0

        # идём по строке, в которой ищем
        for start in range(len(self.string)):

            # если нет совпадения, откатываемся на суффикс
            while (match > 0) & (self.string[start] != self.substring[match]):
                match = self.table[match - 1]

            # если совпали, то идём дальше
            if self.string[start] == self.substring[match]:
                match += 1

            # если дошли  до конца подстроки, то нашли совпадение
            if match == len(self.substring):
                print(f'Match found at position {start - len(self.substring) + 2}')
                return start - len(self.substring) + 2

        print('Match not found')
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
    args = parser.parse_args(args)
    return args


def run(args=None):
    """
    Функция парсит аргументы командной строки и передает их в класс.
    Args:
        args: Аргументы, которые надо парсить. Если None, то парсим аргументы командной строки.
    """
    # Парсим аргументы командной строки или переданную строку с аргументами
    args = parse_args(args)

    print(f"Search started with parameters: {vars(args)}")

    # Создаем экземпляр класса
    # Передаем аргументы командной строки как параметры в метод __init__
    rk = KnuthMorrisPratt(args.string, args.substring)

    # получаем матрицу
    rk.search()


if __name__ == "__main__":
    run()


