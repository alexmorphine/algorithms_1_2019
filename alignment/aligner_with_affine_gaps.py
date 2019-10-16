import numpy as np
import pandas as pd
import argparse


class AffineAligner:

    def __init__(self, seq1, seq2, match=1, gap_start=-1, gap_continued=-0.5,  mismatch=-1, weights=False):
        """
        Инициализация аргументов
        Args:
            seq1 (str): последовательность 1
            seq2 (str): последовательность 2
            match (float): вес совпадения
            gap (float): вес пропуска
            mismatch (float): вес несовпадения
            weights (str): название матрицы весов, может быть 'pam', 'blosum', а также False - весов нет
        """

        # делаем большими и красивыми
        self.seq1 = seq1.upper()
        self.seq2 = seq2.upper()
        self.match = np.float32(match)
        self.gap_start = np.float32(gap_start)
        self.gap_continued = np.float32(gap_continued)
        self.mismatch = np.float32(mismatch)

        # инициализируем матрицы
        self.middle = self.init_middle()
        self.upper = self.init_lower_upper('upper')
        self.lower = self.init_lower_upper('lower')

    def init_middle(self):
        # сначала инициализируем матрицу нулями по размерам посл1 + 1, посл2 + 1; + 1 за пропуск спереди
        matrix = np.zeros([len(self.seq1) + 1, len(self.seq2) + 1])
        matrix[1:, 0] = -np.inf
        matrix[0, 1:] = -np.inf
        return matrix

    def init_lower_upper(self, which):
        matrix = np.zeros([len(self.seq1) + 1, len(self.seq2) + 1])
        if which == 'lower':
            for i in range(matrix.shape[1]):
                matrix[0, i] = self.init_cell(i)
            matrix[1:, 0] = -np.inf
        else:
            for i in range(matrix.shape[0]):
                matrix[i, 0] = self.init_cell(i)
            matrix[0, 1:] = -np.inf
        return matrix

    def get_indeces(self, i, j):
        return (i - 1, j - 1), (i - 1, j), (i, j - 1)

    def init_cell(self, cnt):
        return self.gap_start + cnt * self.gap_continued

    def seq_match(self, i, j):
        return 1 if self.seq1[i - 1] == self.seq2[j - 1] else -1

    def align(self):
        for i in range(1, len(self.seq1) + 1):
            for j in range(1, len(self.seq2) + 1):
                middle, upper, lower = self.get_indeces(i, j)
                match = self.seq_match(i, j)

                self.upper[i, j] = max(self.middle[upper] + self.gap_start + self.gap_continued,
                                       self.upper[upper] + self.gap_continued)
                self.lower[i, j] = max(self.middle[lower] + self.gap_start + self.gap_continued,
                                       self.lower[lower] + self.gap_continued)
                self.middle[i, j] = max(self.middle[middle] + match,
                                        self.upper[middle] + match,
                                        self.lower[middle] + match)
        print('\nUPPER:')
        print(self.upper)

        print('\nMIDDLE:')
        print(self.middle)

        print('\nLOWER:')
        print(self.lower)

    def check_index(self, index):
        index = np.array(index)
        return np.any(index >= np.array([0, 0]))

    def select_alignment(self):
        # значение, из которого пойдём - правый нижний край
        current = tuple(np.array(self.middle.shape) - 1)

        # результат выравнивания
        self.alignment = [[], []]

        # конец (не интересует сочетание двух пропусков)
        end = np.array([1, 1])

        variants = np.array([self.middle[current], self.upper[current], self.lower[current]])
        value = np.where(variants == variants.max())[0]

        while np.any(current >= end):
            # запоминааем текущие индексы по осям
            first_index = current[0]
            second_index = current[1]

            # символы, которые сейчас сравниваем
            if first_index:
                seq1 = self.seq1[first_index - 1]
            else:
                self.alignment[0].append('-')
                self.alignment[1].append(self.seq2[second_index - 1])
                current = np.array([first_index, second_index - 1])
                continue

            if second_index:
                seq2 = self.seq2[second_index - 1]
            else:
                self.alignment[0].append(seq1)
                self.alignment[1].append('-')
                current = np.array([first_index - 1, second_index])
                continue

            middle, upper, lower = self.get_indeces(first_index, second_index)
            nonnegative = [x for x in map(self.check_index, [middle, upper, lower])]

            match = self.seq_match(first_index, second_index)
            for matrix in value:
                if nonnegative[matrix]:
                    if matrix == 0:

                        variants = np.array([self.middle[middle] + match,
                                             self.upper[middle] + match,
                                             self.lower[middle] + match])
                        middle_max = variants.max()
                        if self.middle[current] == middle_max:
                            value = np.where(variants == middle_max)[0]
                            self.alignment[0].append(seq1)
                            self.alignment[1].append(seq2)
                            current = middle
                            break
                    elif matrix == 1:
                        variants = np.array([self.middle[upper] + self.gap_start + self.gap_continued,
                                             self.upper[upper] + self.gap_continued])
                        upper_max = variants.max()
                        if self.upper[current] == upper_max:
                            value = np.where(variants == upper_max)[0]
                            self.alignment[0].append(seq1)
                            self.alignment[1].append('-')
                            current = upper
                            break
                    elif matrix == 2:
                        variants = np.array([self.middle[lower] + self.gap_start + self.gap_continued,
                                             self.lower[lower] + self.gap_continued])
                        lower_max = variants.max()

                        if self.lower[current] == lower_max:
                            value = np.where(variants == lower_max)[0]
                            self.alignment[0].append('-')
                            self.alignment[1].append(seq2)
                            current = lower
                            break

    def print_alignment(self):
        for seq in self.alignment:
            print('\t', *seq[::-1])


def parse_args(args=None):
    """ Функция, которая создает парсер для аргументов командной строки, и потом сама парсит аргументы и результат
    возвращает.
    Args:
        args: Параметры, которые необходимо распарсить вместо аргументов командной строки.
    Returns:
        Namespace: с параметрами, которые смогли распарсить.
    """

    parser = argparse.ArgumentParser(allow_abbrev=False, add_help=True, description="Alignment class")
    parser.add_argument("-seq1", "--seq1", required=True, default=None, type=str,
                        help="First sequence, any case.")
    parser.add_argument("-seq2", "--seq2", required=True, default=None, type=str,
                        help="Second sequence, any case.")
    parser.add_argument("-match", "--match", required=False, type=float, default=1,
                        help="Weight of match. Defaults to 1")

    parser.add_argument("-mismatch", "--mismatch", required=False, type=float, default=-1,
                        help="Weight of mismatch. Defaults to -1")

    parser.add_argument("-gap_start", "--gap_start", required=False, type=float, default=-1,
                        help="Weight of gap start. Defaults to -1")
    parser.add_argument("-gap_continued", "--gap_continued", required=False, type=float, default=-1,
                        help="Weight of gap continuation. Defaults to -1")

    parser.add_argument("-alignment", "--alignment", required=False, action='store_true', dest='alignment',
                        help="Print alignment")

    args = parser.parse_args(args)
    return args


def run(args=None):
    """
    Функция парсит аргументы командной строки и передает их в класс Aligner.
    Args:
        args: Аргументы, которые надо парсить. Если None, то парсим аргументы командной строки.
    """
    # Парсим аргументы командной строки или переданную строку с аргументами
    args = parse_args(args)

    print(f"Alignment started with parameters: {vars(args)}")

    # Создаем экземпляр модели
    # Передаем аргументы командной строки как параметры в метод __init__ модели
    aligner = AffineAligner(args.seq1, args.seq2, match=args.match, mismatch=args.mismatch, gap_start=args.gap_start,
                            gap_continued=args.gap_continued)

    # получаем матрицу
    aligner.align()

    # если надо вывести выравниваение, то выводим
    if args.alignment:
        aligner.select_alignment()
        aligner.print_alignment()


if __name__ == "__main__":
    run()

