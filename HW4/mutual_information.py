from Bio import SeqIO
import pandas as pd
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from functools import lru_cache


class Fasta:
    """ Class for reading fasta files """

    def __init__(self, file, homo, related):
        """
        Args:
            file (str or list): path to file or list with reads
            homo (bool): use only homo sapiens RNA
            related (bool): use only related species RNA
        """

        # если передали строку
        if type(file) == str:

            # создаём словарь вид: RNA
            record_dict = {}

            # чтение из файла
            with open(file, "r") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    record_dict[record.id] = record.seq

            # в датасет, чтобы легче искать
            self.df = pd.DataFrame.from_dict(record_dict, orient='index', columns=['Seq'])

        # если передали список, то сразу в датафрейм
        elif type(file) == list:
            self.df = pd.DataFrame(file, columns=['Seq'])

        # если нужно только РНК человека
        if homo:
            self.homo()

        # если только РНК родственных видов (взяла два с того же сайта, что и РНК)
        if related:
            self.related()

        # максимальная длиина без пропусков
        print(f'To have no NaNs in fasta, please select x and y which are less than {self.min_lenght()}')

    def same_length(self, length=None):
        """
        Брать только записи одной длины

        Args:
            length (int): длина записей, которая нужна
        Returns:
            pd.DataFrame: записи нужной длины
        """

        # если не была задана длина, то берём наибольшую без пропусков
        length = length if length else self.min_lenght()
        return self.df[self.df.Seq.apply(lambda x: len(x) == length)]

    def get_columns(self, i, j):
        """
        Получение колонок i и j

        Args:
            i (int): номер колонки
            j (int): номер колонки
        Returns:
            pd.DataFrame: две колонки
        """
        columns = self.df.Seq.str[i].to_frame().merge(self.df.Seq.str[j], right_index=True, left_index=True)
        if columns.Seq_x.isnull().any():
            print("There are NaNs in column i, consider selecting lower index.")
        if columns.Seq_y.isnull().any():
            print("There are NaNs in column j, consider selecting lower index.")
        return columns

    def min_lenght(self):
        """
        Вычисление минимальной длины без пропусков

        Returns:
            int: длина
        """
        return self.df.Seq.str.len().min()

    def homo(self):
        """
        Использование только РНК человека
        """
        self.df = self.df.loc[self.df.index.str.startswith('Homo')]

    def related(self):
        """
        Использование только РНК родственных видов (два захардкожены)
        """
        self.df = self.df.loc[(self.df.index.str.startswith('Debaryomyces_hansenii')) | (
            self.df.index.str.startswith('Saccharomyces_cerevisiae'))].sort_index()


class MutualInformation:
    """ Класс для вычисления совместной информации """

    def __init__(self, filename, homo=False, related=False):
        """
        Args:
            filename (str or list): path to file or list with reads
            homo (bool): use only homo sapiens RNA
            related (bool): use only related species RNA
        """

        # РНК данные
        self.fasta = Fasta(filename, homo, related)

        # возможные х и у
        self.alphabet = ['A', 'T', 'G', 'C']

        # очистка кэша на всякий
        self.mu.cache_clear()

    def columnwise(self, columns, x, y):
        """
        Вычислениие совместной информации

        Args:
            columns (pd.DataFrame):
            x (str): буква алфавита
            y (str): буква алфавита
        Returns:
            float: результат вычисления
        """
        f_xy = len(columns.loc[(columns.Seq_x == x) & (columns.Seq_y == y)]) / len(columns)
        f_x = len(columns.loc[columns.Seq_x == x]) / len(columns.Seq_x)
        f_y = len(columns.loc[columns.Seq_y == y]) / len(columns.Seq_y)

        # чтобы не ловить ошибку из логарифма
        if f_xy:
            result = f_xy * np.log2(f_xy / (f_y * f_x))
        else:
            result = 0

        return result

    def mutual_information(self, i, j):
        """
        Пользовательский метод для совместной информации двух столбцов

        Args:
            i (int): номер колонки
            j (int): номер колонки
        """
        print(f'Mutual information for columns {i} and {j} is {self.mu(i, j)}')

    @lru_cache(maxsize=None)
    def mu(self, i, j):
        """
        Вычисление совместной информации для всех букв алфавита

        Args:
            i (int): номер колонки
            j (int): номер колонки

        Returns:
            float: совместная информация двух столбцов
        """

        # достаём нужные колонки
        columns = self.fasta.get_columns(i, j)

        # все возможные комбинации алфавита по 2 буквы
        comb = product(self.alphabet, repeat=2)

        # получаем СИ
        mu = 0
        for (x, y) in comb:
            mu += self.columnwise(columns, x, y)
        return mu

    def heatmap(self, title, columns=None):
        """
        Отрисовка хитмэпа

        Args:
            title (str): заголовок графика
            columns (int or list): либо порядковый номер последней используемой колонки: будут использованы от 0 до
            columns - 1 подряд, либо массив из номеров типа [0,  2, 4]. По умолчанию - максимальная длина без пропусков
        """
        plt.figure(figsize=(15, 15))

        # если получили список колонок
        if type(columns) == list:
            hm = np.zeros([len(columns), len(columns)])

            # то вычисляем для каждой пары, вытягивая номер из массива
            for n, i in enumerate(columns):
                j_columns = columns[n:]
                for k, j in enumerate(j_columns):
                    hm[n, k + n] = self.mu(i, j)
                    hm[k + n, n] = hm[n, k + n]

        # если получили int или ничего
        else:

            # берём либо этот int, либо число
            columns = columns if columns else self.fasta.min_lenght()
            min_len = min(self.fasta.min_lenght(), columns)
            hm = np.zeros([min_len, min_len])
            for i in range(min_len):
                for j in range(i, min_len):
                    hm[i, j] = hm[j, i] = self.mu(i, j)
            columns = [i for i in range(min_len)]
        plt.imshow(hm)
        plt.colorbar().set_label('Mutual Information', rotation=270, fontsize='x-large')
        plt.xticks(np.arange(len(hm)), columns)
        plt.yticks(np.arange(len(hm)), columns)
        plt.xlim(0 - 0.5, len(columns) - 0.5)
        plt.ylim(0 - 0.5, len(columns) - 0.5)
        plt.title(title, fontsize='xx-large')
        plt.tight_layout()
        plt.show()
