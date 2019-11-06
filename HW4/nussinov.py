import numpy as np


class Nussinov:
    """ Класс с алгоритмом Нуссинов
    https://www.youtube.com/watch?v=h60raqnvm0s - алгоритм
    http://cgoliver.com/2017/01/15/Nussinov.html - функция для печати"""

    def __init__(self, seq, min_distance):
        """
        Args:
            seq (str): последовательность РНК
            min_distance (int): минимальное расстояние между основаниями для образования пары
        """
        self.seq = seq
        self.n = len(seq)
        self.OPT_matrix = np.zeros([self.n] * 2)
        # self.traceback_matrix = [[0] * self.n] * self.n
        self.min_distance = min_distance
        self.structure = []
        self.__result = None

    @property
    def result(self):
        """
        Property, чтобы не вызывать run. Используется в __repr__

        Returns:
            str: структура
        """
        if self.__result is None:
            self.run()
        return self.structure

    def OPT(self, i, j):
        """
        Функция для построения матрицы DP

        Args:
            i (int): индекс i
            j (int): индекс j

        Returns:
            int: максимальное значение для новой клетки
        """
        # условие на миниимальную дистанцию между основаниями
        if i >= j - self.min_distance:
            return 0

        # если хвосты не с чем соединять
        unpaired = self.OPT(i, j - 1)

        # если есть с чем
        paired = [1 + self.OPT(i, t - 1) + self.OPT(t + 1, j - 1) for t in range(i, j - self.min_distance)
                  if self.match(t, j)]

        # для работы max
        paired.append(0)
        return max(max(paired), unpaired)

    def match(self, t, j):
        """
        Сочетаемость оснований

        Args:
            t (int): индекс t в seq
            j (int): индекс j в seq

        Returns:
            bool: сочетаемы или нет
        """

        tj_set = {self.seq[t], self.seq[j]}
        if tj_set in [{'A', 'U'}, {'C', 'G'}]:
            return True
        return False

    def traceback(self, i, j):
        """
        Обратный проход для получения структуры

        Args:
            i (int): индекс i
            j (int): индекс j
        """
        # если закончили проход
        if j <= i:
            return

        # если скор такой же, как у соседней клетки (не было соединений), то идём в следующую
        elif self.OPT_matrix[i, j] == self.OPT_matrix[i, j - 1]:
            self.traceback(i, j - 1)

        else:
            # если было соединение
            for k in [m for m in range(i, j - self.min_distance) if self.match(m, j)]:
                if k - 1 < 0:

                    # если нашли пару, от которой была получена текущая клетка, шли по диагонали
                    if self.OPT_matrix[i, j] == self.OPT_matrix[k + 1, j - 1] + 1:

                        # добавляем её
                        self.structure.append((k, j))

                        # выясняем, откуда взялась она
                        self.traceback(k + 1, j - 1)
                        break

                # аналогично, из суммы нижних и левых вариантов
                elif self.OPT_matrix[i, j] == self.OPT_matrix[i, k - 1] + self.OPT_matrix[k + 1, j - 1] + 1:

                    self.structure.append((k, j))
                    self.traceback(i, k - 1)
                    self.traceback(k + 1, j - 1)
                    break

    def run(self):
        for k in range(self.min_distance, self.n):
            for i in range(self.n - k):
                j = i + k
                self.OPT_matrix[i, j] = self.OPT_matrix[j, i] = self.OPT(i, j)
        self.traceback(0, self.n - 1)

    def __repr__(self):
        """
        Печать структуры

        Returns:
            str: структура
        """
        dot_bracket = ["." for _ in range(self.n)]

        for s in self.result:
            dot_bracket[min(s)] = "("
            dot_bracket[max(s)] = ")"

        return "".join(dot_bracket)

