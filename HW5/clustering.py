import numpy as np


class Objects:
    """ Класс для хранения информации об объектах и кластерах """

    def __init__(self, objects=None):
        """
        Args:
            objects (int or list): либо количество объектов, либо их названия
        """
        # либо сами генерим - номера
        if type(objects) == int:
            self.names = [_ for _ in range(objects)]

        # либо присваиваем те, что пришли
        else:
            self.names = objects

        # начальное количество объектов
        self.number = len(self.names)

        self.__distances = None
        self.__members = None

    @property
    def distances(self):
        """
        Ленивая инициализация расстояний
        """
        if self.__distances is None:
            self.__distances = [0] * self.number
        return self.__distances

    @property
    def members(self):
        """
        Ленивая инициализация количества объектов
        """
        if self.__members is None:
            self.__members = [1] * self.number
        return self.__members


class Cluster:
    """ Abstract """

    def __init__(self, matrix, objects=None):
        """
        Args:
            matrix (list or np.array): матрица расстояний
            objects (list): имена объектов
        """
        if type(matrix) == list:
            self.matrix = np.array(matrix).astype(float)
        else:
            self.matrix = matrix

        # заполняем диагональные элементы матрицы +бесконечностями, чтобы они не выбирались минимальными
        np.fill_diagonal(self.matrix, np.inf)

        if objects:
            self.objects = Objects(objects)
        else:
            self.objects = Objects(len(matrix))
        self.current_min = None
        self.__linkage = None

    @property
    def linkage(self):
        """
        Newick Standard - красивая строка

        Returns:
            str: linkage
        """
        if self.__linkage is None:
            self.cluster()
        return self.objects.names[0]

    def __repr__(self):
        return self.linkage

    def find_smallest(self, matrix):
        """
        Нахождение минимального индекса в массиве

        Returns:
            tuple: минимальный индекс
        """
        return np.unravel_index(np.argmin(matrix, axis=None), matrix.shape)

    def set_link(self, distance, drop, keep):
        """
        Красивый вывод

        Args:
            keep (int): номер объекта, вместо которого будем писать
            drop (int): номер объекта, который удаляется
            distance (float): расстояние между объектами
        """
        # расстояния равны
        dist1 = dist2 = distance

        # тут получаем имена объектов для записи и расстояния, которые уже имеются до них
        obj_first = self.objects.names[drop]
        obj_first_dist = self.objects.distances[drop]
        obj_second = self.objects.names[keep]
        obj_second_dist = self.objects.distances[keep]

        # переписываем объект
        self.objects.names[keep] = f'({obj_first}: {dist1 - obj_first_dist}, {obj_second}: {dist2 - obj_second_dist})'

        # обновляем его расстояние
        self.objects.distances[keep] = distance

        # удаляем ненужный
        del self.objects.names[drop], self.objects.distances[drop]

    def distance(self):
        """
        Расстояние между объединяемыми объектами

        Returns:
            float: расстояние до каждого
        """
        return self.matrix[self.current_min] / 2

    def new_node(self, *args):
        """
        Обновление узла
        """
        raise NotImplementedError

    def change_matrix(self, distances):
        """
        Обновление матрицы расстояний

        Args:
            distances (float): расстояния до каждого объединяемого узла
        """
        drop, keep = self.current_min

        # объединяем два узла
        self.set_link(distances, drop, keep)

        # обновляем узлы для всех остальных объектов
        idx = set([_ for _ in range(len(self.matrix))]) - set(self.current_min)

        for ind in idx:
            self.matrix[keep, ind] = self.matrix[ind, keep] = self.new_node(ind)

        # удаляем ненужный объект из матрицы
        self.matrix = np.delete(self.matrix, drop, 0)
        self.matrix = np.delete(self.matrix, drop, 1)

        # обновляем кол-во элементов кластера
        self.objects.members[keep] += self.objects.members[drop]
        del self.objects.members[drop]

    def cluster(self):
        """
        Основная функция для кластеризации
        """
        # пока не останется 1 элемент матрицы
        while self.matrix.shape != (1, 1):

            # находим наименьший
            self.current_min = self.find_smallest(self.matrix)

            # получаем расстояния до объединяемых объектов
            distances = self.distance()

            # меняем матрицу расстояний
            self.change_matrix(distances)


class WPGMA(Cluster):
    """ WPGMA """

    def new_node(self, ind):
        """
        Обновление узла

        Args:
            ind (int): номер объекта, расстояние до которого обновляем

        Returns:
            float: новое расстояние
        """
        i, j = self.current_min
        return (self.matrix[i, ind] + self.matrix[j, ind]) / 2


class UPGMA(Cluster):
    """ UPGMA """

    def new_node(self, ind):
        """
        Обновление узла

        Args:
            ind (int): номер объекта, расстояние до которого обновляем

        Returns:
            float: новое расстояние
        """

        drop, keep = self.current_min

        # количества объектов кластеров
        members_keep, members_drop = self.objects.members[keep], self.objects.members[drop]

        return (self.matrix[keep, ind] * members_keep + self.matrix[drop, ind] * members_drop) /\
               (members_keep + members_drop)


class NJ(Cluster):
    """ NJ """

    def __init__(self, matrix, objects=None):
        if type(matrix) == list:
            self.matrix = np.array(matrix).astype(float)
        else:
            self.matrix = matrix
        if objects:
            self.objects = Objects(objects)
        else:
            self.objects = Objects(len(matrix))
        self.current_min = None
        self.__linkage = None
        self.avg_distance = None
        self.min_distance_matrix = np.zeros(self.matrix.shape)

    @property
    def linkage(self):
        """
        Newick Standard

        Returns:
            str: linkage
        """
        if self.__linkage is None:
            self.cluster()
        return self.objects.names[0]

    def distance_to_other(self):
        """
        Расстояния до остальных объектов

        Returns:
            list: список расстояний
        """

        # если достаточно объектов
        if len(self.matrix) > 2:
            return self.matrix.sum(axis=0) / (len(self.matrix) - 2)

        # иначе просто делим на 2
        return self.matrix.sum(axis=0) / 2

    def set_link(self, distance, drop, keep):
        """
        Красивый вывод

        Args:
            keep (int): номер объекта, вместо которого будем писать
            drop (int): номер объекта, который удаляется
            distance (float): расстояние между объектами
        """

        # тут расстояния разные
        dist1, dist2 = distance

        obj_first = self.objects.names[drop]
        obj_second = self.objects.names[keep]

        self.objects.names[keep] = f'({obj_first}: {dist1}, {obj_second}: {dist2})'

        del self.objects.names[drop], self.objects.distances[drop]

    def cluster(self):
        """
        Основная функция для кластеризации
        """

        # пока больше одного объекта в матрице расстояний
        while self.matrix.shape != (1, 1):

            # рассчитываем расстояния до объектов
            self.avg_distance = self.distance_to_other()

            # формируем матрицу минимальных расстояний
            self.min_distance_matrix = np.zeros(self.matrix.shape)

            # заполняем её по формуле
            for i in range(len(self.matrix)):
                for j in range(i + 1, len(self.matrix)):
                    self.min_distance_matrix[i, j] = self.min_distance_matrix[j, i] = self.matrix[i, j] - \
                                                                                      self.avg_distance[i] - \
                                                                                      self.avg_distance[j]

            # заполняем диагонали +бесконечностями, чтобы случайно не выбрать их
            np.fill_diagonal(self.min_distance_matrix, np.inf)

            # обновляем матрицу расстояний
            self.change_matrix()

    def change_matrix(self):
        """
        Обновление матрицы расстоний
        """

        # для удобства - тут минимальное значение из матрицы минимальных расстояний
        drop, keep = self.current_min = self.find_smallest(self.min_distance_matrix)

        # расстояния до объединяемых кластеров
        v_keep = self.matrix[keep, drop] / 2 + (self.avg_distance[keep] - self.avg_distance[drop]) / 2
        v_drop = self.matrix[keep, drop] / 2 + (self.avg_distance[drop] - self.avg_distance[keep]) / 2

        self.set_link((v_drop, v_keep), drop, keep)

        # обновление остальных весов
        idx = set([_ for _ in range(len(self.matrix))]) - set(self.current_min)

        for ind in idx:
            self.matrix[keep, ind] = self.matrix[ind, keep] = self.new_node(ind)

        self.matrix = np.delete(self.matrix, drop, 0)
        self.matrix = np.delete(self.matrix, drop, 1)

    def distance(self, i, j):
        """
        Расстояние между объединяемыми объектами

        Args:
            i (int): номер первого кластера
            j (int): номер второго кластера

        Returns:
            float: расстояние
        """
        return self.matrix[i, j] - self.avg_distance[i] - self.avg_distance[j]

    def new_node(self, ind):
        """
        Обновление узла

        Args:
            ind (int): номер объекта, расстояние до которого обновляем

        Returns:
            float: новое расстояние
        """
        i, j = self.current_min
        return (self.matrix[i, ind] + self.matrix[j, ind] - self.matrix[i, j]) / 2

