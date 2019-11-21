import numpy as np


class Cluster:
    """ Abstract """

    def __init__(self, matrix, objects=None):
        if type(matrix) == list:
            self.matrix = np.array(matrix).astype(float)
        else:
            self.matrix = matrix
        np.fill_diagonal(self.matrix, np.inf)
        self.objects = self.make_objects(objects)
        self.current_min = None
        self.__linkage = None

    @property
    def linkage(self):
        """
        Newick Standard

        Returns:
            str: linkage
        """
        if self.__linkage is None:
            self.cluster()
        return self.objects[0][0]

    def __repr__(self):
        return self.linkage

    def make_objects(self, objects):
        objects = objects or [_ for _ in range(len(self.matrix))]
        return [objects, [0] * len(objects), [1] * len(objects)]

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
            keep:
            drop:
            distance ():

        """
        dist1 = dist2 = distance

        obj_first = self.objects[0][drop]
        obj_first_dist = self.objects[1][drop]
        obj_second = self.objects[0][keep]
        obj_second_dist = self.objects[1][keep]

        self.objects[0][keep] = f'({obj_first}: {dist1 - obj_first_dist}, {obj_second}: {dist2 - obj_second_dist})'
        self.objects[1][keep] = distance

        del self.objects[0][drop], self.objects[1][drop]

    def distance(self):
        """

        Returns:

        """
        return self.matrix[self.current_min] / 2

    def new_node(self, *args):
        """

        Args:
            *args:
        """
        raise NotImplementedError

    def change_matrix(self, distances):

        drop, keep = self.current_min

        self.set_link(distances, drop, keep)

        idx = set([_ for _ in range(len(self.matrix))]) - set(self.current_min)

        for ind in idx:
            self.matrix[keep, ind] = self.matrix[ind, keep] = self.new_node(ind)

        self.matrix = np.delete(self.matrix, drop, 0)
        self.matrix = np.delete(self.matrix, drop, 1)
        self.objects[2][keep] += self.objects[2][drop]
        del self.objects[2][drop]

    def cluster(self):

        while self.matrix.shape != (1, 1):
            self.current_min = self.find_smallest(self.matrix)

            distances = self.distance()

            self.change_matrix(distances)


class WPGMA(Cluster):
    """ WPGMA """

    def new_node(self, ind):
        """

        Args:
            ind:

        Returns:

        """
        i, j = self.current_min
        return (self.matrix[i, ind] + self.matrix[j, ind]) / 2


class UPGMA(Cluster):
    """ UPGMA """

    def new_node(self, ind):
        """

        Args:
            ind:

        Returns:

        """

        drop, keep = self.current_min
        members_keep, members_drop = self.objects[2][keep], self.objects[2][drop]

        return (self.matrix[keep, ind] * members_keep + self.matrix[drop, ind] * members_drop) /\
               (members_keep + members_drop)


class NJ(Cluster):

    def __init__(self, matrix, objects=None):
        if type(matrix) == list:
            self.matrix = np.array(matrix).astype(float)
        else:
            self.matrix = matrix
        self.objects = self.make_objects(objects)
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
        return self.objects[0][0]

    def distance_to_other(self):
        if len(self.matrix) > 2:
            return self.matrix.sum(axis=0) / (len(self.matrix) - 2)
        return self.matrix.sum(axis=0) / 2

    def set_link(self, distance, drop, keep):
        """
        Красивый вывод

        Args:
            keep:
            drop:
            distance ():

        """
        dist1, dist2 = distance

        obj_first = self.objects[0][drop]
        obj_second = self.objects[0][keep]

        self.objects[0][keep] = f'({obj_first}: {dist1}, {obj_second}: {dist2})'

        del self.objects[0][drop], self.objects[1][drop]

    def cluster(self):

        while self.matrix.shape != (1, 1):
            self.avg_distance = self.distance_to_other()
            self.min_distance_matrix = np.zeros(self.matrix.shape)
            for i in range(len(self.matrix)):
                for j in range(i + 1, len(self.matrix)):
                    self.min_distance_matrix[i, j] = self.min_distance_matrix[j, i] = self.matrix[i, j] - \
                                                                                      self.avg_distance[i] - \
                                                                                      self.avg_distance[j]
            np.fill_diagonal(self.min_distance_matrix, np.inf)
            self.change_matrix()

    def change_matrix(self):

        drop, keep = self.current_min = self.find_smallest(self.min_distance_matrix)

        v_keep = self.matrix[keep, drop] / 2 + (self.avg_distance[keep] - self.avg_distance[drop]) / 2
        v_drop = self.matrix[keep, drop] / 2 + (self.avg_distance[drop] - self.avg_distance[keep]) / 2

        self.set_link((v_drop, v_keep), drop, keep)
        idx = set([_ for _ in range(len(self.matrix))]) - set(self.current_min)

        for ind in idx:
            self.matrix[keep, ind] = self.matrix[ind, keep] = self.new_node(ind)

        self.matrix = np.delete(self.matrix, drop, 0)
        self.matrix = np.delete(self.matrix, drop, 1)

    def last_step(self):
        dist = self.matrix[0, 1] / 2
        obj_first = self.objects[0][1]
        obj_second = self.objects[0][0]

        self.objects[0][0] = f'({obj_first}: {dist}, {obj_second}: {dist})'
        del object[0][1]

    def distance(self, i, j):
        return self.matrix[i, j] - self.avg_distance[i] - self.avg_distance[j]

    def new_node(self, ind):
        i, j = self.current_min
        return (self.matrix[i, ind] + self.matrix[j, ind] - self.matrix[i, j]) / 2

