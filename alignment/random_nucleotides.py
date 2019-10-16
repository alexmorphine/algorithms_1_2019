import numpy as np


class RandomNucleotides:
    """ Генератор рандомных последовательностей """

    def __init__(self):

        # желаемые длины
        self.lenght = [100, 1000, 100000]

        # алфавит
        self.alphabeth = np.array(['A', 'T', 'G', 'C'])

    def generate(self):
        """ Запуск генерации """
        result = []

        # по одной последовательности на размер
        for size in self.lenght:
            indices = np.random.choice(self.alphabeth, size=size, p=[0.3, 0.2, 0.4, 0.1])
            result.append(''.join(indices))
        print(*result, sep='\n\n')


def run():
    """
    Запуск генерации
    """

    # Создаем экземпляр класса
    rn = RandomNucleotides()

    # получаем результаты
    rn.generate()


if __name__ == "__main__":
    run()
