import unittest
from alignment.aligner import Aligner
import numpy as np
import pandas as pd


class TestAligner(unittest.TestCase):

    def test_no_weight_1(self):
        """
        Тест для пустой матрицы веса
        проверочная матрица алгоритма взята с
        http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Needleman-Wunsch#
        """
        seq1, seq2 = 'ABCA', 'BDAAAC'

        # создаём экземпляр
        aligner = Aligner(seq1=seq1, seq2=seq2)
        aligner.align()
        aligner.select_global_alignment()

        # нулевые веса нужного размера
        seq = set('-ABCA' + '-BDAAAC')
        weights = pd.DataFrame(np.zeros([len(seq), len(seq)]), index=seq, columns=seq)

        # матрица с сайта
        matrix = pd.read_csv('aligner_no_weights_check.csv', delimiter=';', index_col=0, header=None, skiprows=1)

        # проверяем, что веса есть и они нулевые
        self.assertEqual(np.array_equal(aligner.weights, weights.values), True)

        # матрица с сайта совпала с посчитанной
        self.assertEqual(np.array_equal(aligner.matrix, matrix), True)

        # есть ли пропуски
        has_gaps = ('-' in aligner.alignment[0]) or ('-' in aligner.alignment[1])

        # есть ли несовпадения
        has_mismatches = False
        for first, second in zip(aligner.alignment[0], aligner.alignment[1]):
            if (first != second) and (first != '-') and (second != '-'):
                has_mismatches = True
                break

        # и то, и другое должно быть
        self.assertEqual(has_gaps, True)
        self.assertEqual(has_mismatches, True)

    def test_no_weight_2(self):
        """
        Тест для проверки работоспособности с float. Должно работать
        """
        seq1, seq2 = 'ABCA', 'BDAAAC'

        aligner = Aligner(seq1=seq1, seq2=seq2, gap=-0.499)
        aligner.align()
        self.assertEqual(aligner.select_global_alignment(), None)

    def test_with_weight(self):
        """
        Тест для не пустой матрицы веса
        Должно измениться выравнивание со сменой веса в матрице
        """
        seq1, seq2 = 'ABCA', 'BDAAAC'

        aligner = Aligner(seq1=seq1, seq2=seq2, weights='pam')
        aligner.align()
        aligner.select_global_alignment()
        first_alignment = aligner.alignment

        # меняем веса в матрице сразу в двух местах, она же симметричная
        aligner.weights.loc['B', 'D'] = -10
        aligner.weights.loc['D', 'B'] = -10

        # пересчитываем выравнивание
        aligner.align()
        aligner.select_global_alignment()
        second_alignment = aligner.alignment

        # хотим, чтобы были разные
        self.assertNotEqual(second_alignment, first_alignment)


if __name__ == '__main__':
    unittest.main()
