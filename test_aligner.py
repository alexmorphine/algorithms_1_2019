import unittest
from aligner import Aligner
import numpy as np
import pandas as pd


class TestAligner(unittest.TestCase):

    def test_no_weight(self):
        """
        Тест для пустой матрицы веса
        проверочная матрица алгоритма взята с
        http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Needleman-Wunsch#
        """
        seq1, seq2 = 'ABCA', 'BDAAAC'

        aligner = Aligner(seq1=seq1, seq2=seq2)
        aligner.align()
        aligner.select_alignment()

        seq1, seq2 = set('-ABCA'), set('-BDAAAC')
        weights = pd.DataFrame(np.zeros([len(seq1), len(seq2)]), index=seq1, columns=seq2)
        matrix = pd.read_csv('aligner_no_weights_check.csv', delimiter=';', index_col=0, header=None, skiprows=1)

        self.assertEqual(np.array_equal(aligner.weights, weights.values), True)

        self.assertEqual(np.array_equal(aligner.matrix, matrix), True)
        self.assertEqual(aligner.alignment[0], ['-', 'A', 'C', '-', '-', 'B', 'A'])
        self.assertEqual(aligner.alignment[1], ['C', 'A', 'A', 'A', 'D', 'B', '-'])

    def test_with_weight(self):
        """
        Тест для не пустой матрицы веса
        """
        seq1, seq2 = 'ABCA', 'BDAAAC'

        aligner = Aligner(seq1=seq1, seq2=seq2, weights='pam')
        aligner.align()
        aligner.select_alignment()

        self.assertEqual(aligner.alignment[0], ['-', 'A', 'C', '-', 'B', 'A'])
        self.assertEqual(aligner.alignment[1], ['C', 'A', 'A', 'A', 'D', 'B'])


if __name__ == '__main__':
    unittest.main()
