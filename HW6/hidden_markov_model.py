import pandas as pd
import numpy as np


class HMM:
    """ HMM """

    def __init__(self, pi, A, B, seq, states):
        """
        Args:
            pi (dict): вероятности перехода в скрытое состояние
            A (list): вероятности перехода между скрытыми состояниями
            B (list): вероятности перехода в открытое состояние из скрытого
            seq (str): последовательность
            states (list): возможные скрытые состояния
        """
        self.pi = pi
        self.states = states
        self.alphabet = set(seq)
        self.B = pd.DataFrame(data=B, index=sorted(states), columns=sorted(self.alphabet))
        self.A = pd.DataFrame(data=A, index=sorted(states), columns=sorted(states))
        self.seq = seq

        # alpha
        self.forward_matrix = pd.DataFrame(index=sorted(states),
                                           columns=[y for y in range(len(seq))]).astype(float)

        # beta
        self.backward_matrix = pd.DataFrame(index=sorted(states), columns=[y for y in range(len(seq))]).astype(float)

        self.viterbi_matrix = pd.DataFrame(index=sorted(states), columns=range(len(seq))).astype(float)

    def viterbi(self):
        """
        Алгоритм Витерби
        """
        for i, letter in enumerate(self.seq):
            for state in self.states:
                if i == 0:
                    self.viterbi_matrix.loc[state, i] = self.pi[state] * self.B.loc[state, letter]
                else:
                    greatest = []
                    for other_state in self.states:
                        greatest.append(self.viterbi_matrix.loc[other_state, i - 1] * self.A.loc[state, other_state])

                    self.viterbi_matrix.loc[state, i] = self.B.loc[state, letter] * max(greatest)

    def forward(self):
        """
        Алгоритм вычисления альфа для FB
        """
        for i, letter in enumerate(self.seq):
            result = {}

            for state in self.states:
                if i == 0:

                    result[state] = self.pi[state] * self.B.loc[state, letter]
                else:
                    result[state] = 0
                    for other_state in self.states:
                        result[state] += self.forward_matrix.loc[other_state, i - 1] \
                                        * self.A.loc[other_state, state] * self.B.loc[state, letter]

            self.forward_matrix[i] = pd.Series(result)

    def backward(self):
        """
        Алгоритм вычисления бета для FB
        """
        for i, letter in enumerate(self.seq[::-1]):
            number = len(self.seq) - i - 1
            for state in self.states:
                if not i:
                    self.backward_matrix.loc[state, number] = 1
                else:
                    self.backward_matrix.loc[state, number] = (self.backward_matrix.loc[state, number + 1]
                                                               * self.A.loc[state, :]
                                                               * self.B.loc[:, self.seq[number + 1]]).sum()

    def FB(self):
        """
        Результат алгоритма FB

        Returns:
            pd.DataFrame: результат
        """
        self.forward()
        self.backward()
        return self.forward_matrix * self.backward_matrix / self.forward_matrix.values[:, -1].sum()
