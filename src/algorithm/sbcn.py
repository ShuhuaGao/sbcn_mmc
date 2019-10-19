"""
Switched Boolean Control Network

Author: Gao Shuhua
Date: 2019/10/18
"""


class SBCN:
    """
    Switched BCN
    """
    def __init__(self, n, m, w, Ls, g, Cx=None, Cu=None, Cs=None):
        """
        Initialize a switched Boolean network.

        :param n: number of state variables
        :param m: number of control inputs
        :param w: number of sub-systems
        :param Ls: a list of transition matrices, each for one sub-system
        :param Cx: state constraints, a set of ints or None (no constraints)
        :param Cu: control constraints, a functor: int --> a set of ints, or None
        :param Cs: switching constraints, a functor: int --> a set of ints, or None
        :param g: stage cost function: (int, int, int) --> float
        """
        self.n = n
        self._N = 2 ** n
        self.m = m
        self._M = 2 ** m
        self._all_controls = list(range(1, self._M + 1)) # if no constraints
        self.w = w
        self.g = g
        self._all_subs = list(range(1, self.w + 1))
        self.Ls = Ls
        self.Cx = Cx
        self.Cu = Cu
        self.Cs = Cs

    def compute_successors(self, i):
        """
        Get the succeeding states of i, i.e., R(i, 1), and the associated optimal weight, control, and switch

        :param i: current state
        :return: a list, each element being a tuple (j, weight, control, switch)
        """
        if self.Cx is not None:
            assert i in self.Cx, f"The given state {i} violates the state constraints"
        controls = self.Cu(i) if self.Cu else self._all_controls
        subs = self.Cs(i) if self.Cs else self._all_subs
        successors = {}
        # note that k, l, and i here start from 1
        for k in controls:
            for l in subs:
                L = self.Ls[l - 1]
                blk = L[(k - 1) * self._N: k * self._N ]
                j = blk[i - 1]
                if self.Cx is not None:  # the successor j must satisfy the state constraints
                    if j not in self.Cx:
                        continue
                weight = self.g(i, k, l)
                if j in successors:
                    if weight < successors[j][0]:
                        successors[j] = (weight, k, l)   # a better one
                else:
                    successors[j] = (weight, k, l)
        return [(j, *info) for j, info in successors.items()]