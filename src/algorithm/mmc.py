"""
Infinite-horizon optimal control of SBCNs with the MMC algorithm.

In this program, a logical vector (i.e., a vector with only one entry 1 and all others 0) $\delta_n^i$ is represented
by a single integer i, and a logical matrix is represented by a list of integers, each for one column.

Author: Gao Shuhua
Date: 2019/10/11
"""
from .sbcn import SBCN
import collections, operator
Vertex = collections.namedtuple('Vertex', 'index, weight, control, switch')
Info = collections.namedtuple('Info', 'v_star k_star mmc pt p_star alpha beta')

_INF = float('inf')


class OSTG:
    """
    The optimal state transition graph.

    Refer to Algorithm 1.
    """
    def __init__(self, n, i0: int):
        """
        Create an empty OSTG.

        :param n: number of state variables
        :param i0: the initial state (vertex)
        """
        # it is better to use a `dict` for _adj and _pre_adj to save space
        # self._adj= [{} for _ in range(2 ** n + 1)]  # to get the edge information quickly
        self._adj = collections.defaultdict(dict)
        # self._pre_adj = [[] for _ in range(2 ** n + 1)]
        self._pre_adj = collections.defaultdict(list)
        self._vertices = set()
        self._n = n
        self.i0 = i0

    def add_edge(self, i: int, j: int, weight: float, control: int, switch: int):
        """
        Add an edge for state transition i --> j.

        :param i: a vertex (state)
        :param j: a vertex (state) that *i* can reach in one step
        :param weight: the minimum weight of the transition from *i* to *j*
        :param control: the control to attain the *weight*
        :param switch: the switching signal to attain the *weight*
        """
        self._vertices.add(i)
        self._vertices.add(j)
        v = Vertex(index=j, weight=weight, control=control, switch=switch)
        self._adj[i][j] = v
        self._pre_adj[j].append(i)

    @property
    def adj(self) -> dict:
        """
        The adjacency list representation of this graph. If a vertex (state) i has successors j and k, then
        self.adj[i] is {j: Vertex_j, k: Vertex_k}.

        :return: the adjacency-list representation
        """
        return self._adj

    @property
    def pre_adj(self) -> dict:
        """
        The predecessor list, i.e., a reverse adjacency list. For each edge i --> j, pre_adj[j] contains i.

        :return: pre_adj[j] gives the predecessors (a list) of vertex (state) j
        """
        return self._pre_adj

    @property
    def vertices(self) -> set:
        return self._vertices

    @property
    def n_state_variables(self) -> int:
        return self._n

    @classmethod
    def from_SBCN(cls, sbcn: SBCN, i0: int):
        """
        Build the OSTG for an SBCN

        :param sbcn: a network
        :param i0: the initial state
        :return: the OSTG
        """
        graph = cls(sbcn.n, i0)
        marked = set()
        q = [i0]
        marked.add(i0)

        while q:
            i = q.pop()
            successors = sbcn.compute_successors(i)
            for j, w, k, l in successors:
                graph.add_edge(i, j, w, k, l)
                if j not in marked:
                    marked.add(j)
                    q.append(j)
        return graph


class MMCAlgorithm:
    """
    The MMC based algorithm for control law design

    Refer to Algorithm 2 of the paper.
    """
    def __init__(self, graph: OSTG):
        self.graph = graph

    def solve(self, verbose=False):
        """
        Solve the MMC problem, and get the state-feedback matrices

        :param verbose: if `True`, then return more information about the minimizing vertex, the optimal cycle etc.
        :return: If verbose is `False`, then return a tuple (Ku, Ks, mcm) for the control and switching laws.
            If verbose is `True`, then return a tuple (Ku, Ks, mcm, info), where info is a namedtuple `Info`
                (v_star: int, mmc: list, pt: list, p_star: list, alpha: int, beta: int).
            The state-feedback gain matrices Ku and Ks are given as a dict to save memory. For example, the key-value
            pair <k, z> in Ku means the k-th column of Ku is \delta_M^z. (k starts from 1)
        """
        nv = len(self.graph.vertices)
        N = 2 ** self.graph.n_state_variables
        D = {}  # (int, int) --> number
        B = {}
        D[(0, self.graph.i0)] = 0

        # Task (i)
        for k in range(1, nv + 1):
            for j in self.graph.vertices:
                candidates = ((i, D.get((k - 1, i), _INF) + self.graph.adj[i][j].weight) for i in self.graph.pre_adj[j])
                i_star, d_star = min(candidates, key=lambda ele: ele[1], default=(None, _INF)) # candidates may be empty
                if d_star != _INF:  # j can be reached by a k-edge path
                    D[(k, j)] = d_star
                    B[(k, j)] = i_star
        # get the mcm and its minimizer
        i_star = None
        k_star = None
        mcm = _INF
        for i in self.graph.vertices:
            candidates = [((D.get((nv, i), _INF) - D.get((k, i), _INF)) / (nv - k), k) for k in range(0, nv)]
            m = max(candidates, key=operator.itemgetter(0)) # inf - inf -> nan, but nan will be ignored in `max`
            if m[0] < mcm:
                mcm = m[0]
                i_star = i
                k_star = m[1]

        # Task (ii)
        p = [None] * (nv + 1)
        p[nv] = i_star
        for k in range(nv, 0, -1):
            p[k - 1] = B[(k, p[k])]

        # Task (iii)
        A = [-1] * (N + 1)
        for t in range(0, nv + 1):
            it = p[t]
            if A[it] == -1:
                A[it] = t
            else:
                alpha = A[it]
                beta = t - 1
                break
        # the minimum-mean cycle
        mmc = [p[t] for t in range(alpha, beta + 1)]
        mmc.append(p[alpha])
        # the sub-path
        transient_path = p[:alpha]

        # extract the feedback gain matrix
        Ku = {}
        Ks = {}
        for t in range(0, beta + 1):
            i = p[t]
            if t < beta:
                j = p[t + 1]
            else:
                j = p[alpha]
            uij = self.graph.adj[i][j].control
            sij = self.graph.adj[i][j].switch
            Ku[i] = uij
            Ks[i] = sij
        if verbose:
            info = Info(v_star=i_star, k_star=k_star, mmc=mmc, pt=transient_path, p_star=p, alpha=alpha, beta=beta)
            return Ku, Ks, mcm, info
        return Ku, Ks, mcm