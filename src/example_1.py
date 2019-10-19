"""
Example 1 in the paper: a toy example
"""
from algorithm.mmc import OSTG, MMCAlgorithm
from algorithm.sbcn import SBCN


Qx = [5, 3, 4, 0, 1, 3, 0, 1]
Qu = [3, 1]
Qs = [1, 2]


def Cs(i):
    """
    Switching constraints

    :return: a tuple
    """
    if i in [1, 2, 5]:
        return 1,
    return 1, 2


def g(i, k, l):
    """
    Stage cost

    :param i: state
    :param k: control
    :param l: switching
    :return: the cost
    """
    return Qx[i - 1] + Qu[k - 1] + Qs[l - 1]


if __name__ == '__main__':
    # build the OSTG
    L1 = [7, 6, 8, 6, 3, 5, 7, 1, 3, 5, 7, 1, 7, 6, 8, 6]
    L2 = [4, 2, 4, 2, 3, 5, 3, 2, 3, 1, 3, 2, 4, 6, 4, 2]
    Cx = {1, 2, 3, 5, 6, 7, 8}
    sbcn = SBCN(3, 1, 2, [L1, L2], g, Cx, None, Cs)
    ostg = OSTG.from_SBCN(sbcn, 1)
    # if we want to inspect the graph, print its adjacency-list representation
    print_graph = False
    if print_graph:
        print('................adjacency list of OSTG..............')
        for i in Cx:
            successors = ostg.adj[i]
            if successors:
                print(f'- successors of {i}')
                for v in successors.values():
                    print(v)

        print('................predecessor list..............')
        for i in Cx:
            predecessors = ostg.pre_adj[i]
            if predecessors:
                print(f'----{i}----')
                for v in predecessors:
                    print(v)
    # MMC
    alg = MMCAlgorithm(ostg)
    Ku, Ks, mcm, info = alg.solve(verbose=True)
    print('mu* (minimum cycle mean): ', mcm)
    print('v*: ', info.v_star)
    print('k*: ', info.k_star)
    print('MMC: ', info.mmc)
    print('pt (transient path): ', info.pt)
    print('p* (a path with |V| edges to v*): ', info.p_star)
    print('alpha: ', info.alpha)
    print('beta: ', info.beta)
    print(f'Ku = {Ku}\nKs = {Ks}')