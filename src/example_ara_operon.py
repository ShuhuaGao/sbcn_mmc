"""
The E.coli Ara operon example in [1]

[1] Y. Wu, X.-M. Sun, X. Zhao, and T. Shen, “Optimal control of boolean
control networks with average cost: A policy iteration approach,” Auto-
matica, vol. 100, pp. 378–387, 2019.

"""
import time
from algorithm.sbcn import SBCN
from algorithm.mmc import OSTG, MMCAlgorithm
from algorithm.utils import read_network


def inverse_map(i: int, n: int):
    """
    Accumulative STP of logical variables is bijective.
    Given a result i (\delta_{2^n}^i), find the corresponding logical values.

    :return a list of 0/1
    """
    r = []
    while n > 0:
        if i % 2 == 0:
            r.append(0)
            i = i // 2
        else:
            r.append(1)
            i = (i + 1) // 2
        n = n - 1
    r.reverse()
    return r


def g(i, k, l):
    """
    Stage cost. Please check page 384 and page 385 of [1].

    :param i: state
    :param k: control
    :param l: switching
    :return: the cost
    """
    n = 9
    m = 4
    X = inverse_map(i, n)
    U = inverse_map(k, m)
    A = [-28, -12, 12, 16, 0, 0, 0, 20, 16]
    B = [-8, 40, 20, 40]
    return sum(a * x for a, x in zip(A, X)) + sum(b * u for b, u in zip(B, U))


if __name__ == '__main__':
    # load the ASSR for the E.coli network
    n, m, w, Ls = read_network('./network/ara_operon.txt')
    bcn = SBCN(n, m, w, Ls, g)  # one subsystem, no constraints

    # choose an initial state here
    initial_state = 121

    # build the OSTG: Algorithm 1
    ts = time.time()
    ostg = OSTG.from_SBCN(bcn, initial_state)
    print('Time to build OSTG (s) (Algorithm 1): ', time.time() - ts)
    print(f'Number of vertices in the OSTG from initial state {initial_state}: ', len(ostg.vertices))

    # compute the optimal control law: Algorithm 2
    ts = time.time()
    alg = MMCAlgorithm(ostg)
    Ku, Ks, mcm, info = alg.solve(verbose=True)
    print('Time to solve optimal law (s) (Algorithm 2): ', time.time() - ts)

    print('mu* (minimum cycle mean): ', mcm)
    print('v*: ', info.v_star)
    print('MMC: ', info.mmc)
    print('pt (transient path): ', info.pt)
    print('p* (a path with |V| edges to v*): ', info.p_star)
    print('alpha: ', info.alpha)
    print('beta: ', info.beta)
    print(f'Ku = {Ku}\nKs = {Ks}')
