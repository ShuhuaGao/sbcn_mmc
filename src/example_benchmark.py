"""
The benchmark example: T-LGL network
"""
import time
from algorithm.sbcn import SBCN
from algorithm.mmc import OSTG, MMCAlgorithm
from algorithm.utils import read_network, export_OSTG


def g(i, k, l):
    """
    Stage cost

    :param i: state
    :param k: control
    :param l: switching
    :return: the cost
    """
    if i == 65535:
        return 1
    return 5


if __name__ == '__main__':
    # load the ASSR for the T-LGL network
    n, m, w, Ls = read_network('./network/T_LGL.txt')
    print(n, m, w, len(Ls))
    bcn = SBCN(n, m, w, Ls, g)  # one subsystem, no constraints


    # build the OSTG: Algorithm 1
    ts = time.time()
    initial_state = 58834
    ostg = OSTG.from_SBCN(bcn, initial_state)
    print('Time to build OSTG (s): ', time.time() - ts)
    print(f'Number of vertices in the OSTG from initial state {initial_state}: ', len(ostg.vertices))


    # compute the optimal control law: Algorithm 2
    ts = time.time()
    alg = MMCAlgorithm(ostg)
    Ku, Ks, mcm, info = alg.solve(verbose=True)
    print('Time to solve optimal law (s): ', time.time() - ts)

    print('mu* (minimum cycle mean): ', mcm)
    print(f'v* = {info.v_star}, k* = {info.k_star}')
    print('MMC: ', info.mmc)
    print('pt (transient path): ', info.pt)
    print('p* (a path with |V| edges to v*): ', info.p_star)
    print('alpha: ', info.alpha)
    print('beta: ', info.beta)
    print(f'Ku = {Ku}\nKs = {Ks}')

    # export the OSTG if needed
    # export_OSTG('./network/T-LGL-OSTG.sif', ostg)