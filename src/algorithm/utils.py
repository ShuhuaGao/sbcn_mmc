"""
Some utility functions
"""
from .mmc import OSTG

def read_network(file):
    """
    Read a (switched) Boolean network from a text file:
        Line 1: number of state variables
        Line 2: number of control inputs
        Line 3: number of sub-networks
        Line 4: transition matrix of the first sub-network (linear representation of a logical matrix)
        line 5: transition matrix of the second sub-network (linear representation of a logical matrix)
        ...

    :param file: a text file
    :return: (n, m, w, Ls), where
        n: number of state variables
        m: number of control inputs
        w: number of sub-systems
        Ls: a list of transition matrices, each for one sub-system
    """
    with open(file, 'r') as f:
        n = int(f.readline().strip())
        m = int(f.readline().strip())
        w = int(f.readline().strip())
        N = 2 ** n
        M = 2 ** m
        Ls = []
        for _ in range(w):
            line = f.readline().strip()
            assert line, f'{w} transition matrices must be provided!'
            numbers = line.split()
            assert len(numbers) == M * N, f'The transition matrix must have {M * N} columns'
            Ls.append([int(num) for num in numbers])
        return n, m, w, Ls



def export_OSTG(file: str, g: OSTG):
    """
    Export an OSTG to a .sif file, which can be imported into software like Cytoscape for visualization.

    :param file: specify a text file to be written into
    :param g: a graph
    """
    with open(file, 'w') as f:
        for i in g.vertices:
            for j in g.adj[i]:
                f.write(f'{i} -> {j}\n')
