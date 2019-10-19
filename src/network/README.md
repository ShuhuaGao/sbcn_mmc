You can specify a network in a text file in the following format.

```bash
Line 1: number of state variables
Line 2: number of control inputs
Line 3: number of sub-networks
Line 4: transition matrix of the first sub-network (linear representation of a logical matrix)
Line 5: transition matrix of the second sub-network (linear representation of a logical matrix)
......

```

For example, the linear representation of a logical matrix $[\delta_4^1, \delta_4^3, \delta_4^2]$ is `1 3 2` in the text file.