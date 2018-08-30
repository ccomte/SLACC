# Seize the Longest Available Compatible Chip (SLACC)

Evaluate the performance of the token-based load balancing algorithm proposed in the paper *Dynamic Load Balancing with Tokens* available on [HAL](https://hal.archives-ouvertes.fr/hal-01758912).

Author: Céline Comte

## Files

The following files predict the performance of the dynamic load balancing:
- `dynamic-bipartite-exact.c` and `dynamic-tripartite-exact.c`: Apply the exact formulas proposed in the paper to predict performance metrics like the blocking probability of the jobs of each type or the mean number of jobs within each class. `exact-tripartite.c` considers the general scenario while `exact-bipartite.c` assumes that there is a one-by-one association between classes and servers.
- `dynamic-bipartite-onetype-exact.ipynb`: Performs a state aggregation to simplify the computations when there is a single job type, assuming that there is a one-by-one association between classes and servers.
- `dynamic-bipartite-exp.c` and `dynamic-bipartite-hyp.c`: Compute the performance metrics by simulation, with job sizes that are either exponentially distributed or hyperexponentially distributed, assuming that there is a one-by-one association between classes and servers.

The performance metrics obtained under the static load balancing policies are computed in the Jupyter Notebook `static-bipartite-exact.ipynb`.

We also use the following files:
- `utils.h`: Libraries and random sampling functions used in the the C programs.
- `script.sh`: An example script that builds and runs the above-mentioned C programs. The two scenarios are those described in the paper. The results are saved in the folder *data*.
- `plot.ipynb`: Plot the obtained results.

Each C program can be called with the option `-h` to get more information on its utilization. All the data files are stored in the folder *data*.

## Licence

This project is under [GPL v3.0 license](https://github.com/ccomte/SLACC/blob/master/LICENSE).
