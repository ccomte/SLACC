#!/bin/bash



# compilation

gcc -O3 -o dynamic-bipartite-exact.out dynamic-bipartite-exact.c -lm
gcc -O3 -o dynamic-tripartite-exact.out dynamic-tripartite-exact.c -lm
gcc -O3 -o dynamic-bipartite-exp.out dynamic-bipartite-exp.c -lm
gcc -O3 -o dynamic-bipartite-hyp.out dynamic-bipartite-hyp.c -lm



# folder
mkdir -p data



# a single job type

./dynamic-bipartite-exact.out -f "data/single-dynamic-exact" \
  -k 1 -n 10 -l 6 \
  -c "1 1 1 1 1 1 1 1 1 1" \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .01

./dynamic-tripartite-exact.out -f "data/single-dynamic-exact-3" \
  -k 1 -n 6 -s 10 -l 10 \
  -c "1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 1" \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .005

./dynamic-tripartite-exact.out -f "data/single-dynamic-exact-4" \
  -k 1 -n 4 -s 10 -l 15 \
  -c "1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 1 1 " \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .005

./dynamic-tripartite-exact.out -f "data/single-dynamic-exact-5" \
  -k 1 -n 2 -s 10 -l 30 \
  -c "1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1" \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .005

./dynamic-bipartite-exp.out -f "data/single-dynamic-simu-exp" \
  -k 1 -n 10 -l 6 \
  -c "1 1 1 1 1 1 1 1 1 1" \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -m 3 -d .2 -R 100 -w 6 -t 6

./dynamic-bipartite-hyp.out -f "data/single-dynamic-simu-hyperexp" \
  -k 1 -n 10 -l 6 \
  -c "1 1 1 1 1 1 1 1 1 1" \
  -r "1 1 1 1 1 1 4 4 4 4 4" \
  -s "1 5" \
  -m 3 -d .2 -R 100 -w 6 -t 6



# several job types

./dynamic-bipartite-exact.out -f "data/multi-dynamic-exact" \
  -k 2 -n 10 -l 6 \
  -c "1 0 1 0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 0 1" -r "1 4 1 1 1 1 1 1 1 1 1 1" \
  -m 3 -d .01

./dynamic-bipartite-exp.out -f "data/multi-dynamic-simu-exp" \
  -k 2 -n 10 -l 6 \
  -c "1 0 1 0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 0 1" -r "1 4 1 1 1 1 1 1 1 1 1 1" \
  -m 3 -d .2 -R 100 -w 6 -t 6

./dynamic-bipartite-hyp.out -f "data/multi-dynamic-simu-hyperexp" \
  -k 2 -n 10 -l 6 \
  -c "1 0 1 0 1 0 1 1 1 1 1 1 1 1 0 1 0 1 0 1" -r "1 4 1 1 1 1 1 1 1 1 1 1" \
  -s "1 2 1 5" \
  -m 3 -d .2 -R 100 -w 6 -t 6
