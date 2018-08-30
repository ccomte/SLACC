#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>     /* unsigned integer type uint64_t */
#include <string.h>     /* manipulate intput strings */
#include <unistd.h>     /* getopt */
#include <ctype.h>      /* isprint */

#define ONE 1UL
#define INT_LENGTH 64   /* length of the unsigned integers used to encode the state */
#define MAX_K 10        /* maximum number of job types */
#define MAX_N 10        /* maximum number of token classes */
#define MAX_S 10        /* maximum number of servers */






/*** VARIABLES ***/



/* Inputs */

uint64_t K;       /* number of job types */
uint64_t N;       /* number of token classes */
uint64_t S;       /* number of servers */

uint64_t elli[MAX_N];   /* number of tokens of each class */
uint64_t Ki[MAX_N][MAX_K];          /* compatiblity graph between classes and types,
                                       encoded as an adjacency matrix */
uint64_t notKi[MAX_K][MAX_N + 1];   /* list of classes each type
                                       is **not** compatible with */
uint64_t Si[MAX_N][MAX_S];          /* compatiblity graph between classes and servers,
                                       encoded as an adjacency matrix */
uint64_t notSi[MAX_S][MAX_N + 1];   /* list of classes each server
                                       is **not** compatible with */

double mus[MAX_S];      /* normalized service capacity of each server */
double *muI;            /* normalized service rate function */
double mu;              /* total capacity rate */

double nuk[MAX_K];      /* normalized arrival rate of each job type */
double *nuI;            /* normalized arrival rate function */



/* Variables for the computations */

double Rho;       /* largest load to consider */
double delta;     /* step by which we increase the load */

uint64_t length;  /* number of bits to encode the number of available tokens of each class */
uint64_t mask;    /* mask to detect the number of available tokens of a given class */
uint64_t xymax;   /* maximum aggregate state to consider */
uint64_t Imax;    /* maximum active set of active classes to consider */
uint64_t *ei;     /* N-dimensional unit vectors */

double *Phix;     /* balance function of the departure rates */
double *Phiy;     /* balance function of the departure rates,
                     as a function of y */
double *Lambday;  /* balance function of the arrival rates */
double *piy;      /* stationary measure with piy[0] = 1 */
double pi;        /* normalizating constant of the stationary distribution */



/* Outputs */

char name[40];    /* name of the output file */





/*** UTILITARY ***/

/* print an unsigned integer in binary */
void printb (uint64_t u) {
  uint64_t i;

  for (i = 0 ; i < INT_LENGTH ; ++i) {
    if ( !((INT_LENGTH - i) % length) ) printf(" ");
    printf("%lu", u >> (INT_LENGTH - 1 - i) & ONE);
  }
}

/* verifies whether the bitmap x is valid */
static inline int is_valid (uint64_t x) {
  uint64_t i;

  for (i = 0 ; i < N ; ++i) {
    if ( ((x >> (length * i)) & mask) > elli[i] )
      return 0;
  }

  return 1;
}



/*** COMPUTATIONS ***/

/* basic operations
 * - access the number of available tokens from server i:
 *   (x >> (length * i)) & mask
 * - check whether there is at least one available token from server i:
 *   ((x >> (length * i)) & mask) > 0
 * - remove one available token of server i:
 *   x -= (ONE << (length * i))
 */



void compute_Phix () {
  uint64_t i, x, I;

  /* initialization */
  Phix[0] = 1.;

  /* recursion step */
  for (x = ONE ; x < xymax ; ++x) {
    Phix[x] = 0.;
    if (is_valid(x)) {
      I = 0;
      for (i = 0 ; i < N ; ++i) {
        if ( (x >> (length * i)) & mask ) {
          /* there is at least one job of class i in service */
          Phix[x] += Phix[x - ei[i]];
          I += ONE << i;
        }
      }
      Phix[x] /= muI[I];
    }
  }
}



void compute_Phiy () {
  uint64_t i, x, y;

  for (y = 0 ; y < xymax ; ++y) {
    Phiy[y] = 0.;

    if (is_valid(y)) {
      /* compute x = ell - y */
      x = 0;
      for (i = 0 ; i < N ; ++i)
        x += ( ( elli[i] - ((y >> (length * i)) & mask) ) << (length * i) );

      /* update the array */
      Phiy[y] = Phix[x];
    }
  }
}



void compute_piy (double rho) {
  uint64_t i, y, I;

  /* initialization */
  Lambday[0] = 1.;
  piy[0] = Phiy[0] * Lambday[0];
  pi = piy[0];

  /* recursion step */
  for (y = ONE ; y < xymax ; ++y) {
    Lambday[y] = 0.;
    piy[y] = 0.;

    if (is_valid(y)) {
      /* update Lambday */
      I = 0;
      for (i = 0 ; i < N ; ++i) {
        if ( (y >> (length * i)) & mask ) {
          /* there is at least one available token from class i */
          Lambday[y] += Lambday[y - ei[i]];
          I += ONE << i;
        }
      }
      Lambday[y] /= (rho * nuI[I]);

      /* update piy and pi */
      piy[y] = Phiy[y] * Lambday[y];
      pi += piy[y];
    }
  }
}



double blocking (int k) {
  uint64_t i, t, y, carry;
  double betak;

  betak = 0.;
  y = 0;
  carry = 0;

  while (!carry) {
    /* update the blocking probability */
    betak += piy[y];

    /* update y */
    carry = 1;
    t = 0;
    while ( carry && ((i = notKi[k][t]) != -1) ) {
      if ( ((y >> (i * length)) & mask) < elli[i] ) {
        y += ONE << (i * length);
        carry = 0;
      } else y -= elli[i] << (i * length);
      ++t;
    }
  }

  return betak / pi;
}



double mean (int i) {
  uint64_t y;
  double Li;

  Li = 0.;
  for (y = ONE ; y < xymax ; ++y) {
    if (is_valid(y)) Li += ((y >> (length * i)) & mask) * piy[y];
  }

  return elli[i] - (Li / pi);
}





/*** INPUTS AND OUTPUTS ***/

void read_inputs (int argc, char **argv) {
  int c;
  uint64_t i, k, s, t, current, total;
  char *arg;
  double norm;

  /* default initialization */
  sprintf(name, "toto");
  K = 0; N = 0;
  for (i = 0 ; i < MAX_N ; ++i) {
    elli[i] = 0;
    for (k = 0 ; k < MAX_K ; ++k) Ki[i][k] = 0;
    for (s = 0 ; s < MAX_S ; ++s) Si[i][s] = 0;
  }
  for (s = 0 ; s < MAX_S ; ++s) {
    mus[s] = 0.;
    notSi[s][0] = -1;
  }
  for (k = 0 ; k < MAX_K ; ++k) {
    nuk[k] = 0.;
    notKi[k][0] = -1;
  }
  Rho = 5.;
  delta = .01;

  /* read arguments */
  opterr = 0;
  while ( (c = getopt(argc, argv, "f:k:n:s:l:c:r:m:d:h")) != -1 ) {
    switch(c) {
      case 'f':
        /* output file name */
        strncpy(name, optarg, 40);
        break;
      case 'k':
        /* number of job types */
        K = atoi(optarg);
        if (K > MAX_K) {
          fprintf(stderr, "Option -k: the number of job types"
              "cannot be larger than %d.\n", MAX_K);
          exit(EXIT_FAILURE);
        }
        break;
      case 'n':
        /* number of token classes */
        N = atoi(optarg);
        if (N > MAX_N) {
          fprintf(stderr, "Option -n: the number of token classes"
              "cannot be larger than %d.\n", MAX_N);
          exit(EXIT_FAILURE);
        }
        break;
      case 's':
        /* number of servers */
        S = atoi(optarg);
        if (S > MAX_S) {
          fprintf(stderr, "Option -s: the number of servers"
              "cannot be larger than %d.\n", MAX_S);
          exit(EXIT_FAILURE);
        }
        break;
      case 'l':
        /* number of tokens of each server */
        total = 0;
        arg = strtok(optarg, " ");
        i = 0;
        while (i < N && arg != NULL) {
          elli[i] = atoi(arg);
          total += elli[i];
          ++i;
          arg = strtok(NULL, " ");
        }

        if (i == 0) {
          fprintf(stderr, "Option -l: you should specify the number of tokens.\n");
          exit(EXIT_FAILURE);
        }

        while (i < N) {
          elli[i] = elli[i-1];
          total += elli[i];
          ++i;
        }

        /* number of bits to encode the number of tokens of each server */
        length = 0;
        for (i = 0 ; i < N ; ++i) {
          current = INT_LENGTH;
          while ( (current > length) && !(elli[i] >> (current - 1) & ONE) ) --current;
          if (current > length) length = current;
        }

        if (N * length >= INT_LENGTH) {
          fprintf(stderr, "Option -l: the unsigned integer type does not have enough bits to encode all states.\n");
          exit(EXIT_FAILURE);
        }

        /* other lengths */
        mask = (ONE << length) - ONE;
        xymax = ONE << (N * length);
        Imax = ONE << N;

        break;
      case 'c':
        arg = strtok(optarg, " ");

        /* class-to-type compatibilities,
         * encoded as an adjacency matrix between classes and types */

        for (i = 0 ; i < N ; ++i) {
          for (k = 0 ; k < K ; ++k) {
            if (arg == NULL) {
              fprintf(stderr, "Option -c: you didn't specify all the compatibilities.\n");
              exit(EXIT_FAILURE);
            }

            Ki[i][k] = atoi(arg);
            arg = strtok(NULL, " ");
          }
        }

        /* class-to-server compatibilities,
         * encoded as an adjacency matrix between classes and servers */

        for (i = 0 ; i < N ; ++i) {
          for (s = 0 ; s < S ; ++s) {
            if (arg == NULL) {
              fprintf(stderr, "Option -c: you didn't specify all the compatibilities.\n");
              exit(EXIT_FAILURE);
            }

            Si[i][s] = atoi(arg);
            arg = strtok(NULL, " ");
          }
        }
        break;
      case 'r':
        arg = strtok(optarg, " ");

        /* arrival rates */
        norm = 0.;
        for (k = 0 ; k < K ; ++k) {
          if (arg == NULL) {
            fprintf(stderr, "Option -r: you didn't specify all the arrival rates.\n");
            exit(EXIT_FAILURE);
          }

          nuk[k] = atof(arg);
          norm += nuk[k];
          arg = strtok(NULL, " ");
        }

        /* normalize the arrival rates */
        for (k = 0 ; k < K ; ++k) nuk[k] /= norm;

        /* service rates */
        mu = 0.;
        for (s = 0 ; s < S ; ++s) {
          if (arg == NULL) {
            fprintf(stderr, "Option -r: you didn't specify all the service rates.\n");
            exit(EXIT_FAILURE);
          }

          mus[s] = atof(arg);
          mu += mus[s];
          arg = strtok(NULL, " ");
        }

        /* normalize the service rates */
        for (s = 0 ; s < S ; ++s) mus[s] /= mu;

        break;
      case 'm':
        /* maximum load that is considered */
        Rho = atof(optarg);
        break;
      case 'd':
        /* step by which we increase the load */
        delta = atof(optarg);
        break;
      case 'h':
        /* help */
        printf("\nOptions:\n"
            "-f: output file name\n"
            "-k: number of job types\n"
            "-n: number of servers\n"
            "-l: number of tokens of each class\n"
            "-c: class-to-type compatibilities, encoded as an adjacency matrix,"
            " and class-to-server compatibilities, again encoded as an adjacency matrix\n"
            "-r: per-type arrival rates and per-server capacities\n"
            "-m: maximum load to consider\n"
            "    (default value: 5)\n"
            "-d: step by which the load is increased\n"
            "    (default value: 0.01)\n"
            "\n");
        exit(0);
        break;
      case '?':
        /* error */
        if (optopt == 'c')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        exit(EXIT_FAILURE);
      default:
        abort();
    }
}

  for (i = optind; i < argc; ++i)
    printf ("Non-option argument %s\n", argv[i]);
}

void print_inputs () {
  uint64_t i, k, s;

  printf("\n--- DYNAMIC - TRIPARTITE - EXACT ---\n\n");

  /* output file prefix name */
  printf("Output file prefix name: %s.csv\n\n", name);

  /* computation parameters */
  printf("Maximum load to be considered:      %.3e\n", Rho);
  printf("Step by which we increase the load: %.3e\n\n", delta);

  /* numbers of types and servers */
  printf("%lu job type(s), %lu class(es), and %lu server(s)\n\n", K, N, S);

  /* numbers of tokens of each class */
  for (i = 0 ; i < N ; ++i) printf("Class %lu: %lu token(s)\n", i, elli[i]);
  printf("\n");

  /* class-to-type compatibilities */
  printf("Adjacency matrix of the compatiblity graph between classes and types:\n");
  for (i = 0 ; i < N ; ++i) {
    for (k = 0 ; k < K ; ++k) printf("%lu ", Ki[i][k]);
    printf("\n");
  }
  printf("\n");

  /* class-to-server compatibilities */
  printf("Adjacency matrix of the compatiblity graph between classes and servers:\n");
  for (i = 0 ; i < N ; ++i) {
    for (s = 0 ; s < S ; ++s) printf("%lu ", Si[i][s]);
    printf("\n");
  }
  printf("\n");

  /* arrival rates */
  printf("Normalized arrival rates:\n");
  for (k = 0 ; k < K ; ++k) printf("%.3e  ", nuk[k]);
  printf("\n\n");

  /* service rates */
  printf("Normalized service rates:\n");
  for (s = 0 ; s < S ; ++s) printf("%.3e  ", mus[s]);
  printf("\nTotal service rate: %.3e\n\n", mu);
}

FILE *open_csv_file () {
  char file_name[44];
  FILE *file = NULL;

  sprintf(file_name, "%s.csv", name);
  file = fopen(file_name, "w");

  if (file == NULL) {
    fprintf(stderr, "Couldn't open the output file %s.\n\n", file_name);
    exit(EXIT_FAILURE);
  } else printf("Save data in the file '%s'.\n\n", file_name);

  return file;
}

void initialize_arrays () {
  uint64_t i, k, s, t, I;

  /* stationary measure */
  piy = malloc(xymax * sizeof(double));
  if (piy == NULL) {
    fprintf(stderr, "Couldn't allocate the variable piy.\n\n");
    exit(EXIT_FAILURE);
  }

  /* balance functions */
  Phix = malloc(xymax * sizeof(double));
  if (Phix == NULL) {
    fprintf(stderr, "Couldn't allocate the variable Phix.\n\n");
    exit(EXIT_FAILURE);
  }

  Phiy = malloc(xymax * sizeof(double));
  if (Phiy == NULL) {
    fprintf(stderr, "Couldn't allocate the variable Phiy.\n\n");
    exit(EXIT_FAILURE);
  }

  Lambday = malloc(xymax * sizeof(double));
  if (Lambday == NULL) {
    fprintf(stderr, "Couldn't allocate the variable Lambday.\n\n");
    exit(EXIT_FAILURE);
  }

  /* unit vectors */
  ei = malloc(N * sizeof(uint64_t));
  if (ei == NULL) {
    fprintf(stderr, "Couldn't allocate the variable ei.\n\n");
    exit(EXIT_FAILURE);
  }

  for (i = 0 ; i < N ; ++i) ei[i] = ONE << (length * i);

  /* arrival rate function */
  nuI = malloc(Imax * sizeof(double));
  if (nuI == NULL) {
    fprintf(stderr, "Couldn't allocate the variable nuI.\n\n");
    exit(EXIT_FAILURE);
  }

  nuI[0] = 0.;
  for (I = ONE ; I < Imax ; ++I) {
    nuI[I] = 0.;
    for (k = 0 ; k < K ; ++k) {
      for (i = 0 ; i < N ; ++i) {
        if ( ((I >> i) & ONE) && Ki[i][k] ) {
          nuI[I] += nuk[k];
          break;
        }
      }
    }
  }

  /* service rate function */
  muI = malloc(Imax * sizeof(double));
  if (muI == NULL) {
    fprintf(stderr, "Couldn't allocate the variable muI.\n\n");
    exit(EXIT_FAILURE);
  }

  muI[0] = 0.;
  for (I = ONE ; I < Imax ; ++I) {
    muI[I] = 0.;
    for (s = 0 ; s < S ; ++s) {
      for (i = 0 ; i < N ; ++i) {
        if ( ((I >> i) & ONE) && Si[i][s] ) {
          muI[I] += mus[s];
          break;
        }
      }
    }
  }

  /* list of adjacency lists for the compatibilities between classes and types */
  for (k = 0 ; k < K ; ++k) {
    t = 0;
    for (i = 0 ; i < N ; ++i) {
      if (!Ki[i][k]) {
        notKi[k][t] = i;
        ++t;
      }
    }
    notKi[k][t] = -1;
  }

  /* list of adjacency lists for the compatibilities between classes and servers */
  for (s = 0 ; s < S ; ++s) {
    t = 0;
    for (i = 0 ; i < N ; ++i) {
      if (!Si[i][s]) {
        notSi[s][t] = i;
        ++t;
      }
    }
    notSi[s][t] = -1;
  }
}





/*** MAIN ***/

int main (int argc, char **argv) {
  uint64_t i, k, s;
  double rho, betak, beta, psis, eta, Li, L;
  FILE *file;

  /* input arguments */
  read_inputs(argc, argv);
  print_inputs();
  initialize_arrays();

  /* open the output file */
  file = open_csv_file();

  /* print the column names */
  fprintf(file, "rho");
  for (k = 0 ; k < K ; ++k) fprintf(file, ",betak%lu", k+1);
  for (i = 0 ; i < N ; ++i) fprintf(file, ",Li%lu,gammai%lu", i+1, i+1);
  fprintf(file, ",beta,eta,L,gamma\n");

  /* pre-compute Phi */
  compute_Phix();
  compute_Phiy();

  for (rho = delta ; rho < Rho + .5 * delta ; rho += delta) {
    /* computations */
    printf("rho = %.3e\n", rho);
    compute_piy(rho);

    /* print data in the output file */
    fprintf(file, "%e", rho);

    /* per-type blocking probability */
    beta = 0.;
    for (k = 0 ; k < K ; ++k) {
      betak = blocking(k);
      fprintf(file, ",%e", betak);
      beta += nuk[k] * betak;
    }

    /* per-class mean number of jobs and mean service rate */
    L = 0.;
    for (i = 0 ; i < N ; ++i) {
      Li = mean(i);
      L += Li;

      if (Li > 0.) fprintf(file, ",%e,%e", Li, 0.);
      else fprintf(file, ",%e,%e", 0., 0.);
    }

    /* system-wide metrics */
    fprintf(file, ",%e,%e,%e,%e\n",
        beta, rho * (1. - beta), L, rho * mu * (1. - beta) / L);
  }

  /* close the file and release memory */
  fclose(file);
  free(piy); free(Phix); free(Phiy); free(Lambday); free(ei); free(nuI); free(muI);
}
