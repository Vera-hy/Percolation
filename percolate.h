/*
 *  Main header file for percolation code.
 */

/*
 *  System size L
 */

//#define L 288

//#define P 4
/*
 *  Although overall system is square, i.e. size L x L, we will define
 *  different variables for the first and second dimensions. This is
 *  because, in the parallel code, the local arrays will not be
 *  square. For example, using a simple 1D decomposition over P
 *  processes, then M = L/P and N = L
 */

//#define M L/2
//#define N L/2

//int M,N;

#define TRUE  1
#define FALSE 0
#define ndims 2

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

//void percwrite(char *percfile, int map[L][L], int ncluster);
void percwritedynamic(char *percfile, int **map, int l, int ncluster);
/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);
void *arralloc(size_t size, int ndim, ...);

