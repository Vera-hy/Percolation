/*
 *  Main header file for percolation code.
 */

/*
 *  Although overall system is square, i.e. size L x L, we will define
 *  different variables for the first and second dimensions. This is
 *  because, in the parallel code, the local arrays will not be
 *  square. For example, using a simple 1D decomposition over P
 *  processes, then M = L/P and N = L
 */


#define TRUE  1
#define FALSE 0
#define ndims 2

/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void percwritedynamic(char *percfile, int **map, int l, int ncluster);
/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);

/*
 *  Dynamic allocation
 */
void *arralloc(size_t size, int ndim, ...);

