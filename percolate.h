/*
 *  Main header file for percolation code.
 */

/*
 *  Define different variables.
 */
#define TRUE  1
#define FALSE 0
#define ndims 2

/*
 *  Prototypes for functions shared between parallel version and sequential
 *  version.
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

/*
 *  Test to see if percolation occurred by looking for positive numbers
 *  that appear on both the top and bottom edges.
 */
void test_perc(int l, int** map);

