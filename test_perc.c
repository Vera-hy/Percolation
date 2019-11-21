#include <stdio.h>
#include <stdlib.h>

#include "percolate.h"

void test_perc(int L, int** map){

    int itop, ibot, perc;
    perc = 0;

    for (itop=0; itop < L; itop++)
    {
        if (map[itop][L-1] > 0)
        {
            for (ibot=0; ibot < L; ibot++)
            {
                if (map[itop][L-1] == map[ibot][0])
                {
                    perc = 1;
                }
            }
        }
    }

    if (perc != 0)
    {
        printf("percolate: cluster DOES percolate\n");
    }
    else
    {
        printf("percolate: cluster DOES NOT percolate\n");
    }
}

