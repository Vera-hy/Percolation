#include <stdio.h>
#include <stdlib.h>

#include "percolate.h"

/*
 *  Test to see if percolation occurred by looking for positive numbers
 *  that appear on both the top and bottom edges.
 */
void test_perc(int l, int** map){

    int itop, ibot, perc;
    perc = 0;

    for (itop=0; itop < l; itop++)
    {
        if (map[itop][l-1] > 0)
        {
            for (ibot=0; ibot < l; ibot++)
            {
                if (map[itop][l-1] == map[ibot][0])
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

