#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ldpc_create_pchk.h"

mod2sparse* CreatePchkMatrix (  int nbRows, int nbCols, make_method makeMethod, int leftDegree, int seed, bool no4cycle, SessionType type)
{
	mod2entry *e;
	int added, uneven;
	int i, j, k, t, l;
	int *u;
	mod2sparse *pchkMatrix = NULL;
	int skipCols = 0;		// avoid warning
	int nbDataCols = 0;		// avoid warning

	if (type != TypeLDGM && type != TypeSTAIRS && type != TypeTRIANGLE) {
		return NULL;
	}
	skipCols = nbRows;
	nbDataCols = nbCols-skipCols;

	// Check for some problems.
	if (leftDegree>nbRows)
	{
		return NULL;
	}
	if (no4cycle) { 
		exit(-1);
	}

	ldpc_srand(seed);
	pchkMatrix = mod2sparse_allocate(nbRows, nbCols);

	/* Create the initial version of the parity check matrix. */
	switch (makeMethod)
	{ 
		case Evencol:
			for(j=skipCols; j<nbCols; j++)
			{
				for(k=0; k<leftDegree; k++)
				{
					do
					{
						i = ldpc_rand(nbRows);
					}
					while (mod2sparse_find(pchkMatrix,i,j));
					mod2sparse_insert(pchkMatrix,i,j);
				}
			}
			break;

		case Evenboth:
			u = (int*)calloc(leftDegree*nbDataCols, sizeof(*u));

			for(k = leftDegree*nbDataCols-1; k>=0; k--)
			{
				u[k] = k%nbRows;
			}
			uneven = 0;
			t = 0;	/* left limit within the list of possible choices, u[] */
			for(j = skipCols; j<nbCols; j++)	/* for each source symbol column */
			{
				for(k = 0; k<leftDegree; k++)	/* add left_degree "1s" */
				{ 
					/* check that valid available choices remain */
					for(i = t; i<leftDegree*nbDataCols && mod2sparse_find(pchkMatrix,u[i],j); i++) ;

					if(i < leftDegree*nbDataCols)
					{
						/* choose one index within the list of possible choices */
						do {
							i = t + ldpc_rand(leftDegree*nbDataCols-t);
						} while (mod2sparse_find(pchkMatrix,u[i],j));
						mod2sparse_insert(pchkMatrix,u[i],j);
						/* replace with u[t] which has never been chosen */
						u[i] = u[t];
						t++;
					}
					else
					{
						/* no choice left, choose one randomly */
						uneven += 1;
						do {
							i = ldpc_rand(nbRows);
						} while (mod2sparse_find(pchkMatrix,i,j));
						mod2sparse_insert(pchkMatrix,i,j);
					}
				}
			}

			free(u);	/* VR: added */
			break;

		default: abort();
	}

	/* Add extra bits to avoid rows with less than two checks. */
	added = 0;
	for(i = 0; i<nbRows; i++)
	{
		e = mod2sparse_first_in_row(pchkMatrix,i);
		if(mod2sparse_at_end(e))
		{
			j = (ldpc_rand(nbDataCols))+skipCols;
			e = mod2sparse_insert(pchkMatrix,i,j);
			added ++;
		}
		e = mod2sparse_first_in_row(pchkMatrix,i);
		if(mod2sparse_at_end(mod2sparse_next_in_row(e)) && nbDataCols>1)
		{ 
			do 
			{ 
				j = (ldpc_rand(nbDataCols))+skipCols; 
			} while (j==mod2sparse_col(e));
			mod2sparse_insert(pchkMatrix,i,j);
			added ++;
		}
	}

	/* Add extra bits to try to avoid problems with even column counts. */
	if(leftDegree%2==0 && leftDegree<nbRows && nbDataCols>1 && added<2)
	{
		int a;
		for(a = 0; added+a<2; a++)
		{
			do
			{
				i = ldpc_rand(nbRows);
				j = (ldpc_rand(nbDataCols))+skipCols;
			} while (mod2sparse_find(pchkMatrix,i,j));
			mod2sparse_insert(pchkMatrix,i,j);
		}
	}

	switch (type) {
		case TypeLDGM:
			for (i = 0; i < nbRows; i++) {
				/* identity */
				mod2sparse_insert(pchkMatrix, i, i);
			}
			break;

		case TypeSTAIRS:
			mod2sparse_insert(pchkMatrix, 0, 0);	/* 1st row */
			for (i = 1; i < nbRows; i++) {		/* for all other rows */
				/* identity */
				mod2sparse_insert(pchkMatrix, i, i);
				/* staircase */
				mod2sparse_insert(pchkMatrix, i, i-1);
			}
			break;

		case TypeTRIANGLE:
			mod2sparse_insert(pchkMatrix, 0, 0);	/* 1st row */
			for (i = 1; i < nbRows; i++) {		/* for all other rows */
				/* identity */
				mod2sparse_insert(pchkMatrix, i, i);
				/* staircase */
				mod2sparse_insert(pchkMatrix, i, i-1);
				/* triangle */	
				j = i-1;
				for (l = 0; l < j; l++) { /* limit the # of "1s" added */
					j = ldpc_rand(j);
					mod2sparse_insert(pchkMatrix, i, j);
				}
			}
			break;
	}

	return pchkMatrix;
}

