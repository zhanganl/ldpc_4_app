#ifndef LDPC_MATRIX_SPARSE__
#define LDPC_MATRIX_SPARSE__

/**
 * Structure representing a non-zero entry, or the header for a row or column.
 */
typedef struct mod2entry
{
	/**
	 * Row and column indexes of this entry, starting at 0, and
	 * with -1 for a row or column header
	 */
	int	row;
	int	col;

	/**
	 * Pointers to entries adjacent in row and column, or to headers.
	 * Free entries are linked by 'left'.
	 */
	struct mod2entry	*left,
						*right,
						*down;
	struct mod2entry	*up;
} mod2entry;

#define Mod2sparse_block 10  /* Number of entries to block together for
								memory allocation */

/**
 * Block of entries allocated all at once.
 */
typedef struct mod2block
{
	/** Next block that has been allocated. */
	struct mod2block *next;

	/** Entries in this block. */
	mod2entry	entry[Mod2sparse_block];
} mod2block;

/**
 * Representation of a sparse matrix.
 */
typedef struct mod2sparse
{ 
	int n_rows;		  /* Number of rows in the matrix */
	int n_cols;		  /* Number of columns in the matrix */

	mod2entry *rows;	  /* Pointer to array of row headers */
	mod2entry *cols;	  /* Pointer to array of column headers */

	mod2block *blocks;	  /* Blocks that have been allocated */
	mod2entry *next_free;	  /* Next free entry */
} mod2sparse;

/* MACROS TO GET AT ELEMENTS OF A SPARSE MATRIX.  The 'first', 'last', 'next',
   and 'prev' macros traverse the elements in a row or column.  Moving past
   the first/last element gets one to a header element, which can be identified
   using the 'at_end' macro.  Macros also exist for finding out the row 
   and column of an entry, and for finding out the dimensions of a matrix. */

#define mod2sparse_first_in_row(m,i) ((m)->rows[i].right) /* Find the first   */
#define mod2sparse_first_in_col(m,j) ((m)->cols[j].down)  /* or last entry in */
#define mod2sparse_last_in_row(m,i) ((m)->rows[i].left)   /* a row or column  */
#define mod2sparse_last_in_col(m,j) ((m)->cols[j].up)

#define mod2sparse_next_in_row(e) ((e)->right)  /* Move from one entry to     */
#define mod2sparse_next_in_col(e) ((e)->down)   /* another in any of the three */
#define mod2sparse_prev_in_row(e) ((e)->left)   /* possible directions        */
#define mod2sparse_prev_in_col(e) ((e)->up) 

#define mod2sparse_at_end(e) ((e)->row<0) /* See if we've reached the end     */


#define mod2sparse_row(e) ((e)->row)      /* Find out the row or column index */
#define mod2sparse_col(e) ((e)->col)      /* of an entry (indexes start at 0) */

#define mod2sparse_rows(m) ((m)->n_rows)  /* Get the number of rows or columns*/
#define mod2sparse_cols(m) ((m)->n_cols)  /* in a matrix                      */

/* PROCEDURES TO MANIPULATE SPARSE MATRICES. */
mod2sparse *mod2sparse_allocate (int, int);
void mod2sparse_free            (mod2sparse *);

void mod2sparse_clear    (mod2sparse *);

void mod2sparse_print       (FILE *, mod2sparse *);

mod2entry *mod2sparse_find   (mod2sparse *, int, int);
mod2entry *mod2sparse_insert (mod2sparse *, int, int);
void mod2sparse_delete       (mod2sparse *, mod2entry *);

#endif // #ifndef LDPC_MATRIX_SPARSE__

