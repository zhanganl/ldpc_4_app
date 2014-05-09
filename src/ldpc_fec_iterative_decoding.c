#include "ldpc_fec.h"

/******************************************************************************/
/*
 * Decoder using the Iterative Decoding Algorithm.
 */
/******************************************************************************
 * DecodingStepWithSymbol: Perform a new decoding step with a new (given) symbol.
 * This is the legacy front end to the DecodingStepWithSymbol() method. The actual
 * work will be done in the other DecodingStepWithSymbol() method.
 * => See header file for more informations.
 */ 
	ldpc_error_status
DecodingWithSymbol (
		LDPCFecSession *Session, 
		void*	symbol_canvas[],
		void*	new_symbol,
		int	new_symbol_seqno,
		bool	store_symbol)
{
	void	*new_symbol_dst;	// temp variable used to store symbol
	LDPC_head data_head;

	// Fast path. If store symbol is not set, then call directly
	// the full DecodingStepWithSymbol() method to avoid duplicate processing.

	if (store_symbol == false) {
		return(DecodingStepWithSymbol(Session, symbol_canvas, new_symbol, new_symbol_seqno)); 
	}
	// Step 0: check if this is a fresh symbol, otherwise return
	if ((mod2sparse_last_in_col(Session->m_pchkMatrix, GetMatrixCol(Session, new_symbol_seqno))->row < 0)
			|| (IsSourceSymbol(Session, new_symbol_seqno) && (symbol_canvas[new_symbol_seqno] != NULL))
			|| (IsParitySymbol(Session, new_symbol_seqno) && (Session->m_parity_symbol_canvas[new_symbol_seqno - Session->m_nbSourceSymbols] != NULL))) {
		// Symbol has already been processed, so skip it
		return LDPC_OK;
	}
	// Step 1: Store the symbol in a permanent array if the caller wants it.
	// It concerns only source symbols, since parity symbols are only stored in
	// permanent array if we have a memory gain by doing so, which will
	// be defined later on in the full DecodingStepWithSymbol() method.
	if (IsSourceSymbol(Session, new_symbol_seqno)) {
		// Call any required callback, or allocate memory, and
		// copy the symbol content in it.
		// This is typically something which is done when this
		// function is called recursively, for newly decoded
		// symbols.
		new_symbol_dst = (void *)malloc(Session->m_symbolSize);
		if (new_symbol_dst == NULL) {
			return LDPC_ERROR;
		}
		// Copy data now
		memcpy(GetBufferPtrOnly(new_symbol_dst), GetBuffer(new_symbol), Session->m_symbolSize);
	} else {
		new_symbol_dst = new_symbol;
	}
	/* continue decoding with the full DecodingStepWithSymbol() method */
	return(DecodingStepWithSymbol(Session, symbol_canvas, new_symbol_dst, new_symbol_seqno)); 
}


/******************************************************************************
 * DecodingStepWithSymbol: Perform a new decoding step with a new (given) symbol.
 * => See header file for more informations.
 *
 * This function relies on the following simple algorithm:
 *
 * Given a set of linear equations, if one of them has only one
 * remaining unknown variable, then the value of this variable is
 * that of the constant term.
 * Replace this variable by its value in all remaining linear
 * equations, and reiterate. The value of several variables can
 * therefore be found by this recursive algorithm.
 *
 * In practice, an incoming symbol contains the value of the associated
 * variable, so replace its value in all linear equations in which
 * it is implicated. Then apply the above algorithm and see if decoding
 * can progress by one or more steps.
 *
 * For instance, if {s1, s2} are source symbols, and {f1} a parity symbol:
 *    | s1 + s2 + f1 = A      (eq. 1)
 *    |      s2 + f1 = B      (eq. 2)
 * Now if the node receives symbol s2, he replaces its value in the two
 * equations, he then finds f1, he replaces its value in the first equation
 * and finds s1.
 */
ldpc_error_status
DecodingStepWithSymbol(
		LDPCFecSession *Session,
		void*	symbol_canvas[],
		void*	new_symbol,
		int	new_symbol_seqno)
{
	mod2entry	*e = NULL;	// entry ("1") in parity check matrix
	mod2entry	*delMe;		// temp: entry to delete in row/column
	void		*currChk;	// temp: pointer to Partial sum
	int		row;		// temp: current row value
	int		*CheckOfDeg1 = NULL; // table of check nodes of degree
	// one after the processing of new_symbol
	int		CheckOfDeg1_nb = 0; // number of entries in table
	int		CheckOfDeg1_listSize = 0; // size of the memory block
	LDPC_head data_head;
	// allocated for the table
	bool		keep_symbol;	// true if it's worth to store new_symbol
	// in this function, in case it's a parity
	// symbol, and independantly from the
	// store_symbol argument.

	// Step 0: check if this is a fresh symbol, otherwise return
	if ((mod2sparse_last_in_col(Session->m_pchkMatrix, GetMatrixCol(Session, new_symbol_seqno))->row < 0)
			|| (IsSourceSymbol(Session, new_symbol_seqno) && (symbol_canvas[new_symbol_seqno] != NULL))
			|| (IsParitySymbol(Session, new_symbol_seqno) && (Session->m_parity_symbol_canvas[new_symbol_seqno - Session->m_nbSourceSymbols] != NULL))) {
		// Symbol has already been processed, so skip it
		return LDPC_OK;
	}
	// First, make sure data is available for this new symbol. Must
	// remain valid throughout this function...
	GetBuffer(new_symbol);

	// Step 1: Store the symbol in a permanent array. It concerns only DATA
	// symbols. Parity symbols are only stored in permanent array, if we have
	// a memory gain by doing so (and not creating new partial sums)
	if (IsSourceSymbol(Session, new_symbol_seqno)) {
		// Source symbol
		// There's no need to allocate anything, nor to call
		// anything. It has already been done by the caller...
		symbol_canvas[new_symbol_seqno] = new_symbol;
		if (IsDecodingComplete(Session, symbol_canvas)) {
			// Decoding is now finished, return...
			return LDPC_OK;
		}
		keep_symbol = true;
	} else {
		// Parity symbol
		// Check if parity symbol should be stored or if partial
		// sum should be stored
		if (Session->m_triangleWithSmallFECRatio) {
			// In this case, the symbol will never be stored into
			// permanent array, but directly added to partial sum
			keep_symbol = false;
		} else {
			// The symbol will be stored if more than 1 partial sum
			// is needed
			const int	max_allowed_PS = 1; // threshold before 
			// storing symbol in array
			int		PS_to_create = 0; // nb of partial sums that
			//  would have to be allocated

			// count the number of PS that would be created
			// if we don't keep this parity symbol.
			for (e = mod2sparse_first_in_col(Session->m_pchkMatrix,
						GetMatrixCol(Session, new_symbol_seqno));
					!mod2sparse_at_end(e);
					e = mod2sparse_next_in_col(e))
			{
				if (Session->m_checkValues[e->row] == NULL &&
						Session->m_nb_unknown_symbols[e->row] > 2) {
					PS_to_create++;
					if (PS_to_create > max_allowed_PS) {
						break;	// no need to see further
					}
				}
			}
			// now take a decision...
			if (PS_to_create > max_allowed_PS) {
				// Parity symbol will be stored in a permanent array
				// Alloc the buffer...
				keep_symbol = true;
				Session->m_parity_symbol_canvas[new_symbol_seqno - Session->m_nbSourceSymbols] =
					(void *)malloc(Session->m_symbolSize);
				// copy the content...
				memcpy(GetBufferPtrOnly(Session->m_parity_symbol_canvas[new_symbol_seqno - Session->m_nbSourceSymbols]),
						GetBufferPtrOnly(new_symbol), Session->m_symbolSize);
			} else {
				// Parity symbol will only be added to partial sums
				keep_symbol = false;
			}
		}
	}

	//printf("0: start step 2 for symbol seq=%d\n", new_symbol_seqno);
	// Step 2: Inject the symbol value in each equation it is involved
	// (if partial sum already exists or if partial sum should be created)
	for (e = mod2sparse_first_in_col(Session->m_pchkMatrix, GetMatrixCol(Session, new_symbol_seqno));
			!mod2sparse_at_end(e); ) {
		// for a given row, ie for a given equation where this symbol
		// is implicated, do the following:
		row = e->row;
		Session->m_nb_unknown_symbols[row]--;	// symbol is known
		currChk = Session->m_checkValues[row];	// associated check
		if (currChk == NULL &&
				((keep_symbol == false) || (Session->m_nb_unknown_symbols[row] == 1))
		   ) {
			// we need to allocate a PS (i.e. check node)
			// and add symbol to it, because the parity symbol
			// won't be kept (keep_symbol == false), or it is the
			// last missing symbol of this equation, or because
			// or some particular situation where it is non sense
			// no to allocate a PS (m_triangleWithSmallFECRatio).
			{
				currChk = (void*) calloc(Session->m_symbolSize, 1);
				memset(currChk, 0, Session->m_symbolSize);
			}
			if ((Session->m_checkValues[row] = currChk) == NULL) {
				goto no_mem;
			}
		}
		if (currChk != NULL) {
			// there's a partial sum for this row...
			if (Session->m_nbSymbols_in_equ[row] > 1) {
				// Security: make sure data is available by
				// calling GetBuffer(new_symbol) rather than
				// GetBufferPtrOnly(new_symbol).
				// we can add the symbol content to this PS
				//printf("3': after add to currChk, fromdatabuf=x%x\n", GetBufferPtrOnly(new_symbol));
				AddToSymbol(GetBufferPtrOnly(currChk),
						GetBuffer(new_symbol));
				//printf("3: before add to currChk, to databuf=x%x\n", GetBufferPtrOnly(currChk));
			}
			// else this is useless, since new_symbol is the last
			// symbol of this equation, and its value is necessary
			// equal to the PS. Their sum must be 0 (we don't
			// check it).

			// remove the symbol from the equation (this entry
			// is now useless)
			delMe = e;
			e = mod2sparse_next_in_col(e);
			mod2sparse_delete(Session->m_pchkMatrix, delMe);
			Session->m_nbSymbols_in_equ[row]--;
			if (IsParitySymbol(Session, new_symbol_seqno)) {
				Session->m_nbEqu_for_parity[new_symbol_seqno - Session->m_nbSourceSymbols]--;
			}

			// Inject all permanently stored symbols (DATA and parity)
			// into partial sum
			if ((Session->m_nb_unknown_symbols[row] == 1) ||
					((keep_symbol == false) && (Session->m_triangleWithSmallFECRatio == false)))
			{
				// Inject all permanently stored symbols
				// (DATA and parity) into this partial sum.
				// Requires to scan the equation (ie row).
				mod2entry	*tmp_e;	// curr symbol in this equ
				int		tmp_seqno;// corresponding seq no
				void		*tmp_symbol; // corresponding symbol pointer

				for (tmp_e = mod2sparse_first_in_row(Session->m_pchkMatrix, row);
						!mod2sparse_at_end(tmp_e); ) {

					tmp_seqno = GetSymbolSeqno(Session, tmp_e->col);
					if (IsParitySymbol(Session, tmp_seqno)) {
						tmp_symbol = Session->m_parity_symbol_canvas[tmp_seqno - Session->m_nbSourceSymbols];
					} else {
						// waiting for
						// (m_nb_unknown_symbols[row] == 1)
						// to add source symbols is
						// useless... it's even slower
						// in that case!
						tmp_symbol = symbol_canvas[tmp_seqno];
					}
					if (tmp_symbol != NULL) {
						// add the symbol content now
						//printf("4: ready to add to currChk the tmp_symbol seq=%d, fromdatabuf=x%x\n", tmp_seqno, GetBufferPtrOnly(tmp_symbol));
						AddToSymbol(
								GetBufferPtrOnly(currChk),
								GetBuffer(tmp_symbol));
						//printf("5: add to currChk done, todatabuf=x%x\n", GetBufferPtrOnly(currChk));
						// delete the entry
						delMe = tmp_e;
						tmp_e =  mod2sparse_next_in_row(tmp_e);
						mod2sparse_delete(Session->m_pchkMatrix, delMe);
						Session->m_nbSymbols_in_equ[row]--;
						if (IsParitySymbol(Session, tmp_seqno)) {
							Session->m_nbEqu_for_parity[tmp_seqno - Session->m_nbSourceSymbols]--;
							// check if we can delete
							// parity symbol altogether
							if (Session->m_nbEqu_for_parity[tmp_seqno - Session->m_nbSourceSymbols] == 0) {
								free(tmp_symbol);
								Session->m_parity_symbol_canvas[tmp_seqno - Session->m_nbSourceSymbols] = NULL;
							}
						}
					} else {
						// this symbol not yet known,
						// switch to next one in equ
						tmp_e =  mod2sparse_next_in_row(tmp_e);
					}
				}
			}
		} else {
			// here m_checkValues[row] is NULL, ie. the partial
			// sum has not been allocated
			e = mod2sparse_next_in_col(e);
		}
		if (Session->m_nbSymbols_in_equ[row] == 1) {
			// register this entry for step 3 since the symbol
			// associated to this equation can now be decoded...
			if (CheckOfDeg1 == NULL) {
				// allocate memory for the table first
				CheckOfDeg1_listSize = 4;
				if ((CheckOfDeg1 = (int*)
							calloc(CheckOfDeg1_listSize,
								sizeof(int*))) == NULL) {
					goto no_mem;
				}
			} else if (CheckOfDeg1_nb == CheckOfDeg1_listSize) {
				// not enough size in table, add some more
				CheckOfDeg1_listSize += 4;
				if ((CheckOfDeg1 = (int*)realloc(CheckOfDeg1,
								CheckOfDeg1_listSize * sizeof(int*))) == NULL) {
					goto no_mem;
				}
			}
			CheckOfDeg1[CheckOfDeg1_nb++] = row;
		}
	}

	// Step 3: Check if a new symbol has been decoded and take appropriate
	// measures ...
	int	decoded_symbol_seqno;	// sequence number of decoded symbol
	//for (int i = 0; i < CheckOfDeg1_nb; i++) 
	for (CheckOfDeg1_nb--; CheckOfDeg1_nb >= 0; CheckOfDeg1_nb--) {
		if (IsDecodingComplete(Session, symbol_canvas)) {
			// decoding has just finished, no need to do anything else
			break;
		}
		// get the index (ie row) of the partial sum concerned
		row = CheckOfDeg1[CheckOfDeg1_nb];
		if (Session->m_nbSymbols_in_equ[row] == 1) {
			// A new decoded symbol is available...
			// NB: because of the recursion below, we need to
			// check that all equations mentioned in the
			// CheckOfDeg1 list are __still__ of degree 1.
			e = mod2sparse_first_in_row(Session->m_pchkMatrix, row);
			decoded_symbol_seqno = GetSymbolSeqno(Session, e->col);
			// remove the entry from the matrix
			currChk = Session->m_checkValues[row];	// remember it
			Session->m_checkValues[row] = NULL;
			Session->m_nbSymbols_in_equ[row]--;
			if (IsParitySymbol(Session, decoded_symbol_seqno)) {
				Session->m_nbEqu_for_parity[decoded_symbol_seqno - Session->m_nbSourceSymbols]--;
			}
			mod2sparse_delete(Session->m_pchkMatrix, e);
			if (IsSourceSymbol(Session, decoded_symbol_seqno)) {
				// source symbol.
				void	*decoded_symbol_dst;// temp variable used to store symbol

				// First copy it into a permanent symbol.
				// Call any required callback, or allocate memory, and
				// copy the symbol content in it.
				decoded_symbol_dst =
					(void *)malloc(Session->m_symbolSize);
				if (decoded_symbol_dst == NULL) {
					goto no_mem;
				}
				memcpy(GetBufferPtrOnly(decoded_symbol_dst),
						GetBuffer(currChk), Session->m_symbolSize);
				// Free partial sum which is no longer used.
				// It's important to free it before calling
				// DecodingStepWithSymbol recursively to reduce max
				// memory requirements.
				free(currChk);
				//printf("get data buf seqno:%d\n", decoded_symbol_seqno);
				memcpy(&data_head, decoded_symbol_dst, sizeof(data_head));
				data_head.type_flag = 0;
				data_head.longest_length = Session->m_symbolSize;
				memcpy(decoded_symbol_dst, &data_head, sizeof(data_head));

				// And finally call this method recursively
				DecodingStepWithSymbol(Session, symbol_canvas, decoded_symbol_dst,
						decoded_symbol_seqno);

			} else {
				//printf("get fec buf seqno:%d\n", decoded_symbol_seqno);
				memcpy(&data_head, currChk, sizeof(data_head));
				data_head.type_flag = 1;
				data_head.longest_length = Session->m_symbolSize;
				memcpy(currChk, &data_head, sizeof(data_head));

				// Parity symbol.
				// Call this method recursively first...
				DecodingStepWithSymbol(Session, symbol_canvas, currChk,
						decoded_symbol_seqno);
				// Then free the partial sum which is no longer needed.
				free(currChk);	
			}
		}
	}
	if (CheckOfDeg1 != NULL) {
		free(CheckOfDeg1);
	}
	return LDPC_OK;

no_mem:
	return LDPC_ERROR;
}

