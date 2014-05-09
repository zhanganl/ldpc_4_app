#include "ldpc_fec.h"


/******************************************************************************
 * InitSession : Initializes the LDPC session.
 * => See header file for more informations.
 */
ldpc_error_status
InitSession (
		LDPCFecSession *Session,
		int nbSourceSymbols,
		int nbParitySymbols,
		int symbolSize,
		int flags,
		int seed,
		SessionType codecType,
		int leftDegree)
{
	int row, seq;
	mod2entry	*e;

	Session->m_initialized	= false;
	Session->m_sessionFlags	= flags;
	Session->m_sessionType	= codecType;
	Session->m_symbolSize	= symbolSize;

	if (nbSourceSymbols + nbParitySymbols > GetMaxN()) 
		return LDPC_ERROR;

	Session->m_nbSourceSymbols	= nbSourceSymbols;
	Session->m_nbParitySymbols	= nbParitySymbols;

	// sanity check: setting a left degree != 3 is only meaningfull with
	// LDGM codes, not with LDPC-Staircase/Triangle codes.
	if ((Session->m_sessionType != TypeLDGM) && (leftDegree != 3)) {
		fprintf(stderr, "LDPCFecSession::InitSession: ERROR: setting a leftDegree!=3 (here %d) is only meaningfull with LDGM codes, not with LDPC-Staircase/Triangle codes!\n", leftDegree);
		return LDPC_ERROR;
	}
	Session->m_leftDegree	= leftDegree;
	Session->m_firstNonDecoded = 0;

	Session->m_pchkMatrix = CreatePchkMatrix(Session->m_nbParitySymbols, Session->m_nbSourceSymbols + Session->m_nbParitySymbols, Evenboth, Session->m_leftDegree, seed, false, Session->m_sessionType);
	if (Session->m_pchkMatrix == NULL) 
		return LDPC_ERROR;

	if (Session->m_sessionFlags & FLAG_CODER) {
		Session->m_nb_unknown_symbols_encoder = (int*)calloc(Session->m_nbParitySymbols, sizeof(int));
		if (Session->m_nb_unknown_symbols_encoder == NULL) 
			return LDPC_ERROR;

		for (row = 0; row < Session->m_nbParitySymbols; row++) {
			mod2entry *e;
			for (e = mod2sparse_first_in_row(Session->m_pchkMatrix, row);
					!mod2sparse_at_end(e);
					e = mod2sparse_next_in_row(e)) {
				Session->m_nb_unknown_symbols_encoder[row]++;
			}
		}
	} else {
		Session->m_nb_unknown_symbols_encoder = NULL;
	}

	if (Session->m_sessionFlags & FLAG_DECODER) {
		if (((Session->m_checkValues	= (void**)calloc(Session->m_nbParitySymbols, sizeof(void*))) == NULL) ||
				((Session->m_nbSymbols_in_equ = (int*)calloc(Session->m_nbParitySymbols, sizeof(int))) == NULL) ||
				((Session->m_nb_unknown_symbols = (int*)calloc(Session->m_nbParitySymbols, sizeof(int))) == NULL) ||
				((Session->m_nbEqu_for_parity = (int*)calloc(Session->m_nbParitySymbols, sizeof(int))) == NULL) ||
				((Session->m_parity_symbol_canvas = (void**)calloc(Session->m_nbParitySymbols, sizeof(void*))) == NULL)) {
			return LDPC_ERROR;
		}
		// and update the various tables now
		for (row = 0; row < Session->m_nbParitySymbols; row++) {
			for (e = mod2sparse_first_in_row(Session->m_pchkMatrix, row);
					!mod2sparse_at_end(e);
					e = mod2sparse_next_in_row(e))
			{
				Session->m_nbSymbols_in_equ[row]++;
				Session->m_nb_unknown_symbols[row]++;
			}
		}
		for (seq = Session->m_nbSourceSymbols; seq < (Session->m_nbParitySymbols+Session->m_nbSourceSymbols); seq++) {
			for (e = mod2sparse_first_in_col(Session->m_pchkMatrix,
						GetMatrixCol(Session, seq));
					!mod2sparse_at_end(e);
					e = mod2sparse_next_in_col(e))
			{
				Session->m_nbEqu_for_parity[seq - Session->m_nbSourceSymbols]++;
			}
		}
	} else {
		// CODER session
		Session->m_checkValues = NULL;
		Session->m_nbSymbols_in_equ = NULL;
		Session->m_nb_unknown_symbols = NULL;
		Session->m_nbEqu_for_parity = NULL;
		Session->m_parity_symbol_canvas = NULL;
	}
	if ((Session->m_sessionType == TypeTRIANGLE) && (((Session->m_nbParitySymbols+Session->m_nbSourceSymbols)/Session->m_nbSourceSymbols) < 2.0)) {
		Session->m_triangleWithSmallFECRatio = true;
	} else {
		Session->m_triangleWithSmallFECRatio = false;
	}
	Session->m_initialized = true;
	return LDPC_OK;
}


	void
EndSession(LDPCFecSession *Session)
{
	int i;

	if (Session->m_initialized == true) {
		Session->m_initialized = false;
		mod2sparse_free(Session->m_pchkMatrix);
		free(Session->m_pchkMatrix);	/* mod2sparse_free does not free it! */

		if (Session->m_checkValues != NULL) {
			for (i = 0; i < Session->m_nbParitySymbols; i++) {
				if (Session->m_checkValues[i] != NULL) {
					{
						free(Session->m_checkValues[i]);
					}
				}
			}
			free(Session->m_checkValues);
		}
		if (Session->m_parity_symbol_canvas != NULL) {
			for (i = 0; i < Session->m_nbParitySymbols; i++) {
				if (Session->m_parity_symbol_canvas[i] != NULL) {
					free(Session->m_parity_symbol_canvas[i]);
				}
			}
			free(Session->m_parity_symbol_canvas);
		}
		if (Session->m_nbSymbols_in_equ != NULL) {
			free(Session->m_nbSymbols_in_equ);
		}
		if (Session->m_nbEqu_for_parity != NULL) {
			free(Session->m_nbEqu_for_parity);
		}
		if (Session->m_nb_unknown_symbols != NULL) {
			free(Session->m_nb_unknown_symbols);
		}
		if (Session->m_nb_unknown_symbols_encoder != NULL) {
			free(Session->m_nb_unknown_symbols_encoder);
		}
	}
}

/******************************************************************************
 * Calculates the XOR sum of two symbols: to = to + from.
 * => See header file for more informations.
 */
	void
AddToSymbol	(void	*to,
		void	*from)
{
	unsigned int		i;
	LDPC_head data_head;
	unsigned int 	data_size;

	memcpy(&data_head, from, sizeof(data_head));

	if(data_head.type_flag)
		data_size = data_head.longest_length - sizeof(data_head) + sizeof(unsigned short);
	else
		data_size = data_head.current_length - sizeof(data_head) + sizeof(unsigned short);

	UINT8		*t = (UINT8*)to;	// to pointer to 32-bit integers
	UINT8		*f = (UINT8*)from;	// from pointer	to 32-bit integers
	t += sizeof(data_head) - sizeof(unsigned short);
	f += sizeof(data_head) - sizeof(unsigned short);

	for (i = data_size; i > 0; i--) {
		*t ^= *f;
		t++;
		f++;
	}
}


/******************************************************************************
 * BuildParitySymbol: Builds a new parity symbol.
 * => See header file for more informations.
 */
	ldpc_error_status
BuildParitySymbol (
		LDPCFecSession *Session, 
		void* symbol_canvas[],
		int paritySymbol_index,
		void* paritySymbol)
{
	uintptr_t	*fec_buf;	// buffer for this parity symbol
	uintptr_t	*to_add_buf;	// buffer for the  source.parity symbol to add
	mod2entry	*e;
	int seqno, k=0;
	static int n=0;

	fec_buf = (uintptr_t*)GetBufferPtrOnly(paritySymbol);

	e = mod2sparse_first_in_row(Session->m_pchkMatrix, paritySymbol_index);

	while (!mod2sparse_at_end(e)) {
		// paritySymbol_index in {0.. n-k-1} range, so this test is ok
		if (e->col != paritySymbol_index) {
			// don't add paritySymbol to itself
			seqno = GetSymbolSeqno(Session, e->col);
			//fprintf(stderr, "[%s:%d] total:%d now:%d seqno:%d, paritySymbol_index:%d, col:%d\n", __FILE__, __LINE__, n++, k++, seqno, paritySymbol_index, e->col);
			to_add_buf = (uintptr_t*)
				GetBuffer(symbol_canvas[seqno]);
			if (to_add_buf == NULL) {
				return LDPC_ERROR;
			}
			AddToSymbol(fec_buf, to_add_buf);
		}
		e = mod2sparse_next_in_row(e);
	}
	return LDPC_OK;
}

/******************************************************************************
 * GetMatrixCol:
 * => See header file for more informations.
 */
	int
GetMatrixCol(LDPCFecSession *Session, int symbolSeqno)
{
	if (symbolSeqno < Session->m_nbSourceSymbols) {
		/* source symbol */
		return (symbolSeqno + Session->m_nbParitySymbols);
	} else {
		/* parity symbol */
		return (symbolSeqno - Session->m_nbSourceSymbols);
	}
}


/******************************************************************************
 * GetSymbolSeqno:
 * => See header file for more informations.
 */
	int
GetSymbolSeqno(LDPCFecSession *Session, int matrixCol)
{
	int colInOrder;

	colInOrder = matrixCol;
	if (colInOrder < Session->m_nbParitySymbols) {
		/* parity symbol */
		return (colInOrder + Session->m_nbSourceSymbols);
	} else {
		/* source symbol */
		return (colInOrder - Session->m_nbParitySymbols);
	}
}


/******************************************************************************
 * IsDecodingComplete: Checks if all DATA symbols have been received/rebuilt.
 * => See header file for more informations.
 */
	bool
IsDecodingComplete(LDPCFecSession *Session, void*	symbol_canvas[])
{
	int i;

	if (!Session->m_initialized) {
		fprintf(stderr, "LDPCFecSession::IsDecodingComplete: ERROR: LDPC Session is NOT initialized!\n");
		return false;
	}

	for (i = Session->m_firstNonDecoded; i < Session->m_nbSourceSymbols; i++) {
		if (symbol_canvas[i] == NULL) {
			/* not yet decoded! */
			Session->m_firstNonDecoded = i;	/* remember for next time */
			return false;
		}
	}
	return true;
}

//------------------------------------------------------------------------------
// Inlines for all classes follow
//------------------------------------------------------------------------------

	bool
IsInitialized( LDPCFecSession *Session )
{
	return Session->m_initialized;
}

	bool
IsSourceSymbol	(LDPCFecSession *Session, int	symbolSeqno)
{
	return ((symbolSeqno < Session->m_nbSourceSymbols) ? true : false);
}

	bool	
IsParitySymbol	(LDPCFecSession *Session, int	symbolSeqno)
{
	return ((symbolSeqno < Session->m_nbSourceSymbols) ? false : true);
}
