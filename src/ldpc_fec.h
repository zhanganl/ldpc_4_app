#ifndef LDPC_FEC_H /* { */
#define LDPC_FEC_H

#include <math.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "ldpc_create_pchk.h"
#include "ldpc_types.h"
#include "ldpc_matrix_sparse.h"

/****** CONSTANT AND CLASS DEFINITION *****************************************/

/**             
 * Error status returned by functions.
 */             
typedef enum {
	LDPC_OK = 0,
	LDPC_ERROR = 1
}ldpc_error_status;


/**
 * Is the session a coding or decoding session, or both.
 */
#define FLAG_CODER	0x00000001
#define FLAG_DECODER	0x00000002
#define FLAG_BOTH (FLAG_DECODER|FLAG_CODER)

typedef struct {
	unsigned int 	type_flag:1;
	unsigned int 	group_id:31;
	unsigned int 	sequence_no;
	unsigned short 	total_data;
	unsigned short 	total_fec;
	unsigned short 	longest_length;
	unsigned short 	current_length;
}LDPC_head;

typedef struct {
	bool	m_initialized;	// is TRUE if session has been initialized
	int	m_sessionFlags;	// Mask containing session flags
	// (FLAG_CODER, FLAG_DECODER, ...)
	SessionType	m_sessionType;	// Type of the session. Can be one of
	// LDGM, LDPC STAIRS, and
	// LDPC TRIANGLE.

	unsigned int	m_symbolSize;	// Size of symbols in BYTES

	int	m_nbSourceSymbols;	// number fo source symbol (K)
	int	m_nbParitySymbols;	// number of parity symbol (=m_nbCheck)

	mod2sparse*	m_pchkMatrix;	// Parity Check matrix in sparse mode 
	// format. This matrix is also used as
	// a generator matrix in LDGM-* modes.

	int		m_leftDegree;	// Number of equations per data symbol

	// Encoder specific...
	int*		m_nb_unknown_symbols_encoder; // Array: nb unknown symbols
	// per check node. Used during per column
	// encoding.

	// Decoder specific...
	void**		m_checkValues;	// Array: current check-nodes value.
	// Each entry is the sum (XOR) of some
	// or all of the known symbols in this
	// equation.
	int*		m_nbSymbols_in_equ;// Array: nb of variables per check
	// node, ie. per equation
	int		m_firstNonDecoded; // index of first symbol not decoded.
	// Used to know whether decoding is
	// finished or not.
	int*		m_nb_unknown_symbols; // Array: nb unknown symbols per check node
	int*		m_nbEqu_for_parity; // Array: nb of equations where
	// each parity symbol is included
	void**		m_parity_symbol_canvas; //Canvas of stored parity symbols.

#if 0
	uintptr_t* 	m_builtSymbol; 	// symbol built by decoder, used for 
	// recursive calls of DecodingStepWithSymbol
#endif
	bool		m_triangleWithSmallFECRatio;
	// with LDGM Triangle and a small FEC
	// ratio (ie. < 2), some specific
	// behaviors are needed...

	void*		m_context_4_callback; // used by callback functions
}LDPCFecSession;


/**
 * InitSession: Initializes the LDPC session.
 * @param nbSourceSymbols	(IN) number of source symbols (i.e. k).
 * @param nbParitySymbols	(IN) number of parity symbols (i.e. n-k).
 *			Be careful that n-k cannot be less than the left
 *			degree (i.e. 3 by default), otherwise an error is
 *			returned.
 * @param symbolSize	(IN) symbol size in bytes. It does NOT need to be
 *			multiple of 4, any value is accepted.
 * @param flags		(IN) session flags (FLAG_CODER, FLAG_DECODER, ...).
 * @param seed		(IN) seed used to build the parity check matrix (H).
 * @param codecType	(IN) Type of codec algorithm and matrix to use.
 *			Can be on of TypeLDGM, TypeSTAIRS, or TypeTRIANGLE.
 * @param leftDegree	(IN) number of equations in which a symbol is involved.
 *			3 (default) is the optimal value for TypeSTAIRS
 *			and TypeTRIANGLE codes, DO NOT change.
 *			With TypeLDGM, higher values are usually preferable
 *			(see INRIA Research Report 5225, June 2004).
 * @return		Completion status (LDPC_OK or LDPC_ERROR).
 */
ldpc_error_status InitSession  (
		LDPCFecSession *Session, 
		int	nbSourceSymbols,
		int	nbParitySymbols,
		int	symbolSize,
		int	flags,
		int	seed,
		SessionType	codecType,
		int	leftDegree);


/**
 * EndSession: Ends the LDPC session, cleans up everything.
 */
void EndSession (LDPCFecSession *Session);


/**
 * IsInitialized: Check if the LDPC session has been initialized.
 * @return	  TRUE if the session is ready and initialized, FALSE if not.
 */
bool IsInitialized (LDPCFecSession *Session);


/**
 * Build a new parity symbol.
 * @param symbol_canvas	(IN)	Array of source and parity symbols.
 *				This is a table of n pointers to buffers
 *				containing the source and parity symbols.
 * @param paritySymbol_index	(IN)	Index of parity symbol to build in {0.. n-k-1}
 *				range (!)
 * @param paritySymbol	(IN-OUT) Pointer to the parity symbol buffer that will
 *				be built. This buffer MUST BE allocated
 *				before,	but NOT cleared (memset(0)) since
 *				this function will do it.
 * @return			Completion status (LDPC_OK or LDPC_ERROR).
 */
ldpc_error_status BuildParitySymbol (
		LDPCFecSession *Session, 
		void*		symbol_canvas[],
		int		paritySymbol_index,
		void*		paritySymbol); 


/**
 * Perform a new decoding step thanks to the newly received symbol.
 * @param symbol_canvas	(IN-OUT) Global array of received or rebuilt source
 * 				symbols (parity symbols need not be stored here).
 *				This is a table of k pointers to buffers.
 * 				This array must be cleared (memset(0)) upon
 * 				the first call to this function. It will be
 * 				automatically updated, with pointers to
 * 				symbols received or decoded, by this function.
 * @param new_symbol	(IN)	Pointer to the buffer containing the new symbol.
 * @param new_symbol_seqno	(IN)	New symbol's sequence number in {0.. n-1} range.
 * @param store_symbol	(IN)	true if the function needs to allocate memory,
 *				copy the symbol content in it, and call
 *				any required callback.
 *				This is typically done when this function is
 *				called recursively, for newly decoded symbols,
 *				or under special circunstances (e.g. perftool).
 * @return			Completion status (LDPC_OK or LDPC_ERROR).
 */
ldpc_error_status DecodingWithSymbol (
		LDPCFecSession *Session, 
		void*	symbol_canvas[],
		void*	new_symbol,
		int	new_symbol_seqno,
		bool	store_symbol);


/**
 * Perform a new decoding step thanks to the newly received symbol.
 * Same as the other DecodingStepWithSymbol method, without the store_symbol argument
 * (prefered solution).
 * @param symbol_canvas	(IN-OUT) Global array of received or rebuilt source
 * 				symbols (parity symbols need not be stored here).
 *				This is a table of k pointers to buffers.
 * 				This array must be cleared (memset(0)) upon
 * 				the first call to this function. It will be
 * 				automatically updated, with pointers to
 * 				symbols received or decoded, by this function.
 * @param new_symbol	(IN)	Pointer to the buffer containing the new symbol.
 * @param new_symbol_seqno	(IN)	New symbol's sequence number in {0.. n-1} range.
 * @return			Completion status (LDPC_OK or LDPC_ERROR).
 */
ldpc_error_status DecodingStepWithSymbol (
		LDPCFecSession *Session, 
		void*	symbol_canvas[],
		void*	new_symbol,
		int	new_symbol_seqno);


/**
 * Checks if all DATA symbols have been received/rebuilt.
 * @param symbol_canvas	(IN)	Array of received/rebuilt source symbols.
 * @return			TRUE if all DATA symbols have been received
 * 				or decoded.
 */
bool IsDecodingComplete (LDPCFecSession *Session, void*	symbol_canvas[] );


/**
 * Return true if this is a DATA source symbol.
 * @param symbolSeqno	(IN) symbol sequence number in {O; n-1} range
 * @return		true if source symbol, false if parity symbol
 */
bool	IsSourceSymbol	(LDPCFecSession *Session, int	symbolSeqno);

/**
 * Return true if this is a parity (AKA FEC) symbol.
 * @param symbolSeqno	(IN) symbol sequence number in {O; n-1} range
 * @return		true if parity symbol, false if source symbol
 */
bool	IsParitySymbol	(LDPCFecSession *Session, int	symbolSeqno);

/**
 * symbol sequence number to column index translation.
 * @param symbolSeqno	(IN) symbol sequence number in {O; n-1} range
 * @return		corresponding column number in matrix
 */
int	GetMatrixCol	(LDPCFecSession *Session, int symbolSeqno);

/**
 * Internal column index to symbol sequence number translation.
 * @param matrixCol	(IN) column number in matrix
 * @return		corresponding symbol sequence number in
 *			{O; n-1} range
 */
int	GetSymbolSeqno	(LDPCFecSession *Session, int matrixCol);

/**
 * Calculates the XOR sum of two symbols: to = to + from.
 * @param to		(IN/OUT) source symbol
 * @param from		(IN/OUT) symbol added to the source symbol
 */
void	AddToSymbol	(void	*to,
		void	*from);


/**
 * Returns the maximum encoding block length (n parameter).
 * This limit is not LDPC-* specific that are nature large bloc FEC codes,
 * meaning that (k,n) can both be very very large. This is a codec specific
 * limit, due to the way the codec is implemented.
 * See ldpc_profile.h:
 *	If SPARSE_MATRIX_OPT_SMALL_INDEX is defined, then
 *		k <= n < 2^15;
 *	Else
 *		k <= n < 2^31
 * The limits are essentially over the n parameter, but given the
 * desired FEC Expansion ratio n/k (or the code rate, k/n), it will also
 * limit the source block length (k parameter).
 * @return			Maximum n value (A.K.A. encoding block length).
 */
	static inline int
GetMaxN ()
{
	return 0x7FFFFFFF;
}

/**
 * Get the data buffer associated to a symbol stored in the
 * symbol_canvas[] / m_parity_symbol_canvas[] / m_checkValues[] tables.
 * This function is usefull in EXTERNAL_MEMORY_MGMT_SUPPORT mode
 * when a the Alloc/Get/Store/Free callbacks are used, but it does
 * nothing in other mode. This is due to the fact that with these
 * callbacks, the various canvas do not point to data buffers but
 * to intermediate structures, and therefore accessing the associated
 * buffer needs extra processing.
 * @param symbol	(IN) pointer stored in the various canvas
 * @return		associated buffer
 */
	static inline void*
GetBuffer	(void	*symbol)
{
	return symbol;		// nothing to do here
}

/**
 * Same as GetBuffer, except that this call does not use the
 * GetData_callback but GetDataPtrOnly_callback instead.
 * For instance, in EXTERNAL_MEMORY_MGMT_SUPPORT, it will not
 * make sure that data is actually available and up-to-date,
 * perhaps because this is a destination buffer in a memcpy
 * that has just been allocated!
 * @param symbol	(IN) pointer stored in the various canvas
 * @return		associated buffer
 */
	static inline void*
GetBufferPtrOnly	(void	*symbol)
{
	return symbol;		// nothing to do here
}

#endif /* } LDPC_FEC_H */
