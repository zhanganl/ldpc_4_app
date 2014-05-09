#ifndef LDPC_TYPES_H
#define LDPC_TYPES_H

/****** type specifications ******/
#ifndef UINT32
#define INT8    char
#define INT16   short
#define UINT8   unsigned char
#define UINT16  unsigned short
#if defined(__LP64__) || (__WORDSIZE == 64) /* 64 bit architectures */
#define INT32   int		// yes, it's also true in LP64!
#define UINT32  unsigned int	// yes, it's also true in LP64!
#else  /* 32 bit architectures */
#define INT32   int		// int creates less compilations pbs than long
#define UINT32  unsigned int	// int creates less compilations pbs than long
#endif /* 32/64 architectures */
#endif /* !UINT32 */

#ifndef UINT64
#define INT64   long long
#define UINT64  unsigned long long
#endif /* !UINT64 */

#endif /* LDPC_TYPES_H */
