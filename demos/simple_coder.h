#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/time.h>	/* for gettimeofday */

#include "../src/ldpc_fec.h"
#include "../src/ldpc_group.h"

/*
 * OS dependant definitions
 */
#define SOCKET		int
#define INVALID_SOCKET	(-1)
#define SOCKET_ERROR	(-1)
#define closesocket	close
#define MSLEEP(t)	usleep(t*1000)


/*
 * Simulation parameters...
 * Change as required
 */
#define GROUP 	5
#define PKTSZ	1024		// Packets size, in bytes (multiple of 4).
#define NBDATA	10		// Number of original DATA packets to send.
#define NBFEC	10		// Number of FEC packets to build.
#define NBDROP	3		// NBPKT/NBDROP is Drop percent.
#define NBPKT	(NBDATA+NBFEC)	// Total number of packets to send.
#define LEFT_DEGREE	3	// Left degree of data nodes in the checks graph

/*
 * The Session Type.
 * Can be one of: TypeTRIANGLE, TypeLDGM, TypeSTAIRS, TypeLDPC
 */
//#define SESSION_TYPE	TypeLDPC
//#define SESSION_TYPE	TypeLDGM
//#define SESSION_TYPE	TypeSTAIRS
#define SESSION_TYPE  TypeTRIANGLE

#define SEED		2003	// Seed used to initialize LDPCFecSession

#define DEST_IP		"127.0.0.1"	// Destination IP
//#define DEST_IP		"194.188.224.103"	// Destination IP
#define DEST_PORT	10978		// Destination port (UDP)


