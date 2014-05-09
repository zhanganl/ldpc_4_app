#include <sys/types.h>
#include <sys/socket.h>
#include <mcheck.h>
#include "simple_coder.h"

/* Prototypes */
SOCKET	initSocket( void );
void	randomizeArray( int**, int );
void	DumpBuffer( char*, int );

int main(int argc, char* argv[])
{
	char**	packetsArray	= NULL;

	SOCKET	mySock		= INVALID_SOCKET;
	int*	randOrder1	= NULL, *randOrder2 = NULL;
	int	ret		= -1;
	int	i		= 0, j=0;
	int	data_size;
	LDPC_head data_head;
	char **all_data[GROUP];
	LDPCFecSession Session[GROUP];

	mtrace();

	for(j=0; j<GROUP; j++)
	{
		memset(&Session[j], 0, sizeof(Session[i]));

		// Initialize the LDPC session
		if(InitSession(&Session[j], NBDATA, NBFEC, PKTSZ+sizeof(LDPC_head), FLAG_CODER, SEED, SESSION_TYPE, LEFT_DEGREE ) == LDPC_ERROR)
		{
			printf("Error: Unable to initialize LDPC Session\n");
			ret = -1; goto cleanup;
		}
		packetsArray = (char**)calloc( NBPKT, sizeof(char*) );
		if( packetsArray == NULL ) {
			printf("Error: insufficient memory (calloc failed for packetsArray)\n");
			ret = -1; goto cleanup;
		}

		srand((unsigned int)time(NULL));
		for( i=0; i<NBDATA; i++ )
		{	// First packet filled with 0x1111..., second with 0x2222..., etc.
			data_size = PKTSZ - (rand()%10)+sizeof(LDPC_head);
			packetsArray[i] = (char*) malloc(data_size);
			if( packetsArray[i] == NULL ) {
				ret = -1; goto cleanup;
			}
			printf("data_size:%d\n", data_size);
			data_head.type_flag 		= 	0;
			data_head.group_id 			= 	j+1;
			data_head.total_data 		= 	(unsigned short)NBDATA;
			data_head.total_fec 		= 	(unsigned short)NBFEC;
			data_head.sequence_no 		= 	(unsigned int)i;
			data_head.current_length	= 	(unsigned short)data_size;
			data_head.longest_length 	= 	(unsigned short)(PKTSZ + sizeof(LDPC_head));
			memcpy(packetsArray[i], &data_head, sizeof(data_head));
			memset( packetsArray[i]+sizeof(data_head), (char)(i+((j+1)<<4)), data_size-sizeof(data_head));
			printf( "DATA[%03d]= \n", i );
			DumpBuffer( packetsArray[i], data_size );
		}

		// Now creating some FEC Packets...
		printf("\nCreating FEC Packets...\n");

		for( i=0; i < NBFEC; i++ )
		{
			packetsArray[i+NBDATA] = (char*)malloc(PKTSZ+sizeof(data_head));
			if( packetsArray[i+NBDATA] == NULL ) {
				ret = -1; 
				goto cleanup;
			}
			data_head.type_flag 		= 	1;
			data_head.group_id 			= 	j+1;
			data_head.total_data 		= 	(unsigned short)NBDATA;
			data_head.total_fec 		= 	(unsigned short)NBFEC;
			data_head.sequence_no 		= 	(unsigned int)(i+NBDATA);
			data_head.current_length	= 	0;
			data_head.longest_length 	= 	(unsigned short)(PKTSZ + sizeof(LDPC_head));
			memcpy(packetsArray[i+NBDATA], &data_head, sizeof(data_head));
			memset( packetsArray[i+NBDATA]+sizeof(data_head), 0, PKTSZ);
			BuildParitySymbol(&Session[j],  (void**)packetsArray, i, packetsArray[i+NBDATA] );
			//printf( "DATA[%03d]= \n", i+NBDATA );
			//DumpBuffer( packetsArray[i+NBDATA], PKTSZ+sizeof(data_head) );  // dump packet to screen
		}
		all_data[j]=packetsArray;
	}

	// Randomize packets order...
	printf("\nRandomizing transmit order...\n");
	randOrder1 = (int*)malloc(NBPKT*sizeof(int));
	randOrder2 = (int*)malloc(GROUP*sizeof(int));
	if( (randOrder1==NULL) ||(randOrder2 == NULL) ) {
		ret = -1; 
		goto cleanup;
	}
	randomizeArray( &randOrder1, NBPKT );
	randomizeArray( &randOrder2, GROUP );

	// ... and finally, throw our packets to space using UDP socket :-)
	struct sockaddr_in destHost;
	destHost.sin_family = AF_INET;
	destHost.sin_port = htons((short)DEST_PORT);
	destHost.sin_addr.s_addr = inet_addr(DEST_IP);
	mySock = initSocket();
	if( mySock == INVALID_SOCKET )
	{
		ret = -1; 
		goto cleanup;
	}
	printf( "Sending packets (DATA&FEC) to %s/%d\n", DEST_IP, DEST_PORT );
	for(j=0; j<GROUP; j++)
	{
		for( i=0; i < NBPKT; i++ )
		{
			if(i%NBDROP == 1)
				continue;

			memcpy(&data_head, all_data[randOrder2[j]][randOrder1[i]], sizeof(data_head));
			if(data_head.type_flag)
				data_size = data_head.longest_length;
			else
				data_size = data_head.current_length;

			//printf("Sending packet %-5d %-5d (%s)\n", randOrder2[j], randOrder1[i], randOrder1[i]<NBDATA ? "DATA" : "FEC");
			ret = sendto(mySock, all_data[randOrder2[j]][randOrder1[i]], data_size, 0, (struct sockaddr *)&destHost, sizeof(destHost));
			if (ret == SOCKET_ERROR) {
				printf( "main: Error! sendto() failed!\n" );
				ret = -1;
				break;
			}
			MSLEEP(10); // SLEEP avoid UDP flood (in milliseconds)
		}
	}
	if( i==NBPKT ) {
		printf( "\nComplete! %d packets sent successfully.\n", (i-NBPKT/NBDROP-1)*GROUP);
		ret = 1;
	}


cleanup:
	// Cleanup...
	if( mySock!= INVALID_SOCKET ) closesocket(mySock);
	for(j=0; j<GROUP; j++)
	{
		if( IsInitialized(&Session[j]) ) 
			EndSession(&Session[j] );
	}
	if( randOrder1 ) { free(randOrder1); }
	if( randOrder2 ) { free(randOrder2); }

	for(j=0; j<GROUP; j++)
	{
		packetsArray = all_data[j];
		if( packetsArray ) {
			for( i=0; i<NBPKT; i++ ) {
				free(packetsArray[i]);
			}
			free(packetsArray);
		}
	}

	// Bye bye! :-)
	return ret;
}



// Randomize an array of integers
void randomizeArray( int** array, int arrayLen )
{
	int backup=0,randInd=0;
	int	seed, i;	/* random seed for the srand() function */

	struct timeval	tv;
	if (gettimeofday(&tv, NULL) < 0) {
		perror("randomizeArray: gettimeofday() failed:");
		exit(-1);
	}
	seed = (int)tv.tv_usec;
	srand(seed);
	for( i=0; i<arrayLen; i++ )
		(*array)[i]=i;

	for( i=0; i<arrayLen; i++ )
	{
		backup = (*array)[i];
		randInd = rand()%arrayLen;
		(*array)[i] = (*array)[randInd];
		(*array)[randInd] = backup;
	}
}



/* Initialize Winsock engine and our UDP Socket */
SOCKET initSocket()
{
	SOCKET s;

	s = socket(AF_INET, SOCK_DGRAM, 0);
	if (s == INVALID_SOCKET)
	{
		printf("Error: call to socket() failed\n");
		return INVALID_SOCKET;
	}
	return s;
}



void DumpBuffer( char* buf, int len )
{
	char *ptr; int j = 0;

	//printf("0x");
	for (ptr = buf; len > 0; len--, ptr++) {
		/* convert to big endian format to be sure of byte order */
		printf( "%02X", *ptr);
		if (++j == 40)
		{
			j = 0;
			printf("\n");
		}
	}
	printf("\n");
}

