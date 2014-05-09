#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <mcheck.h>
#include "simple_coder.h"

/* Prototypes */
SOCKET initSocket( );
void DumpBuffer( char*, int );

int main(int argc, char* argv[])
{
	int i, data_size;
	char is_new=0;
	LDPC_group_list *group_list=NULL, *group_head=NULL;

	// Received (and rebuilt) packets (DATA and FEC) are stored in a 
	// packets array where each packet is an array of bytes (char).

	SOCKET	mySock	= INVALID_SOCKET;
	//char*	recvPkt	= NULL;
	char*	buff	= NULL;
	int	ret	= -1;
	int	decodeSteps = 0;
	int 	total = 0;
	LDPC_head data_head;

	struct timeval timeout;
	fd_set readfds;
	timeout.tv_sec = 5;
	timeout.tv_usec = 0;

	mtrace();

	// Initialize our UDP socket
	mySock = initSocket();
	if( mySock == INVALID_SOCKET ) {
		printf("Error initializing socket\n");
		ret = -1;
		goto cleanup;
	}

	buff	= (char*) malloc(PKTSZ+sizeof(LDPC_head));
	if ( buff==NULL) {
		printf("Error: insufficient memory (calloc failed in main())\n");
		ret = -1; goto cleanup;
	}

	//printf( "Decoding in progress...\nWaiting for new packets...\n" );
	while(1)
	{
		decodeSteps++;
		FD_ZERO(&readfds);
		FD_SET(mySock, &readfds);

		if(select(mySock+1,&readfds,NULL, NULL, &timeout) > 0)
		{
			ret = recvfrom( mySock, buff, PKTSZ+sizeof(data_head), 0, NULL, NULL );
			if(ret < 0)
			{
				ret = -1;
				goto cleanup;	
			}
			// OK, new packet received...
			memcpy(&data_head, buff, sizeof(data_head));
			//printf("------------------------------------\n");
			//printf("--- Step %d : new packet received: %02d, size:%d %d, group_id:%d, buffer:0x%x\n", decodeSteps, data_head.sequence_no, data_head.current_length, data_head.longest_length, data_head.group_id, buff);
			if(NULL == group_head)
			{
				group_head = group_list_init(data_head.group_id, data_head.total_data + data_head.total_fec);
				if(NULL == group_head)
				{
					ret = -1;
					goto cleanup;
				}

				group_list = group_head;

				// Initialize the LDPC session
				if(InitSession(group_list->Session, NBDATA, NBFEC, data_head.longest_length, FLAG_DECODER, SEED, SESSION_TYPE, LEFT_DEGREE ) == LDPC_ERROR)
				{
					printf("Error: Unable to initialize LDPC Session\n");
					ret = -1; goto cleanup;
				}

				DecodingWithSymbol(group_list->Session, (void**)(group_list->packet), buff, data_head.sequence_no, true);
			}
			else
			{
				group_list = group_list_search(group_head, &is_new, data_head.group_id, data_head.total_data + data_head.total_fec);
				if(NULL == group_list)
				{
					ret = -1;
					goto cleanup;
				}

				if(is_new)
				{
					// Initialize the LDPC session
					if(InitSession(group_list->Session, NBDATA, NBFEC, data_head.longest_length, FLAG_DECODER, SEED, SESSION_TYPE, LEFT_DEGREE ) == LDPC_ERROR)
					{
						printf("Error: Unable to initialize LDPC Session\n");
						ret = -1; goto cleanup;
					}
				}
				DecodingWithSymbol(group_list->Session, (void**)(group_list->packet), buff, data_head.sequence_no, true);
			}
			if(IsDecodingComplete(group_list->Session, (void**)(group_list->packet)))
			{
				printf("group:%d\n", group_list->group_id);
				for(i=0;i<group_list->total_pkt;i++)
				{
					if(group_list->packet[i] == NULL)
						continue;
					memcpy(&data_head, group_list->packet[i], sizeof(data_head));
					data_size = data_head.current_length;
					printf("group_id:\t\t%d\ntotal_data:\t\t%d\ntotal_fec:\t\t%d\nsequence_no:\t\t%d\ncurrent_length:\t%d\nlongest_length:\t%d\n",
							data_head.group_id,
							data_head.total_data,
							data_head.total_fec,
							data_head.sequence_no,
							data_head.current_length,
							data_head.longest_length);
					printf("size=%d, DATA[%d]= \n",data_size, i);
					DumpBuffer(group_list->packet[i], data_size);
					total++;
				}
				if( IsInitialized(group_list->Session) ) 
					EndSession(group_list->Session);
				group_list_delete(&group_head, group_list->group_id);
			}
		}
		else
			break;
	}

	printf("===============================================================================================\n");
	printf("=======Above is wonderful packet decoder, next is decoder some packet and rest=================\n");
	printf("===============================================================================================\n");

	group_list = group_head;
	total=0;
	while(group_list != NULL)
	{
		printf("group:%d\n", group_list->group_id);
		for(i=0;i<group_list->total_pkt;i++)
		{
			if(group_list->packet[i] == NULL)
				continue;
			memcpy(&data_head, group_list->packet[i], sizeof(data_head));
			data_size = data_head.current_length;
			printf("group_id:\t\t%d\ntotal_data:\t\t%d\ntotal_fec:\t\t%d\nsequence_no:\t\t%d\ncurrent_length:\t%d\nlongest_length:\t%d\n",
					data_head.group_id,
					data_head.total_data,
					data_head.total_fec,
					data_head.sequence_no,
					data_head.current_length,
					data_head.longest_length);
			printf("size=%d, DATA[%d]= \n",data_size, i);
			DumpBuffer(group_list->packet[i], data_size);
			total++;
		}
		group_list = group_list->next;
	}

	printf("%d packets rebuilt\n", total);
	printf("Done! All DATA packets rebuilt in %d decoding steps [%d-%d]\n", decodeSteps, NBDATA, NBPKT);

	// Cleanup...
cleanup:
	if( mySock!= INVALID_SOCKET ) closesocket(mySock);
	group_list = group_head;
	while(group_list != NULL){
		if( IsInitialized(group_list->Session) ) 
			EndSession(group_list->Session);
		group_list = group_list->next;
	}
	group_list_deinit(group_head);

	if(buff) free(buff);

	// Bye bye! :-)
	return ret;
}


/* Initialize Winsock engine and our UDP Socket */
SOCKET initSocket()
{
	SOCKET s = INVALID_SOCKET;
	int err  = SOCKET_ERROR;

	s = socket(AF_INET, SOCK_DGRAM, 0);
	if (s == INVALID_SOCKET)
	{
		printf("Error: call to socket() failed\n");
		return INVALID_SOCKET;
	}

	if (-1 == fcntl(s, F_SETFL, O_NONBLOCK))
	{
		printf("fcntl socket error!\n");
	}

	struct sockaddr_in bindAddr;
	bindAddr.sin_family = AF_INET;
	bindAddr.sin_port = htons((short)DEST_PORT);
	bindAddr.sin_addr.s_addr = INADDR_ANY;

	err = bind( s, (struct sockaddr *)&bindAddr, sizeof(bindAddr));
	if( err == SOCKET_ERROR)
	{
		printf("initSocket: bind() failed. Port %d may be already in use\n", DEST_PORT);
		return INVALID_SOCKET;
	}
	return s;
}

void DumpBuffer( char* buf, int len )
{
	int j = 0;
	char *ptr;

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

