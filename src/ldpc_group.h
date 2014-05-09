#ifndef LDPC_GROUP_H
#define LDPC_GROUP_H

#include "ldpc_fec.h"

typedef struct group_list {
	struct group_list *next;
	unsigned int group_id;
	unsigned int total_pkt;
	char** 	packet;
	LDPCFecSession *Session;
}LDPC_group_list;

LDPC_group_list* group_list_init(unsigned int group_id, unsigned int total_pkt);

LDPC_group_list* group_list_search(LDPC_group_list *head, char *is_new, unsigned int group_id, unsigned int total_pkt);

void group_list_delete(LDPC_group_list **head, unsigned int group_id);
#endif
