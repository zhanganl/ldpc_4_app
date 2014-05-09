#include <stdio.h>
#include <stdlib.h>
#include "ldpc_group.h"

LDPC_group_list* group_list_init(unsigned int group_id, unsigned int total_pkt)
{
	LDPC_group_list *head;
	int i;

	head = (LDPC_group_list *)malloc(sizeof(LDPC_group_list));
	if(NULL != head)
	{
		head->next = NULL;
		head->group_id = group_id;
		head->total_pkt = total_pkt;
		head->packet = (char **)malloc(total_pkt*sizeof(char*));
		if(NULL == head->packet)
		{
			printf("[%s:%d] malloc err!\n", __FILE__, __LINE__);
			free(head);
			head = NULL;
		}
		head->Session = (LDPCFecSession *)malloc(sizeof(LDPCFecSession));
		if(NULL == head->Session)
		{
			printf("[%s:%d] malloc err!\n", __FILE__, __LINE__);
			free(head->packet);
			free(head);
			head = NULL;
		}
		memset(head->packet, 0, total_pkt*sizeof(char *));
		//for(i=0; i<total_pkt; i++)
		//	head->packet[i] = NULL;
		memset(head->Session, 0, sizeof(LDPCFecSession));
	}
	return head;
}

LDPC_group_list* group_list_search(LDPC_group_list *head, char *is_new, unsigned int group_id, unsigned int total_pkt)
{
	LDPC_group_list *P, *N=head;
	int i;

	*is_new = 0;

	while(N->next != NULL)
	{
		if(N->group_id != group_id)
			N = N->next;
		else
			return N;
	}
	if(N->group_id == group_id)
		return N;

	*is_new = 1;

	P = (LDPC_group_list *)malloc(sizeof(LDPC_group_list));   //申请头结点空间  
	if(NULL != P)
	{
		P->next = NULL;
		P->group_id = group_id;
		P->total_pkt = total_pkt;
		P->packet = (char **)malloc(total_pkt*sizeof(char*));
		if(NULL == P->packet)
		{
			printf("[%s:%d] malloc err!\n", __FILE__, __LINE__);
			free(P);
			P = NULL;
		}
		P->Session = (LDPCFecSession *)malloc(sizeof(LDPCFecSession));
		if(NULL == P->Session)
		{
			printf("[%s:%d] malloc err!\n", __FILE__, __LINE__);
			free(P->packet);
			free(P);
			P = NULL;
		}
		//for(i=0; i<total_pkt; i++)
		//	P->packet[i] = NULL;
		memset(P->packet, 0, total_pkt*sizeof(char *));
		memset(P->Session, 0, sizeof(LDPCFecSession));
	}
	else
		printf("[%s:%d] malloc err!\n", __FILE__, __LINE__);
	N->next = P;
	return P;
}

////////////////////////////////////////////   
////单链表的删除，在链表中删除值为x的元素  
void group_list_delete(LDPC_group_list **head, unsigned int group_id)  
{  
	LDPC_group_list *p,*pre, *tmp;                   //pre为前驱结点，p为查找的结点。   
	int i;
	
	p = *head;  
	while(NULL != p)              //查找值为x的元素   
	{    
		if(p->group_id == group_id)
		{
			break;
		}
		pre = p;   
		p = p->next;  
	}

	if(p == NULL)
	{
		printf("[%s:%d] can't find the group id!\n", __FILE__, __LINE__);
		return;
	}
	if(p==*head)
		*head = p->next;
	else
		pre->next = p->next;          //删除操作，将其前驱next指向其后继。  
	for(i=0; i<p->total_pkt; i++)
	{
		if(p->packet[i] != NULL)
		{
			free(p->packet[i]);
			p->packet[i] = NULL;
		}
	}
	if(IsInitialized(p->Session))
		EndSession(p->Session);
	free(p->packet);
	free(p->Session);
	free(p);
	p=NULL;
	
	tmp = *head;  
	while(NULL != tmp)              //查找值为x的元素   
	{    
		tmp = tmp->next;  
	}

	return;
}

//didn't think about free p->packet[i] like group_list_delete
void group_list_deinit(LDPC_group_list *head)  
{  
	LDPC_group_list *p,*pre;                   //pre为前驱结点，p为查找的结点。   
	p = head;  
	while(p!=NULL)              //查找值为x的元素   
	{     
		if( IsInitialized(p->Session) ) 
			EndSession(p->Session);
		pre = p;   
		p = p->next;  
		free(pre->packet);
		free(pre->Session);
		free(pre);
		pre=NULL;
	}
	
	return;
}
