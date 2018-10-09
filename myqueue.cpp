#include  <iostream>
#include  <stdlib.h>
#include "myqueue.h"

/****************************************************************/

myqueue::myqueue()
{
   end=new q_element;
   end->next=NULL;
   end->pre=NULL;
   end->value=-1;
   start=end;
}

myqueue::~myqueue(){
  q_element *q=start;
  q_element *p;
  while (q!=NULL) {
    p=q;
    q=q->next;
    delete p;
  }
}

void myqueue::enqueue(int value)
{
   q_element *temp;

   temp=new q_element;
   temp->value=value;
   temp->next=start;
   start->pre=temp;
   temp->pre=NULL;
   start=temp;

}

int myqueue::dequeue(void)
{
   q_element *temp1;
   q_element *temp2;
   int value;

   //cout<<"in dequeue..."<<endl;
   if(is_empty()==1) return -1;
   else
   {
     //cout<<"temp1=end->pre"<<endl;
      temp1=end->pre;
      //cout<<"temp2=temp1->pre"<<endl;
      temp2=temp1->pre;
      //cout<<"value=temp1->value"<<endl;
      value=temp1->value;
      //cout<<"value to be dequeued: "<<value<<endl;
      //cout<<"temp2->next=end"<<endl;
      if(temp2!=NULL)
	temp2->next=end;
      //cout<<"end->pre=temp2"<<endl;
      end->pre=temp2;

      if(temp2==NULL) start=end;
  
      //cout<<"delete temp1"<<endl;
      //cout<<"temp1: "<<temp1<<endl;
      delete temp1;
      //cout<<"return value"<<endl;
      return value;
   }
}

int myqueue::is_empty(void)
{
   if(start==end) return 1;
   else return 0;
}

int myqueue::have_this(int value)
{
   q_element *current;

   current=start;
   while(current!=NULL)
   {
      if(current->value==value) return 1;
      current=current->next;
   }
   return 0;
}

void myqueue::clear(void) {
  q_element *q=start;
  q_element *p;
  while (q->next!=NULL) {
    p=q;
    q=q->next;
    delete p;
  }
  start=end;
}
