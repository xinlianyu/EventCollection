#pragma comment(lib,"ws2_32.lib")
#ifndef MYQUEUE
#define MYQUEUE

struct q_element
{
   int value;
   q_element *next;
   q_element *pre;
};

class myqueue
{
   public:
  q_element *start,*end; //end does not store any member of the queue. it's just a symbol of the ending
   public:
      myqueue();
      ~myqueue();
      void enqueue(int value);
      int dequeue(void);
	  int is_empty(void);
      int have_this(int value);
      void clear(void);
};

#endif