typedef int QUEUEITEM;
typedef struct {
	int front;
	int rear;
	int size;
	QUEUEITEM item[1];
} QUEUETYPE;
typedef QUEUETYPE *QUEUE;
