struct elem{
	int node_num;
	double val;
	struct elem * next; 
	struct elem * prev;
};

typedef struct elem node;

class linear_tree{
	public:
//		node *tree_head=NULL;
		node *tree_head;
		node *node_pointer[100000];
		double node_val[100000];
	// CONSTRUCT PUBLIC FUNCTIONS TO MODIFY THE TREE
	public:
		void change_parent(int node_num, node *new_parent)
		{
			node *temp=node_pointer[node_num];
			node_val[node_num]=((new_parent->next->val)-(new_parent->val))/2;
			temp->val=(new_parent->val)+node_val[node_num];
			temp->prev->next=temp->next;
			temp->next->prev=temp->prev;
			temp->next=new_parent->next;
			temp->prev=new_parent;
			new_parent->next=temp;
			temp->next->prev=temp;
		}
		void linearize()
		{
			
		}
};

