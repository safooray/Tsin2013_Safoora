#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <list>
//#include <ext/hash_map>
#include <map>
#include <iostream>
using namespace std;


#define max_length 10
/*Since we used "int" variables to store vertices, the lergest positive number 
can be represented by a 32bit "int" data type is 2147483647. 
the number 2147483647 has 10 digits.
*/

//Boundary Vertices
struct compBndryInfoStruct
{
	int** compBndryVrtxDists;
	vector<int> compBndrys;
};
compBndryInfoStruct* compBndryInfo;

//Boundary Vertices
struct adjacent_with_v_in_G { 
   /*is used to define adjacency lists*/
   int w;
   /*to represent an edge (v,w) in the list of vertex v;
     note that this edge is repeated in the list of vertex w as well.
   */			
   char status;
   /*to determine if the edge (v,w) in the list of vertex v, is an incoming backedge of v (Y).
   */
//DeleteCutPairs  
	 char deleted;
     //to determine if the edge (v,w) in the list of vertex v, is Deleted (Y).
//DeleteCutPairs	 

   struct adjacent_with_v_in_G* more;
};
typedef struct adjacent_with_v_in_G* adjacentG;

adjacentG* LG;
//DeleteCutPairs
adjacentG* outgoing_tree_edge;
   //in order to have a list of pointers. since we have |V| outgoing tree edges
   //each pointer of this list refers to one of those outgoing tree edges. 
//DeleteCutPairs

/*//OnlyFindCutPairs
char* outgoing_tree_edge;
*///OnlyFindCutPairs

adjacentG edge, edge2;

//*****************************Added for simple 3edge***************************
int Vnum;
int comIndex = 1;
int edgeNum = 0; /*initilizing the number of edges in G*/
int dfs = 1;
int compNum = 0;
int bridgeNum = 0;
bool* isBoundry;
int orignCC = 0;
int cutPairNum = 0;
char *visited;
char *v_to_low_bedge;
int *pre;
int *nd;
int *lowpt;
int *second_lowpt;
int *low;
int *second_low;
int *to_low;
//3edge
int *newlink;
int *components;
//3edge
//DeleteCutPairs
int *gen_backedge;
//DeleteCutPairs
//------------ONE STACK (ELEMENTS)----------------------------------------------
int wtop = 0;
int wbot = 0;
char *Sbedge; //is (x,y) a back edge or not?
int *Sx;
int *Sy;
int *Sp;
int *Sq;
//DeleteCutPairs
adjacentG* Sgen_pointer;
adjacentG* Sgen_other_pointer;
//DeleteCutPairs
//******************************************************************************

void abrt(char ch[50]) {
   printf("\n%s\n", ch);
   exit(1);
}

void freeMem(int v) {
   int r;
   for (r=1;r<v;r++) {
        while (LG[r] != NULL) {
		     edge = LG[r]->more;
			 free(LG[r]);
			 LG[r] = edge;
		}
   }
   
   free(LG);
   free(nd);
   free(Sx);
   free(Sy);
   free(Sp);
   free(Sq);
   free(components);
   free(Sbedge);
   free(pre);
   free(low);
   free(lowpt);
   free(to_low);  
   free(second_lowpt); 
   free(second_low);
   free(visited);
   free(isBoundry);
   free(v_to_low_bedge);
   free(outgoing_tree_edge);
   
//3edge
   free(newlink);
//3edge

//DeleteCutPairs   
   free(gen_backedge);
   free(Sgen_pointer);
   free(Sgen_other_pointer);
//DeleteCutPairs
}

//3edge
void DFS(int v) {
   int w;
   adjacentG edge;
//PRINT-3edge
   printf("%d,",v);
   components[v] = comIndex;
//PRINT-3edge
   visited[v] = 'Y';
   edge = LG[v];
   while (edge != NULL) {
      w = edge->w;
      if ((edge->deleted == 'N') && (visited[w] == 'N'))
         DFS(w);
      edge = edge->more;
   }
}
//3edge



void Orig_DFS(int v, bool* visited) 
{
   int w;
   adjacentG edge;
   visited[v] = true;
   edge = LG[v];
   while (edge != NULL) {
      w = edge->w;
      if ((edge->deleted == 'N') && (visited[w] == false))
		  Orig_DFS(w, visited);
      edge = edge->more;
   }
}

int nConComp(adjacentG* Graph)
{
	bool* visited = new bool[Vnum];
	int nCC = 0;
	for(int i = 0; i < Vnum; i++)
	{
		visited[i] = false;
	}

	for(int i = 0; i < Vnum; i++)
	{
		if(!visited[i])
		{
			nCC++;
			Orig_DFS(i, visited);
		}
	}
	delete[] visited;
	return nCC++;
}

int* BFS(adjacentG* graph, int s, bool considerDeleted)
{
	int* dists = new int[Vnum];
	adjacent_with_v_in_G* edge;
	for (int i = 1; i < Vnum; i++)
	{
		dists[i] = -1;
	}
	std::list<int> q;
	q.push_back(s);
	dists[s] = 0;
	while (q.size() != 0)
	{
		s = q.front();
		q.pop_front();
		int dist = dists[s];
		edge = graph[s];
		while (edge)
		{
			if(considerDeleted && edge->deleted == 'Y')
			{
				edge = edge->more;
				continue;
			}
			int curChild = edge->w;
			if (dists[curChild] < 0)
			{
				dists[curChild] = dist + 1;
				q.push_back(curChild);
			}
			edge = edge->more;
		}
	}
	return dists;
}

int* farness(adjacentG* graph, bool considerDeleted)
{
	int* fars = new int[Vnum];
	for(int i = 1; i < Vnum; i++)
	{
		int* dists = BFS(graph, i, considerDeleted);
		int far = 0;
		for(int j = 1; j < Vnum; j++)
		{
			if(dists[j] > 0)
			{
				far += dists[j];
			}
		}
		fars[i] = far;
		delete[] dists;
	}
	return fars;
}

void calcCompBndryInfo(int nComp)
{
	//if (( compBndryInfo = (compBndryInfoStruct*)malloc(nComp * sizeof(compBndryInfoStruct)) ) == NULL)
	//	abrt("Not enough memory to allocate buffer");
	compBndryInfo = new compBndryInfoStruct[nComp];
	//for (int r)
	for (int r = 1; r < Vnum; r++)
	{
		if(isBoundry[r])
			compBndryInfo[components[r]].compBndrys.push_back(r);
	}
	for (int i = 1; i <nComp; i++)
	{
		compBndryInfo[i].compBndryVrtxDists = new int*[compBndryInfo[i].compBndrys.size()];
		for(unsigned int j = 0; j < compBndryInfo[i].compBndrys.size(); j++)
		{
			compBndryInfo[i].compBndryVrtxDists[j] = BFS(LG, compBndryInfo[i].compBndrys[j], false);
		}
	}
}
void printArray(int* arr, int n)
{
	for(int i = 1; i < n; i++)
	{
		//printf("%s\n", arr[i] ? "B" : "NB");
		printf("%d\n", arr[i]);
	}
	printf("\n");
}
void printArrayOArrays(int** arr, int n, int m)
{
	for(int i = 0; i < n; i++)
	{
		for (int j = 1; j < m; j++)
		{
			//printf("%s\n", arr[i] ? "B" : "NB");
			printf("%d ", arr[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void readLinkFromat(FILE* in)
{
	typedef pair<int, int> Int_Pair;
	map<int, int> hash;
	map<Int_Pair, int> edges;
	const int READ_COMMENT = 0, READ_W = 1, READ_V = 2;
	int w,v,indx = 0,state = READ_W, cVnum = 1;
	char nst[20];
	char ch;
	while ( (ch = fgetc(in)) != EOF) 
	{
		switch(ch)
		{
		case '#':
			state = READ_COMMENT;
			break;
		case '\t':
		case 32:
			if(state == READ_W)
			{
				nst[indx] = 0;
				indx = 0;
				w = atoi(nst);
				state = READ_V;
			}
			else if(state != READ_COMMENT)
				abrt("Problem reading file: state != READ_COMMENT");
			break;
		case '\n':
			if(state == READ_V)
			{
				nst[indx] = 0;
				indx = 0;
				v = atoi(nst);
				
				
				
				if( hash.count(w) )
					w = hash[w];
				else
				{
					hash.insert(Int_Pair (w, cVnum));
					w = cVnum++;
				}				

				if( hash.count(v) )
					v = hash[v];
				else
				{
					hash.insert(Int_Pair (v, cVnum));
					v = cVnum++;
				}

				
				//int addEdges = 0;
				if( edges.count(Int_Pair (w,v)) == 0 )
				{
					//add link (w, v) to the graph here!
					if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_v_in_G)) ) == NULL)
					   abrt("Not enough memory to allocate buffer");
					edge->more = NULL;
					edge->w = v;
					edge->status = 'N';
					//DeleteCutPairs			  
					edge->deleted = 'N';						
					//DeleteCutPairs
					edge->more = LG[w];
					LG[w] = edge;
					edgeNum++;
					

					// add (v,w)
					if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_v_in_G)) ) == NULL)
					   abrt("Not enough memory to allocate buffer");
					edge->more = NULL;
					edge->w = w;
					edge->status = 'N';
					//DeleteCutPairs			  
					edge->deleted = 'N';						
					//DeleteCutPairs
					edge->more = LG[v];
					LG[v] = edge;
					
					edges[Int_Pair(w,v)] = 0;
					edges[Int_Pair(v,w)] = 0;
					edgeNum = edgeNum + 1;
				}

				state = READ_W;
			}
			else if(state == READ_COMMENT)
				state = READ_V;
			else
				abrt("Reading file Problem: state == READ_V");
			break;
		default:
			if((state == READ_V || state == READ_W) && ch >= '0' && ch <= '9')
				nst[indx++] = ch;	
			break;
		}	
	}
	edgeNum = edgeNum / 2;	 
	fclose(in);
	printf("Reading file finished!\n");
}
void find_cut_pair(int v,int parent) {

   adjacentG edge;
//DeleteCutPairs
   adjacentG v_to_low_pointer;
   adjacentG low_to_v_pointer;
//DeleteCutPairs
   int w,u,p;
   int vtop;
   int vbot;
   
//3edge
     int y;
     char bedge;
//3edge   
   
//DeleteCutPairs
     outgoing_tree_edge[v] = NULL;
//DeleteCutPairs

/*//OnlyFindCutPairs
   outgoing_tree_edge[v] = 'N';
*///OnlyFindCutPairs
   
   vtop = wtop; //intialize top indicator for stack v (wtop & wbot are global)
   vbot = wtop; //intialize bottom indicator for stack v (vtop & vbot are local)
   visited[v] = 'Y';
   pre[v] = dfs;
   dfs = dfs + 1;
   nd[v] = 1;
   lowpt[v] = pre[v];
   low[v] = v;
   second_lowpt[v] = pre[v];
   second_low[v] = v;
	 
// Step 1.
   edge = LG[v];
   while (edge != NULL) { /*Do the followings for every edge e=(v,w) which is represented in the adjacency list of v*/
      w = edge->w;
// Step 1.1
      if (visited[w]=='N') {
         find_cut_pair(w,v);

// Step 1.1.1
      if (lowpt[w] == pre[w]) { //the edge (v,w) is a bridge, so we delete it.
                                //note that at this time stack[w] is empty. because
                                //there is no vertex from which w recieves lowpt and
                                //the last possible cut-pair on stack[w] was (w,u),(x,y)
                                //where u is the child of w and x,y are descendants of w.
/*//YesOrNo
   printf("It's a NO instance!\n");
   freeMem(Vnum); 
         exit(1);
*///YesOrNo

//PRINT-bridge
            printf("\n(%d,%d) is a bridge!",v,w);
			bridgeNum = bridgeNum + 1;
//PRINT-bridge

         compNum = compNum + 1;
		 
//DeleteCutPairs
         edge->deleted = 'Y';
         //mark the edge (v,w) as deleted in the list of v.
         outgoing_tree_edge[w]->deleted = 'Y';
         //mark the edge (w,v) as deleted in the list of w.
//DeleteCutPairs

      }
      else
         if ( wtop!=wbot && w==Sq[wtop] ) { //if stack[w] is not empty, [(x,y),p~>q] denotes the top of the stack.
/*//YesOrNo	 
   printf("It's a NO instance!\n");
      freeMem(Vnum);
	 exit(1);
*///YesOrNo	 
	 				 
//PRINT-CutPair			
            printf("\n(%d,%d) and (%d,%d) is a cut pair!",Sx[wtop],Sy[wtop],v,w);
			isBoundry[Sx[wtop]] = true;
			isBoundry[Sy[wtop]] = true;
			isBoundry[v] = true;
			isBoundry[w] = true;

//PRINT-CutPair

            compNum = compNum + 1;
			cutPairNum = cutPairNum + 1;	
			
//DeleteCutPairs			
            edge->deleted = 'Y';
            //mark the edge (v,w) as deleted in the list of v.

            outgoing_tree_edge[w]->deleted = 'Y';
            //mark the edge (w,v) as deleted in the list of w.

            Sgen_pointer[wtop]->deleted = 'Y';
            //mark the edge (x,y) as deleted in the list of x.

            if (Sbedge[wtop]=='Y')
               gen_backedge[Sx[wtop]]=Sy[wtop];
            else
               Sgen_other_pointer[wtop]->deleted = 'Y';
               //mark the edge (y,x) as deleted in the list of y.
//DeleteCutPairs

//3edge
            bedge=Sbedge[wtop];
            y=Sy[wtop];
//3edge			

            p=Sp[wtop];
            wtop--;
            /*pop STACK[w];
            */
   	        if ( v != p ) {
               wtop++;
               Sq[wtop]=v;
               /*push [(x,y),p~>v] onto STACK[w]
               */
            }
//3edge			
		   		  if (bedge=='N') {
              newlink[y] = v;
            }
//3edge
              /*newlink[y] is continuously updated and finally 
			  denotes to the right additional edge (y,v). Later, after the DFS 
			  is finished, the list is scanned in order to add correct
			  additional edges.
			  It is not correct that we change this 'if' to an else if' (i.e. v==p)
			  inorder to add additinal edges right away at this stage. Because 
			  sometimes in section 3 we encounter an incoming back edge of v 
			  while v is located in <p,q> path. So we will miss the actual p (=v)
			  since	we pop the candidate in section 3. Hence we are not able to
			  add the (y,v) edge.   
              */

         }
         nd[v] = nd[v] + nd[w];
// Step 1.1.2
         if (lowpt[w] <= lowpt[v]) {
            second_lowpt[v] = lowpt[v];
            lowpt[v] = lowpt[w];
            second_low[v] = low[v];
            low[v] = low[w];

            vtop = wtop;
            vbot = wbot;
            /*STACK[v] = STACK[w];
            */
            to_low[v] = w;
            v_to_low_bedge[v] = 'N';
            /*(v,to_low[v]) is not a back edge, so v_to_low_bedge is set to 'N'.
            */
//DeleteCutPairs
            v_to_low_pointer = edge;
            //store the address of the edge(v,w)(in the list of v) in v_to_low_pointer.

            low_to_v_pointer = outgoing_tree_edge[w];
            //outgoing_tree_edge[w] is the pointer to the edge(w,v)(in the list of w), 
            //such that v is the parent of w.

//DeleteCutPairs			
         }
// Step 1.1.3
         else {
            if (lowpt[w] < second_lowpt[v]) {
               second_lowpt[v] = lowpt[w];
               second_low[v] = low[w];
            }
            wtop = wbot;
            /*empty(STACK[w]);
            */
         }
      }
// Step 1.2
      else {
//DeleteCutPairs	  
         if (w == parent && !outgoing_tree_edge[v]) {
            outgoing_tree_edge[v] = edge;         
//DeleteCutPairs		 
	
/*//OnlyFindCutPairs
		 if (w == parent && outgoing_tree_edge[v] == 'N') {
            outgoing_tree_edge[v] = 'Y';
*///OnlyFindCutPairs			
            /*At the beginning of this iteration outgoing_tree_edge has been set
              to NULL or 'N'. So once we encounter to an (v,w) edge, such that w is the
              parent of v, we consider this edge as the incoming tree edge of w
              or outgoing_tree_edge of v. we store the address of this edge in
              the outgoing_tree_edge[v] because after we come back to the parent, it
              might turn out that this edge is a cut edge so we also want to mark it as
              deleted in the list of current v. Therefore, at that time this edge can
              be marked as deleted in the list of current v(it would be w after we
              come back to the parent; so we will be able to set outgoing_tree_edge[w]->deleted='Y'

              From now, no more tree edge connecting current v and its parent w must be encountered! 
			  We then consider other (v,w) edges(if any exists) connecting current v and its parent w
              , as backedges (parallel edges).
            */
         }
         else if (pre[v] > pre[w]) { //the edge (v,w) is an outgoing backedge of v
            if (pre[w] <= lowpt[v]) {
               second_lowpt[v] = lowpt[v];
               lowpt[v] = pre[w];
               second_low[v] = low[v];
               low[v] = w;
               vtop = vbot;
               /*empty (STACK[v]);
               */
               to_low[v] = w;
               v_to_low_bedge[v] = 'Y';
               /*(v,to_low[v]) is a backedge, so set v_to_low_bedge to 'Y'.
               */
//DeleteCutPairs		   
               v_to_low_pointer = edge;
               low_to_v_pointer = NULL;
//DeleteCutPairs
               /*because we have not yet visited the incoming backedge (w,v); so
                 we don't know the address of that edge in the list of w. Later when we encounter 
				 the incoming backedge in the list of w, we use gen_backedge[v] in order to determine 
				 if it is a cut-edge and must be marked as deleted.
               */
            }
            else if (pre[w] < second_lowpt[v]){
               second_lowpt[v] = pre[w];
               second_low[v] = w;
            }
         }
         else {
            edge->status = 'B';
            /*mark edge (v,w) as an incoming backedge in the list of vertex v.
            */
//DeleteCutPairs			
            if (gen_backedge[w]==v)
               edge->deleted = 'Y';
//DeleteCutPairs
         }
      }
      edge = edge->more;
   }
// Step 2.
// Step 2.1
   if (vtop == vbot) {
   /*if STACK[v] is empty
   */
// Step 2.2.1
      if (second_lowpt[v] > lowpt[v]) {
//         push((v, to_low[v]), low[v]~>second_low[v]) onto STACK[v]
         vtop++;
         Sbedge[vtop]=v_to_low_bedge[v];
//DeleteCutPairs		 
         Sgen_pointer[vtop]=v_to_low_pointer;
         //This is a pointer to the potential generator (v,to_low[v]) in
         //the list of v, but not to (to_low[v],v) in the list of to_low[v].

         Sgen_other_pointer[vtop]=low_to_v_pointer;
         //This is a pointer to the potential generator (to_low[v],v). If
         //it is NULL, that means we don't know its address yet because this 
		 //is a backedge generator and we have not yet visited the incoming 
		 //backedge (to_low[v],v) in the list of to_low[v].
		 //We later take care of this situation by using gen_backedge[v].
//DeleteCutPairs		 
         Sx[vtop]=v;
         Sy[vtop]=to_low[v];
         Sp[vtop]=low[v];
         Sq[vtop]=second_low[v];
      }
   }
// Step 2.2
   else {
// Step 2.2.1
      if (second_lowpt[v] > pre[Sq[vtop]]) {
			   vtop++;
       	 Sbedge[vtop]=v_to_low_bedge[v];
//DeleteCutPairs	 
         Sgen_pointer[vtop]=v_to_low_pointer;
         Sgen_other_pointer[vtop]=low_to_v_pointer;
//DeleteCutPairs
				 Sx[vtop]=v;
				 Sy[vtop]=to_low[v];
				 Sp[vtop]=Sq[vtop-1];
				 Sq[vtop]=second_low[v];
         /*push((v, to_low[v]), q~>second_low[v]) onto STACK[v]
         */
   }
// Step 2.2.2
		else {
      	while (vtop != vbot && second_lowpt[v] <= pre[Sp[vtop]])  {
				vtop--;
            /*pop STACK[v]
            */
         }
         if (vtop != vbot && second_lowpt[v] < pre[Sq[vtop]]) {
				Sq[vtop]=second_low[v];
            /*pop STACK[v]
      	   push ((x,y), p~>second_low[v]) onto STACK[v]
            */

         }
      }
	}
// Step 3.
   edge = LG[v];
   while ((vtop != vbot) && (edge != NULL)) {
   /*Do the followings for every edge e=(v,w) in the list of v
   */
      u = edge->w;
      if ( (edge->status == 'B') && (pre[v] < pre[u]) ) {
      	while ( (vtop != vbot) && (Sbedge[vtop]=='N') && (pre[u]>=pre[Sy[vtop]]) && pre[u]<=(pre[Sy[vtop]]+nd[Sy[vtop]]-1) ) {
//            printf("\n(%d,%d), (%d,%d)-(%d,%d) is popped up!",v,u,Sx[vtop],Sy[vtop],Sp[vtop],Sq[vtop]);//PRINT
//         	pop STACK[v]
         	vtop--;
         }


      }
      edge = edge->more;
   }
   wtop = vtop;
   wbot = vbot;
}
//******************************************************************************

int main(int argc, char **argv) 
{

   if (argc < 2) 
   {
      printf("No input!\n");
      abort();
   }
   int r;	 /*determines the starting vertex for DFS*/
   FILE *in;//, *out;
   const char* in_filename = argv[1];
   int n, ch, v, indx;
   int next_list = 1;
   char ch2[max_length];

   printf("\nFinding cut-pairs using\nWOR algorithm...\n");	 	 
	 	 
   if ((in = fopen(in_filename, "rt")) == NULL)
      abrt("Cannot open input file.");
	 indx = 0;		
   while ( (ch = fgetc(in)) != 10) { //reading the first number which is the # of vertices; the maximum allowed number of digits is max_length
      ch2[indx] = (char)ch;
      indx = indx + 1;
   }
   ch2[indx] = '\0';	 
   Vnum = atoi(ch2) + 1;
    
//*********************************Memory allocation	 						
   if (( isBoundry = (bool*)malloc(Vnum * sizeof(bool)) ) == NULL)
      abrt("Not enough memory to allocate buffer");
   for (r = 1; r <Vnum; r++)
   {
	   isBoundry[r] = false;
   }

   if (( pre = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( nd = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( lowpt = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( second_lowpt = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( low = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( second_low = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( to_low = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( components = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

//DeleteCutPairs
   if (( gen_backedge = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( outgoing_tree_edge = (adjacentG*)malloc(Vnum * sizeof(struct adjacent_with_v_in_G)) ) == NULL)
      abrt("Not enough memory to allocate buffer");
//DeleteCutPairs

//3edge			
   if (( newlink = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");
//3edge

/*//OnlyFindCutPairs
   if (( outgoing_tree_edge = (char*)malloc(Vnum * sizeof(char)) ) == NULL)
	   abrt("Not enough memory to allocate buffer");
*///OnlyFindCutPairs
			
   if (( LG = (adjacentG*)malloc(Vnum * sizeof(struct adjacent_with_v_in_G)) ) == NULL)
      abrt("Not enough memory to allocate buffer");
			
   if (( v_to_low_bedge = (char*)malloc(Vnum * sizeof(char)) ) == NULL)
      abrt("Not enough memory to allocate buffer");			
			
   if (( visited = (char*)malloc(Vnum * sizeof(char)) ) == NULL)
      abrt("Not enough memory to allocate buffer");							 			

   for (indx = 0; indx < Vnum; indx++) {
      LG[indx] = NULL;
      visited[indx] = 'N';
//3edge	  
      newlink[indx] = 0;
//3edge
   }	 
//	 indx=0;
	readLinkFromat(in);
	//printf("# of original Graph Connected Componentn: %d\n", nConComp(LG));
	orignCC = nConComp(LG);
//   while ( (ch = fgetc(in)) != EOF) {
//      if (indx == 10 && (ch != 62 || ch != 10)) 
//			/*10 is a carriage return and 62 is a '>', after we reach to the maximum
//			length for a number we expect new line or '>'. 
//			*/
//   		   abrt("Your input file has an error!");
//      else if (ch == 62) {
//         ch2[indx] = '\0';
//         indx = 0;
//			   n = atoi(ch2);
//         if (next_list) {
//				      v = n;
//							/*determine the next list
//							*/
//            next_list = 0;
//         }
//         else {
//            if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_v_in_G)) ) == NULL)
//               abrt("Not enough memory to allocate buffer");
//            edge->more = NULL;
//   	        edge->w = n;
//         	  edge->status = 'N';
////DeleteCutPairs			  
//         	  edge->deleted = 'N';						
////DeleteCutPairs
//      	    edge->more = LG[v];
//         	  LG[v] = edge;
//            edgeNum = edgeNum + 1;						
//         }
//		  }
//      else if (ch == 10) {
//         ch2[indx] = '\0';
//         indx = 0;
//				 if (!next_list) {
//   			    n = atoi(ch2);
//            if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_v_in_G)) ) == NULL)
//               abrt("Not enough memory to allocate buffer");
//            edge->more = NULL;
//            edge->w = n;
//            edge->status = 'N';
////DeleteCutPairs		
//         	  edge->deleted = 'N';
////DeleteCutPairs
//            edge->more = LG[v];
//            LG[v] = edge;
//            edgeNum = edgeNum + 1;						
//            next_list = 1;						
//         }
//      }
//      else {
//         ch2[indx] = (char)ch;
//         indx = indx + 1;
//      }
//   }
//	 if (!next_list) {
//      ch2[indx] = '\0';
//      n = atoi(ch2);
//      if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_v_in_G)) ) == NULL)
//         abrt("Not enough memory to allocate buffer");
//      edge->more = NULL;
//      edge->w = n;
//      edge->status = 'N';
////DeleteCutPairs
//      edge->deleted = 'N';			
////DeleteCutPairs
//      edge->more = LG[v];
//      LG[v] = edge;
//      edgeNum = edgeNum + 1;						
//   }
//   edgeNum = edgeNum / 2;	 
//   fclose(in);
	 
	 printf("\nComplexity of the input graph:\n|V| + |E| = %d + %d = %d\n",Vnum-1,edgeNum,Vnum+edgeNum-1);
	 
//DeleteCutPairs
   if (( Sgen_pointer = (adjacentG*)malloc(Vnum * sizeof(struct adjacent_with_v_in_G)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( Sgen_other_pointer = (adjacentG*)malloc(Vnum * sizeof(struct adjacent_with_v_in_G)) ) == NULL)
      abrt("Not enough memory to allocate buffer");
//DeleteCutPairs

   if (( Sbedge = (char*)malloc(Vnum * sizeof(char)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( Sx = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( Sy = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( Sp = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( Sq = (int*)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");
	clock_t start_time = clock(); 			
//   r = 1;
   for (r=1;r<Vnum;r++)
   {
      if (visited[r] == 'N')
	  {

         if (r > 1)
		 {
            compNum = compNum + 1;
         }

/*//YesOrNo		 
         if (r>1) {
            printf("It's a NO instance!\n");
            freeMem(Vnum); 
            exit(1);
         }
*///YesOrNo		 

         find_cut_pair(r,0);
      }
   }
   clock_t finish_time = clock() - start_time;

/*//YesOrNo
   freeMem(Vnum); 
   printf("It's a YES instance!\n");	 
   exit(1);
*///YesOrNo	 

//3edge
   for (v=1;v<Vnum;v++) {
      if (newlink[v] != 0) {
//PRINT			
         //printf("\n(%d,%d) is an additinal edge!",v,newlink[v]);
//PRINT
         if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_v_in_G)) ) == NULL)
            abrt("Not enough memory to allocate buffer");
         edge->w = newlink[v];
         edge->status = 'N';
         edge->deleted = 'N';				 
         edge->more = LG[v];
         LG[v] = edge;

         if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_v_in_G)) ) == NULL)
            abrt("Not enough memory to allocate buffer");
         edge->w = v;
         edge->status = 'N';
         edge->deleted = 'N';				 
         edge->more = LG[newlink[v]];
         LG[newlink[v]] = edge;
      }
   }
   for (v=1;v<Vnum;v++)
   visited[v] = 'N';
   for (v=1;v<Vnum;v++)
      if (visited[v] == 'N')
	  {
//PRINT-3edge
         printf("\nA 3edge-connected component: ");
//PRINT-3edge
         DFS(v);
		 comIndex++;
      }
	 
//3edge
	  
//component size analysis
	int* compSize =  (int*)malloc(sizeof(int) * comIndex); 
	for (int i = 1; i <Vnum; i++)
	{
		if(compSize[components[i] - 1] >= 1)
			compSize[components[i] - 1]++;
		else compSize[components[i] - 1] = 1;
	}
	FILE *f = fopen("compSizes.txt", "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}
	for (int i = 0; i < comIndex - 1; i++)
	{
		fprintf(f, "%d\n", compSize[i]);
	}
	fclose(f);
	free(compSize);
//component size analysis

//Boundary Vertices
	calcCompBndryInfo(comIndex);
//Boundary Vertices

//Closeness
	int* fars = farness(LG, true);
	for(int i = 1; i < Vnum; i++)
	{
		int cci = components[i];
		for(int j = 1; j < Vnum; j++)
		{
			if(components[j] != cci)
			{
				int minFar = INT_MAX;
				int** bDists = compBndryInfo[cci].compBndryVrtxDists;
				for(int bi = 0; bi < compBndryInfo[cci].compBndrys.size(); bi++)
				{
					if(bDists[bi][i] < 0 || bDists[bi][j] < 0)
					{
						continue;
					}
					int far = bDists[bi][i] + bDists[bi][j];
					if(far < minFar)
					{
						minFar = far;
					}
				}
				if(minFar != INT_MAX)
				{
					fars[i] += minFar;
				}
				
			}
		}
	}
	printf("\n");
	printArray(farness(LG, false), Vnum);
	printArray(fars, Vnum);
	delete[] fars;
//Closeness
	compNum = compNum + 1;
	printf("\n# of nodes: %d",Vnum);
	printf("\n# of edges: %d",edgeNum);
	printf("\nRunning Time: %2.3f seconds.", finish_time / (double)CLOCKS_PER_SEC);
	printf("\n# of cut-pairs: %d",cutPairNum);
	printf("\n# of 3edge-connected components: %d\n",compNum);   
	printf("\n# of 2edge-connected components: %d\n", bridgeNum + orignCC);
//	printArray(isBoundry, Vnum);
//	printArray(dists, Vnum);
//	printArray(farness(LG, true), Vnum);
//	printArray(farness(LG, false), Vnum);
	freeMem(Vnum);	
	return 0;
}