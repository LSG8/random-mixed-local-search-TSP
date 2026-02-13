// [random3mixed] - [implementation of mixed local search (node swap,2opt,link swap) algorithm for open-loop tsp]
// input: Location of tsps with TSPLIB format
// output: result_RMLS.txt (final order with path length),swap_cum_info.txt (contains detail information of the operations)
// owner: [Lahari]
// version 1.0.0 - [1.3.2018]
//
// [LS]: [1.3.2018]: [implemented the algorithm in c]

//............................................................................................
//
// [unlocked] 
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<sys/time.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<dirent.h>

/*------------------End of header----------------------------*/
int getXnY(char[],float*, float*);
int localsearch_random_mixed(int,float*,float*,int*,int*,int*,int*, int*, int*, char[], char[], char[],float*,float*);
int executeSwap(float*,float*,int*,int,float*);
int executeRC(float*,float*,int*,int,float*);
int executeSwap_link(float*,float*,int*,int,float*);
int swap(int,int,int,int*,int);
float pathLength(int,float*,float*,int*);
int findMin(float*);
int reverse(int*, int, int);
float addlink(int,int,int,int,float*,float*,int*);
float remlink(int,int,int,int,float*,float*,int*);
// int result(int,float*,float*,int*,char[],int,int*,int*,int*, int*, char[], float[], int[],float*,float*);
int print_tsp(int*,int);
//int printDM(int,float[][100]);
int writeToSteps(char[]);
int writeOrders(int*,int);
float distanceval(int,int,float*,float*);
float findHvsDist(float,float,float,float);
float findEucDist(float,float,float,float);
int getRandomStart(int, int*, int*);
int allFunctions(char[],char[],char[]);
/*---------------End of function declarations----------------*/
	//get distance value
float distanceval(int row, int column, float *x, float *y)
{	float x1=0.0,x2=0.0,y1=0.0,y2=0.0,dist=0.0;
	//printf("row and col is:%d and %d\n",row,column);
	if(row!=0 && column!=0)
	{	x1=x[row-1];
		y1=y[row-1];
		x2=x[column-1];
		y2=y[column-1];
		dist=findHvsDist(x1,x2,y1,y2);
		// dist = findEucDist(x1,x2,y1,y2);
	}
	//printf("row and col: %d %d and %d %d\n",x1,y1,x2,y2);
	return dist;
}	
/*------------------End of distanceval--------------------------*/
	//get x and y array and no. of goals
int getXnY(char path[], float *r, float *c)
{	int stop=0,i=0,index;
	FILE *fp;
	char buff[255], id[100],*end;
	//strcat(path,"/game.txt");
	strcpy(buff,"start");
	strcpy(id,"start");
	fp = fopen(path, "r");
	while(strcmp(buff,"EOF")!=0 && stop==0)
	{	fgets(buff, 255, fp);
		if(strncmp(buff,"NODE_COORD_SECTION",5)==0)
			stop=1;
	}
	fgets(id,100,fp);
	while(strncmp(id,"EOF",3)!=0)
	{	index = strtof(id, &end);
		r[i] = strtof(end, &end);
		c[i] = strtof(end, &end);
		//printf("row is %f and col is %f \n",r[i],c[i]);
		fgets(id,100,fp);
		i++;
	}
	fseek(fp,0,SEEK_SET);
	fclose(fp);
	//printf("dimension is %d\n", i);
	return i;
}
/*-------------------End of getXnY------------------------*/
float findHvsDist(float lon1,float lon2,float lat1,float lat2)
{	// L = R \cdot \arccos ( \sin \theta_1 \cdot \sin \theta_2 + \cos \theta_1 \cdot \cos \theta_2 \cdot \cos (\phi_1 - \phi_2) ).
	float r= 6356.752; //km
	float pi = 3.141593,diffLat, diffLon,a,coef,len;
    
    lat1 = (pi/180)*(lat1);
    lat2 = (pi/180)*(lat2);
	
    diffLon = (pi/180)*(lon1 - lon2);
    diffLat = lat1 - lat2;
	
    a = sin(diffLat/2)*sin(diffLat/2)+cos(lat1)*cos(lat2)*sin(diffLon/2)*sin(diffLon/2);
    coef = 2 * asin(sqrt(a));
    len = r * coef;
		
	return len;
}
/*-------------------End of findHvsDist------------------------*/
	//find euc distance
float findEucDist(float x1,float x2,float y1,float y2)
{	float dist;
	dist= sqrt(pow((x1-x2),2)+pow((y1-y2),2));
	return dist;
}
/*-------------------End of euc dist finding------------------------*/
//start the algorithm
int localsearch_random_mixed(int N,float *x, float *y,int *tsp_c,int *suc,int *trial, int *ns, int *rc, int *ls, char *steps, char *lengths, char *trials, float *msec_all,float *msec)
{	int i,p1,p2,pos=0,pos2=0,success=0,op,tr=0,n=0,r=0,l=0;
	float cost=0.00,timing_algo=0.00, timing_all=0.0,newLen, preLen,imp=0.0;
	char algo[10],ln[10],tn[5],next[5];
	struct timeval start_algo, end_algo, end_all, t_algo, t_all;
	//printf("\nNow list: ");
	//print_tsp(tsp_c,N);
	/**** for step by step improved length, trial calculation****/
	/*preLen=pathLength(Noftarget,cur,tsp_c);
	newLen=preLen;
	strcpy(steps,"0");
	sprintf(ln,"%0.3f",preLen);
	strcpy(lengths,ln);
	strcpy(trials,"0");*/
	/****start working****/
	gettimeofday(&start_algo, NULL);
	srand(time(NULL));
	for(i=0;i<10000;i++)
	{	op=(1+rand()%(3));
		//op=1;
		switch(op){
			case 1:
				//executeSwap_random
				strcpy(algo,"swap");
				//tsp_c=executeSwap(cur,tsp_c,Noftarget,&newLen);//for step by step length calculation; slower
				executeSwap(x,y,tsp_c,N,&imp);//only final length calculation; faster
				break;
			case 2:
				//executeRC_first
				strcpy(algo,"RC");
				//tsp_c=executeRC(cur,tsp_c,Noftarget,&newLen);//for step by step length calculation; slower
				executeRC(x,y,tsp_c,N,&imp);//only final length calculation; faster
				break;
			case 3:
				//executeSwapLink
				strcpy(algo,"swapLink");
				//tsp_c=executeSwap_link(cur,tsp_c,Noftarget,&newLen);//for step by step length calculation; slower
				executeSwap_link(x,y,tsp_c,N,&imp);//only final length calculation; faster
				break;
		}
		//printf("start=%ld and end=%ld",start_t,end_t);
		/**** if step by step length calculation is needed****/
		/*if(newLen<preLen)
		{	print_tsp(tsp_c,N);
			writeOrders(tsp_c,N);
			//printf(" %f ",newLen);
			preLen=newLen;
			strcat(steps,";");
			strcat(lengths,";");
			strcat(trials,";");
			success++;
			tr=i+1;
			strcat(steps,algo);
			sprintf(ln,"%0.3f",preLen);
			strcat(lengths,ln);
			sprintf(tn,"%d",tr);
			strcat(trials,tn);
			if(strcmp(algo,"swap")==0)
				n++;
			else if(strcmp(algo,"RC")==0)
				r++;
			else if(strcmp(algo,"swapLink")==0)
				l++;
			timing=((float)(end_t-start_t)/CLOCKS_PER_SEC);
		}*/
		/****for only final calculation;faster****/
		if(imp==1.0)
		{	imp=0.0;
			success++;
			tr=i+1;
			gettimeofday(&end_algo, NULL);
			// timersub(&end_algo, &start, &t_algo);
			// timing_algo = (((int)t_algo.tv_sec)*1000)+(((int)t_algo.tv_usec)*0.001);
			timing_algo = end_algo.tv_sec - start_algo.tv_sec;
		}
		/**** if iteration by iteration result is needed****/
		/*if(i==2)
			writeToSteps(ln);
		if(i==4)
			writeToSteps(ln);
		if(i==5)
			writeToSteps(ln);
		if(i==6)
			writeToSteps(ln);
		if(i==7)
			writeToSteps(ln);
		if(i==9)
			writeToSteps(ln);
		if(i==19)
			writeToSteps(ln);
		if(i==29)
			writeToSteps(ln);
		if(i==39)
			writeToSteps(ln);
		if(i==49)
			writeToSteps(ln);
		if(i==69)
			writeToSteps(ln);
		if(i==89)
			writeToSteps(ln);
		if(i==99)
			writeToSteps(ln);
		if(i==149)
			writeToSteps(ln);
		if(i==199)
			writeToSteps(ln);
		if(i==299)
			writeToSteps(ln);
		if(i==399)
			writeToSteps(ln);
		if(i==499)
			writeToSteps(ln);
		if(i==699)
			writeToSteps(ln);
		if(i==999)
			writeToSteps(ln);
		if(i==4999)
			writeToSteps(ln);
		if(i==9999)
			writeToSteps(ln);*/
	}
	// gettimeofday(&end_all, NULL);
	// timersub(&end_all, &start_algo, &t_all);
	// timing_all = (((int)t_all.tv_sec)*1000)+(((int)t_all.tv_usec)*0.001);
	
	//printf("sec val is %d and microsec val is %d",(int)t_base.tv_sec,(int)t_base.tv_usec);
	//printf("Timing all: %f ",timing_all);
	//printf("Timing algo: %f",timing_algo);
	/*late = time(0);
	printf("current time is %f", (float(late)));*/
	//printf("\nsteps:%s and timing:%f\n",steps,timing);
	/**** if iteration by iteration result is needed****/
	/*strcpy(next,"\n");
	writeToSteps(next);*/
	*trial=tr;
	*suc=success;
	*msec_all=timing_all;
	*msec=timing_algo;
	//*ns=n;*rc=r;*ls=l;//calculation of step by step number of operations used
	return 0;
}
/*-----------End of localsearch_random_first function----------------*/
//write order
int writeOrders(int *tsp_c, int N)
{	char cwd[200];
	int i;
	char filename1[]="/steps_mst.txt";
	getcwd(cwd, sizeof(cwd));
	FILE *outfile;
	strcat(cwd,filename1);
	outfile=fopen(cwd, "a");
	fprintf(outfile,"\nResult for LS: ");
	for(i=0;i<N;i++)
		fprintf(outfile,"%d ",tsp_c[i]);
	fclose(outfile);
	return 0;
}
/*-----------End of writeOrders function----------------*/
//write
int writeToSteps(char len[])
{	char cwd[200];
	char filename1[]="/steps.txt";
	getcwd(cwd, sizeof(cwd));
	FILE *outfile;
	strcat(cwd,filename1);
	outfile=fopen(cwd, "a");
	fprintf(outfile,"%s\t",len);
	fclose(outfile);
	return 0;
}
/*-----------End of writetosteps function----------------*/
//swap
int executeSwap(float *x, float *y,int *tsp_c,int N,float *len) 
{	int T,D,pos,posT;
	float cost=0.0;
	posT=(rand()%(N-1));
	//printf("position of T=%d",posT);
	T=tsp_c[posT];
	D=(rand()%(N-1));
	//printf("\nT and D:%d and %d\n",T,D);
	
	pos=D;
	if(posT>D)
		D--;
	cost= addlink(T,D,posT,N,x,y,tsp_c)-remlink(T,D,posT,N,x,y,tsp_c);
	//printf("cost: %f\n",cost);
	if(cost<0)
	{	//printf("\nNow swap: ");
		swap(T,pos,posT,tsp_c,N);
		*len = 1.0;//for only final calculation;faster
		//print_tsp(tsp_c,N);
	}	
	//*len=pathLength(N,x,y,tsp_c);//for all step length calculation; slower
	return 0;
}
/*-----------End of executeswap function----------------*/
	//swapping
int swap(int target,int new_pos,int old_pos, int *tsp_h,int N)
{	int i;
	if(new_pos<old_pos)
	{	//--
		for(i=old_pos-1;i>=new_pos;i--)
			tsp_h[i+1]=tsp_h[i];
	}
	else
	{	//++
		for(i=old_pos+1;i<=new_pos;i++)
			tsp_h[i-1]=tsp_h[i];
	}
	//printf("\nnewpos=%d,t at newpos=%d,target=%d",new_pos,tsp_h[new_pos],target);
	tsp_h[new_pos]=target;
	return 0;
}
/*------------------End of swap------------------------------*/
	//start counting adding links
float addlink(int T,int D,int posT,int N,float *x, float *y,int *tsp_c)
{	float a=0.0;
	int link1_1,link1_2,link2_1,link2_2,link3_1,link3_2;
	if(posT!=D)
	{	//print_tsp(tsp_c,N);
		if(posT>0 && posT<N){
			link1_1=tsp_c[posT-1];
			link1_2=tsp_c[posT+1];
			a=a+distanceval(link1_1,link1_2,x,y);
		} else 
			a=0;
		if(D>=0){
			link2_1=tsp_c[D];
			link2_2=T;
			a=a+distanceval(link2_1,link2_2,x,y);
		}
		//printf("\nafter 2 a=%f ",a);
		if(D<N){
			link3_1=T;
			link3_2=tsp_c[D+1];
			a=a+distanceval(link3_1,link3_2,x,y);
		}
		/*printf("\nafter 3 a=%f ",a);
		printf("\nADD={(%d,%d)+(%d,%d)+(%d,%d)}=%f",link1_1,link1_2,link2_1,link2_2,link3_1,link3_2,a);*/
	}
	return a;	
}
/*-------------------End of addlink---------------------------*/
	//start counting removing links
float remlink(int T,int posD,int posT,int N,float *x, float *y,int *tsp_c)
{	float r=0.0;
	int link1_1,link1_2,link2_1,link2_2,link3_1,link3_2;
	if(posT!=posD)
	{	//print_tsp(tsp_c,N);
		if(posT>0){
			link1_1=tsp_c[posT-1];
			link1_2=T;
			r=r+distanceval(link1_1,link1_2,x,y);
		}else
			r=0;
		//printf("\nr=%f ",r);
		if(posT<N){
			link2_1=T;
			link2_2=tsp_c[posT+1];
			r=r+distanceval(link2_1,link2_2,x,y);
		}
		//printf("\nr=%f ",r);
		if(posD>=0 && posD<N){
			link3_1=tsp_c[posD];
			link3_2=tsp_c[posD+1];
			r=r+distanceval(link3_1,link3_2,x,y);
		}
		/*printf("\nr=%f ",r);
		printf("\nREM={(%d,%d)+(%d,%d)+(%d,%d)}=%f",link1_1,link1_2,link2_1,link2_2,link3_1,link3_2,r);*/
	}
	return r;
}
/*-------------------End of remlink---------------------------*/
	//RC
int executeRC(float *x, float *y,int *tsp_h,int N,float *len) 
{	int a,b,c,d,success=0,start=0,end=0,pair1_2,pair2_1;
	float cost=0.00,link1_p,link2_p,link1_n,link2_n,new_link,cur_link;
	pair1_2=rand()%(N-2);
	pair2_1=(1+rand()%(N-1));
	if(pair2_1>pair1_2)
	{	start=pair1_2;
		end=pair2_1;
	}	
	else if(pair2_1<pair1_2)
	{	start=pair2_1;
		end=pair1_2;
	}
	else
		start=end;
	b=tsp_h[start];
	c=tsp_h[end];
	if(start==0)
	{	link1_p=0.0;a=0;
		link1_n=0.0;
	}
	else
	{	a=tsp_h[start-1];
		link1_p=distanceval(a,b,x,y);
		link1_n=distanceval(a,c,x,y);
	}
	if(end==(N-1))
	{	d=0;
		link2_p=0.0;
		link2_n=0.0;
	}	
	else
	{	d=tsp_h[end+1];
		link2_p=distanceval(c,d,x,y);
		link2_n=distanceval(b,d,x,y);
	}	
	new_link=link1_n+link2_n;
	cur_link=link1_p+link2_p;;
	cost=(new_link-cur_link);
	//printf("\na=%d,b=%d,c=%d,d=%d,startV=%d,endV=%d,link1p=%f,link2p=%f,link1n=%f,link2n=%f,cost=%f\n",a,b,c,d,tsp_h[start],tsp_h[end],link1_p,link2_p,link1_n,link2_n,cost);	
	if((start!=end) && cost<0)
	{	//printf("\nNow 2opt: ");
		//reverse(address of b, address of d)
		//printf("\na=%d,b=%d,c=%d,d=%d \n",a,b,c,d);
		reverse(tsp_h,start,end);
		//printf("\nreversed list");
		//print_tsp(tsp_h,N);
		*len = 1.0;//for only final calculation;faster
	}
	//*len=pathLength(N,x,y,tsp_h);//for all step length calculation; slower
	return 0;
}
/*-----------End of executeRC function----------------*/
//swapLink
int executeSwap_link(float *x, float *y,int *tsp_c,int N,float *len)
{	int i,j,minIndex,temp,change=0;;
	float linkLen[4], len1=0.0, len2=0.0;
	//print_tsp(tsp_h,N);
	i=(1+rand()%(N-1));
	//printf("i is: %d",i);
		
	if(i!=1 && i!=(N-1))
	{	linkLen[0]=distanceval(tsp_c[i-1],tsp_c[i],x,y);
		linkLen[1]=distanceval(tsp_c[i-1],tsp_c[N-1],x,y);
		linkLen[2]=distanceval(tsp_c[0],tsp_c[i],x,y);
		linkLen[3]=distanceval(tsp_c[0],tsp_c[N-1],x,y);
		//printf("\nlen1, len2, len3, and len4 are: %f %f %f %f\n",linkLen[0],linkLen[1],linkLen[2],linkLen[3]);
		minIndex=findMin(linkLen);
		//printf("\nminimum len index is: %d\n",minIndex);
		switch(minIndex) {
			case 1:
			//printf("No change\n");
			//print_tsp(tsp_c,N);
			break;
			case 2:
			//printf("change1\n");
			reverse(tsp_c,i,N-1);
			change=1;
			//print_tsp(tsp_c,N);
			break;
			case 3:
			//printf("change2\n");
			reverse(tsp_c,0,(i-1));
			change=1;
			//print_tsp(tsp_c,N);
			break;
			case 4:
			//printf("change3\n");
			reverse(tsp_c,0,(i-1));
			reverse(tsp_c,i,N-1);
			change=1;
			//print_tsp(tsp_c,N);
			break;
	   }	
	}
	else
	{	len1=distanceval(tsp_c[i-1],tsp_c[i],x,y);
		len2=distanceval(tsp_c[0],tsp_c[N-1],x,y);
		//printf("\nlen1 and len2 is: %f %f\n",len1,len2);
		if(i==1)
		{	if(len2<len1)
			{	//printf("now change 0 to N\n" );
				change=1;
				temp=tsp_c[0];
				for(j=1;j<N;j++)
				{	tsp_c[j-1]=tsp_c[j];
				}
				tsp_c[N-1]=temp;
			}
			//print_tsp(tsp_c,N);
		}
		else
		{	if(len2<len1)
			{	//printf("now change N to 0\n" );
				change=1;
				temp=tsp_c[N-1];
				for(j=(N-2);j>=0;j--)
				{	tsp_c[j+1]=tsp_c[j];
				}
				tsp_c[0]=temp;
			}
			//print_tsp(tsp_c,N);
		}
	}
	
	if(change)
		*len = 1.0;//for only final calculation;faster
	//*len=pathLength(N,dist,t);//for all step length calculation; slower
	return 0;
}
/*-----------End of executeSwapLink function----------------*/
int findMin(float *len)
{	int i, minIndex=0;
	float min=len[0];
	for (i =1; i <4; i++)
    {   if (min > len[i])
		{	min = len[i];
			minIndex=i;
		}
    }
	return (minIndex+1);
}
/*-------------------End of length------------------------*/
	//Function to reverse the array
int reverse(int *tsp_h, int start, int end)
{	int temp;
	while (start < end)
	{	temp = tsp_h[start];
		tsp_h[start] = tsp_h[end];
		tsp_h[end] = temp;
		start++;             
		end--;
	}
	return 0;
}
   
/*--------------------End of reverse---------------------------*/
	//get path length
float pathLength(int Noftarget,float *x, float *y,int *tsp_c)
{	float len=0.0, val=0.0;
	int i;
	for(i=0;i<(Noftarget-1);i++)
	{   if(tsp_c[i]!=0 && tsp_c[i+1]!=0)
			val=distanceval(tsp_c[i],tsp_c[i+1],x,y);
		else
			val=0.0;
		//printf("\nThe distance value between %d and %d is %.17f ",tsp_c[i],tsp_c[i+1],val);
		len+=val;
		//printf("\nNow length is %.17f",len);
	}
	return len;
}
/*-------------------End of length------------------------*/
	//final result
// int result(int Noftarget,float *x, float *y,int *tsp_h,char path[],
// 	int szd_n,int *suc,int *trial,int *ns, int *rc, int *ls, char *steps, 
// 	char *lengths, char *trials, float *msec_all, float *msec)
// int result(int,float*,float*,int*,char[],int,int*,int*,int*, int*, char[], float[], int[],float*,float*);
// {	char filename1[]="/result.txt";
// 	char filename2[]="/swap_cum_info.txt";
// 	char path2[200];
// 	FILE *outfile,*infofile;
// 	int i;
// 	strcpy(path2,path);
// 	strcat(path,filename1);
// 	strcat(path2,filename2);
// 	//printf("\nPath in output%s",path3);
// 	outfile=fopen(path, "w");
// 	infofile=fopen(path2,"w+");
// 	fprintf(outfile,"\n");
// 	for(i=0;i<Noftarget;i++)
// 	{	fprintf(outfile," %d ",tsp_h[i]);
// 		//fprintf(tempfile,"%d",tsp_h[i]);
// 	}
// 	fprintf(outfile,"\n\n&%f",pathLength(Noftarget,x,y,tsp_h));
// 	//fprintf(infofile,"{\"suc\":%d,\"trial\":%d,\"nswaps\":\"%d\",\"rcs\":\"%d\",\"lswaps\":\"%d\",\"steps\":\"%s\",\"lengths\":\"%s\",\"trials\":\"%s\",\"time\":%0.9f,\"timeAll\":%0.9f}",*suc,*trial,*ns,*rc,*ls,steps,lengths,trials,*msec,*msec_all);
// 	fprintf(infofile,"{\"suc\":%d,\"trial\":%d,\"time\":%0.9f,\"timeAll\":%0.9f}",*suc,*trial,*msec,*msec_all);
	
// 	fclose(outfile);
// 	fclose(infofile);
	
// 	return 0;
// }
/*------------------End of result-------------------------*/
	//print order
int print_tsp(int *tsp_c, int N)
{	//printf("Current tsp:");
		int i;
		for(i=0;i<N;i++)
		{	printf("%d ",tsp_c[i]);
		}
		printf("\n");
	return 0;
}
/*------------------End of print_tsp-------------------------*/

	//find random start
int getRandomStart(int N, int *tsp_prev, int *tsp_h)
{	int num, i=0, j, duplicate=0,t;
	for (i = 0; i<N; i++) 
    	tsp_h[i]=tsp_prev[i];
	srand(time(NULL));
	for (i = 0; i<N-1; i++) 
    {	size_t j = i + rand() / (RAND_MAX / (N - i) + 1);
		//j = rand() % (i+1);
		t = tsp_h[j];
		tsp_h[j] = tsp_h[i];
		tsp_h[i] = t;
	}
	for (i = 0; i<N; i++) 
    	tsp_prev[i]=tsp_h[i];
	/*printf("random start:");
	print_tsp(tsp_h,N);*/
	return 0;
}
/*-------------------End of random start------------------------*/
	//execute all functions
int allFunctions(char fileName[],char name[],char w_file[])
{	int N,tsp_h[100],tsp_prev[100],i,j,suc=0,trial=0,ns=0,rc=0,ls=0;
	float x[100],y[100],pre_len=1000.0, new_len,msec=0.00,msec_all=0.00;
	char type[10],steps[1], lengths[1], trials[1];
	N=getXnY(fileName,x,y);
	for (i=0;i<N;i++)
	{	tsp_h[i]=i+1;
		tsp_prev[i]=i+1;
	}
	srand(time(NULL));
	for(i=0;i<25;i++)
	{	getRandomStart(N,tsp_prev,tsp_h);
		localsearch_random_mixed(N,x,y,tsp_h,&suc,&trial,&ns,&rc,&ls,steps,lengths,trials,&msec_all,&msec);
		//printf("Path_length: %f",pathLength(N,x,y,tsp_h));
		new_len=pathLength(N,x,y,tsp_h);
		if(new_len<pre_len)
		{	pre_len=new_len;
			
		}
	}
	// char w_file[]="/usr/local/www_root/o-mopsi/tsp/exp/O-MopsiGames/O_mopsi_result_LS.txt";
	FILE *outfile;
	outfile=fopen(w_file, "a");
	fprintf(outfile,"%s\t%f\n",name,pre_len);
	fclose(outfile);
	
	
	return 0;
}
/*-------------------End of all functions------------------------*/
	//main function
int main()
{	char cwd[1024],dataLoc[2500], filename[500],path[500], name[100];
	FILE *outfile;
	char w_file[2500];
	getcwd(cwd, sizeof(cwd));
	DIR *dir;
	struct dirent *files;
	
	strcpy(dataLoc,cwd);
	strcat(dataLoc,"/data/");
	
	strcpy(w_file,cwd);
	strcat(w_file,"/output/result_RMLS.txt");
	outfile=fopen(w_file, "w");
	fprintf(outfile,"Name\tFinalLen\n");
	fclose(outfile);
	
	if ((dir = opendir (dataLoc)) != NULL) {
	// print all the files and directories within directory 
		while ((files = readdir (dir)) != NULL) {
			if(strcmp(files->d_name, ".")!=0 && strcmp(files->d_name, "..")!=0)
			{	strcpy(filename, dataLoc);
				strcpy(name, files->d_name);
				strcat(filename,files->d_name);
				printf ("%s\n", filename);
				allFunctions(filename,name,w_file);
				//break;
			}	
			//printf ("%s\n", filename);
		}
		closedir (dir);
	}
	else
		printf("can't open dir");
	
	return 0;
}