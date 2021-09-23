#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <fstream>
#define n 150  //the number of agents
#define m 300  //twice of times of the number of agents
#define T 1000 //the simulation periods
#define wgt 1  // herding factor w
#define alpha 0.01 //parameter alpha in the amount calculation
#define ESP 0.00000001 //control the accuracy of bisection method
#define F(x) x*(exp(x))- y //transforming equation to solve the price which determined by the cash position
using namespace std;

int main()

{
    int r=0,t,k,j,i,h;
    double k1,p; //p-reference price observed by all agents
    int *a=new int[n];//initial outgoing links
    int *a1 = new int[n];//newly random formed links
    int *a2 = new int[n];//final links selected by the probability
    double *c= new double[n]();//agent's cash
    double *s = new double[n]();//agent's stock
    double *pt=new double [T]();//all reference price at time t
    double *w = new double[n]();//agent's wealth
    double *f = new double[n]();//agent's fitness
    double *proi = new double[n]();//agent's probability calculated by fitness
    double *volatility = new double[n]();//agent's volatility
    double *ret=new double [n]();//agent's expected return
    double *p1=new double [n]();//agent's expectation price
    double *pcash=new double [n]();//agent's price calculated by the cash position
    double *p2=new double [n]();//agent's final submitted price in the order market
    double *amount=new double [n]();//agent's final submitted amount in the order market
    int *ran=new int[n];//includes n numbers
    int *sequence=new int[n];//random sequence
    double order[n][6]={};//order book information before the transaction
    double order1[m][6]={};//order book information after the transaction
    double buy[n][6]={};//all buy orders
    double sell[n][6]={};//all sell orders
    double sum,sum1,average,variance,hold;
    int tp;

    h=200;//agent's time horizon
    tp=0;
    double *ptemp=new double [h]();
    sum=0;sum1=0;average=0;variance=0;

    int flag = 0;//control variable
    unsigned int seed=10; //set the seed equals to 1
    srand(seed); //using the seed to generate random number
    //srand((unsigned)time(NULL));//using the time seed to generate the random number
    for (i=0;i<n;i++){a[i]=1000;a1[i]=1000;a2[i]=1000;}//initialization for the outgoing links
    for(i=0;i <n;i++) {c[i] = 100000;s[i]=100;}//initialization for cash and stock position
    p = 1000;//initialization for reference price

while (r < n){

for(j=0;j<=n;j++)

{ k=rand()% n;//k is a random number from 0 to n-1
if(r==k)//if r equals to k
{
flag = 0;
break; //then break out of the loop

}
if(r!=k)//if k not equals to r (the agents not copy themselves
{


flag = 1;// go to next loop
}
}

if(flag==1)//if k not equals to r
{
    a [r]= k;//assign k into the matrix
    /*cout << r << setw(2) << a [r] << endl;*/
   r++;//go to next element
}
}
r=0;//after these loops r becomes zero
	ofstream outfile("material.txt",ios::out|ios::app); //create or open the txt file named material
	outfile<<("initial matrix\n");// print the character initial matrix and goes to the next line
	for(i=0;i<n;i++)
	{
        outfile<< a[i] <<"  ";//record the matrix a (initial outgoing links) into the material file
	}
	outfile<<endl;//goes to the next line
 for (t=0;t<T;t++){cout<<t<<endl;//repeat for t times print the time t each time
     if (t>0){for(i=0;i < n;i++) {a[i]=a2[i];/*c[i]=order1[i][4];s[i]=order1[i][3];/*cout << i << setw(2)<< a[i] << endl;*/}}
     /*srand((unsigned)time(NULL));*/pt[t]=p; cout << pt[t] << endl;//after t gets bigger than 0, the old outgoing matrix
     //is the previously final matrix, record all price in matrix pt accordingly
    outfile << ' ' << t << endl;
    outfile <<' '<< pt[t] << endl;//record the t and pt in material file
	outfile<<("existing matrix\n");//write a line with existing matrix
	for(i=0;i<n;i++) //for i from 0 to n-1
	{
        outfile<< a[i] <<"  "; //record all old outgoing matrix called a
	}
	outfile<<endl; //end
	//outfile.close();
r=0; //initialize r to 0
while (r < n){ // when r is smaller than n

for(j=0;j<=n;j++)//for j from 0 to n-1

{ k=rand()% n;//k is a random number from 0 to n-1
if(r==k||k==a[r])//one constrain is k can not be r or a[r]
{
flag = 0;
break;

}
if(r!=k&&k!=a[r])
{


flag = 1;
}
}

if(flag==1)
{
    a1 [r]= k;//if k is not equal to r or a[r] (the old outgoing link), the new outgoing is k in a1[r]
    /*cout << r << setw(2) << a1 [r] << endl;*/
   r++;//next one
}
}

outfile<<("new matrix\n"); // print the new line - new matrix in material file
for(i=0;i<n;i++)
	{
        outfile<< a1[i] <<"  ";
	}
outfile<<endl;//print all elements in a1

r=0;//initialize r to 0

    int n1,n2;
    double betai;
    double wmax=0;

outfile<<("wealth and max wealth\n");//print wealth and max wealth
   for(i=0;i<n;i++){
             w[i] = c[i] + p * s[i];//formulation for wealth
             if (w[i] > wmax)//find the max wealth
             {wmax = w[i];}
            outfile << w [i] << ' ';}
            outfile << endl;
            outfile << wmax << endl;// print the wealth for each agent and the maximum wealth
            outfile<<("fitness,beta and probability\n");// print the fitness, beta and probability in material file
    for (i=0;i<n;i++) {f[i] = w[i] / wmax;//the fitness formula
     outfile << f[i]<< ' ';}
     outfile <<endl;// print fitness for each agent
     for (i=0;i<n;i++) {n1=a[i];n2=a1[i];betai = (double) rand()/RAND_MAX * 40 + 5;proi[i] = 1 / (1 + (exp (-betai*(f[n2] - f[n1]))));
    outfile <<betai<< ' '<<proi[i]<< endl;}//choose a random number distributed in uniform for beta and then calculate the probability and print
     outfile << endl;
     outfile<<("final matrix\n");// print final matrix
     for(i=0;i < n;i++) {
       k1 =((double)rand())/RAND_MAX; //produce random number k1
      outfile <<"random number"<< k1 << endl;//print k1
       if (k1 >= proi[i]) {a2 [i]= a [i];}
       if (k1 < proi[i]){ a2[i]=a1[i];}
       outfile << i << setw(5) << a2 [i] << endl;}//according to the rules, the final links formed
    double *frequency=new double[n]();

    for (i=0;i < n;i++) {++frequency[a2[i]];}
   outfile << "Rating" << setw (17) << "Frequency" <<endl;
    for (int rating=0; rating<n;rating++)
        {outfile <<setw(6)<< rating <<setw(17) <<frequency [rating] <<endl;}//number of links for each agent is calculated

    static double V1, V2, S;
    static int phase = 0;
    double X,pmax=0;
    int n3,n4,total=0;

    double std0;
    //srand(time(NULL));
    outfile<<("number, std0, volatility, normal noise and idiosyncratic return\n");
    for (i=0;i < n;i++) {
    std0=(double) rand()/RAND_MAX * 0.01; volatility[i]=std0 * (1+((frequency[i]/n )*(1-wgt)));
    if(phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
            } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        {X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;}//random number from normal distribution X is produced
    ret[i]=volatility[i]*X;//formula for return
    outfile << i << ' '<< std0 <<' '<< volatility [i] << ' '<< X << ' '<< ret[i]<<endl;}//print in material file


    int index,x=0;
    for(i=0;i<n;i++)
    ran [i]=i+1;
    //srand((unsigned)time(0));
    for(i=0;i<n;)
    {
    index=rand()%n;
    if(ran[index]!=0)
    {
    sequence[i]=ran[index]-1;
 //printf("%d\n",sequence[i]);
    ran[index]=0;
    ++i;

 }
 }//produce random sequence
 outfile<<("sequence, final matrix with sequence and revised return\n");//print the words
    for (i=0;i < n;i++) {n4=sequence [i];
    n3=a2[n4];
   // ret1[n4]=(wgt*ret[n4]) + ((1-wgt)*ret[n3]);
     ret[n4]=(wgt*ret[n4]) + ((1-wgt)*ret[n3]);
   // cout << n4<< ' '<<n3 << ' '<< ret[n4]<<' '<< ret[n3] << ' '<<ret1[n4] <<endl;//}
    //for (i=0;i < n;i++) {p1[i]=p[i]*exp(ret1[i]);
    outfile << n4<< ' '<<n3 <<' '<<ret[n4] <<endl;}//revised their return and print in the material file

 //for (r=0;r<n;r++)   {cout<<ret[i];order1[i][0]=sequence[i];order1[i][1]=ret[i];}

 outfile<<("revised return,expected price, max price and total\n");//print
    for (i=0;i < n;i++) {p1[i]=p*exp(ret[i]*sqrt(h));if (p1[i]>pmax){pmax=p1[i];} total+=p1[i];
    outfile << ret[i]<<' '<< p1[i]<<endl;}
    outfile << endl;
    outfile <<pmax<<' '<<total <<endl;
 outfile<<("time period, temp period, sum of return, average, sum of distance from average, unconditional variance \n");
//calculate the future price and find the maximum one
if (tp==h) {tp=0;}//variance based on the latest h returns
pt[t]=p;pt[-1]=999;
ptemp[tp]=log(pt[t]/pt[t-1]);for(i=0;i<h;i++){sum+=ptemp[i];} if (t<h) {average=sum/(t+1);}
else {average=sum/h;}
if (t<h) {for(i=0;i<=t;i++){outfile <<pt[t]<<' '<<ptemp[i]<< ' ';sum1+=(average-ptemp[i])*(average-ptemp[i]);variance=sum1/h;if(variance==0){variance=variance+0.0001;}}}
else {for(i=0;i<h;i++){outfile <<pt[t]<<' '<<' '<<ptemp[i]<< ' ';sum1+=(average-ptemp[i])*(average-ptemp[i]);variance=sum1/h;}}
outfile << t <<' '<<tp <<' '<<sum <<' '<< average<<' '<< sum1<<' '<< variance << endl;
tp++;sum=0;sum1=0;//calculate the variance
 outfile<<("variance for each agent\n");
    double *var=new double [n]();
    for (i=0;i < n;i++) {var[i]=variance*(1-((1-wgt)*(frequency[i]/n)));
    outfile << var[i]<<' ';}
    outfile<<endl;//calculate the variance for each agent and print
    // outfile<<("current portfolio\n");

    double x0,x1,x2;
    double f1,f2,f0,y,y1,y2,y3,y4,y5,y6,y7,y8;
    /*for (j=0;j < n;j++) {if (s[j]==0) {y1=alpha*var[j];y2=y1*p1[j];y=y2*s[j];pcur[j]=p1[j];}else{i=1;y1=alpha*var[j];y2=y1*p1[j];y=y2*s[j]; x0=0; x1=p1[j];
  do
  {
  x2=(x0+x1)/2;
  f0=F(x0);
  f1=F(x1);
  f2=F(x2);
  if(f0*f2<0)
   {
    x1=x2;
   }
   else
   {
    x0=x2;
   }
   i++;
  }while(fabs(f2)>ESP); pcur[j]=x2/y*p1[j];}outfile << j <<' '<<y1<<' ' << y2 <<' '<< y<< ' '<< x2<< ' '<< pcur[j]<< endl;
}*/
 outfile<<("cash position\n");
for (j=0;j < n;j++) {if (s[j]==0) {y1=alpha*var[j];y2=y1*p1[j];y3=exp(y1*c[j]);y=y2*s[j]/y3;pcash[j]=y2/y1/y3;}else {i=1;y1=alpha*var[j];y2=y1*p1[j];y3=exp(y1*c[j]);y=y2*s[j]/y3; x0=0; x1=y2*s[j];
  do
  {
  x2=(x0+x1)/2;
  f0=F(x0);
  f1=F(x1);
  f2=F(x2);
  if(f0*f2<0)
   {
    x1=x2;
   }
   else
   {
    x0=x2;
   }
   i++;
  }while(fabs(f2)>ESP); pcash[j]=x2/y1/s[j];} y4=p1[j];y5=log(y4/(pcash[j]));y6=y5/y1/pcash[j];y7=y6-s[j];y8=pcash[j]*y7;
  outfile << j <<' '<<y1<<' ' << y2 <<' '<< y<< ' '<< x2<< ' '<< pcash[j] <<' '<<y8<< endl;//bisection method to find the price for cash position
}
 outfile<<("price and amount in the order book\n");

    for (i=0;i < n;i++) {y1=pcash[i];y2=p1[i];p2[i] = (double) rand()/RAND_MAX * (y2-y1) + y1;y3=p2[i];amount[i]=(log(p1[i]/y3))/(alpha*var[i]*y3);y=amount [i];amount[i]=floor(y);
    outfile << p2[i]<<' '<< amount[i]<< endl;}//find the summited price and amount

    double ask,bid,ask1=0,bid1=0,ask2=0,bid2=0,p3,p4;
    int sumbuy=0,sumsell=0,j1,n5,bn=0,sn=0,mbn=0,msn=0;

    for (r=0;r<n;r++) { for (i=0;i<6;i++) buy[r][i]=-1;}
    for (r=0;r<n;r++) { for (i=0;i<6;i++) sell[r][i]=pmax+301;;}//initialize the buy and sell order
    //in case that the blank one will not sort in the front of the order book


     bid=0;ask=pmax;
     outfile<<("order book:sequence,price,amount,stock,cash,buy/sell\n");

        for(i=0;i<n;i++)
    ran [i]=i+1;
    //srand((unsigned)time(0));
    for(i=0;i<n;)
    {
    index=rand()%n;
    if(ran[index]!=0)
    {
    sequence[i]=ran[index]-1;
 //printf("%d\n",sequence[i]);
    ran[index]=0;
    ++i;}}//produce the random sequence
    for(r=0;r < n;r++){
                       order[r][0] = sequence[r];n3=sequence[r];y=p2[n3];y1=amount[n3];if(y1>s[n3])//if (y<pcur[n3])
            {order[r][5] = 0;order[r][1]=p2[n3];} if (y1<s[n3]) {order[r][5] = 1;order[r][1]=p2[n3];}
            if(y1==s[n3]){order[r][5]=-1;order[r][2]=order[r][1]=0;}
            if (y1>s[n3]){order[r][2]= y1-s[n3];} if (y1<s[n3]) {order[r][2]=s[n3]-y1;}
            if((order[r][5]==0)&&(order[r][1]>bid)){bid=order[r][1];}
            if((order[r][5]==1)&&(order[r][1]<ask)){ask=order[r][1];}
            if(order[r][5]==0){sumbuy++;}
            if(order[r][5]==1){sumsell++;}
            order[r][3]=s[n3];order[r][4]=c[n3];
            for (i=0;i<6;i++) outfile << order[r][i]<<' ';outfile<<endl;}
            //outfile <<sumbuy<<' '<<sumsell<<' '<<mbn<<' '<<msn<<endl;}
//record sequence, price, amount, stock, cash, buy or sell position
    for (r=0;r<n;r++) {for (i=0;i<6;i++){if (order[r][5]==0||order[r][5]==-1) buy[r][i]=order[r][i];
                     if (order [r][5]==1) sell[r][i]=order[r][i];}
                     for (n5=0;n5<=r-1;n5++){//repeat for r-1 times
                     for (j=0;j<=r-1;j++) {//compare two of them bubble sort
                       if (buy[j][1] < buy[j+1][1] )
                        {hold=buy[j][0]; buy[j][0]=buy[j+1][0];buy[j+1][0]=hold;
                        hold=buy[j][1]; buy[j][1]=buy[j+1][1];buy[j+1][1]=hold;
                        hold=buy[j][2]; buy[j][2]=buy[j+1][2];buy[j+1][2]=hold;
                        hold=buy[j][3]; buy[j][3]=buy[j+1][3];buy[j+1][3]=hold;
                        hold=buy[j][4]; buy[j][4]=buy[j+1][4];buy[j+1][4]=hold;
                        hold=buy[j][5]; buy[j][5]=buy[j+1][5];buy[j+1][5]=hold;}
                        if (sell[j][1]> sell[j+1][1] )
                        {hold=sell[j][0]; sell[j][0]=sell[j+1][0];sell[j+1][0]=hold;
                        hold=sell[j][1]; sell[j][1]=sell[j+1][1];sell[j+1][1]=hold;
                        hold=sell[j][2]; sell[j][2]=sell[j+1][2];sell[j+1][2]=hold;
                        hold=sell[j][3]; sell[j][3]=sell[j+1][3];sell[j+1][3]=hold;
                        hold=sell[j][4]; sell[j][4]=sell[j+1][4];sell[j+1][4]=hold;
                        hold=sell[j][5]; sell[j][5]=sell[j+1][5];sell[j+1][5]=hold;}}}
                        bid2=buy[bn][1];ask2=sell[sn][1];
                      for(j1=0;j1<=r;j1++){
                        if((buy[bn][1]>=sell[sn][1])&&(buy[bn][2]<sell[sn][2]))
                            {buy[bn][3]=buy[bn][3]+buy[bn][2];buy[bn][4]=buy[bn][4]-(buy[bn][1]*buy[bn][2]);
                            sell[sn][3]=sell[sn][3]-buy[bn][2];sell[sn][4]=sell[sn][4]+(sell[sn][1]*buy[bn][2]);
                            sell[sn][2]=sell[sn][2]-buy[bn][2];
                            bid1=buy[bn][1];ask1=sell[sn][1];
                            buy[bn][1]=-1;buy[bn][2]=0;bn=bn+1;sn=sn;mbn=mbn+1;msn=msn+1;
                            bid2=buy[bn][1];ask2=sell[sn][1];}

                        if((buy[bn][1]>=sell[sn][1])&&(buy[bn][2]==sell[sn][2]))
                            {buy[bn][3]=buy[bn][3]+buy[bn][2];buy[bn][4]=buy[bn][4]-(buy[bn][1]*buy[bn][2]);
                            sell[sn][3]=sell[sn][3]-buy[bn][2];sell[sn][4]=sell[sn][4]+(sell[sn][1]*buy[bn][2]);
                            bid1=buy[bn][1];ask1=sell[sn][1];
                            //if(buy[bn][1]<ask1){ask1=buy[bn][1];}if(sell[sn][1]>bid1){bid1=sell[sn][1];}
                            buy[bn][1]=-1;sell[sn][1]=pmax+301;buy[bn][2]=sell[sn][2]=0;bn=bn+1;sn=sn+1;mbn=mbn+1;msn=msn+1;
                            bid2=buy[bn][1];ask2=sell[sn][1];}

                        if((buy[bn][1]>=sell[sn][1])&&(buy[bn][2]>sell[sn][2]))
                            {buy[bn][3]=buy[bn][3]+sell[sn][2];buy[bn][4]=buy[bn][4]-(buy[bn][1]*sell[sn][2]);
                            sell[sn][3]=sell[sn][3]-sell[sn][2];sell[sn][4]=sell[sn][4]+(sell[sn][1]*sell[sn][2]);
                            buy[bn][2]=buy[bn][2]-sell[sn][2];
                            bid1=buy[bn][1];ask1=sell[sn][1];
                            //if(buy[bn][1]<ask1){ask1=buy[bn][1];}if(sell[sn][1]>bid1){bid1=sell[sn][1];}
                            sell[sn][1]=pmax+301;sell[sn][2]=0;sn=sn+1;bn=bn;mbn=mbn+1;msn=msn+1;
                            bid2=buy[bn][1];ask2=sell[sn][1];}}
                           sn=bn=0;}
                          outfile<<("buyers,sellers and market orders(buyer/seller)\n");
                          outfile <<sumbuy<<' '<<sumsell<<' '<<mbn<<' '<<msn<<endl;
//once one new order entered, sorted it. once transaction happened, the smaller amount one(buy or sell) go to the next line to compare
                      //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << buy[r][i]<<' '; cout<<endl;}
                     //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << sell[r][i]<<' '; cout<<endl;}
                     outfile<<("after trading order book:sequence,price,amount,stock,cash,buy/sell\n");
                        for (r=0;r<m;r++){for (i=0;i<6;i++){if(r<n){order1[r][i]=buy[r][i];if(order1[r][0]==-1){order1[r][0]=10000;}}//if(buy[0][1]>ask) {ask=buy[0][1];}
                        if(r>=n){order1[r][i]=sell[r-n][i];if(order1[r][0]==-1){order1[r][0]=10000;}/*if(sell[0][1]<bid){bid=sell[0][1];}*/}}}

                        for (r=0;r<m;r++) {for (n5=0;n5<=r-1;n5++){for (j=0;j<=r-1;j++) {
                        if (order1[j][0] > order1[j+1][0])
                        {hold=order1[j][0]; order1[j][0]=order1[j+1][0];order1[j+1][0]=hold;
                        hold=order1[j][1]; order1[j][1]=order1[j+1][1];order1[j+1][1]=hold;
                        hold=order1[j][2]; order1[j][2]=order1[j+1][2];order1[j+1][2]=hold;
                        hold=order1[j][3]; order1[j][3]=order1[j+1][3];order1[j+1][3]=hold;
                        hold=order1[j][4]; order1[j][4]=order1[j+1][4];order1[j+1][4]=hold;
                        hold=order1[j][5]; order1[j][5]=order1[j+1][5];order1[j+1][5]=hold;}
                        }}}//have all buy and sell order in 300 matrix and sort by the sequence
                        for (r=0;r<n;r++){for (i=0;i<6;i++) outfile <<order1[r][i]<<' ';outfile<< endl;}//leave first n lines in the order1 matrix
                         // for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << buy[r][i]<<' '; cout<<endl;}
                       //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << sell[r][i]<<' '; cout<<endl;}
                         outfile<<("cash, stock and bid/ask price\n");
                       for (i=0;i<n;i++)  {c[i]=order1[i][4];s[i]=order1[i][3];outfile<<c[i]<<setw(15)<<s[i]<<setw(15);}
                       outfile<<bid<<' '<<ask <<' '<<endl;
                       if(sumbuy==0||sumsell==0){p=p;} else { p=(bid+ask)/2;}//if (mbn==0||msn==0) {p3= (bid+ask)/2;} else {p3=(bid1+ask1)/2;} }
                       //if (sumbuy==0) {p=ask;} if (sumsell==0) {p=bid;} if(bid2==-1||ask2==pmax+301||bid2==0||ask2==0) {p=(bid1+ask1)/2;} else {p=(bid2+ask2)/2;}
                        cout <<bid<<' '<<bid1<<' '<<bid2<<' '<<ask<<' '<<ask1<<' '<<ask2<<endl;
                       ofstream outfile;//updated cash, stock and price position accordingly. print them.
                       outfile.open("price.xls",ios::out|ios::app);
                       outfile <<p <<endl;
                       outfile.close();
                       outfile.open("expected price.xls",ios::out|ios::app);
                       outfile <<total<<endl;
                      outfile.close();
                       outfile.open("total buyers.xls",ios::out|ios::app);
                       outfile <<sumbuy<<endl;
                     outfile.close();
                    outfile.open("total sellers.xls",ios::out|ios::app);
                       outfile <<sumsell <<endl;
                     outfile.close();
                     outfile.open("market buyers.xls",ios::out|ios::app);
                       outfile<<mbn<<endl;
                     outfile.close();
                     outfile.open("market sellers.xls",ios::out|ios::app);
                       outfile <<msn<<endl;
                     outfile.close();
                      outfile.open("variance.xls",ios::out|ios::app);
                       outfile <<variance<<endl;
                     outfile.close();

                      outfile.open("matrix.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{
        outfile<< a2[i] <<endl;
	}
	outfile.close();
	                  outfile.open("frequency.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{
        outfile<< frequency[i] <<endl;
	}
	outfile.close();
		                  outfile.open("fitness.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{
        outfile<< f[i] <<endl;
	}
	outfile.close();
	 outfile.open("wealth.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{
        outfile<< w[i] <<endl;
	}
	outfile.close();
 outfile.open("stock.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{
        outfile<< s[i] <<endl;
	}
	outfile.close();
	 outfile.open("order size.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{
        outfile<< order[i][2] <<endl;
	}
	outfile.close();//print the most important information separately in the named file

                        }//t1++;}}
                        }


   /*for (r=0;r<n;r++) {for (i=0;i<6;i++){if (order[r][5]==0) buy[r][i]=order[r][i]; //cout << buy[r][i]<<' ';
                     if (order [r][5]==1) sell[r][i]=order[r][i]; }}//cout << sell[r][i]<<' ';}cout <<endl;}
                     //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << buy[r][i]<<' '; cout << endl;}
                     //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout <<' '<< sell[r][i]<< ' '; cout << endl;}
 for (r=1;r<n;r++) {
     for (j=0;j<n-1;j++)
                if (buy[j][1] < buy[j+1][1] )
                        {hold=buy[j][1]; buy[j][1]=buy[j+1][1];buy[j+1][1]=hold;
                        hold=buy[j][0]; buy[j][0]=buy[j+1][0];buy[j+1][0]=hold;
                        hold=buy[j][2]; buy[j][2]=buy[j+1][2];buy[j+1][2]=hold;
                        hold=buy[j][3]; buy[j][3]=buy[j+1][3];buy[j+1][3]=hold;
                        hold=buy[j][4]; buy[j][4]=buy[j+1][4];buy[j+1][4]=hold;
                        hold=buy[j][5]; buy[j][5]=buy[j+1][5];buy[j+1][5]=hold;}}
for (r=0;r<n;r++) { for (i=0;i<6;i++)cout << buy[r][i]<<' '; cout<<endl;}
 for (r=1;r<n;r++) {
     for (j=0;j<n-1;j++)
                if (sell[j][1]> sell[j+1][1] )
                        {hold=sell[j][1]; sell[j][1]=sell[j+1][1];sell[j+1][1]=hold;
                        hold=sell[j][0]; sell[j][0]=sell[j+1][0];sell[j+1][0]=hold;
                        hold=sell[j][2]; sell[j][2]=sell[j+1][2];sell[j+1][2]=hold;
                        hold=sell[j][3]; sell[j][3]=sell[j+1][3];sell[j+1][3]=hold;
                        hold=sell[j][4]; sell[j][4]=sell[j+1][4];sell[j+1][4]=hold;
                        hold=sell[j][5]; sell[j][5]=sell[j+1][5];sell[j+1][5]=hold;}}
for (r=0;r<n;r++) { for (i=0;i<6;i++)cout << sell[r][i]<<' '; cout<<endl;}*/
  //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << buy[r][i]<<' '; cout << endl;}
  //for (r=1;r<=n;r++) {for (j=0;j<r;j++)if (sell[j][1]>buy[j+1][1] ) {hold[i]=sell[j][i]; sell[j][i]=sell[j+1][i];sell[j+1][i]=hold[i];}for (i=0;i<6;i++) cout <<' '<< sell[r][i]<< ' '; cout << endl;}
   //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout <<' '<< sell[r][i]<< ' '; cout << endl;}
   //for (r=0;r<n;r++)  {cout << buy[r][i]<< endl;}

            //cout << order[r][0]<<' '<<order[r][1]<<' '<<order[r][2]<<' '<<order[r][3]<<' '<<order[r][4]<<' '
            //<<order[r][5]<<' '<<ask<<' '<<bid<< endl;}


      //for(r=0;r < n;r++){cout <<order1[r]<<endl;}



