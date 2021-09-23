#include <iostream>
#include <cstdlib>
//#include <ctime>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <fstream>
#define n 150 //the number of agents
#define m 300//150*2
#define T 1000//the simulation periods
#define wgt 1 // herding factor w
#define alpha 0.01 //parameter alpha in the amount calculation
#define run 100
#define detail 1
//#define ESP 0.00000001 //control the accuracy of bisection method
//#define F(x) x*(exp(x))- y //transforming equation to solve the price which determined by the cash position
using namespace std;

int main()

{
    int simulation; double *totalprice=new double [T](); double *squareprice=new double [T]();
    for(simulation=0;simulation<run;simulation++) { cout<<simulation<<endl;
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
    double *pricet=new double [n]();
    double *order_volume = new double[n]();
    //double *un_trade = new double[n]();
    //double *percent= new double[n](); //percentage of excedding profit
    //double *pcash=new double [n]();//agent's price calculated by the cash position
    double *p2=new double [n]();//agent's final submitted price in the order market
    double *amount=new double [n]();//agent's final submitted amount in the order market
    int *ran=new int[n];//includes n numbers
    int *sequence=new int[n];//random sequence
    int *trade_state=new int[n];
    int *guru_lifetime=new int[T];
    double order[m][4]={};//order book information before the transaction
    double new_entrance [4]= {0};
    double best_ask [5] = {0,100000,0,0,0};
    double best_bid [5] = {0};
    //double order[n][6]={};//order book information before the transaction
    //double order1[m][6]={};//order book information after the transaction
    //double buy[n][6]={};//all buy orders
    //double sell[n][6]={};//all sell orders
    double sum,sum1,average,variance,hold,pricetemp,amount_temp,stock_temp,guru_link,guru_fitness,average_expected,total_guru,total_imitator,total_noise,total_herding,sim_herding=0,average_gurulife;
    int tp,i1,i2,ordertemp,guru,imitator,noise,un_trade,guru_life=1,old_guru;
    int tp1=0;
    static double V1, V2, S;
    static int phase = 0;
    double X,pmax=0,percent=0;
    int n3,n4,n5;
    double total=0,trades=0,volumebuy=0,volumesell=0,tradevolume=0,sumsell=0,sumbuy=0,sumguru=0,number_guru=0,std0,maxpercent,average_guru,average_imitator,average_noise,herding,mkt_buy,mkt_sell,average_guru_order,average_rest_order;
    double dif=50,profit=0;
    double sumguruw=0,sumimitatorw=0,sumnoisew=0,sgo=0,sro=0;

    h=1;//agent's time horizon
    tp=0;
    double *ptemp=new double [h]();
    sum=0;sum1=0;average=0;variance=0;

    //for(simulation=0;simulation<run;simulation++) {
    int flag = 0;//control variable
    unsigned int seed=simulation; //set the seed equals to 1
    srand(seed); //using the seed to generate random number
    //srand((unsigned)time(NULL));//using the time seed to generate the random number
    //for(simulation=0;simulation<run;simulation++) {
    for (i=0;i<n;i++){a[i]=1000;a1[i]=1000;a2[i]=1000;}//initialization for the outgoing links
    for(i=0;i <n;i++) {c[i] = 100000;s[i]=100;}//initialization for cash and stock position
    for (i=0;i<n;i++){
      for(j=0;j<5;j++){
        order[i][j]=0;
      }
    }
    p = 1000;//initialization for reference price
    //cout<<"simu"<<trades<<endl;
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
/*	 while (r<h){
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
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;
    ptemp[r]=0.03333*X;
    r++;
    }
    r=0;
     for(r=0;r<h;r++){
      outfile << ptemp[r] << endl;
    }*/
 for (t=0;t<T;t++){//cout<<t<<endl; //repeat for t times print the time t each time
     if (t>0){for(i=0;i < n;i++) {a[i]=a2[i];/*c[i]=order1[i][4];s[i]=order1[i][3];/*cout << i << setw(2)<< a[i] << endl;*/}}
     /*srand((unsigned)time(NULL));*///pt[t]=p; totalprice[t]=totalprice[t]+ pt[t]; squareprice[t]=squareprice[t]+(pt[t]*pt[t]);//cout << best_bid[1]<<' '<<best_ask[1]<<' '<<pt[t] << endl;//after t gets bigger than 0, the old outgoing matrix
     //is the previously final matrix, record all price in matrix pt accordingly
    //outfile << ' ' << t << endl;
    //outfile <<' '<< pt[t] << endl;//record the t and pt in material file
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

   ofstream outfile("material.txt",ios::out|ios::app); outfile<<("wealth and max wealth\n");//print wealth and max wealth
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
       outfile << i << setw(5) << a2 [i] << endl;}
    double *frequency=new double[n]();

    for (i=0;i < n;i++) {++frequency[a2[i]];}
   outfile << "Rating" << setw (17) << "Frequency" <<endl;
   if (t>0) {old_guru=guru;} guru_link=0;
    for (int rating=0; rating<n;rating++)
        {   if(frequency[rating]>guru_link) {guru=rating; guru_link=frequency[rating];guru_fitness=f[rating];/*cout<<"guru"<<guru<<' '<<guru_link;*/}
            outfile <<setw(6)<< rating <<setw(17) <<frequency [rating] <<endl;}
    guru_link=guru_link/n; total_imitator=0; total_noise=0; imitator=0;noise=0;

    if (old_guru==guru) {guru_life=guru_life+1; guru_lifetime[t]=0;} if(old_guru!=guru) {guru_lifetime[t]=guru_life; guru_life=1;}
    //cout<<guru_life<<' '<<t;
    if ((t==(T-1))&&(old_guru==guru)){guru_lifetime[t]=guru_life;}
    if ((guru_lifetime[t]!=0)&&(t!=0)) {number_guru=number_guru+1;/*sumguru=sumguru+guru_lifetime[t];*/}
    if (guru_lifetime[T-1]==1) {number_guru=number_guru+1;}
    //cout << sumguru<<endl;

    //static double V1, V2, S;
    //static int phase = 0;
    //double X,pmax=0;
    //int n3,n4,total=0;

    //double std0;
    //srand(time(NULL));
    outfile<<("number, std0, volatility, normal noise and idiosyncratic return\n");
    for (i=0;i < n;i++) {
    std0=(double) rand()/RAND_MAX * 0.07; volatility[i]=std0 * (1+((frequency[i]/n )*(1-wgt)));
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

    phase = 1 - phase;}
    ret[i]=volatility[i]*X;
    outfile << i << ' '<< std0 <<' '<< volatility [i] << ' '<< X << ' '<< ret[i]<<endl;
    }


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
 }
 outfile<<("sequence, final matrix with sequence and revised return\n");
    for (i=0;i < n;i++) {n4=sequence [i];
    n3=a2[n4];
   // ret1[n4]=(wgt*ret[n4]) + ((1-wgt)*ret[n3]);
     ret[n4]=(wgt*ret[n4]) + ((1-wgt)*ret[n3]);
   // cout << n4<< ' '<<n3 << ' '<< ret[n4]<<' '<< ret[n3] << ' '<<ret1[n4] <<endl;//}
    //for (i=0;i < n;i++) {p1[i]=p[i]*exp(ret1[i]);
    outfile << n4<< ' '<<n3 <<' '<<ret[n4] <<endl;
    }

 //for (r=0;r<n;r++)   {cout<<ret[i];order1[i][0]=sequence[i];order1[i][1]=ret[i];}

 /*outfile<<("revised return,expected price, max price and total\n");
    for (i=0;i < n;i++) {p1[i]=p*exp(ret[i]*sqrt(h));if (p1[i]>pmax){pmax=p1[i];} total+=p1[i];
    outfile << ret[i]<<' '<< p1[i]<<endl;}
    outfile << endl;
    outfile <<pmax<<' '<<total <<endl;*/

     //cout << best_bid[1] <<' ' << best_ask[1] <<endl;


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
    ++i;}}

    //j=tp1*n;
    j=tp1*n;  trades=0; sumbuy=0; sumsell=0; total=0; //cout<<trades<<endl;
    for (i=0;i<n;i++){
        for (i2=0;i2<m;i2++) { if (t>0&&sequence[i]==order[i2][0]) {order[i2][1]=order[i2][2]=order[i2][3]=0;ordertemp=i2;
           outfile<<"the replace order is "<<i2<<' '<<order[i2][0]<<' '<< order[i2][1]<<' '<<order[i2][2]<<' '<<order[i2][3]<<endl;
          }
        //cout<<i2<<' '<<order[i2][0]<<' '<< order[i2][1]<<' '<<order[i2][2]<<' '<<order[i2][3]<<endl;
        }
        best_bid[1]=0;best_ask[1]=100000;
        for (i1=0;i1<m;i1++) {
        if((order[i1][3]==1) && (best_bid[1]<order[i1][1])) {
        best_bid[0]=order[i1][0]; best_bid[1]=order[i1][1]; best_bid[2]=order[i1][2]; best_bid[3]=order[i1][3];best_bid[4]=i1;
        }// best ask updated,need to compare them again
        }

         for (i1=0;i1<m;i1++) {
         if((order[i1][3]==-1) && (best_ask[1]>order[i1][1])) {
         best_ask[0]=order[i1][0]; best_ask[1]=order[i1][1]; best_ask[2]=order[i1][2]; best_ask[3]=order[i1][3]; best_ask[4]=i1;
         }
         }
         outfile<<"best bid and ask" << best_bid[1]<<' '<<best_bid[4]<<' '<<best_ask[1]<<' '<<best_ask[4]<<endl;
        if (t==0) {pricetemp=1000;}
        else {
        if ( best_bid[1]==0&&best_ask[1]==100000) {pricetemp=p;}
        if(best_bid[1]!=0&&best_ask[1]==100000) {pricetemp=best_bid[1];}
        if(best_bid[1]==0&&best_ask[1]!=100000) {pricetemp=best_ask[1];}
        if(best_bid[1]!=0&&best_ask[1]!=100000) {pricetemp=(best_bid[1]+best_ask[1])/2;}}
        //outfile <<"the current price is "<<pricetemp<<endl;
        n3=sequence[i]; new_entrance[0]=n3; pricet[n3]=pricetemp;

        p1[n3]=pricetemp*(1+ret[n3]);if (p1[n3]>pmax){pmax=p1[n3];} total+=p1[n3]; //cout<<n3<<' '<<total<<endl;
        outfile <<"the expected price is "<<n3<<' '<<p1[n3]<<endl;

        if (p1[n3]>pricetemp) {new_entrance[3]=1;trade_state[n3]=1;p2[n3]=pricetemp+dif; profit=p1[n3]-pricetemp; if (profit<dif) {p2[n3]=p1[n3]-dif;} //cout<<"buy"<<profit<<' '<<p2[n3]<<endl;
        new_entrance[1]=p2[n3];
        percent=(p1[n3]-p2[n3])/p2[n3]; //cout <<n3<<"percent"<<percent<<endl;
        maxpercent=(double) rand()/RAND_MAX * (1-0.5) + 0.5;
        //cout<<n3<<"maxpercent"<<maxpercent<<endl;
        if (percent>=maxpercent) {amount_temp=c[n3]/p2[n3];}
        if (percent<maxpercent) {amount_temp=c[n3]/maxpercent*percent/p2[n3];}
        amount[n3]=s[n3] + amount_temp;
        order_volume[n3]=amount_temp;
        //cout<<"the amount willint to hold is"<<amount_temp<<' '<<amount[n3]<<endl;
        }

        if (p1[n3]<pricetemp) {new_entrance[3]=-1;trade_state[n3]=-1;p2[n3]=pricetemp-dif; profit=pricetemp-p1[n3]; if (profit<dif) {p2[n3]=p1[n3]+dif;} //cout<<"sell"<<profit<<' '<<p2[n3]<<endl;
        new_entrance[1]=p2[n3];
        percent=(p2[n3]-p1[n3])/p2[n3]; //cout<<"percent" <<percent<<endl;
        maxpercent=(double) rand()/RAND_MAX * (1-0.5) + 0.5; //cout<<"maxpercent"<<maxpercent<<endl;
        stock_temp=s[n3]; amount_temp=c[n3]/maxpercent*percent/p2[n3];
        if (stock_temp>amount_temp) {amount[n3]=s[n3]-amount_temp; order_volume[n3]=amount_temp;}
        if (stock_temp<=amount_temp) {amount[n3]=0; order_volume[n3]=s[n3];}
        //cout<<"the amount willint to hold is"<<amount_temp<<' '<<amount[n3]<<endl;
        }
        //cout <<n3<<"percent"<<percent<<endl;
        //cout<<n3<<"maxpercent"<<maxpercent<<endl;

        new_entrance[2]= fabs(amount[n3]-s[n3]);
        if (p1[n3]==pricetemp) {new_entrance[1]=new_entrance[2]=new_entrance[3]=0;}



        outfile << "new entrance comes" << new_entrance[0] <<' '<< new_entrance[1]<<' '<<new_entrance[2]<<' '<<new_entrance[3]<<endl;
           if(new_entrance[3]==1){sumbuy++;}  if(new_entrance[3]==1&&t!=0){volumebuy=volumebuy+new_entrance[2]; /*cout<<" t>0 volumebuy "<<new_entrance[0]<<' '<<volumebuy<<endl;*/}
          if(new_entrance[3]==-1){sumsell++;}  if(new_entrance[3]==-1&&t!=0){volumesell=volumesell+new_entrance[2];/*cout<<" t>0 volumesell"<<new_entrance[0]<<' '<<volumesell<<endl;*/}
        /*if ((best_ask[1]==100000)&&(best_bid[1]==0)) {
          order[j][0]=new_entrance[0]; order[j][1]=new_entrance[1]; order[j][2]=new_entrance[2]; order[j][3]=new_entrance[3];
          //order[j][4]=tp1;
          cout<<order[j][0]<<order[j][1]<<order[j][2]<<order[j][3]<<j<<endl; j=j+1;
          //s[n3]= s[n3]+ (new_entrance[3]*new_entrance[2]); c[n3]=c[n3]-(new_entrance[3]*new_entrance[2]*new_entrance[1]);
            if (new_entrance[3]==1) {// best bid/ask updated
            best_bid[0]=new_entrance[0];best_bid[1]=new_entrance[1];best_bid[2]=new_entrance[2];best_bid[3]=new_entrance[3];best_bid[4]=j-1;
            cout<<best_bid[0]<<best_bid[1]<<best_bid[2]<<best_bid[3]<<best_bid[4]<<endl;
            }
            if (new_entrance[3]==-1) {// best bid/ask updated
            best_ask[0]=new_entrance[0];best_ask[1]=new_entrance[1];best_ask[2]=new_entrance[2];best_ask[3]=new_entrance[3];best_ask[4]=j-1;
            cout<< best_ask[0]<<best_ask[1]<<best_ask[2]<<best_ask[3]<<best_ask[4]<<endl;
            }
        }*/
        //best_ask[1]=100000;best_bid[1]=0;
        //cout<<"start position"<<trades<<endl;
        if ((best_ask[1]!=100000) || (best_bid[1]!=0)) {
          if (new_entrance[3]==1) {
            //if (new_entrance[1]>=best_ask[1]) {
              while ((new_entrance[1]>=best_ask[1])&&(new_entrance[2]>=best_ask[2])) { outfile<<"trading update ask"<<endl; trades=trades+1; //cout<<trades<<endl;
              if(t!=0) {tradevolume=tradevolume+best_ask[2]; //cout<<"t>0tradevolume" <<tradevolume<<endl;
                       }
              //if (new_entrance[2]>=best_ask[2]) {
              if (t==0) { n3= new_entrance[0];n4=best_ask[4]; n5=best_ask[0]; }
               if (t!=0) {n3= new_entrance[0]; s[n3]=s[n3]+best_ask[2]; c[n3]=c[n3]-(best_ask[2]*best_ask[1]); n4=best_ask[4]; n5=best_ask[0];
                s[n5]=s[n5]- best_ask[2]; c[n5]=c[n5]+(best_ask[2]*best_ask[1]);}
                //if(order[n4][4]!=tp1) {c[n5]=c[n5]+(best_ask[2]*best_ask[1]);}
                order[n4][1]=100000; order[n4][2]=0; order[n4][3]=0;
                new_entrance[2]=new_entrance[2]-best_ask[2];best_ask[2]=0; best_ask[1]=100000;
                outfile<<"stock and cash"<<s[n3]<<' '<<c[n3]<<' '<<s[n5]<<' '<<c[n5]<<endl;
                  for (i1=0;i1<m;i1++) {
                   if((order[i1][3]==-1) && (best_ask[1]>order[i1][1])) {
                     best_ask[0]=order[i1][0]; best_ask[1]=order[i1][1]; best_ask[2]=order[i1][2]; best_ask[3]=order[i1][3]; best_ask[4]=i1;
                     }// best ask updated,need to compare them again
                   }
                 outfile<<" best ask"<< best_ask[0]<<' '<<best_ask[1]<<' '<<best_ask[2]<<' '<<best_ask[3]<<' '<<best_ask[4]<<endl;
                 }
               if ((new_entrance[1]>=best_ask[1])&&(new_entrance[2]<best_ask[2])) { outfile<< "trading update ask amount"<< best_ask[2]<<new_entrance[2]<<endl; trades=trades+1; //cout<<trades<<endl;
                if (t==0) { n3= new_entrance[0];n4=best_ask[4]; n5=best_ask[0]; }
                if(t!=0) {n3= new_entrance[0]; s[n3]=s[n3]+new_entrance[2]; c[n3]=c[n3]-(new_entrance[2]*best_ask[1]); n4=best_ask[4]; n5=best_ask[0];
                s[n5]=s[n5]- new_entrance[2]; c[n5]=c[n5]+(new_entrance[2]*best_ask[1]);}
                 if(t!=0) {tradevolume=tradevolume+new_entrance[2]; //cout<<"t>0tradevolume" <<tradevolume<<endl;
                       }
                //if(order[n4][4]!=tp1) {c[n5]=c[n5]+(new_entrance[2]*best_ask[1]);}
                best_ask[2]=best_ask[2]-new_entrance[2];order[n4][2]=best_ask[2];new_entrance[2]=0;
                outfile<<"order updated"<<order[n4][2]<<"stock and cash"<<s[n3]<<' '<<c[n3]<<' '<<s[n5]<<' '<<c[n5]<<endl;
               }
            if (new_entrance[1]<best_ask[1]) { outfile << "no trading update bid"<<endl;
              order[j][0]=new_entrance[0]; order[j][1]=new_entrance[1]; order[j][2]=new_entrance[2]; order[j][3]=new_entrance[3];
              //order[j][4]=tp1;
              if (new_entrance[1]>best_bid[1]) {
                best_bid[0]=new_entrance[0];best_bid[1]=new_entrance[1];best_bid[2]=new_entrance[2];best_bid[3]=new_entrance[3];best_bid[4]=j;
                 outfile<<"best bid"<<best_bid[0]<<' '<<best_bid[1]<<' '<<best_bid[2]<<' '<<best_bid[3]<<' '<<best_bid[4]<<endl;
              }
              outfile<<"order book updating"<<order[j][0]<<' '<<order[j][1]<<' '<<order[j][2]<<' '<<order[j][3]<<' '<<j<<
              ' '<<endl;
              j=j+1;
            }
           }

           if (new_entrance[3]==-1) {
            //if (new_entrance[1]<=best_bid[1]) {
             while ((new_entrance[1]<=best_bid[1])&&(new_entrance[2]>=best_bid[2])) { outfile<<"trading update bid"<<endl; trades=trades+1; //cout<<trades<<endl;
              if(t!=0) {tradevolume=tradevolume+best_bid[2]; //cout<<"t>0tradevolume" <<tradevolume<<endl;
                       }
              //if (new_entrance[2]>=best_ask[2]) {
               // new_entrance[2]=new_entrance[2]-best_bid[2]; best_bid[2]=0; best_bid[1]=0;
                if (t==0) { n3= new_entrance[0];n4=best_bid[4];n5=best_bid[0]; }
                if (t!=0) {n3= new_entrance[0]; s[n3]=s[n3]-best_bid[2]; c[n3]=c[n3]+(best_bid[2]*best_bid[1]);n4=best_bid[4];n5=best_bid[0];
                s[n5]=s[n5]+best_bid[2]; c[n5]=c[n5]-(best_bid[2]*best_bid[1]);}
                //if(order[n4][4]!=tp1) {s[n5]=s[n5]+best_bid[2];}
                order[n4][1]=0; order[n4][2]=0; order[n4][3]=0;
                new_entrance[2]=new_entrance[2]-best_bid[2]; best_bid[2]=0; best_bid[1]=0;
                outfile<<"stock and cash"<<s[n3]<<' '<<c[n3]<<' '<<s[n5]<<' '<<c[n5]<<endl;
                  for (i1=0;i1<m;i1++) {
                   if((order[i1][3]==1) && (best_bid[1]<order[i1][1])) {
                     best_bid[0]=order[i1][0]; best_bid[1]=order[i1][1]; best_bid[2]=order[i1][2]; best_bid[3]=order[i1][3];best_bid[4]=i1;
                     }// best ask updated,need to compare them again
                   }
                   outfile<<"best bid"<<best_bid[0]<<' '<<best_bid[1]<<' '<<best_bid[2]<<' '<<best_bid[3]<<' '<<best_bid[4]<<endl;
                 }
               if ((new_entrance[1]<=best_bid[1])&&(new_entrance[2]<best_bid[2])) { outfile<< "trading update bid amount"<<endl; trades=trades+1; //cout<<trades<<endl;
                if(t!=0) {tradevolume=tradevolume+new_entrance[2]; //cout<<"t>0tradevolume" <<tradevolume<<endl;
                       }
                if (t==0) { n3= new_entrance[0];n4=best_bid[4];n5=best_bid[0]; }
                if (t!=0) {n3= new_entrance[0]; s[n3]=s[n3]-new_entrance[2]; c[n3]=c[n3]+(new_entrance[2]*best_bid[1]);n4=best_bid[4];n5=best_bid[0];
                s[n5]=s[n5]+new_entrance[2]; c[n5]=c[n5]-(new_entrance[2]*best_bid[1]);}
                //if(order[n4][4]!=tp1) {s[n5]=s[n5]+new_entrance[2];}
                best_bid[2]=best_bid[2]-new_entrance[2]; order[n4][2]=best_bid[2]; new_entrance[2]=0;
                outfile<<"order updating"<<order[n4][2]<<"stock and cash"<<s[n3]<<' '<<c[n3]<<' '<<s[n5]<<' '<<c[n5]<<' '<<endl;
               }
            if (new_entrance[1]>best_bid[1]) { outfile<<"no trading update ask"<<endl;
              order[j][0]=new_entrance[0]; order[j][1]=new_entrance[1]; order[j][2]=new_entrance[2]; order[j][3]=new_entrance[3];
              //order[j][4]=tp1;
              if (new_entrance[1]<best_ask[1]) {
                best_ask[0]=new_entrance[0];best_ask[1]=new_entrance[1];best_ask[2]=new_entrance[2];best_ask[3]=new_entrance[3];best_ask[4]=j;
                outfile<<"best ask"<< best_ask[0]<<' '<<best_ask[1]<<' '<<best_ask[2]<<' '<<best_ask[3]<<' '<<best_ask[4]<<endl;
              }
              outfile<<"order book updating"<<order[j][0]<<' '<<order[j][1]<<' '<<order[j][2]<<' '<<order[j][3]<<' '<<j<<endl;
              j=j+1;
            }
           }
          }
          if ((best_ask[1]==100000)&&(best_bid[1]==0)) {
          order[j][0]=new_entrance[0]; order[j][1]=new_entrance[1]; order[j][2]=new_entrance[2]; order[j][3]=new_entrance[3];
          outfile<<"order book updating"<<order[j][0]<<' '<<order[j][1]<<' '<<order[j][2]<<' '<<order[j][3]<<' '<<j<<' '<<endl; j=j+1;
          //s[n3]= s[n3]+ (new_entrance[3]*new_entrance[2]); c[n3]=c[n3]-(new_entrance[3]*new_entrance[2]*new_entrance[1]);
            if (new_entrance[3]==1) {// best bid/ask updated
            best_bid[0]=new_entrance[0];best_bid[1]=new_entrance[1];best_bid[2]=new_entrance[2];best_bid[3]=new_entrance[3];best_bid[4]=j-1;
            outfile<<"best bid"<<best_bid[0]<<' '<<best_bid[1]<<' '<<best_bid[2]<<' '<<best_bid[3]<<' '<<best_bid[4]<<endl;
            }
            if (new_entrance[3]==-1) {// best bid/ask updated
            best_ask[0]=new_entrance[0];best_ask[1]=new_entrance[1];best_ask[2]=new_entrance[2];best_ask[3]=new_entrance[3];best_ask[4]=j-1;
            outfile<<"best ask"<< best_ask[0]<<' '<<best_ask[1]<<' '<<best_ask[2]<<' '<<best_ask[3]<<' '<<best_ask[4]<<endl;
            }
        } //cout << best_bid[1]<<' '<<best_ask[1] <<' '<<endl;
        }
        //for (i=0;i<((tp1+1)*n);i++) {
         // cout<<"order book updating"<<order[i][0]<<' '<<order[i][1]<<' '<<order[i][2]<<' '<<order[i][3]<<' '<<order[i][4]<<endl;
        //}
      /*for (i=(tp1*n);i<((tp1+1)*n);i++) {n3=order[i][0];
        if(order[i][3]==-1){s[n3]=s[n3]+(order[i][2]*order[i][3]);}
        if (order[i][3]==1) {c[n3]=c[n3]-(order[i][1]*order[i][2]*order[i][3]);}
        //cout<<i<<n3<<' '<<s[i]<<' '<<c[i]<<endl;
      }*/
      un_trade=0;
      if(best_bid[1]==0||best_ask[1]==0||best_ask[1]==100000){p=p;} else { p=(best_bid[1]+best_ask[1])/2;}
      for (i=0;i<m;i++)
        { if (order[i][3]!=0) un_trade=un_trade+1;
        }
        for(i=0;i<n;i++){
             w[i] = c[i] + p * s[i];//formulation for wealth
             if (w[i] > wmax)//find the max wealth
             {wmax = w[i];}
            outfile << w [i] << ' ';}
        average_guru=w[guru];
   for (i=0;i<n;i++) {
        if(a2[i]==guru) {total_imitator=total_imitator+w[i];imitator=imitator+1;/*cout<<"imitator"<<total_imitator<<' '<<imitator<<endl;*/}
        else {total_noise=total_noise+w[i];noise=noise+1;/*cout<<"noise"<<total_noise<<' '<<noise<<endl;*/}
    }
average_imitator=total_imitator/imitator; average_noise=(total_noise-average_guru)/(noise-1);


      for (i=0;i<m;i++)
      {   if ((order[i][3]==1)&&(t==0)) volumebuy=volumebuy+order[i][2]; //cout<<"t=0 volumebuy" <<volumebuy<<endl;
          if ((order[i][3]==-1)&&(t==0)) volumesell=volumesell+order[i][2]; //cout<<"t=0 volumesell" <<volumesell<<endl;
      }
      //cout<<"t=0 volumebuy" <<volumebuy<<endl;
      //cout<<"t=0 volumesell" <<volumesell<<endl;
      average_expected=total/n; mkt_buy=trades/sumbuy; mkt_sell=trades/sumsell; average_guru_order=order_volume[guru]; //cout<<guru<<endl;cout<<order_volume[guru]<<endl;
      average_rest_order=0;
      for (i=0;i<n;i++)
      {
          average_rest_order=average_rest_order+order_volume[i]; //cout<<average_rest_order<<' ';
      }
      average_rest_order=(average_rest_order-average_guru_order)/(n-1);

      if (t==1) {
        for (i=0;i<n;i++){for(j=0;j<4;j++)
                            {order[i][j]=0;/*cout<<order[i][j];*/}
        /*cout<<endl;  */               }
         for (i1=n;i1<m;i1++) {
        if((order[i1][3]==1) && (best_bid[1]<order[i1][1])) {
        best_bid[0]=order[i1][0]; best_bid[1]=order[i1][1]; best_bid[2]=order[i1][2]; best_bid[3]=order[i1][3];best_bid[4]=i1;
        }// best ask updated,need to compare them again
        }
         for (i1=n;i1<m;i1++) {
         if((order[i1][3]==-1) && (best_ask[1]>order[i1][1])) {
         best_ask[0]=order[i1][0]; best_ask[1]=order[i1][1]; best_ask[2]=order[i1][2]; best_ask[3]=order[i1][3]; best_ask[4]=i1;
         }
         }
         /*cout<<"best bid and ask" << best_bid[1]<<' '<<best_bid[4]<<' '<<best_ask[1]<<' '<<best_ask[4]<<endl;*/
      }
      tp1=tp1+1; if (tp1==2) {tp1=0;}
      outfile<< "stock and cash updating after the trading"<<endl;
      for (i=0;i<n;i++) {
        outfile<<' '<<s[i]<<' '<<c[i]<<endl;;
      }
      sumguruw=sumguruw+average_guru; sumimitatorw=sumimitatorw+average_imitator; sumnoisew=sumnoisew+average_noise; sgo=sgo+average_guru_order;
      sro=sro+average_rest_order;
      //cout << sumguruw<<' '<<sumimitatorw<<' '<<sumnoisew<<' '<<sgo<<' '<<sro<<endl;
     //tp1=tp1+1; //if (tp1==h) {tp1=0;}
     /* if (t>=(h-1)) {
         for (i=(tp1*n);i<((tp1+1)*n);i++) {n3=order[i][0];
           if(order[i][3]==-1){s[n3]=s[n3]-(order[i][2]*order[i][3]);}
           if (order[i][3]==1) {c[n3]=c[n3]+(order[i][1]*order[i][2]*order[i][3]);}
           order[i][0]=order[i][1]=order[i][2]=order[i][3]=order[i][4]=0;
      }
      }
      outfile<<"stock and cash updating after the releasing"<<endl;
      for (i=0;i<n;i++) {
        outfile<<' '<<i<<' '<<s[i]<<' '<<c[i]<<endl;
      }
    /*for(r=0;r < n;r++){
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

    for (r=0;r<n;r++) {for (i=0;i<6;i++){if (order[r][5]==0||order[r][5]==-1) buy[r][i]=order[r][i];
                     if (order [r][5]==1) sell[r][i]=order[r][i];}
                     for (n5=0;n5<=r-1;n5++){
                     for (j=0;j<=r-1;j++) {
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
                      //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << buy[r][i]<<' '; cout<<endl;}
                     //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << sell[r][i]<<' '; cout<<endl;}
                     outfile<<("after trading order book:sequence,price,amount,stock,cash,buy/sell\n");
                        for (r=0;r<m;r++){for (i=0;i<6;i++){if(r<n){order1[r][i]=buy[r][i];if(order1[r][0]==-1){order1[r][0]=10000;}}//if(buy[0][1]>ask) {ask=buy[0][1];}
                        if(r>=n){order1[r][i]=sell[r-n][i];if(order1[r][0]==-1){order1[r][0]=10000;}/*if(sell[0][1]<bid){bid=sell[0][1];}*///}}}

                       /* for (r=0;r<m;r++) {for (n5=0;n5<=r-1;n5++){for (j=0;j<=r-1;j++) {
                        if (order1[j][0] > order1[j+1][0])
                        {hold=order1[j][0]; order1[j][0]=order1[j+1][0];order1[j+1][0]=hold;
                        hold=order1[j][1]; order1[j][1]=order1[j+1][1];order1[j+1][1]=hold;
                        hold=order1[j][2]; order1[j][2]=order1[j+1][2];order1[j+1][2]=hold;
                        hold=order1[j][3]; order1[j][3]=order1[j+1][3];order1[j+1][3]=hold;
                        hold=order1[j][4]; order1[j][4]=order1[j+1][4];order1[j+1][4]=hold;
                        hold=order1[j][5]; order1[j][5]=order1[j+1][5];order1[j+1][5]=hold;}
                        }}}
                        for (r=0;r<n;r++){for (i=0;i<6;i++) outfile <<order1[r][i]<<' ';outfile<< endl;}
                         // for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << buy[r][i]<<' '; cout<<endl;}
                       //for (r=0;r<n;r++) {for (i=0;i<6;i++) cout << sell[r][i]<<' '; cout<<endl;}
                         outfile<<("cash, stock and bid/ask price\n");
                       for (i=0;i<n;i++)  {c[i]=order1[i][4];s[i]=order1[i][3];outfile<<c[i]<<setw(15)<<s[i]<<setw(15);}
                       outfile<<bid<<' '<<ask <<' '<<endl;*/
                       //if (mbn==0||msn==0) {p3= (bid+ask)/2;} else {p3=(bid1+ask1)/2;} }
                       //if (sumbuy==0) {p=ask;} if (sumsell==0) {p=bid;} if(bid2==-1||ask2==pmax+301||bid2==0||ask2==0) {p=(bid1+ask1)/2;} else {p=(bid2+ask2)/2;}
                        //cout <<best_bid[1]<<' '<<best_ask[1]<<endl;
                        pt[t]=p; totalprice[t]=totalprice[t]+ pt[t]; squareprice[t]=squareprice[t]+(pt[t]*pt[t]);//cout << best_bid[1]<<' '<<best_ask[1]<<' '<<pt[t] << endl;//after t gets bigger than 0, the old outgoing matrix
                       //is the previously final matrix, record all price in matrix pt accordingly
                        outfile << ' ' << t << endl;
                       outfile <<' '<< pt[t] << endl;//record the t and pt in material file
                       if(detail==0) {ofstream outfile;
                       outfile.open("guru_life.xls",ios::out|ios::app);
                       outfile <<guru_life<<' '<<guru_lifetime[t] <<endl;
                       outfile.close();
                       outfile.open("price.xls",ios::out|ios::app);
                       outfile <<p <<endl;
                       outfile.close();
                       outfile.open("expected price.xls",ios::out|ios::app);
                       outfile <<total<<endl;
                      outfile.close();
                      outfile.open("trades.xls",ios::out|ios::app);
                       outfile <<trades <<endl;
                       outfile.close();
                       outfile.open("total buyers.xls",ios::out|ios::app);
                       outfile <<sumbuy<<endl;
                     outfile.close();
                    outfile.open("total sellers.xls",ios::out|ios::app);
                       outfile <<sumsell <<endl;
                     outfile.close();
                     outfile.open("best bid.xls",ios::out|ios::app);
                       outfile <<best_bid[1]<<endl;
                     outfile.close();
                      outfile.open("best ask.xls",ios::out|ios::app);
                       outfile <<best_ask[1]<<endl;
                     outfile.close();
                    /* outfile.open("market buyers.xls",ios::out|ios::app);
                       outfile<<mbn<<endl;
                     outfile.close();
                     outfile.open("market sellers.xls",ios::out|ios::app);
                       outfile <<msn<<endl;
                     outfile.close();*/
                     // outfile.open("variance.xls",ios::out|ios::app);
                     //  outfile <<variance<<endl;
                    // outfile.close();

            //outfile.open("guru_life.xls",ios::out|ios::app);
            //for(i=0;i<T;i++)
            //{
            //outfile <<guru_life<<' '<<guru_lifetime[i] <<endl;
            //}
            //outfile.close();

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
		 outfile.open("sequence.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{
        outfile<< sequence[i] <<endl;
	}
	outfile.close();
		 outfile.open("reference price.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{n3=sequence[i];
        outfile<< pricet[n3] <<endl;
	}
	outfile.close();
		 outfile.open("return.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{n3=sequence[i];
        outfile<< ret[n3] <<endl;
	}
	outfile.close();
		 outfile.open("expected price for each agent.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{   n3=sequence[i];
        outfile<< p1[n3] <<endl;
	}
	outfile.close();
	 outfile.open("order size.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{n3=sequence[i];
        outfile<< amount[n3] <<endl;
	}
	outfile.close();
		 outfile.open("order price.xls",ios::out|ios::app);
                     	for(i=0;i<n;i++)
	{ n3=sequence[i];
        outfile<< p2[n3] <<endl;
	}
	outfile.close();
    outfile.open("trading.txt",ios::out|ios::app);
        for(i=0;i<n;i++)
	{   n3=sequence[i];
        outfile<<setprecision(4)<<setw(6)<<n3<<setprecision(6)<<setw(8)<<pricet[n3] <<setprecision(6)<<setw(8)<<p1[n3]<<setprecision(6)<<setw(8)<<p2[n3]<<setprecision(6)<<setw(8)<<amount[n3]
        <<setprecision(6)<<setw(8)<<trade_state[n3]<<setprecision(4)<<setw(8)<<order_volume[n3]<<setprecision(4)<<setw(8)<<trades<<setw(8)<<un_trade<<endl;
	}
	//cout<<sumbuy<<' '<<n<<endl;

	//cout<<herding<<endl;
	outfile.close();}
	herding=sumbuy/n; sim_herding=sim_herding+herding; //cout<<herding<<' '<<sim_herding<<endl;
	if(detail==1){ofstream outfile; outfile.open("result.txt",ios::out|ios::app);


        outfile<<setprecision(0)<<setw(6)<<guru<<setprecision(2)<<setw(8)<<guru_link <<setprecision(3)<<setw(8)<<guru_fitness<<setprecision(8)<<setw(11)<<average_guru<<setprecision(8)<<setw(11)<<average_imitator<<setprecision(8)<<setw(11)<<average_noise
        <<setprecision(6)<<setw(10)<<p<<setprecision(6)<<setw(10)<<average_expected<<setprecision(2)<<setw(8)<<herding<<setprecision(3)<<setw(11)<<average_guru_order<<setprecision(3)<<setw(11)<<average_rest_order
        <<endl;

	outfile.close();}

                        }//t1++;}}

       average_gurulife=T/number_guru; //cout << number_guru<<' '<<average_gurulife<<endl;
       sumguruw=sumguruw/T; sumimitatorw=sumimitatorw/T; sumnoisew=sumnoisew/T; sgo=sgo/T; sro=sro/T;

      number_guru=0; for (i=0;i<T;i++) {//cout <<i<<' '<<totalprice[i]<<' '<<squareprice[i]<<endl;
      }
      volumebuy=tradevolume/volumebuy; volumesell=tradevolume/volumesell; //cout<<simulation<<' '<< volumebuy<<' '<<volumesell<<endl;
      sim_herding=sim_herding/T; //if (simulation==(run-1)) {sim_herding=sim_herding/run;} //cout<<simulation<<' '<<sim_herding<<endl;
      if(detail==1) {ofstream outfile; outfile.open("result(each run).txt",ios::out|ios::app);

        outfile<<setprecision(4)<<setw(8)<<average_gurulife<<setprecision(3)<<setw(8)<<volumebuy<<setprecision(3)<<setw(8)<<volumesell<<setprecision(3)<<setw(8)<<sim_herding
        <<setprecision(8)<<setw(11)<<sumguruw<<setprecision(8)<<setw(11)<<sumimitatorw<<setprecision(8)<<setw(11)<<sumnoisew<<setprecision(3)<<setw(8)<<sgo
        <<setprecision(3)<<setw(8)<<sro<<endl;

      outfile.close();}

      if(simulation==(run-1)) { for(i=0;i<T;i++){
      totalprice[i]=totalprice[i]/run; squareprice[i]=sqrt((squareprice[i]/run)-(totalprice[i]*totalprice[i]));
      //cout <<totalprice[i]<<' '<<squareprice[i]<<endl;
      if (detail==1) {ofstream outfile;
       outfile.open(" mean and standard deviation.txt",ios::out|ios::app);

        outfile<< setprecision(6)<<setw(8)<< totalprice[i]<<' '<<setprecision(6)<<setw(8)<<squareprice[i]<<endl;

	outfile.close();}
      }}
      }
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



