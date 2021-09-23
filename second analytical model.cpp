#include <iostream>
#include <cstdlib>
//#include <ctime>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <fstream>
//#define n 500
#define n1 20
#define n2 100
//#define run 1
//#define detail 2
using namespace std;

int main()
{
    ifstream inFile ("parameter.txt", ios::in);
    int count=0;
    string parameterName;
    double parameterValue;
    double *par=new double [n1]();

    if(!inFile)
    {
        cerr<<"Error: File could not be opened" <<endl;

    }
    cout<<"parameter name, value" <<endl;

    while(!inFile.eof())
    {

        inFile>>parameterName>>parameterValue;
        cout<<setw(3)<<count<<setw(15)<<parameterName<<setw(6)<<parameterValue<<endl;
        par[count]=parameterValue;
        count++;
    }
    inFile.close();
    int n,run,T,detail;
    n=par[0]; run=par[15]; T=par[14]; detail=par[16];
    double beta,tc,pf,tf,alpha1,alpha2,v1;
    double alpha3,r,v2,R,s,mu,deviation;
    double opinion,price_trend,expu1,exp_u1,profit_cha, profit_fun,expu21,exp_u21,expu22,exp_u22,pro_p,pro_n,pro_f,current;
    v1=par[1]; v2=par[2]; beta=par[3]; tc=par[4]; tf=par[5]; alpha1=par[6]; alpha2=par[7]; alpha3=par[8];
    pf=par[9]; r=par[10]; R=par[11]; s=par[12]; deviation=par[13];
    int simulation,temp_agent;
    for(simulation=(run-1);simulation<run;simulation++) { cout<<"simulation"<<simulation<<endl;
    int *a=new int[n];
    int i,k,mi=4,j=1,timeinterval=100;
    int t1,treal=0,t100=0;
    double tinterval=0;
    double k1,nf=0,nc=0,np=0,nn=0,z,x;
    double t=0.01;
    //double alpha3=0.5,r=0.004,v2=2,R=0.0004,s=0.75;
    double stc,gamma;
    double dp,u1,p,piepn,pienp,p1,dif=0;
    double u21,u22,piepf,piefp,pienf,piefn,pief,pien,piep;
    double pexp_oppe=0, pexp_opf=0, pexp_peop=0, pexp_pef=0, pexp_fop=0, pexp_fpe=0, preal_oppe=0;
    int preal_opf=0, preal_peop=0, preal_pef=0, preal_fop=0, preal_fpe=0;
    double X,edc,edf,ed,pieup=0, piedown=0;
    double *last_t=new double [100]();
    double *last_p=new double [100]();
    static double V1, V2, S;
    static int phase = 0;
    unsigned int seed=simulation; //set the seed equals to 1
    srand(seed); //using the seed to generate random number
    pien=0.1; piep=0.1; pief=pien+piep;
    stc=tc/n; gamma=tf/n;
    //cout<<k<<" "<<k1<<endl;
    for (i=0;i<100;i++){

      last_p[i]=0; last_t[i]=-1;
    }

    for (i=0;i<n;i++)
    {   k1=((double)rand())/RAND_MAX;
        //cout<<"k1"<<k1<<" "<<i<<" ";
        if (k1<piep)
        {
           a[i]=1;
           //cout<<"kp"<<" "<<endl;
        }
        if (k1>=piep && k1<pief)
         {
           a[i]=-1;
           //cout<<"kn"<<" "<<endl;
         }
        if (k1>=pief)
        {
           a[i]=0;
           //cout<<"kf"<<" "<<endl;
        }
    }

    //k =((double)rand())/RAND_MAX;
    //if (k >= proi[i]) {a2 [i]= a [i];}
    //if (k < proi[i]){ a2[i]=a1[i];}
    for (i=0;i<n;i++)
    {
        if(a[i]==-1)
        {
           nn=nn+1; nc=nc+1;
        }
          if(a[i]==1)
        {
           np=np+1; nc=nc+1;
        }
          if(a[i]==0)
        {
           nf=nf+1;
        }
    }
    p=pf; dp=0;
    if (detail==0||detail==1) {ofstream outfile("material.txt",ios::out|ios::app);
    outfile<<("agents' category\n");
    outfile<<nf<<" "<<nc<<" "<<np<<" "<<nn<<" "<<endl;}

    for (t1=0;t1<T;t1++){ //cout<<"t"<<" "<<t1<<" "<<treal<<endl;
    //cout<< dp<<endl;
    //if (treal>=17236&&treal<17237) {detail=3;} if (treal>=17237) {detail=2;} //change the display of time period
    if (detail==0||detail==1) {ofstream outfile("material.txt",ios::out|ios::app);
    outfile<<("time\n");
    outfile<<"t"<<" "<<t1<<" "<<treal<<endl;}
    z=nc/n; x=(np-nn)/nc;
    if (detail==0||detail==1) {ofstream outfile("material.txt",ios::out|ios::app);
    outfile<<("agents' category\n");
    outfile<<nf<<" "<<nc<<" "<<np<<" "<<nn<<" "<<endl;}
    //cout<<nf<<" "<<nc<<" "<<np<<" "<<nn<<" "<<z<<" "<<x<<endl;
    /*dp1=0; dp5=0;
     for(i=0;i<n1;i++)
     {
        dp1=dp1+pdif[i]; //cout<<i<<" "<<pdif[i]<<" "<<dp<<" ";
     }
     for(i=0;i<n2;i++)
     {
        dp5=dp5+pdif5[i]; //cout<<i<<" "<<pdif[i]<<" "<<dp<<" ";
     }*/
     //cout<<endl;
     //dp=0;
    opinion=alpha1*x; price_trend=alpha2*dp/v1;
     u1=opinion+price_trend;
     //u1=alpha1*x+(alpha2*dp/v1);
     expu1=exp(u1); exp_u1=exp(-u1);
     //piepn=v1*(nc/n*exp(u1))*t; pienp=v1*(nc/n*exp(-u1))*t;
     piepn=0.0001*np+0.0002; pienp=0.0001*nn+0.0002;
    if (detail==0||detail==1) {ofstream outfile("material.txt",ios::out|ios::app);
    outfile<<("parameters\n");
    outfile<<"z"<<" "<<z<<" "<<"x"<<" "<<x<<" "<<"dp"<<" "<<dp<<" "<<"u1"<<" "<<u1<<"piepn"<<" "<<piepn<<" "<<"pienp"<<" "<<pienp<<" "<<endl;}
    current=p; //dp=0;
    profit_cha=(( r+(dp/v2) )/ p) - R; profit_fun=s* abs((pf-p)/p);
    u21=alpha3*(profit_cha-profit_fun);
    //u21=alpha3*((( r+(dp/v2) )/ p) - R - (s* abs((pf-p)/p)));
    //u21=(( r+(dp/v2) )/ p); u21=u21-R; cout<<u21<<" "<<R<<endl;
    //u21=((( r+(dp/v2) )/ p) - R); cout<<u21<<endl;
    //u21=(r+dp/v2)/p; cout<< u21; u21=u21-R; cout<< u21<<endl; u21=u21-s*(abs(pf-p)/p); cout<< u21<<endl; u21=u21*alpha3; cout<< u21<<endl;
    u22=alpha3*(-profit_cha-profit_fun);
    //u22=alpha3*(R-((r+(dp/v2))/p)-(s*abs((pf-p)/p)));
    //u21=0; u22=0;
    expu21=exp(u21); exp_u21=exp(-u21); expu22=exp(u22); exp_u22=exp(-u22); pro_p=np/n; pro_n=nn/n; pro_f=nf/n;
    //piepf=v2*(np/n*exp(u21))*t;
    //piefp=v2*(nf/n*exp(-u21))*t;
    //pienf=v2*(nn/n*exp(u22))*t;
    //piefn=v2*(nf/n*exp(-u22))*t;
    piepf=0.00004*np+0.000005;
    piefp=0.00004*nf+0.0001;
    pienf=0.00004*nn+0.000005;
    piefn=0.00004*nf+0.0001;
    pief=piepf+pienf;
    pien=piepn+piefn;
    piep=pienp+piefp;

    if (detail==0||detail==1) {ofstream outfile("material.txt",ios::out|ios::app);
            outfile<<("parameters(u21,u22,piepf,piefp,pienf,piefn,pief)\n");
            outfile<<u21<<" "<<u22<<" "<<piepf<<" "<<piefp<<" "<<pienf<<" "<<piefn<<" "<<pief<<endl;}
    //cout<<endl<<u21<<" "<<u22<<" "<<piepf<<" "<<piefp<<" "<<pienf<<" "<<piefn<<endl;
    i=0;
    for (i=0;i<n;i++)
    {   temp_agent=a[i];
        if(temp_agent==1&&np>=mi)
        {
            //piefp=v2*(np/n*exp(-u21))*t;
            k1=((double)rand())/RAND_MAX; //cout<<k1<<" k1 "<<piefp<<" ";
            /*if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
            outfile<<("parameters\n");
            outfile<<"k1"<<" "<<k1<<" "<<"piefp"<<" "<<piefp<<endl;}*/
             if (detail==3) {ofstream outfile("to_pe.txt",ios::out|ios::app);
                 outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<pienp<<setprecision(4)<<setw(20)<<"np"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<endl;}
            if (detail==3) {ofstream outfile("to_nf.txt",ios::out|ios::app);
              outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<piep<<setprecision(4)<<setw(20)<<"np"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<setprecision(4)<<setw(20)<<pienp<<endl;}
              pexp_oppe=pexp_oppe+pienp; pexp_opf=pexp_opf+piefp;
            if (k1 < pienp)
            {   a[i]=-1; np=np-1; nn=nn+1; z=nc/n; x=(np-nn)/nc; //u1=alpha1*x+(alpha2*dp/v1);
            preal_oppe=preal_oppe+1;
            if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
                outfile<<("parameters\n");
                outfile<<"k1"<<" "<<k1<<" "<<"pienp"<<" "<<pienp<<endl;
                outfile<<("update\n");
                outfile<<"np"<<" "<<np<<" "<<"nn"<<" "<<nn<<" "<<"z"<<" "<<z<<" "<<"x"<<" "<<x<<" "<<"dp"<<" "<<dp<<" "<<"u1"<<" "<<u1<<endl;}
                /*if (detail==3) {ofstream outfile("equation1-random.txt",ios::out|ios::app);
                 outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<pienp<<setprecision(4)<<setw(20)<<"np"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<endl;}*/
            //cout<<a[i]<<" "<<np<<" "<<nn<<" "<<dp<<" "<<u1<<" "<<endl;
            }
            if (k1 >= pienp&& k1<piep)
            {   a[i]=0; np=np-1; nf=nf+1; nc=nc-1; z=nc/n; x=(np-nn)/nc; //u21=alpha3*((r+(dp/v2))/p-R-(s*abs((pf-p)/p))); u22=alpha3*(R-((r+(dp/v2))/p)-(s*abs((pf-p)/p)));
                //cout<<a[i]<<" "<<nf<<" "<<np<<" "<<dp<<" "<<u21<<" "<<u22<<" ";
                 preal_opf=preal_opf+1;
                if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
                outfile<<("parameters\n");
                outfile<<"k1"<<" "<<k1<<" "<<"piefp"<<" "<<piefp<<endl;
                outfile<<("update\n");
                outfile<<"np"<<" "<<np<<" "<<"nf"<<" "<<nf<<" "<<"z"<<" "<<z<<" "<<"x"<<" "<<x<<" "<<"dp"<<" "<<dp<<" "<<"u21"<<" "<<u21<<" "<<"u22"<<" "<<u22<<endl;}
                /*if (detail==3) {ofstream outfile("equation2-random.txt",ios::out|ios::app);
              outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<piefp<<setprecision(4)<<setw(20)<<"np"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<endl;}*/
            }
        }

       if(temp_agent==0&&nf>=mi)
        {
            //piepf=v2*(np/n*exp(u21))*t;
            k1=((double)rand())/RAND_MAX;
            /*if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
            outfile<<("parameters\n");
            outfile<<"k1"<<" "<<k1<<" "<<"piepf"<<" "<<piepf<<endl;} //cout<<k1<<" k1 "<<piepf<<" ";*/
            if (detail==3) {ofstream outfile("to_op.txt",ios::out|ios::app);
                outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<piepf<<setprecision(4)<<setw(20)<<"nf"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<endl;}
            if (detail==3) {ofstream outfile("to_pe.txt",ios::out|ios::app);
                outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<pief<<setprecision(4)<<setw(20)<<"nf"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<setprecision(4)<<setw(20)<<piepf<<endl;}
              pexp_fop=pexp_fop+piepf; pexp_fpe=pexp_fpe+pienf;
            if (k1 < piepf)
            {   a[i]=1; np=np+1; nf=nf-1; nc=nc+1; z=nc/n; x=(np-nn)/nc; //u21=alpha3*((r+(dp/v2))/p-R-(s*abs((pf-p)/p))); u22=alpha3*(R-((r+(dp/v2))/p)-(s*abs((pf-p)/p)));
                //cout<<a[i]<<" "<<np<<" "<<nf<<" "<<dp<<" "<<u21<<" "<<u22<<" ";
                preal_fop=preal_fop+1;
                if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
                outfile<<("parameters\n");
                outfile<<"k1"<<" "<<k1<<" "<<"piepf"<<" "<<piepf<<endl;
                outfile<<("update\n");
                outfile<<"np"<<" "<<np<<" "<<"nf"<<" "<<nf<<" "<<"z"<<" "<<z<<" "<<"x"<<" "<<x<<" "<<"dp"<<" "<<dp<<" "<<"u21"<<" "<<u21<<" "<<"u22"<<" "<<u22<<endl;}
                /*if (detail==3) {ofstream outfile("equation2-random.txt",ios::out|ios::app);
                outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<pief<<setprecision(4)<<setw(20)<<"nf"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<setprecision(4)<<setw(20)<<piepf<<endl;}*/
            }
            if (k1>= piepf && k1 <pief)
            {  a[i]=-1; nf=nf-1; nn=nn+1; nc=nc+1; z=nc/n; x=(np-nn)/nc; //u21=alpha3*((r+(dp/v2))/p-R-(s*abs((pf-p)/p))); u22=alpha3*(R-((r+(dp/v2))/p)-(s*abs((pf-p)/p)));
                //cout<<a[i]<<" "<<nf<<" "<<nn<<" "<<dp<<" "<<u21<<" "<<u22<<" "<<endl;
                preal_fpe=preal_fpe+1;
                if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
                outfile<<("parameters\n");
                outfile<<"k1"<<" "<<k1<<" "<<"piepf"<<" "<<piepf<<" "<<"pief"<<" "<<pief<<endl;
                outfile<<("update\n");
                outfile<<"nn"<<" "<<nn<<" "<<"nf"<<" "<<nf<<" "<<"z"<<" "<<z<<" "<<"x"<<" "<<x<<" "<<"dp"<<" "<<dp<<" "<<"u21"<<" "<<u21<<" "<<"u22"<<" "<<u22<<endl;}
                /* if (detail==3) {ofstream outfile("equation2-random.txt",ios::out|ios::app);
                outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<pief<<setprecision(4)<<setw(20)<<"nf"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<setprecision(4)<<setw(20)<<piepf<<endl;}*/
            }

           //cout<<endl<<i<<" ";
        }


        if(temp_agent==-1&&nn>=mi)
        {
            //piefn=v2*(nn/n*exp(-u22))*t;
            k1=((double)rand())/RAND_MAX; //cout<<k1<<" k1 "<<piefn<<" ";
            /*if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
            outfile<<("parameters\n");
            outfile<<"k1"<<" "<<k1<<" "<<"piefn"<<" "<<piefn<<endl;}*/
              if (detail==3) {ofstream outfile("to_op.txt",ios::out|ios::app);
                outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<piepn<<setprecision(4)<<setw(20)<<"nn"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<endl;}
            if (detail==3) {ofstream outfile("to_nf.txt",ios::out|ios::app);
                 outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<pien<<setprecision(4)<<setw(20)<<"nn"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<setprecision(4)<<setw(20)<<piepn<<endl;}
              pexp_peop=pexp_peop+piepn; pexp_pef=pexp_pef+piefn;
            if (k1 < piepn)
            {   a[i]=1; np=np+1; nn=nn-1; z=nc/n; x=(np-nn)/nc; //u1=alpha1*x+(alpha2*dp/v1);
                preal_peop=preal_peop+1;
                if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
                outfile<<("parameters\n");
                outfile<<"k1"<<" "<<k1<<" "<<"piepn"<<" "<<piepn<<endl;
                outfile<<("update\n");
                outfile<<"np"<<" "<<np<<" "<<"nn"<<" "<<nn<<" "<<"z"<<" "<<z<<" "<<"x"<<" "<<x<<" "<<"dp"<<" "<<dp<<" "<<"u1"<<" "<<u1<<endl;}
               /* if (detail==3) {ofstream outfile("equation1-random.txt",ios::out|ios::app);
                outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<piepn<<setprecision(4)<<setw(20)<<"nn"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<endl;}*/
            //cout<<a[i]<<" "<<np<<" "<<nn<<" "<<dp<<" "<<u1<<" ";
            }
            if (k1 >=piepn&&k1 < pien)
            {   a[i]=0; nf=nf+1; nn=nn-1; nc=nc-1; z=nc/n; x=(np-nn)/nc; //u21=alpha3*((r+(dp/v2))/p-R-(s*abs((pf-p)/p))); u22=alpha3*(R-((r+(dp/v2))/p)-(s*abs((pf-p)/p)));
                //cout<<a[i]<<" "<<nn<<" "<<nf<<" "<<dp<<" "<<u21<<" "<<u22<<" "<<endl;
                preal_pef=preal_pef+1;
                if (detail==0) {ofstream outfile("material.txt",ios::out|ios::app);
                outfile<<("parameters\n");
                outfile<<"k1"<<" "<<k1<<" "<<"piefn"<<" "<<piefn<<endl;
                outfile<<("update\n");
                outfile<<"nn"<<" "<<nn<<" "<<"nf"<<" "<<nf<<" "<<"z"<<" "<<z<<" "<<"x"<<" "<<x<<" "<<"dp"<<" "<<dp<<" "<<"u21"<<" "<<u21<<" "<<"u22"<<" "<<u22<<endl;}
                /*if (detail==3) {ofstream outfile("equation2-random.txt",ios::out|ios::app);
                 outfile<<setprecision(4)<<setw(20)<<i<<setprecision(4)<<setw(20)<<k1
              <<setprecision(4)<<setw(20)<<piefn<<setprecision(4)<<setw(20)<<"nn"
              <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<endl;}*/
            }

             //cout<<endl<<i<<" ";
        }
    }
    for (i=0;i<n;i++)
    {
        //cout<<a[i]<<" ";
    }
    edc=(np-nn)*stc;
    edf=nf*gamma*(pf-p);
    ed=edc+edf;
    if(phase == 0)
    {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
            } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    }
    else
    {
        X = V2 * sqrt(-2 * log(S) / S);
        phase = 1 - phase;
    }
    mu=deviation*X;
    if (detail==0||detail==1) {ofstream outfile("material.txt",ios::out|ios::app);
            outfile<<("parameters\n");
            outfile<<"mu"<<" "<<mu<<endl;}
    pieup=beta*(ed+mu);
    if(pieup<0)
    {
        pieup=0;
    }

    piedown=-(beta*(ed+mu));
    if(piedown<0)
    {
        piedown=0;
    }
    k1=((double)rand())/RAND_MAX; //cout<<k1<<" k1 "<<" ";
    if (detail==0||detail==1) {ofstream outfile("material.txt",ios::out|ios::app);
            outfile<<("parameters\n");
            outfile<<"k1"<<" "<<k1<<" "<<"pieup"<<" "<<pieup<<" "<<"piedown"<<" "<<piedown<<endl;}
    //cout<<"t100"<<t100<<endl;
    if (k1 <= pieup)
    {
        p=p+0.01; last_p[t100]=0.01; dif=dif+0.01;
    }
    if (k1 <= piedown)
    {
        p=p-0.01; last_p[t100]=-0.01; dif=dif+0.01;
    }
    if (pieup!=0&&k1>pieup)
    {
        last_p[t100]=0; dif=dif;
    }
    if (piedown!=0&&k1>piedown)
    {
        last_p[t100]=0; dif=dif;
    }

    //cout<<endl<<edc<<" "<<edf<<" "<<ed<<" "<<mu<<" "<<pieup<<" "<<piedown<<" "<<p<<endl;
    if(t1==0){p1=p;}
    if(detail==0||detail==1)
        {ofstream outfile;
         outfile.open("detailed price.xls",ios::out|ios::app);
         outfile<< p <<endl;
         outfile.close();

         /*outfile.open("detailed loop.txt",ios::out|ios::app);
         outfile<<setprecision(4)<<setw(20)<< dp <<setprecision(4)<<setw(20)<< dp1<<setprecision(4)<<setw(20)<< dp5  <<setprecision(4)<<setw(20)<<dif<<setprecision(4)<<setw(20)<<j<<setprecision(4)<<setw(20)<<j1<<setprecision(4)<<setw(20)<<j5<<endl;
         outfile.close();*/
         //outfile.open("detailed chartists.xls",ios::out|ios::app);
         //outfile<< z <<endl;
         //outfile.close();
         //outfile.open("price1.xls",ios::out|ios::app);
         //outfile<< p1 <<endl;
         //outfile.close();
         /*outfile.open("j.xls",ios::out|ios::app);
         outfile<< j <<endl;
          outfile.open("detailed time.xls",ios::out|ios::app);
         outfile<< t <<endl;
         outfile.close();
          outfile.open("detailed fundamentalists.xls",ios::out|ios::app);
      outfile<< nf <<endl;
      outfile.close();
      outfile.open("detailed optimistic chartists.xls",ios::out|ios::app);
      outfile<< np <<endl;
      outfile.close();
      outfile.open("detailed pessimistic chartists.xls",ios::out|ios::app);
      outfile<< nn <<endl;
      outfile.close();*/
      outfile.open("detailed probability.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(20)<<piepn
      <<setprecision(4)<<setw(20)<<pienp<<setprecision(4)<<setw(20)<<piepf
      <<setprecision(4)<<setw(20)<<piefp<<setprecision(4)<<setw(20)<<pienf
      <<setprecision(4)<<setw(20)<<piefn<<setprecision(4)<<setw(20)<<pieup
      <<setprecision(4)<<setw(20)<<piedown<<endl;
      outfile.close();
      outfile.open("detailed equation1.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(20)<<x
      <<setprecision(4)<<setw(20)<<dp<<setprecision(4)<<setw(20)<<opinion
      <<setprecision(4)<<setw(20)<<price_trend<<setprecision(4)<<setw(20)<<u1
      <<setprecision(4)<<setw(20)<<z<<setprecision(4)<<setw(20)<<expu1
     <<setprecision(4)<<setw(20)<<exp_u1<<setprecision(4)<<setw(20)<<piepn<<setprecision(4)<<setw(20)<<pienp<<endl;
      outfile.close();
      outfile.open("detailed equation2.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(20)<<dp
      <<setprecision(4)<<setw(20)<<current<<setprecision(4)<<setw(20)<<profit_cha
      <<setprecision(4)<<setw(20)<<profit_fun<<setprecision(4)<<setw(20)<<u21
      <<setprecision(4)<<setw(20)<<u22<<setprecision(4)<<setw(20)<<expu21
      <<setprecision(4)<<setw(20)<<exp_u21<<setprecision(4)<<setw(20)<<expu22
      <<setprecision(4)<<setw(20)<<exp_u22<<setprecision(4)<<setw(20)<<pro_p
      <<setprecision(4)<<setw(20)<<pro_n<<setprecision(4)<<setw(20)<<pro_f
      <<setprecision(4)<<setw(20)<<piepf<<setprecision(4)<<setw(20)<<piefp
      <<setprecision(4)<<setw(20)<<pienf<<setprecision(4)<<setw(20)<<piefn<<endl;
      outfile.close();
       outfile.open("detailed equation3.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(12)<<edc
      <<setprecision(4)<<setw(12)<<edf<<setprecision(4)<<setw(12)<<ed
      <<setprecision(4)<<setw(12)<<mu<<setprecision(4)<<setw(12)<<pieup
      <<setprecision(4)<<setw(12)<<piedown<<endl;
      outfile.close();
      outfile.open("detailed information.txt",ios::out|ios::app);
      outfile<<setw(8)<<treal<<setw(12)<<p
      <<setw(12)<<z <<setw(12)<<nf
      <<setw(12)<<nc<<setw(12)<<np
      <<setw(12)<<nn<<setw(12)<<j<<endl;
      outfile.close();
        }
       if(detail==3)
        {ofstream outfile;
         outfile.open("detailed change.txt",ios::out|ios::app);
         outfile<<setprecision(4)<<setw(20)<< a[0] <<setprecision(4)<<setw(20)<<t1<<setprecision(4)<<setw(20)<<treal<<endl;
         outfile.close();}

    if(j==timeinterval)
     {//j=0;
      ofstream outfile;
      outfile.open("price.xls",ios::out|ios::app);
      outfile<< p <<endl;
      outfile.close();
      outfile.open("chartists.xls",ios::out|ios::app);
      outfile<< z <<endl;
      outfile.close();
       outfile.open("fundamentalists.xls",ios::out|ios::app);
      outfile<< nf <<endl;
      outfile.close();
      outfile.open("optimistic chartists.xls",ios::out|ios::app);
      outfile<< np <<endl;
      outfile.close();
      outfile.open("pessimistic chartists.xls",ios::out|ios::app);
      outfile<< nn <<endl;
      outfile.close();
       outfile.open("price fluctuation.xls",ios::out|ios::app);
      outfile<< dif <<endl;
      outfile.close();
      outfile.open("probability.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(20)<<piepn
      <<setprecision(4)<<setw(20)<<pienp<<setprecision(4)<<setw(20)<<piepf
      <<setprecision(4)<<setw(20)<<piefp<<setprecision(4)<<setw(20)<<pienf
      <<setprecision(4)<<setw(20)<<piefn<<setprecision(4)<<setw(20)<<pieup
      <<setprecision(4)<<setw(20)<<piedown<<endl;
      outfile.close();
      outfile.open("equation1.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(20)<<x
      <<setprecision(4)<<setw(20)<<dp<<setprecision(4)<<setw(20)<<opinion
      <<setprecision(4)<<setw(20)<<price_trend<<setprecision(4)<<setw(20)<<u1
      <<setprecision(4)<<setw(20)<<z<<setprecision(4)<<setw(20)<<expu1
     <<setprecision(4)<<setw(20)<<exp_u1<<setprecision(4)<<setw(20)<<piepn<<setprecision(4)<<setw(20)<<pienp<<endl;
      outfile.close();
      outfile.open("equation2.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(20)<<dp
      <<setprecision(4)<<setw(20)<<current<<setprecision(4)<<setw(20)<<profit_cha
      <<setprecision(4)<<setw(20)<<profit_fun<<setprecision(4)<<setw(20)<<u21
      <<setprecision(4)<<setw(20)<<u22<<setprecision(4)<<setw(20)<<expu21
      <<setprecision(4)<<setw(20)<<exp_u21<<setprecision(4)<<setw(20)<<expu22
      <<setprecision(4)<<setw(20)<<exp_u22<<setprecision(4)<<setw(20)<<pro_p
      <<setprecision(4)<<setw(20)<<pro_n<<setprecision(4)<<setw(20)<<pro_f
      <<setprecision(4)<<setw(20)<<piepf<<setprecision(4)<<setw(20)<<piefp
      <<setprecision(4)<<setw(20)<<pienf<<setprecision(4)<<setw(20)<<piefn<<endl;
      outfile.close();
       outfile.open("equation3.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(12)<<edc
      <<setprecision(4)<<setw(12)<<edf<<setprecision(4)<<setw(12)<<ed
      <<setprecision(4)<<setw(12)<<mu<<setprecision(4)<<setw(12)<<pieup
      <<setprecision(4)<<setw(12)<<piedown<<endl;
      outfile.close();
      outfile.open("information.txt",ios::out|ios::app);
      outfile<<setw(8)<<treal<<setw(12)<<p
      <<setw(12)<<z <<setw(12)<<nf
      <<setw(12)<<nc<<setw(12)<<np
      <<setw(12)<<nn<<setw(12)<<j<<endl;
      outfile.close();
      outfile.open("accumulative pro.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(20)<<pexp_oppe<<setprecision(4)<<setw(20)<<pexp_opf
      <<setprecision(4)<<setw(20)<<pexp_peop<<setprecision(4)<<setw(20)<<pexp_pef
      <<setprecision(4)<<setw(20)<<pexp_fop<<setprecision(4)<<setw(20)<<pexp_fpe
      <<setprecision(4)<<setw(20)<<preal_oppe<<setprecision(4)<<setw(20)<<preal_opf
      <<setprecision(4)<<setw(20)<<preal_peop<<setprecision(4)<<setw(20)<<preal_pef
      <<setprecision(4)<<setw(20)<<preal_fop<<setprecision(4)<<setw(20)<<preal_fpe<<endl;
      outfile.close();
     }
     /*dp1=0; dp5=0;
     for(i=0;i<n1;i++)
     {
        dp1=dp1+pdif[i]; //cout<<i<<" "<<pdif[i]<<" "<<dp<<" ";
     }
     for(i=0;i<n2;i++)
     {
        dp5=dp5+pdif5[i]; //cout<<i<<" "<<pdif[i]<<" "<<dp<<" ";
     }
     if (j==500)
     {
         dif=dif/5;
     }*/

     /*ofstream outfile;
      outfile.open("detail dif.txt",ios::out|ios::app);
      outfile<<setprecision(4)<<setw(20)<<dif<<setprecision(4)<<setw(20)<<timeinterval<<endl;
      outfile.close();*/
      last_t[t100]=tinterval;
      t100=t100+1;
      if (t100==100)
      {
          t100=0;
      }

      dp=0; k=0;
     for (i=0;i<100;i++){
       if  (last_t[i]>=(tinterval-0.2))
       {
           dp=dp+last_p[i]; k=k+1;
       }
    }


     if (k==0) {dp=0;}
     if (k!=0) {dp=dp/k;}

    if (j==timeinterval)
    {
        //dif=abs(p-p1);
        cout<<treal<<endl;
        pexp_oppe=0; pexp_opf=0; pexp_peop=0; pexp_pef=0; pexp_fop=0; pexp_fpe=0; preal_oppe=0; preal_opf=0; preal_peop=0; preal_pef=0; preal_fop=0; preal_fpe=0;
        if (dif>10)
        {
            t=0.002; timeinterval=500;
        }
        if(dif<=10)
        {
            t=0.01; timeinterval=100;
        }
        p1=p; j=0; treal=treal+1; dif=0;
    }
    tinterval=tinterval+t;
    //timeinterval=100; t=0.01;
    /*if (timeinterval==100)
    {
        dp=dp1;
    }
    if (timeinterval==500)
    {
        dp=dp5;
    }*/
    //cout<<j1<<" "<<j5<<" "<<dif<<endl;
     j=j+1;
     //j1=j1+1; j5=j5+1;
     //if(j1==20){j1=0;}
     //if(j5==100) {j5=0;}
     //cout<<j<<" "<<j1<<" ";


    //outfile.close();
    if(detail==0||detail==1)
        {ofstream outfile;
        outfile.open("agents.xls",ios::out|ios::app);
        for(i=0;i<n;i++)
	        {
             outfile<< a[i] <<endl;
	        }
        outfile.close();
        outfile.open("last time.xls",ios::out|ios::app);
        for(i=0;i<100;i++)
	        {
             outfile<< last_t[i]<<endl;
	        }
        outfile.close();
        outfile.open("last price.xls",ios::out|ios::app);
        for(i=0;i<100;i++)
	        {
             outfile<< last_p[i]<<endl;
	        }
        outfile.close();
      outfile.open("k.xls",ios::out|ios::app);
      outfile<< k <<endl;
      outfile.close();
        outfile.open("dp.xls",ios::out|ios::app);
      outfile<< dp <<endl;
      outfile.close();


        /*outfile.open("timeinterval.xls",ios::out|ios::app);
        outfile<< timeinterval <<endl;
        outfile.close();*/
        }
	//outfile.close();
    }
    }
}
