class Integral{
  public:
    //Construtores
    Integral(double aF(double), double x1, double x2, int n){
      iF1=aF;
      gaulegGer(x1,x2,n);
    }
    Integral(double aF(double,double,double), double x1, double x2, int n){
      iF3=aF;
      gaulegGer(x1,x2,n);
    }
    Integral(double aF(double,double,double,double), double x1, double x2, int n){
      iF4=aF;
      gaulegGer(x1,x2,n);
    }
    Integral(double aF(double,double,double,double,char [],double,double,int), double x1, double x2, int n){
      iF8=aF;
      gaulegGer(x1,x2,n);
    }
    Integral(double aF(double,double,double,double,double,char [],double,double,int), double x1, double x2, int n){
      iF9=aF;
      gaulegGer(x1,x2,n);
    }
    Integral(){};

    //Metodos membro
    double calcIntegral() const;
    double calcIntegral(double, double, int) const;
    double calcIntegral(double, double, double) const;
    double calcIntegral(double, double, double, char [],double, double, int, int) const;
    double calcIntegral(double, double, double, double, char [],double, double, int, int) const;
    void gaulegGer(double x_inf, double x_sup, int N);
    //x_inf,x_sup - limite inferior e superior da integral
    //N - divisoes da integral

    //membros


    private:
    int nDivisoes;
    double X_int[201],W_int[201];
    double (*iF1)(double);
    double (*iF3)(double,double,double);
    double (*iF4)(double,double,double,double);
    double (*iF8)(double,double,double,double,char [],double,double,int);
    double (*iF9)(double,double,double,double,double,char [],double,double,int);

};
/*******************************************************************/
  double Integral::calcIntegral() const{
    double INT=0.0;
    for(int i=1;i<=nDivisoes;i++)
    {
    INT += W_int[i]*iF1(X_int[i]);
    }

    return INT;
  }
  double Integral::calcIntegral(double X1, double X2, int i = 3) const{
    double INT=0.0;

    if (i==1){
      for(int i=1;i<=nDivisoes;i++){
        INT += W_int[i]*iF3(X_int[i],X1,X2);
      }
    }

    if (i==2){
      for(int i=1;i<=nDivisoes;i++){
        INT += W_int[i]*iF3(X1,X_int[i],X2);
      }
    }

    if (i==3){
      for(int i=1;i<=nDivisoes;i++){
        INT += W_int[i]*iF3(X1,X2,X_int[i]);
      }
    }
    return INT;
  }
  double Integral::calcIntegral(double X1, double X2, double X3) const{
    double INT=0.0;
    for(int i=1;i<=nDivisoes;i++)
    {
    INT += W_int[i]*iF4(X1,X2,X3,X_int[i]);
    }

    return INT;
  }
  double Integral::calcIntegral(double X1, double X2, double X3, char X4[], double X5, double X6, int X7, int i=3) const{
    double INT=0.0;

    //int size = sizeof(X4);
    //char IN[size];
    //for(int j=0; j < size ;j++) IN[j] = X4[j];

    if(i==2){
      for(int i=1;i<=nDivisoes;i++)
      {
        INT += W_int[i]*iF8(X1,X_int[i],X2,X3, X4 ,X5,X6,X7);
      }
    }else{printf("ERROR");}
    return INT;
  }
  double Integral::calcIntegral(double X1, double X2, double X3, double X4, char X5[], double X6, double X7, int X8, int i=3) const{
    double INT=0.0;

    //int size = sizeof(X4);
    //char IN[size];
    //for(int j=0; j < size ;j++) IN[j] = X4[j];

    if(i==5){
      for(int i=1;i<=nDivisoes;i++)
      {
        INT += W_int[i]*iF9(X1,X2,X3,X4,X_int[i], X5 ,X6,X7,X8);
      }
    }else{printf("ERROR");}
    return INT;
  }
/*******************************************************************/
  void Integral::gaulegGer(double x1, double x2, int N){

    double EPS = 3E-15;     //precisao relativa
    double PI = 3.1415926535897932384626433832795;
    double z1, z, xm, xl, pp, p3, p2, p1;
    int m;

    if(N>200){std::cout<<"O numero de divisoes deve ser <= 200; fazendo N="<< int(N=200) <<std::endl;}

    m = (N+1)/2;
    xm = 0.5*(x2+x1);
    xl = 0.5*(x2-x1);

    for(int i=1;i<=m;i++)
    {
        z = cos(PI*(i-0.25)/(N+0.5));

        do{
            p1=1.0;
            p2=0.0;
            for(int i=1;i<=N;i++)
            {
                p3=p2;
                p2=p1;
                p1=((2.0*i-1.0)*z*p2-(i-1.0)*p3)/i;
            }
            pp=N*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        }while (fabs(z-z1) > EPS);

        X_int[i] = xm-xl*z;
        X_int[N+1-i] = xm+xl*z;
        W_int[i] = 2.0*xl/((1.0-z*z)*pp*pp);
        W_int[N+1-i]=W_int[i];

    }
    nDivisoes = N;
    //for(int i=1;i<=N;i++){printf("x[%i]=%.15E;\tw[%i]=%.15E;\n",i,X_int[i],i,W_int[i]);}//teste
  }




class NunberArcs{
  public:
    //constructors
    NunberArcs(){}

    //methods
    double Zm(double mlim){
      return 0.23*(mlim - 20.6);
    }

    double Zstar(double zm){
      return zm/1.412;
    }

    double N(double zs, double mlim){//total surface density of sources, in units of arcmin^(-2)
      //double zstar = Zstar(Zm(mlim));
      //double zstar3 = zstar*zstar*zstar;
      //double zs2 = zs*zs;

      //return 45.0*pow(Zm(mlim),3.4)*zs2/zstar3*exp( -sqrt(zs2*zs/zstar3) );
      double E=2.718281828;
     // return (1.8316785562615905e9*pow(-20.6 + mlim,0.3999999999999999)*pow(zs,2))/pow(E,15.211088427357856*pow(zs/(-20.6 + mlim),1.5));

      return (8.316785562615905e8*pow(-20.6 + mlim,0.3999999999999999)*pow(zs,2))/pow(E,15.211088427357856*pow(zs/(-20.6 + mlim),1.5));

    }

    double Neff(double zs, double mlim, double mu){//total surface density of sources, in units of arcmin^(-2), considering the magnification efect
    if(std::isnan(N(zs,1.25*log10(mu)+mlim))) return 0.0;
    return N(zs,1.25*log10(mu)+mlim);
    }

};

class LensParameters{
  public:

    LensParameters(){
      Om0=0.4;
      Ol0=1.0-Om0;
    }

    double invE(double z){//inverse of function E;E=H/H_0
      double zm1 =z+1.0;
      double zm12=zm1*zm1;
      double zm13=zm1*zm12;
      return 1.0/sqrt(Om0*zm13 + Ol0);
    }

    double Om(double z){//Density of Matter
      double zm1 =z+1.0;
      double zm12=zm1*zm1;
      double zm13=zm1*zm12;
      double invEz = invE(z);
      double invEz2 = invEz*invEz;
      return Om0*zm13*invEz2;
    }

    double DELTAvir(double z){//virial density
      double cons=177.652879266;
      double Omzm1 = Om(z)-1.0;
      double Omzm12 = Omzm1*Omzm1;
      return cons + 82.0*Omzm1 - 39.0*Omzm12;
    }

    double conc(double M, double z){
      double c_0 = 5.97;
      double alc = -0.102;
      double M_st = 8.0;
      return c_0/(1.0+z)*pow(M/M_st,alc);
    }

    double g(double M, double z){
      double C = conc(M,z);
      return C*C/( log(1.0+C) - C/(1.0+C) );
    }

    double Kappa_s_incomplete(double zl, double M){
      return 7.35667E-4*g(M,zl)*pow(DELTAvir(zl)*Om0,2.0/3.0)*(1.0+zl)*(1.0+zl)*pow(invE(zl),-2.0/3.0)*pow(M,1.0/3.0);
    }

  private:
    double Om0,Ol0;


};

/*class SigmaInt{
  public:

    SigmaInt(int nlinhas, char file[]){

      for(int i=0;i<=500;i++) veck_s[i] = &k_s[i];
      for(int i=0;i<=500;i++) vecsig[i] = &sig[i];

        c2d(file,90,nlinhas,1,16,veck_s);
        c2d(file,90,nlinhas,37,53,vecsig);

    }

  private:

    double *veck_s[500];
    double k_s[500];

    double sig[500];
    double *vecsig[500];



    int c2d(char file[],int ncolunas,int nlinhas, int colunainicial, int colunafinal, double *vecout[])
    {

      FILE *input = fopen(file,"r");

      if (input==NULL){
        std::cout<<"Erro->Arquivo "<<file<<" nao existe"<<std::endl;
        return 1;
      }

      ncolunas+=2;

      char buffer[ncolunas];

      for(int i=1; i<=nlinhas ;i++){

        fgets (buffer , ncolunas , input);

        if(feof(input)){
          std::cout<<"Arquivo "<<file<<" possui "<<i-1<<" linhas"<<std::endl;
          break;
        }

        char X[colunafinal-colunainicial];
        for(int j=colunainicial; j<=colunafinal ;j++){
          X[j-colunainicial] = buffer[j-1];
        }

        //printf("%i-%E\n",i,atof(X));
        *vecout[i-1] = atof(X);
      }

      fclose(input);
      return 0;
    }



};*/
