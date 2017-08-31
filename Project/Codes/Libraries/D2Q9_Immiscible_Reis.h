//#include <sstream>
#include <string>
//#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
//const int Lx=200;
//const int Ly=200;

enum TipoCelda{Fluido,Pared,FuenteVelocidad,FuentePresion};

class LatticeBoltzmann{
    private:
    int Lx,Ly; //Dimensiones de la cuadricula
    double w[9]; //Pesos de la colisión BGK
    double B[9];  //Pesos del colisionador de color de Reis y Phillips
    int V[2][9];  //V[x=0,y=1][i]  , i=el vector velocidad
    //Variables del fluido marcado como R (RED)
    double ***fR; //funciones de distribución (old) de las partículas marcadas como R (RED)
    double ***fRnew; //funciones de distribución (new) de las partículas marcadas como R (RED)
    double TauR; //Tiempo de relajación del fluido marcado como R (RED)
    double RHOR; //Densidad del fluido marcado como R (RED)
    double OmegaR; //Inverso del tiempo de relajación del fluido marcado como R (RED)
    double AlphaR; //Contiene los parámetros que determinan la velocidad del sonido del fluido marcado como R (RED)
    double SurfaceTensionR; //Contiene el parámetro que determinan el valor de la tensión superficial del fluido marcado como R (RED)
    double DensidadR; //Contiene la densidad total del fluido R (RED)
    //Variables del fluido marcado como B (BLUE)
    double ***fB; //funciones de distribución (old) de las partículas marcadas como B (BLUE)
    double ***fBnew; //funciones de distribución (new) de las partículas marcadas como B (BLUE)
    double TauB; //Tiempo de relajación del fluido marcado como B (BLUE)
    double OmegaB; //Inverso del tiempo de relajación del fluido marcado como B (BLUE)
    double RHOB; //Densidad del fluido marcado como B (BLUE)
    double AlphaB; //Contiene los parámetros que determinan la velocidad del sonido del fluido marcado como B (BLUE)    
    double SurfaceTensionB; //Contiene el parámetro que determinan el valor de la tensión superficial del fluido marcado como B (BLUE)
    double DensidadB; //Contiene la densidad total del fluido tipo B (BLUE)
    //Variables que determinan la interacción entre los fluidos
    double Gamma; //Las razones entre las densidades de los fluidos
    double Delta; //Los parametros que determinan el grosor de las interfaces entre los fluidos
    double Phi; //Variable indicadora del tipo de fluido en cada celda
    double O1,O2,O3,O4,O5; //Parametros de la función que hace suave el tiempo de relajación en la interfaz
    double A1,A2,A3,A4,A5; //Parametros de la función que hace suave el parámetro AlphaR(B) en la interfaz
    double *ftotal; //Contiene la función densidad total en cada dirección
    double *FLujoColor; //Vector que indica la diferencia entre los flujos de cada color
    double *GradienteColor; //Vector que indica el gradiente del cambio de color
    double GradienteColor2, GradienteColorNorma, GradienteColorDotCi,GradienteColorDotCi2, Colisionador; //Variables utilizadas para implementar
                                                                                                         //el colisionador de Reis y el Recoloring
    double Presicion;
    TipoCelda **Celda;
    double ***Frontera;

    public:
    LatticeBoltzmann(void);
    void CrearVariables(int Lx0, int Ly0, double TauR0, double TauB0, double RHOR0, double RHOB0, double AlphaR0, double AlphaB0, double Delta0, double SurfaceTensionR0, double SurfaceTensionB0, double Presicion0);
    void CondicionesDeFrontera(string tipo);
    void CondicionesInicialesCuadradoEstatico(double razon);
    void CondicionesInicialesPoiseuille(double presion);
    double Omega(int ix, int iy, int t, bool CalculeConNew);
    double Alpha(int ix, int iy, int t, bool CalculeConNew);
    //Se utilizará index=0 para fluido tipo R, e index=1 para fluido tipo B
    double Rho(int ix,int iy,int t,bool CalculeConNew,int index);
    double Ux(int ix,int iy,int t,bool CalculeConNew,int index);
    double Uy(int ix,int iy,int t,bool CalculeConNew,int index);
    double feq(double Rho0,double Ux0,double Uy0,int ii,int ix, int iy, int index, bool AlphaContinuo , bool CalculeConNew, int t);
    void ColisioneBGK(int t, bool OmegaContinuo, bool AlphaContinuo, bool CalculeConNew);
    //Nuevo colisionador y tipo varios tipos de Recoloring
    void ColisioneReis(int t, bool CalculeConNew);
    void RecoloringReis(int t, bool CalculeConNew);
    void RecoloringLatvaKokko(int t, bool CalculeConNew, double Beta);
    void RecoloringTolke(int t), bool CalculeConNew;//Incompleto
    void Adveccione(void);
    void ImprimaseVelocidad(int t,char const *NombreArchivo);
    void ImprimaseIndicador(int t,char const *NombreArchivo);
    void Borrar(void);
};

// Se cargan los vectores del lattice D2Q9 y los pesos respectivos
LatticeBoltzmann::LatticeBoltzmann(void){
    //D2Q9
    //Cargar los pesos del colisionador BGK
    w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;

    //Cargar los pesos del colisionador de color de Reis y Phillips
    B[0]=-4.0/27.0; B[1]=B[2]=B[3]=B[4]=2.0/27.0; B[5]=B[6]=B[7]=B[8]=5.0/108.0;

    //Cargar los vectores velocidad
    V[0][0]=0; 
    V[1][0]=0;

    V[0][1]=1;   V[0][2]=0;   V[0][3]=-1;  V[0][4]=0; 
    V[1][1]=0;   V[1][2]=1;   V[1][3]=0;   V[1][4]=-1;

    V[0][5]=1;   V[0][6]=-1;  V[0][7]=-1;  V[0][8]=1; 
    V[1][5]=1;   V[1][6]=1;   V[1][7]=-1;  V[1][8]=-1;
}

//Se crean las variables necesarias para realizar la simulación
void LatticeBoltzmann::CrearVariables(int Lx0, int Ly0, double TauR0, double TauB0, double RHOR0, double RHOB0, double AlphaR0, 
                                      double AlphaB0, double Delta0, double SurfaceTension0, double Presicion0){
    Lx=Lx0; Ly=Ly0; //Dimensiones del lattice
    Cantidad=Cantidad0; //Cantidad de fluidos inmiscibles
    
    //Tiempo de relajación, densidades de los fluidos, velocidad del sonido, Delta de la interfaz y tensión superficial
    TauR = TauR0;
    TauB = TauB0;
    OmegaR=1.0/TauR;
    OmegaB=1.0/TauB;
    RHOR = RHOR0;
    RHOB = RHOB0;
    AlphaR = AlphaR0;
    AlphaB = AlphaB0;
    Delta = Delta0;
    SurfaceTensionR = SurfaceTensionR0;
    SurfaceTensionB = SurfaceTensionB0;
    Presicion = Presicion0;
    //Parametros para hacer suave el tiempo de relajación en la interfaz
    O1=(2.0*OmegaB*OmegaR)/(OmegaB+OmegaR);
    O2=2.0*(OmegaR-O1)/Delta;
    O3=-(O2/(2.0*Delta));
    O4=2.0*(O1-OmegaB)/Delta;
    O5=A4/(2.0*Delta);
    //Parametros para hacer suave el tiempo de relajación en la interfaz
    A1=(2.0*AlphaB*AlphaR)/(AlphaB+AlphaR);
    A2=2.0*(AlphaR-A1)/Delta;
    A3=-(A2/(2.0*Delta));
    A4=2.0*(A1-AlphaB)/Delta;
    A5=A4/(2.0*Delta);
    //Se guardan las razones de las densidades de los fluidos y los parámetros que determinan el grosor de la interfaz
    Gamma = RHOR/RHOB;
    
    //Función de distribución de las partículas (old)
    
   fR = new double**[Lx];
   fB = new double**[Lx];
   for(int ii=0; ii<Lx; ii++){
       fR[ii]=new double*[Ly]; 
       fB[ii]=new double*[Ly];
       for(int jj=0;jj<Ly;jj++){
           fR[ii][jj]= new double[9]; 
           fB[ii][jj]= new double[9];
       }
   }
   
    
    //Función de distribución de las partículas (new)
    fRnew = new double**[Lx];
    fBnew = new double**[Lx];
    for(int ii=0; ii<Lx; ii++){
        fRnew[ii]=new double*[Ly]; 
        fBnew[ii]=new double*[Ly];
        for(int jj=0;jj<Ly;jj++){
            fRnew[ii][jj]= new double[9]; 
            fBnew[ii][jj]= new double[9];
        }
    }

    //Se crean los vectores GradienteColor y FLujoColor que serán utilizados en el colisionador de Reis
    GradienteColor = new double[2];
    FLujoColor = new double [2];
    
    //Se crea el vector ftotal con la densidad total del fluido en cada dirección i.e. desidad de R más densidad de B e
    ftotal = new double[9];
}

//Se construye la geometría de la simulación
void LatticeBoltzmann::CondicionesDeFrontera(string tipo){
    //Se crean las variables que guardan la información del tipo de celda a utilizar
    Celda = new TipoCelda*[Lx];
    for(int ii=0; ii<Lx; ii++){
        Celda[ii]=new TipoCelda[Ly]; 
    }
    
    //Se crean las variables que guradan la información sobre la velocidad y presión en las fronteras
    Frontera = new double**[Lx];
    for(int kk=0;kk<Lx;kk++){
        Frontera[kk]=new double*[Ly];
            for(int jj=0; jj< Ly; jj++){
                Frontera[kk][jj] = new double[3]; //Las primeras dos entradas guardan información sobre la velocidad
                                                //y la tercera entrada da información sobre la presión
            }
    }
    
    //Se crea la geometría de la simulación
    if(tipo=="caja"){
        //En esta caso se tienen 4 paredes estáticas
        for(int kk=0;kk<Lx;kk++){
            for(int jj=0; ll<Ly; jj++){
                if((kk==0)||(kk==Lx-1)||(jj==0)||(jj==Ly-1)){
                    Celda[kk][jj]=Pared; 
                    Frontera[kk][jj][0]=0.0;//Velocidad nula en la dirección x
                    Frontera[kk][jj][1]=0.0;//Velocidad nula en la dirección y
                    Frontera[kk][jj][2]=0.0;//El valor de la presión no es importante, se utiliza un valor cualquiera ya que no va a ser utilizada
                }
                else{
                    Celda[kk][jj]=Fluido;
                }
            }
        }
    }
    
    if(tipo=="periodicas"){
        for(int kk=0;kk<Lx;kk++){
            for(int jj=0; ll<Ly; jj++){
                Celda[kk][jj]=Fluido;
            }
        }    
    }
    
    
}

//Se dan las condiciones iniciales a la simulación
void LatticeBoltzmann::CondicionesInicialesCuadradoEstatico(double razon){
   int Ladox = int(Lx*razon/2.0);
   int Ladoy = int(Ly*razon/2.0);
   int Centrox = Lx/2;
   int Centroy = Ly/2;
    
   for(int ix=0;ix<Lx;ix++){
       for(int iy=0;iy<Ly;iy++){
           for(int ii=0;ii<9;ii++){
                if((abs(Centrox-ix)<Ladox)||(abs(Centroy-iy)<Ladoy)){
                    fR[ix][iy][ii]=fRnew[ix][iy]ii]=feq(RHOR,0.0,0.0,ii,0);
                    fB[ix][iy][ii]=fBnew[ix][iy]ii]=0.0;    
                }
                else{
                    fR[ix][iy][ii]=fRnew[ix][iy]ii]=0.0;
                    fB[ix][iy][ii]=fBnew[ix][iy]ii]=feq(RHOB,0.0,0.0,ii,1);   
                }
            }
        }
    }
   
}

//Se implementa la función que vuelve continuo el factor de relajación
double LatticeBoltzmann::Omega(int ix, int iy,int t, bool CalculeConNew){
    
    Phi=(Rho(ix,iy,t,CalculeConNew,0)-Rho(ix,iy,t,CalculeConNew,1))/(Rho(ix,iy,t,CalculeConNew,0)+Rho(ix,iy,t,CalculeConNew,1));
       
    if(Phi>Delta){
        return OmegaR;
    }
    else if(0<Phi<=Delta){
        return O1+O2*Phi+O3*Phi*Phi;
    }
    else if(-Delta<Phi<=0){
        return O1+O4*Phi+O5*Pḧi*Phi;
    }
    else{
    return OmegaB;
    }
    
}

//Se implementa la función que vuelve continuo el factor Alpha que controla la velocidad del sonido en la interfaz
double LatticeBoltzmann::Alpha(int ix, int iy,int t, bool CalculeConNew){
    
    Phi=(Rho(ix,iy,t,CalculeConNew,0)-Rho(ix,iy,t,CalculeConNew,1))/(Rho(ix,iy,t,CalculeConNew,0)+Rho(ix,iy,t,CalculeConNew,1));
       
    if(Phi>Delta){
        return AlphaR;
    }
    else if(0<Phi<=Delta){
        return A1+A2*Phi+A3*Phi*Phi;
    }
    else if(-Delta<Phi<=0){
        return A1+A4*Phi+A5*Pḧi*Phi;
    }
    else{
    return AlphaB;
    }
    
}

//Se calcula la densidad correspondiente a cada tipo de fluido
double LatticeBoltzmann::Rho(int ix,int iy,int t,bool CalculeConNew,int index){
    double suma=0.0;
    if(index==0){//Fluido tipo R (RED)
        for(int ii=0;ii<9;ii++){
            if(CalculeConNew){suma+=fRnew[ix][iy][ii];}
            else{suma+=fR[ix][iy][ii];}
        }
    }
    else{//Fluido tipo B (BLUE)
        for(int ii=0;ii<9;ii++){
            if(CalculeConNew){suma+=fBnew[ix][iy][ii];}
            else{suma+=fB[ix][iy][ii];}
        }  
    }
    return suma;
}

// Se calcula la velocidad en la dirección x del correspondiente tipo de fluido
double LatticeBoltzmann::Ux(int ix,int iy,int t,bool CalculeConNew,int index){
    
    if(Celda[ix][iy]==Pared){return Frontera[ix][iy][0];}
    else{
        if(index==0){//Fluido tipo R (RED)
            double suma=0.0;
            for(int ii=0;ii<9;ii++){
                if(CalculeConNew){suma+=fRnew[ix][iy][ii]*V[0][ii];}
                else{suma+=fR[ix][iy][ii]*V[0][ii];}
            }
            return suma/Rho(ix,iy,t,CalculeConNew,index);
        }
        else{//Fluido tipo B (BLUE)
            double suma=0.0;
            for(int ii=0;ii<9;ii++){
                if(CalculeConNew){suma+=fBnew[ix][iy][ii]*V[0][ii];}
                else{suma+=fB[ix][iy][ii]*V[0][ii];}
            }
            return suma/Rho(ix,iy,t,CalculeConNew,index);
                
        }
    }
}

// Se calcula la velociad en la dirección y del correspondiente tipo de fluido
double LatticeBoltzmann::Uy(int ix,int iy,int t,bool CalculeConNew,int index){
    
    if(Celda[ix][iy]==Pared){return Frontera[ix][iy][1];}
    else{
        if(index==0){//Fluido tipo R (RED)
            double suma=0.0;
            for(int ii=0;ii<9;ii++){
                if(CalculeConNew){suma+=fRnew[ix][iy][ii]*V[1][ii];}
                else{suma+=fR[ix][iy][ii]*V[1][ii];}
            }
            return suma/Rho(ix,iy,t,CalculeConNew,index);
        }
        else{//Fluido tipo B (BLUE)
            double suma=0.0;
            for(int ii=0;ii<9;ii++){
                if(CalculeConNew){suma+=fBnew[ix][iy][ii]*V[1][ii];}
                else{suma+=fB[ix][iy][ii]*V[1][ii];}
            }
            return suma/Rho(ix,iy,t,CalculeConNew,index);
                
        }
    }
  
}

// Se calcula la función de equilibrio para cada tipo de fluido
double LatticeBoltzmann::feq(double Rho0,double Ux0,double Uy0,int ii,int ix, int iy, int index, bool AlphaContinuo , bool CalculeConNew, int t){
    
    double UdotVi=Ux0*V[0][ii]+Uy0*V[1][ii];
    double U2=Ux0*Ux0+Uy0*Uy0;
    double Normal=w[ii]*(1.0+3.0*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
    
    if(AlphaContinuo){
        if(ii==0){
            return=Rho0*(Alpha(ix,iy,t,CalculeConNew)-U2*(2.0/3.0)); //v0
        }
        else if((1<=ii)&(ii<=4)){
            return Rho0*(((1.0+Alpha(ix,iy,t,CalculeConNew))/5.0) + Normal);//v1,...,v4
        }
        else{
            return Rho0*(((1.0+Alpha(ix,iy,t,CalculeConNew))/20.0) + Normal);//v5,...,v8
        }
    }
    else{
        if(index==0){
            if(ii==0){//Fluido tipo R (RED)
                return=Rho0*(AlphaR-U2*(2.0/3.0)); //v0
            }
            else if((1<=ii)&(ii<=4)){
                return Rho0*(((1.0+AlphaR)/5.0) + Normal);//v1,...,v4
            }
            else{
                return Rho0*(((1.0+AlphaR)/20.0) + Normal);//v5,...,v8
            }
        }
        else{
            if(ii==0){//Fluido tipo B (BLUE)
                return=Rho0*(AlphaB-U2*(2.0/3.0)); //v0
            }
            else if((1<=ii)&(ii<=4)){
                return Rho0*(((1.0+AlphaB)/5.0) + Normal);//v1,...,v4
            }
            else{
                return Rho0*(((1.0+AlphaB)/20.0) + Normal);//v5,...,v8
            }   
        }
    }
    
}

//Se implementa la colisión de BGK
void LatticeBoltzmann::ColisioneBGK(int t, bool OmegaContinuo, bool AlphaContinuo, bool CalculeConNew){
    int ix,iy,ii,kk; double Rho0,Ux0,Uy0;
    
    if(OmegaContinuo){
        for(kk=0;kk<2;kk++){
            for(ix=0;ix<Lx;ix++){ //Para Cada Celda
                for(iy=0;iy<Ly;iy++){//Para cada celda
                    Rho0=Rho(ix,iy,t,CalculeConNew,kk); Ux0=Ux(ix,iy,t,CalculeConNew,kk);  Uy0=Uy(ix,iy,t,CalculeConNew,kk);
                    if(kk==0){
                        for(ii=0;ii<9;i++){ //en cada direccion      
                            if(Celda[ix][iy]==Pared){
                                fRnew[ix][iy][ii]=feq(Rho0,Frontera[ix][iy][0],Frontera[ix][iy][1],ii,kk);
                            }
                            else{
                            fRnew[ix][iy][ii]=fR[ix][iy][ii]-Omega(ix.iy,t,CalculeConNew)*(fR[ix][iy][ii]-feq(Rho0,Ux0,Uy0,ii,ix,iy,kk,AlphaContinuo,CalculeConNew,t));
                            }
                        }
                    }
                    else{
                        for(ii=0;ii<9;i++){ //en cada direccion      
                            if(Celda[ix][iy]==Pared){
                                fBnew[ix][iy][ii]=feq(Rho0,Frontera[ix][iy][0],Frontera[ix][iy][1],ii,kk);
                            }
                            else{
                            fBnew[ix][iy][ii]=fB[ix][iy][ii]-Omega(ix,iy,t,CalculeConNew)*(fB[ix][iy][ii]-feq(Rho0,Ux0,Uy0,ii,ix,iy,kk,AlphaContinuo,CalculeConNew,t));
                            }
                        }
                    }
                }
            }
        }
    }
    else{
        for(int kk=0;kk<2;kk++){
            for(ix=0;ix<Lx;ix++){ //Para Cada Celda
                for(iy=0;iy<Ly;iy++){//Para cada celda
                    Rho0=Rho(ix,iy,t,false,kk); Ux0=Ux(ix,iy,t,false,kk);  Uy0=Uy(ix,iy,t,false,kk);
                    if(kk==0){
                        for(ii=0;ii<9;i++){ //en cada direccion      
                            if(Celda[ix][iy]==Pared){
                                fRnew[ix][iy][ii]=feq(Rho0,Frontera[ix][iy][0],Frontera[ix][iy][1],ii,kk);
                            }
                            else{
                            fRnew[ix][iy][ii]=fR[ix][iy][ii]-OmegaR*(fR[ix][iy][ii]-feq(Rho0,Ux0,Uy0,ii,ix,iy,kk,AlphaContinuo,CalculeConNew,t));
                            }
                        }
                    }
                    else{
                        for(ii=0;ii<9;i++){ //en cada direccion      
                            if(Celda[ix][iy]==Pared){
                                fBnew[ix][iy][ii]=feq(Rho0,Frontera[ix][iy][0],Frontera[ix][iy][1],ii,kk);
                            }
                            else{
                            fBnew[ix][iy][ii]=fB[ix][iy][ii]-OmegaB*(fB[ix][iy][ii]-feq(Rho0,Ux0,Uy0,ii,ix,iy,kk,AlphaContinuo,CalculeConNew,t));
                            }
                        }
                    }
                }
            }
        }    
    }
}

//Se implementa la colisión de Color de Reis y Phillips
void LatticeBoltzmann::ColisioneReis(int t, bool CalculeConNew){
    int ix,iy,ii,kk; double Rho0,Ux0,Uy0;
    
    for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            //Los vectores GradienteColor y FLujoColor son reiniciados en 0.0
            GradienteColor[0]=GradienteColor[1]=0.0;
                        
            //Se calcula el valor del gradiente de color en cada posición entre cada par de fluidos inmiscibles
            for(int ii=0;ii<9;ii++){
                GradienteColor[0] += v[0][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,0)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,1)); 
                GradienteColor[1] += v[1][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,0)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,1));
            }
            
            //Se aplica el operador de colisión de color si el vector Gradiente Color tiene un valor mayor a un valor dado,
            //esto se realiza para evitar errores numéricos al operar con números muy pequeños
            GradienteColor2 = GradienteColor[0]*GradienteColor[0] + GradienteColor[1]*GradienteColor[1];
            GradienteColorNorma = sqrt(GradienteColor2); 
            if(GradienteColor2>Presicion){
                //Se aplica el colisionador de Reis a cada componente de la distribución de los fluidos
                for(int ii=0; ii<9;ii++){
                    GradienteColorDotCi = 0.0; //Se reinicia la variable GradienteColorDotCi en 0.0
                    GradienteColorDotCi = GradienteColor[0]*v[0][ii] + GradienteColor[1]*v[1][ii];
                    GradienteColorDotCi2 = GradienteColorDotCi*GradienteColorDotCi;
                    Colisionador = (GradienteColorNorma/2.0)*(w[ii]*(GradienteColorDotCi2/GradienteColor2)-B[ii])
                    fRnew[ix][iy][ii] += SurfaceTensionR*Colisionador;
                    fBnew[ix][iy][ii] += SurfaceTensionB*Colisionador;
                }
            }
        }
    }
}

//Se realiza el Recoloring de Reis
void LatticeBoltzmann::RecoloringReis(int t, bool CalculeConNew){
    int ix,iy,kk,indice; double temp;
    
    for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            //Los vectores GradienteColor y FLujoColor son reiniciados en 0.0
            GradienteColor[0]=GradienteColor[1]=0.0;
            FLujoColor[0]=FLujoColor[1]=0.0;
            
            //Se calcula el valor del gradiente de color en cada posición entre cada par de fluidos inmiscibles
            for(int ii=0;ii<9;ii++){
                GradienteColor[0] += v[0][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,0)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,1)); 
                GradienteColor[1] += v[1][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,0)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,1));
                FLujoColor[0] += v[0][ii]*(Rho(ix,iy,t,CalculeConNew,0)-Rho(ix,iy,t,CalculeConNew,1));
                FLujoColor[1] += v[1][ii]*(Rho(ix,iy,t,CalculeConNew,0)-Rho(ix,iy,t,CalculeConNew,1));
            }
            
            //Se guarda la información sobre la densidad completa de los tipos de fluidos
            for(int kk=0;kk<9;kk++){
                ftotal[kk]=fRnew[ix][iy][kk]+fBnew[ix][iy][kk];
            }
            
            //Se guarda la información sobre la cantidad de fluido R en la posición a evaluar
            DensidadR=0.0;
            for(int kk=0; kk<9;kk++){
                DensidadR+=fRnew[ix][iy][kk];
            }
            
            //Se aplica el operador de colisión de color si el vector Gradiente Color tiene un valor mayor a un valor dado,
            //esto se realiza para evitar errores numéricos al operar con números muy pequeños
            
            GradienteColor2 = GradienteColor[0]*GradienteColor[0] + GradienteColor[1]*GradienteColor[1];
            GradienteColorNorma = sqrt(GradienteColor2); 
            //Aplicar el método de Recoloring de Reis & Phillips
            if(GradienteColor2>Presicion){
                //Se calcula el producto entre el gradiente de color y los vectores v_{kk} y se divide por la norma del vector v_{kk}
                for(int kk=0;kk<8;kk++){
                    GradienteDotCi[kk][0]=kk+1;
                    GradienteDotCi[kk][1]=v[0][kk+1]*GradienteColor[0]+v[1][kk+1]*GradienteColor[1];
                    GradienteDotCi[kk][1]/=sqrt(v[0][kk+1]*v[0][kk+1]+v[1][kk+1]*v[1][kk+1]);//Se divide por la norma del vector velocidad
                }
                //Se organizan los elementos del arreglo GradienteDotCi
                for (int cc = 0 ; cc < ( 8 - 1 ); cc++)
                {
                    for (int dd = 0 ; dd < 8 - cc - 1; dd++)
                    {
                        if (GradienteDotCi[dd][1] < GradienteDotCi[dd+1][1]) /* For decreasing order use < */
                        {
                            temp = GradienteDotCi[dd][1];
                            GradienteDotCi[dd][1] = GradienteDotCi[dd+1][1];
                            GradienteDotCi[dd+1][1] = temp;
                            temp = GradienteDotCi[dd][0];
                            GradienteDotCi[dd][0] = GradienteDotCi[dd+1][0];
                            GradienteDotCi[dd+1][0] = temp;
                        }
                    }
                }
                //Se distribuyen las funciones densidad para maximizar el flujo de partículas tipo R en la dirección del gradiente de color
                for(int cc=0;cc<8;cc++){
                    indice=(int)GradienteDotCi[cc][0];
                    if(DensidadR>=ftotal[indice]){
                        fRnew[ix][iy][indice]=ftotal[indice];
                        fBnew[ix][iy][indice]=0.0;
                        DensidadR-=ftotal[indice];
                    }
                    else if(DensidadR>0.0){
                        fRnew[ix][iy][indice]=DensidadR;
                        fBnew[ix][iy][indice]=ftotal[indice]-DensidadR;
                        DensidadR-=DensidadR;
                    }
                    else{
                        fBnew[ix][iy][indice]=ftotal[indice];
                        fRnew[ix][iy][indice]=0.0;
                    }
                }
            }
        }
    }
    
}

//Se implementa el método de Recoloring de Latva-Kokko
void RecoloringLatvaKokko(int t, bool CalculeConNew, double Beta, double AlphaContinuo){
     int ix,iy,kk,indice; double temp,Rho0;
    
    for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            //Los vectores GradienteColor y FLujoColor son reiniciados en 0.0
            GradienteColor[0]=GradienteColor[1]=0.0;
                        
            //Se calcula el valor del gradiente de color en cada posición entre cada par de fluidos inmiscibles
            for(int ii=0;ii<9;ii++){
                GradienteColor[0] += v[0][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,0)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,1)); 
                GradienteColor[1] += v[1][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,0)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,1));
            }
            
            //Se guarda la información sobre la densidad completa de los tipos de fluidos
            for(int kk=0;kk<9;kk++){
                ftotal[kk]=fRnew[ix][iy][kk]+fBnew[ix][iy][kk];
            }
            
            //Se guarda la información sobre la cantidad de fluido R y fluido B en la posición a evaluar
            DensidadR=0.0;
            DensidadB=0.0;
            for(int kk=0; kk<9;kk++){
                DensidadR+=fRnew[ix][iy][kk];
                DensidadB+=fBnew[ix][iy][kk];
            }
            
            //Se aplica el operador de colisión de color si el vector Gradiente Color tiene un valor mayor a un valor dado,
            //esto se realiza para evitar errores numéricos al operar con números muy pequeños
            
            if(GradienteColor2>Presicion){
                //Se calcula el producto entre el gradiente de color y los vectores v_{kk} y se divide por la norma del vector v_{kk}
                for(int kk=0;kk<8;kk++){
                    GradienteDotCi[kk][0]=kk+1;
                    GradienteDotCi[kk][1]=v[0][kk+1]*GradienteColor[0]+v[1][kk+1]*GradienteColor[1];
                    GradienteDotCi[kk][1]/=sqrt(v[0][kk+1]*v[0][kk+1]+v[1][kk+1]*v[1][kk+1]);//Se divide por la norma del vector velocidad
                    GradienteDotCi[kk][1]/=sqrt(GradienteColor[0]*GradienteColor[0]+GradienteColor[1]*GradienteColor[1]);//Se divide por la norma del vector GradienteColor
                }
                //Se realiza el Recoloring de Latva-Kokko
                Rho0=Rho(ix,iy,t,false,0)+Rho(ix,iy,t,false,1);
                for(int kk=1;kk<9;kk++){
                    fRnew[ix][iy][kk]=(DensidadR/(DensidadB+DensidadR)*ftotal[kk]);
                    fRnew+=Beta*DensidadB*DensidadR*GradienteDotCi[kk-1][1]*(feq(Rho0,0.0,0.0,kk,ix,iy,0,AlphaContinuo,false,t))/((DensidadB+DensidadR)*(DensidadB+DensidadR));
                    fBnew[ix][iy][kk]=(DensidadB/(DensidadB+DensidadR)*ftotal[kk]);
                    fBnew-=Beta*DensidadB*DensidadR*GradienteDotCi[kk-1][1]*(feq(Rho0,0.0,0.0,kk,ix,iy,1,AlphaContinuo,false,t))/((DensidadB+DensidadR)*(DensidadB+DensidadR));
                }
            }
        }
    }
    
}

//Se implementa el método de Recoloring de Tolke, incompleto
void RecoloringTolke(int t){
    int ix,iy,kk,indice; double temp,Rho0;
    
    for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            //Los vectores GradienteColor y FLujoColor son reiniciados en 0.0
            GradienteColor[0]=GradienteColor[1]=0.0;
                        
            //Se calcula el valor del gradiente de color en cada posición entre cada par de fluidos inmiscibles
            for(int ii=0;ii<9;ii++){
                GradienteColor[0] += v[0][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,0)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,1)); 
                GradienteColor[1] += v[1][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,0)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,1));
            }
            
            //Se guarda la información sobre la densidad completa de los tipos de fluidos
            for(int kk=0;kk<9;kk++){
                ftotal[kk]=fRnew[ix][iy][kk]+fBnew[ix][iy][kk];
            }
            
            //Se guarda la información sobre la cantidad de fluido R y fluido B en la posición a evaluar
            DensidadR=0.0;
            DensidadB=0.0;
            for(int kk=0; kk<9;kk++){
                DensidadR+=fRnew[ix][iy][kk];
                DensidadB+=fBnew[ix][iy][kk];
            }
            
            //Se aplica el operador de colisión de color si el vector Gradiente Color tiene un valor mayor a un valor dado,
            //esto se realiza para evitar errores numéricos al operar con números muy pequeños
            
            if(GradienteColor2>Presicion){
                //Se calcula el producto entre el gradiente de color y los vectores v_{kk} y se divide por la norma del vector v_{kk}
                for(int kk=0;kk<8;kk++){
                    GradienteDotCi[kk][0]=kk+1;
                    GradienteDotCi[kk][1]=v[0][kk+1]*GradienteColor[0]+v[1][kk+1]*GradienteColor[1];
                    GradienteDotCi[kk][1]/=sqrt(v[0][kk+1]*v[0][kk+1]+v[1][kk+1]*v[1][kk+1]);//Se divide por la norma del vector velocidad
                    GradienteDotCi[kk][1]/=sqrt(GradienteColor[0]*GradienteColor[0]+GradienteColor[1]*GradienteColor[1]);//Se divide por la norma del vector GradienteColor
                }
                //Se realiza el Recoloring de Latva-Kokko
                Rho0=Rho(ix,iy,t,false,0)+Rho(ix,iy,t,false,1);
                for(int kk=1;kk<9;kk++){
                    fRnew[ix][iy][kk]=(DensidadR/(DensidadB+DensidadR)*ftotal[kk]);
                    fRnew+=Beta*DensidadB*DensidadR*GradienteDotCi[kk-1][1]*(feq(Rho0,0.0,0.0,kk,ix,iy,0,AlphaContinuo,false,t))/((DensidadB+DensidadR)*(DensidadB+DensidadR));
                    fBnew[ix][iy][kk]=(DensidadB/(DensidadB+DensidadR)*ftotal[kk]);
                    fBnew-=Beta*DensidadB*DensidadR*GradienteDotCi[kk-1][1]*(feq(Rho0,0.0,0.0,kk,ix,iy,1,AlphaContinuo,false,t))/((DensidadB+DensidadR)*(DensidadB+DensidadR));
                }
            }
        }
    }
    
}

//Se realiza la advección
void LatticeBoltzmann::Adveccione(void){
    
    for(int kk=0;kk<Cantidad;kk++){  
        for(int ix=0;ix<Lx;ix++){
            for(int iy=0;iy<Ly;iy++){
                for(int ii=0;ii<9;ii++){
                f[kk][(ix+V[0][ii]+Lx)%Lx][(iy+V[1][ii]+Ly)%Ly][ii]=fnew[kk][ix][iy][ii]; //Condiciones de frontera periodicas
                }  
            }
        }
    }
}

//Se imprime la información de la velocidad de los fluidos
void LatticeBoltzmann::ImprimaseVelocidad(int t, char const * NombreArchivo){
  ofstream Imprimir;
  Imprimir.open(NombreArchivo);
  //Imprimir<<"x, y, vx, vy"<<endl;
    for(int ix=0;ix<Lx;ix+=4){
       for(int iy=0;iy<Ly;iy+=4){
            Imprimir<<ix<<", "<<iy<<", "<<4.0*Ux(ix,iy,t,true)<<", "<<4.0*Uy(ix,iy,t,true)<<endl;
       }
    }
    Imprimir.close();
}

//Se imprime la información de la función marcador Phi
void LatticeBoltzmann::ImprimaseIndicador(int t, char const * NombreArchivo){
    double Phi;
    Phi=(Rho(ix,iy,t,false,0)-Rho(ix,iy,t,false,1))/(Rho(ix,iy,t,false,0)+Rho(ix,iy,t,false,1));
    ofstream Imprimir;
    Imprimir.open(NombreArchivo);
    //Imprimir<<"x, y, Phi"<<endl;
    for(int ix=0;ix<Lx; ix++){
        for(int iy=0;iy<Ly;iy++){
            Imprimir<<ix<<" "<<iy<<" "<<Phi<<endl; 
        }
    }
    Imprimir.close();
}

//Se borran los punteros creados para guardar la información de la función de distribución de los fluidos
//y la geometría de la simulación
void LatticeBoltzmann::Borrar(){
    //Se borra la información de las funciones de distribución 
    for(int ii = 0; ii < Lx; ++ii) {
        for (int jj = 0; jj < Ly; ++jj){
            delete [] fR[ii][jj]; delete [] fB[ii][jj];
            delete [] fRnew[ii][jj]; delete [] fBnew[ii][jj];
        }
        delete [] fR[ii]; delete [] fB[ii];
        delete [] fRnew[ii]; delete [] fBnew[ii];
    }
    delete [] fR; delete [] fB;  
    delete [] fRnew; delete fBnew;
    
    //Se borra la información de la geomatría de la simulación 
    for(int ii=0; ii<Lx;ii++){
        delete [] Celda[ii];            
    }
    delete [] Celda;

    //Se borra la información sobre las condiciones e frontera
    for(int ii = 0; ii < Lx; ++ii) {
        for (int jj = 0; jj < Ly; ++jj){
            delete [] Frontera[ii][jj]; 
        }
        delete [] Frontera[ii];
    }
    delete [] Frontera;
    
    //Se borra la información sobre el GradienteColor y FLujoColor
    delete [] FLujoColor; delete [] GradienteColor;
    
    //Se borra la información sobre la función densidad total
    delete [] ftotal;
}    

