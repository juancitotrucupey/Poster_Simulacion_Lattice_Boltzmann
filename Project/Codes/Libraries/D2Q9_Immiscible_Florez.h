//#include <sstream>
//#include <string>
//#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
//const int Lx=200;
//const int Ly=200;

const double Tau=0.53;
const double RHO0=1.0, UX0=0.06, UY0=0;

enum TipoCelda{aire,obstaculo,ventilador};

class LatticeBoltzmann{
    private:
    int Cantidad;
    int Lx,Ly; //Dimensiones de la cuadricula
    double w[9]; //Pesos de la colisión BGK
    double B[9];  //Pesos del colisionador de color de Reis y Phillips
    int V[2][9];  //V[x=0,y=1][i]  , i=el vector velocidad
    double ****f; //funciones de distribución (old)
    double ****fnew; //funciones de distribución (new)
    double *Tau; //Tiempo de relajación de los fluidos
    double *Omega; //Inverso del tiempo de relajación de los fluidos
    double **Gamma; //Las razones entre las densidades de los fluidos
    double **Delta; //Los parametros que determinan el grosor de las interfaces entre los fluidos
    double *RHO; //Densidades de los fluidos
    double *Alpha; //Contiene los parámetros que determinan la velocidad del sonido
    double **SurfaceTension; //Contiene los parámetros que determinan el valor de la tensión superficial entre fluidos
    double ***FLujoColor; //Vector que indica la diferencia entre los flujos de cada color
    double ***GradienteColor; //Vector que indica el gradiente del cambio de color
    
    /*
    double ***f0;
    double ***f0new; //f[ix][iy][i]  ix,iy=la celda  ,  i=el vector velocidad
    double ***f1;
    double ***f1new;
    */
    TipoCelda **Celda;

    public:
    LatticeBoltzmann(void);
    void CrearVariables(int Lx0, int Ly0, int Cantidad0, double *Tau0[], double *RHO0[]);
    void ConstruyaLaGeometria(void);
    void Inicie(void);
    double Rho(int ix,int iy,int t,bool CalculeConNew,int index);
    double Ux(int ix,int iy,int t,bool CalculeConNew,int index);
    double Uy(int ix,int iy,int t,bool CalculeConNew,int index);
    double feq(double Rho0,double Ux0,double Uy0,int i,int index);
    void ColisioneBGK(int t);
    void ColisioneReis(int t);
    void Recoloring(int t);
    void Adveccione(void);
    void Imprimase(int t,char const *NombreArchivo);
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
void LatticeBoltzmann::CrearVariables(int Lx0, int Ly0, int Cantidad0, double *Tau0[], double *RHO0[], double *Alpha0[], 
                                      double *Delta0[][], double *SurfaceTension0[][]){
    Lx=Lx0; Ly=Ly0; //Dimensiones del lattice
    Cantidad=Cantidad0; //Cantidad de fluidos inmiscibles
    
    //Tiempo de relajación y densidades de los fluidos
    Tau = new double[Cantidad];
    RHO = new double[Cantidad];
    Omega = new double[Cantidad];
    //Se guardan los valores de los tiempos de relajación de los fluidos y su densidad 
    for(int kk=0;kk<Cantidad;kk++){
        Tau[kk]=Tau0[kk]; RHO[kk]=RHO0[kk];
        Omega[kk]=1.0/Tau[kk];
    }
    
    //Se guardan los parámetros que determinan la velocidad del sonido en cada fluido
    Alpha = new double[Cantidad];
    for(int kk=0;kk<Cantidad;kk++){
        Alpha[kk]=Alpha0[kk];
    }
    
    //Se guardan las razones de las densidades de los fluidos y los parámetros que determinan el grosor de la interfaz
    Gamma = new double*[Cantidad];
    Delta = new double*[Cantidad];
    SurfaceTension = new double*[Cantidad];
    for(int kk=0;kk<Cantidad;kk++){
        Gamma[kk] = new double[Cantidad];
        Delta[kk] = new double[Cantidad];
        SurfaceTension[kk] = new double[Cantidad];
    }
    
    for(int ii=0;ii<Cantidad;ii++){
        for(int jj=0;jj<Cantidad;jj++){
            Gamma[ii][jj] = RHO[ii]/RHO[jj]; //La matriz Gamma guarda la información sobre las razones entre las densidades de los fluidos
            SurfaceTension[ii]{jj] = SurfaceTension0[ii][jj]; //La matriz SurfaceTension guarda a información de los parámetros
                                                              //que determinan el valor de la tensión superficial entre cada par de fluidos ii, jj
        }
    }
    
    for(int ii=0;ii<Cantidad;ii++){
        for(int jj=0;jj<Cantidad;jj++){
            Delta[ii][jj] = Delta0[ii][jj]; //La matriz Delta guarda el grosor de las interfaces entre los fluidos ii, jj
        }
    }

    
    //Función de distribución de las partículas (old)
    f = new double***[Cantidad]
    for(int kk=0;kk<Cantidad;kk++){
        //Función de distribución de las partículas marcadas como kk (old) 
        f[kk] = new double**[Lx];
        for(int ii=0; ii<Lx; ii++){
            f[kk][ii]=new double*[Ly]; 
            for(int jj=0;jj<Ly;jj++){
                f0[kk][ii][jj]= new double[9];  
            }
        }
    }
    
    //Función de distribución de las partículas (new)
    fnew = new double***[Cantidad]
    for(int kk=0; kk<Cantidad;kk++){
        // Función de distribución de las partículas marcadas como kk (new)
        fnew[kk] = new double**[Lx];
        for(int ii=0; ii<Lx; ii++){
            fnew[kk][ii]=new double*[Ly]; 
            for(int jj=0;jj<Ly;jj++){
                fnew[kk][ii][jj]= new double[9];  
            }
        }
    }
    
    //Vectores que guardan la información sobre la los GradienteColor y FLujoColor
    GradienteColor = new double**[Lx];
    FLujoColor = new double**[Lx];
    for(int ii=0;ii<Lx;ii++){
        GradienteColor = new double*[Ly];
        FLujoColor = new double*[Lx];
        for(int jj=0;jj<Ly;jj++){
            GradienteColor = new double[2];
            FLujoColor = new double[2];
        }
    }
    
}

//Se construye la geometría de la simulación
void LatticeBoltzmann::ConstruyaLaGeometria(void){
    //Se implementan las reglas correspondientes i.e. fuente, pared, etc
    Celda = new TipoCelda*[Lx];
    for(int ii=0; ii<Lx; ii++){
        Celda[ii]=new TipoCelda[Ly]; 
    }
    
    //Se crea la geometría de la simulación
    int ixc=Lx/4, iyc=Ly/2, R=Ly/5;
    for(int ix=0;ix<Lx;ix++){
        for(int iy=0;iy<Ly;iy++){
            if(ix==0){Celda[ix][iy]=ventilador;}
            else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R*R){Celda[ix][iy]=obstaculo;}
            else{Celda[ix][iy]=aire;}
        }
    }
    Celda[ixc][iyc+R+1]=obstaculo;
}

//Se dan las condiciones iniciales a la simulación
void LatticeBoltzmann::Inicie(void){
    
    for(int kk=0;kk<Cantidad;kk++){  
        for(int ix=0;ix<Lx;ix++){
            for(int iy=0;iy<Ly;iy++){
                for(int ii=0;ii<9;ii++){
                    if(Celda[ix][iy]==obstaculo){f[kk][ix][iy][ii]=fnew[kk][ix][iy][ii]=feq(RHO0,0,0,ii,kk);}
                    else{f[kk][ix][iy][ii]=fnew[kk][ix][iy][ii]=feq(RHO0,UX0,UY0,ii,kk);}
                }
            }
        }
    }
}

//Se calcula la densidad correspondiente a cada tipo de fluido
double LatticeBoltzmann::Rho(int ix,int iy,int t,bool CalculeConNew,int index){
    double suma=0.0;
    for(int ii=0;ii<9;ii++){
        if(CalculeConNew){suma+=fnew[index][ix][iy][ii];}
        else{suma+=f[index][ix][iy][ii];}
    }
    return suma;
}

// Se calcula la velocidad en la dirección x del correspondiente tipo de fluido
double LatticeBoltzmann::Ux(int ix,int iy,int t,bool CalculeConNew,int index){
    if(Celda[ix][iy]==ventilador){return UX0;}
    else if(Celda[ix][iy]==obstaculo){return 0.0;}
    else{
        double suma=0.0;
        for(int ii=0;ii<9;ii++){
            if(CalculeConNew){suma+=fnew[index][ix][iy][ii]*V[0][ii];}
            else{suma+=f[index][ix][iy][ii]*V[0][ii];}
        }
        return suma/Rho(ix,iy,t,CalculeConNew,index);
    }
}

// Se calcula la velociad en la dirección y del correspondiente tipo de fluido
double LatticeBoltzmann::Uy(int ix,int iy,int t,bool CalculeConNew,int index){
    if(Celda[ix][iy]==ventilador){return UY0;}
    else if(Celda[ix][iy]==obstaculo){return 0.0;}
    else{
        double suma=0.0;
        for(int ii=0;ii<9;ii++){
            if(CalculeConNew){suma+=fnew[index][ix][iy][ii]*V[1][ii];}
            else{suma+=f[index][ix][iy][ii]*V[1][ii];}
        }
        return suma/Rho(ix,iy,t,CalculeConNew,index);
    }
}

// Se calcula la función de equilibrio para cada tipo de fluido
double LatticeBoltzmann::feq(double Rho0,double Ux0,double Uy0,int ii,int index){
    double UdotVi=Ux0*V[0][ii]+Uy0*V[1][ii];
    double U2=Ux0*Ux0+Uy0*Uy0;
    double Normal=w[ii]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
        
    if(ii==0){
        return=Rho0*(Alpha[index]-U2*(2.0/3.0)); //v0
    }
    else if((1<=ii)&(ii<=4)){
        return Rho0*(((1.0+Alpha[index])/5.0) + Normal);//v1,...,v4
    }
    else{
        return Rho0*(((1.0+Alpha[index])/20.0) + Normal);//v5,...,v8
    }
}

//Se implementa la colisión de BGK
void LatticeBoltzmann::ColisioneBGK(int t){
    int ix,iy,ii.kk; double Rho0,Ux0,Uy0;
    
    for(kk=0;kk<Cantidad;kk++){ //Para cada tipo de fluido
        for(ix=0;ix<Lx;ix++){ //Para Cada Celda
            for(iy=0;iy<Ly;iy++){
                Rho0=Rho(ix,iy,t,false,kk); Ux0=Ux(ix,iy,t,false,kk);  Uy0=Uy(ix,iy,t,false,kk);
                for(ii=0;ii<9;i++){ //en cada direccion      
                    if(Celda[ix][iy]==ventilador){
                        fnew[kk][ix][iy][ii]=feq(Rho0,UX0,UY0,ii,kk);
                    }
                    else if(Celda[ix][iy]==obstaculo){
                    fnew[kk][ix][iy][ii]=feq(Rho0,0,0,ii,kk);
                    }
                    else{
                    fnew[kk][ix][iy][ii]=f[kk][ix][iy][ii]-Omega[kk]*(f[kk][ix][iy][ii]-feq(Rho0,Ux0,Uy0,ii,kk));
                    }
                }
            }
        }
    }
}

//Se implementa la colisión de Color de Reis y Phillips
void LatticeBoltzmann::ColisioneReis(int t){
    int ix,iy,ii.kk; double Rho0,Ux0,Uy0;
    
    //Los vecotres GradienteColor y FLujoColor son reiniciados en 0.0
    for(int kk=0;kk<Lx;kk++){
         for(int ll=0;ll<Ly;ll++){
             GradienteColor[kk][ll][0]=GradienteColor[kk][ll][1]=0.0;
             FLujoColor[kk][ll][0]=FLujoColor[kk][ll][1]=0.0;
        }
    }
    
    //Se calcula el valor del gradiente de color en cada posición entre cada par de fluidos inmiscibles
    for(int ix=0;ix<Lx;ix++){
         for(int iy=0;iy<Ly;iy++){
              for(int ii=0;ii<9;ii++){
                   GradienteColor[ix][iy][0] += v[0][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,Rojo)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,Azul)); 
                   GradienteColor[ix][iy][1] += v[1][ii]*(Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,Rojo)-Rho(ix+v[0][ii],iy+v[1][ii],t,CalculeConNew,Azul));
              }
        }
    }
      
    //Se aplica el operador de colisión de color si el vector Gradiente Color tiene un valor mayor a un valor dado,
    //esto se realiza para evitar errores numéricos al operar con números muy pequeños
    double GradienteColor2, GradienteColorNorma, GradienteColorDotCi,GradienteColorDotCi2, Colisionador;
    for(int ll=0; ll< Cantidad; ll++){
        for(int mm = ll+1; mm < Cantidad; mm++){
            for(int ix=0;ix<Lx;ix++){
                for(int iy=0;iy<Ly;iy++){
                    GradienteColor2 = GradienteColor[ix][iy][0]*GradienteColor[ix][iy][0] + GradienteColor[ix][iy][1]*GradienteColor[ix][iy][1];
                    GradienteColorNorma = sqrt(GradienteColor2); 
                    //Se aplica el colisionador de Reis a cada componente de la distribución de los fluidos
                    for(int kk=0; kk<9;kk++){
                        GradienteColorDotCi = 0.0; //Se reinicia la variable GradienteColorDotCi en 0.0
                        GradienteColorDotCi = GradienteColor[ix][iy][0]*v[0][kk] + GradienteColor[ix][iy][1]*v[1][kk];
                        GradienteColorDotCi2 = GradienteColorDotCi*GradienteColorDotCi;
                        Colisionador = (SurfaceTension[ll][mm]*GradienteColorNorma/2.0)*(w[kk]*(GradienteColorDotCi2/GradienteColor2)-B[kk])
                        fnew[ll][ix][iy][kk]+=Colisionador;
                        fnew[mm][ix][iy][kk]+=Colisionador;
                    }
                }
            } 
        }
    }
    
}

//Se realiza el Recoloring de Reis
void LatticeBoltzmann::Recoloring(int t){
    
    
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

//Se imprime la información de los fluidos
void LatticeBoltzmann::Imprimase(int t, char const * NombreArchivo){
  ofstream Imprimir;
  Imprimir.open(NombreArchivo);
  Imprimir<<"x, y, vx, vy"<<endl;
    for(int ix=0;ix<Lx;ix+=4)
       for(int iy=0;iy<Ly;iy+=4)
            Imprimir<<ix<<", "<<iy<<", "<<4.0/UX0*Ux(ix,iy,t,true)<<", "<<4.0/UX0*Uy(ix,iy,t,true)<<endl;
  Imprimir.close();
}

//Se borran los punteros creados para guardar la información de la función de distribución de los fluidos
//y la geometría de la simulación
void LatticeBoltzmann::Borrar(){
    //Se borra la información de las funciones de distribución 
    for(int kk=0;kk<Cantidad;kk++){
        for(int ii = 0; ii < Lx; ++ii) {
            for (int jj = 0; jj < Ly; ++jj){
                delete [] f[kk][ii][jj];
                delete [] fnew[kk][ii][jj];
            }
            delete [] f[kk][ii];
            delete [] fnew[kk][ii];
        }
        delete [] f[kk];   
        delete [] fnew[kk];
    }
    delete [] f;
    delete [] fnew;
    
    //Se borra la información de la geomatría de la simulación 
    for(int ii=0; ii<Lx;ii++){
        delete [] Celda[ii];            
    }
    delete [] Celda;
}    

 
