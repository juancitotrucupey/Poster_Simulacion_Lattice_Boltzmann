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
    double w[9];
    int V[2][9];  //V[x=0,y=1][i]  , i=el vector velocidad
    double ****f; //funciones de distribución (old)
    double ****fnew; //funciones de distribución (new)
    double *Tau; //Tiempo de relajación de los fluidos
    double *RHOInicial; //Densidades de los fluidos
    double **FLujoColor; //Vector que indica la diferencia entre los flujos de cada color
    double **GradienteColor; //Vector que indica el gradiente del cambio de color
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
    void Colisione(int t);
    void Adveccione(void);
    void Imprimase(int t,char const *NombreArchivo);
    void Borrar(void);
};

// Se cargan los vectores del lattice D2Q9 y los pesos respectivos
LatticeBoltzmann::LatticeBoltzmann(void){
    //D2Q9
    //Cargar los pesos
    w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;
    //Cargar los vectores velocidad
    V[0][0]=0; 
    V[1][0]=0;

    V[0][1]=1;   V[0][2]=0;   V[0][3]=-1;  V[0][4]=0; 
    V[1][1]=0;   V[1][2]=1;   V[1][3]=0;   V[1][4]=-1;

    V[0][5]=1;   V[0][6]=-1;  V[0][7]=-1;  V[0][8]=1; 
    V[1][5]=1;   V[1][6]=1;   V[1][7]=-1;  V[1][8]=-1;
}

//Se crean las variables necesarias para realizar la simulación
void LatticeBoltzmann::CrearVariables(int Lx0, int Ly0, int Cantidad0, double *Tau0[], double *RHO0[]){
    Lx=Lx0; Ly=Ly0; //Dimensiones del lattice
    Cantidad=Cantidad0; //Cantidad de fluidos inmiscibles
    
    Tau =  new double[Cantidad];
    RHOInicial = new double[Cantidad];
    //Se guardan los valores de los tiempos de relajación de los fluidos y su densidad 
    for(int kk=0;kk<Cantidad;kk++){
        Tau[kk]=Tau0[kk]; RHOInicial[kk]=RHO0[kk];
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
    for(int ix=0;ix<Lx;ix++)
        for(int iy=0;iy<Ly;iy++)
        if(ix==0)
        Celda[ix][iy]=ventilador;
        else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R*R)
        Celda[ix][iy]=obstaculo;
        else
        Celda[ix][iy]=aire;
    Celda[ixc][iyc+R+1]=obstaculo;
}

//Se dan las condiciones iniciales a la simulación
void LatticeBoltzmann::Inicie(void){
    for(int kk=0;kk<Cantidad;kk++){  
        for(int ix=0;ix<Lx;ix++){
            for(int iy=0;iy<Ly;iy++){
                for(int i=0;i<9;i++){
                    if(Celda[ix][iy]==obstaculo)
                    f[kk][ix][iy][i]=fnew[kk][ix][iy][i]=feq(RHO0,0,0,i,kk);
                    else
                    f[kk][ix][iy][i]=fnew[kk][ix][iy][i]=feq(RHO0,UX0,UY0,i,kk);
                }
            }
        }
    }
}

//Se calcula la densidad correspondiente a cada tipo de fluido
double LatticeBoltzmann::Rho(int ix,int iy,int t,bool CalculeConNew,int index){
  double suma=0;
  for(int i=0;i<9;i++)
    if(CalculeConNew)
      suma+=fnew[index][ix][iy][i];
    else
      suma+=f[index][ix][iy][i];
  return suma;
}

// Se calcula la velocidad en la dirección x del correspondiente tipo de fluido
double LatticeBoltzmann::Ux(int ix,int iy,int t,bool CalculeConNew,int index){
  if(Celda[ix][iy]==ventilador)
    return UX0;
  else if(Celda[ix][iy]==obstaculo)
    return 0;
  else{
    double suma=0.0;
    for(int i=0;i<9;i++)
      if(CalculeConNew)
	suma+=fnew[index][ix][iy][i]*V[0][i];
      else
	suma+=f[index][ix][iy][i]*V[0][i];
    return suma/Rho(ix,iy,t,CalculeConNew,index);
  }
}

// Se calcula la velociad en la dirección y del correspondiente tipo de fluido
double LatticeBoltzmann::Uy(int ix,int iy,int t,bool CalculeConNew,int index){
  if(Celda[ix][iy]==ventilador)
    return UY0;
  else if(Celda[ix][iy]==obstaculo)
    return 0;
  else{
    double suma=0.0;
    for(int i=0;i<9;i++)
      if(CalculeConNew)
	suma+=fnew[index][ix][iy][i]*V[1][i];
      else
	suma+=f[index][ix][iy][i]*V[1][i];
    return suma/Rho(ix,iy,t,CalculeConNew);
  }
}

// Se calcula la función de equilibrio para cada tipo de fluido
double LatticeBoltzmann::feq(double Rho0,double Ux0,double Uy0,int i,int index){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i];
  double U2=Ux0*Ux0+Uy0*Uy0;
  return w[i]*Rho0*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

//Se implementa la colisión 
void LatticeBoltzmann::Colisione(int t){
  int ix,iy,i; double Rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++) //Para Cada Celda
    for(iy=0;iy<Ly;iy++){
      Rho0=Rho(ix,iy,t,false); Ux0=Ux(ix,iy,t,false);  Uy0=Uy(ix,iy,t,false);
      for(i=0;i<9;i++) //en cada direccion      
	if(Celda[ix][iy]==ventilador)
	  fnew[ix][iy][i]=feq(Rho0,UX0,UY0,i);
	else if(Celda[ix][iy]==obstaculo)
	  fnew[ix][iy][i]=feq(Rho0,0,0,i);
	else
	  fnew[ix][iy][i]=f[ix][iy][i]-1.0/Tau*(f[ix][iy][i]-feq(Rho0,Ux0,Uy0,i));
    }
}

//Se realiza la advección
void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<9;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
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

