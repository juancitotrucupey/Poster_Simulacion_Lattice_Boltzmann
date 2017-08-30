#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
 
const int Lx=128;
const int Ly=128;

const double Taur=0.53;
const double Taub=0.53;
const double Ar=1.0;
const double Ab=1.0;
const double Alphar=0.5;
const double Alphab=0.5;
const double Betta=1.0;
const double Err=pow(10,-9);


const double Grosor=0.5;
const double RhoRInicial=1.0;
const double RhoBInicial=(1.0-Alphar)*RhoRInicial/(1.0-Alphab);

enum TipoCelda{fluidor,fluidob,frontera,fuente,obstaculo};

class LatticeBoltzmann{
private:
  double w[9];
  double Col2[9];
  int V[2][9];  //V[x=0,y=1][i]  , i=el vector velocidad
  double fr[Lx][Ly][9],frnew[Lx][Ly][9]; //f[ix][iy][i]  ix,iy=la celda  ,  i=el vector velocidad
  double fb[Lx][Ly][9],fbnew[Lx][Ly][9]; //f[ix][iy][i]  ix,iy=la celda  ,  i=el vector velocidad

  TipoCelda Celda[Lx][Ly];
public:
  LatticeBoltzmann(void);
  void ConstruyaLaGeometria(void);
  void Inicie(void);
  double Rho1(int ix,int iy,int t,bool CalculeConNew,bool Rojo);
  double Ux1(int ix,int iy,int t,bool CalculeConNew);
  double Uy1(int ix,int iy,int t,bool CalculeConNew);
  double feq(double Rho0,double Ux0,double Uy0,int i,bool Rojo);
  double Tau(int ix, int iy,int t);
  void Colisione1(int t);
  double Colisionador2(int ix, int iy,int t,int j,bool Rojo);
  void ColorNuevo(int ix,int iy,int t);
  void Adveccione1(int t);
  void ImprimaseVelocidad(char const * NombreArchivo,int t,bool Rojo);
  void ImprimaseDensidad(char const * NombreArchivo,int t,bool Rojo);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //D2Q9
  //Cargar los pesos
  w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Cargar los pesos del colisionador
  Col2[0]=-4/27; Col2[1]=Col2[2]=Col2[3]=Col2[4]=2/27; Col2[8]=Col2[7]=Col2[6]=Col2[5]=5/108;
  //Cargar los vectores velocidad
  V[0][0]=0; 
  V[1][0]=0;

  V[0][1]=1;   V[0][2]=0;   V[0][3]=-1;  V[0][4]=0; 
  V[1][1]=0;   V[1][2]=1;   V[1][3]=0;   V[1][4]=-1;

  V[0][5]=1;   V[0][6]=-1;  V[0][7]=-1;  V[0][8]=1; 
  V[1][5]=1;   V[1][6]=1;   V[1][7]=-1;  V[1][8]=-1;
}
void LatticeBoltzmann::ConstruyaLaGeometria(void){
  int ix,iy;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      if(((ix-64)*(ix-64)+(iy-84)*(iy-84)>20*20)&&((ix-64)*(ix-64)+(iy-44)*(iy-44)>20*20))
	Celda[ix][iy]=fluidob;
            else  Celda[ix][iy]=fluidor;
    }
}
void LatticeBoltzmann::Inicie(void){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<9;i++){
	if(Celda[ix][iy]==fluidor) {
	  fr[ix][iy][i]=frnew[ix][iy][i]=feq(RhoRInicial,0,0,i,true);
	   fb[ix][iy][i]=fbnew[ix][iy][i]=0;
	}
	else{
	   fb[ix][iy][i]=fbnew[ix][iy][i]=feq(RhoBInicial,0,0,i,true);
	   fr[ix][iy][i]=frnew[ix][iy][i]=0;
	}
}
}
double LatticeBoltzmann::Rho1(int ix,int iy,int t,bool CalculeConNew,bool Rojo){
  double suma=0;
  if(Rojo){
  for(int i=0;i<9;i++)
    if(CalculeConNew)
      suma+=frnew[ix][iy][i];
    else
      suma+=fr[ix][iy][i];
  return suma;
}
  else{
    for(int i=0;i<9;i++)
    if(CalculeConNew)
      suma+=fbnew[ix][iy][i];
    else
      suma+=fb[ix][iy][i];
  return suma;
  }
}  
double LatticeBoltzmann::Ux1(int ix,int iy,int t,bool CalculeConNew){
  double suma=0,Rho;
  Rho=Rho1(ix,iy,t,CalculeConNew,true)+Rho1(ix,iy,t,CalculeConNew,false);
     for(int i=0;i<9;i++)
      if(CalculeConNew)
	suma+=(frnew[ix][iy][i]*V[0][i]+fbnew[ix][iy][i]*V[0][i]);
      else
	suma+=(fr[ix][iy][i]*V[0][i]+fb[ix][iy][i]*V[0][i]);
    return suma/Rho;
}
double LatticeBoltzmann::Uy1(int ix,int iy,int t,bool CalculeConNew){
  double suma=0,Rho;
  Rho=Rho1(ix,iy,t,CalculeConNew,true)+Rho1(ix,iy,t,CalculeConNew,false);
     for(int i=0;i<9;i++)
      if(CalculeConNew)
	suma+=(frnew[ix][iy][i]*V[1][i]+fbnew[ix][iy][i]*V[1][i]);
      else
	suma+=(fr[ix][iy][i]*V[1][i]+fb[ix][iy][i]*V[1][i]);
    return suma/Rho;
 }
double LatticeBoltzmann::feq(double Rho0,double Ux0,double Uy0,int i,bool Rojo){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i];
  double U2=Ux0*Ux0+Uy0*Uy0;
  if(Rojo){
    if(i==0) return Rho0*(Alphar-2*U2/3);
    if((0<i)&&(i<5)) return Rho0*((1-Alphar)/5.0+w[i]*(3*UdotVi+9*UdotVi*UdotVi/2-3*U2/2));
    if((4<i)&&(i<9)) return Rho0*((1-Alphar)/20.0+w[i]*(3*UdotVi+9*UdotVi*UdotVi/2-3*U2/2)); 
}
  else{
    if(i==0) return Rho0*(Alphab-2*U2/3);
    if((0<i)&&(i<5)) return Rho0*((1-Alphab)/5.0+w[i]*(3*UdotVi+9*UdotVi*UdotVi/2-3*U2/2));
    if((4<i)&&(i<9)) return Rho0*((1-Alphab)/20.0+w[i]*(3*UdotVi+9*UdotVi*UdotVi/2-3*U2/2)); 
}
}
double LatticeBoltzmann::Tau(int ix,int iy,int t){
  double Rho; Rho=Rho1(ix,iy,t,false,true)+Rho1(ix,iy,t,false,false);
  double Phy; Phy=(Rho1(ix,iy,t,false,true)-Rho1(ix,iy,t,false,false))/Rho;
  double Alpha; Alpha=2*Taur*Taub/(Taur+Taub);
  double Gamma; Gamma=2*(Taur-Alpha)/Grosor;
  double Epsilon; Epsilon=(-1)*Gamma/(2*Grosor);
  double Etta; Etta=2*(Alpha-Taub)/Grosor;
  double Omega; Omega=Etta/(2*Grosor);
		if(Phy>Grosor) return Taur;
                else if((0<=Phy)&&(Phy<=Grosor)) return Alpha+Gamma*Phy+Epsilon*Phy*Phy;
                else if((0>Phy)&&(Phy>=(-1)*Grosor)) return Alpha+Etta*Phy+Omega*Phy*Phy;
		else if(Phy<(-1)*Grosor) return Taub;
}
void LatticeBoltzmann::Colisione1(int t){
  int ix,iy,i; double Rho0,Ux0,Uy0;
  //Rojo
  for(ix=0;ix<Lx;ix++) {//Para Cada Celda
    for(iy=0;iy<Ly;iy++){
      Rho0=Rho1(ix,iy,t,false,true); Ux0=Ux1(ix,iy,t,false);  Uy0=Uy1(ix,iy,t,false);
      for(i=0;i<9;i++) {//en cada direccion      
	frnew[ix][iy][i]=fr[ix][iy][i]-1.0/Tau(ix,iy,t)*(fr[ix][iy][i]-feq(Rho0,Ux0,Uy0,i,true))+Colisionador2(ix,iy,t,i,true);
    }
  }
  }
      //Blue
  for(ix=0;ix<Lx;ix++){ //Para Cada Celda
    for(iy=0;iy<Ly;iy++){
      Rho0=Rho1(ix,iy,t,false,false); Ux0=Ux1(ix,iy,t,false);  Uy0=Uy1(ix,iy,t,false);
      for(i=0;i<9;i++){ //en cada direccion      
	fbnew[ix][iy][i]=fb[ix][iy][i]-1.0/Tau(ix,iy,t)*(fb[ix][iy][i]-feq(Rho0,Ux0,Uy0,i,false))+Colisionador2(ix,iy,t,i,false);
    }
  }
  }
  for(ix=0;ix<Lx;ix++){
      for(iy=0;iy<Ly;iy++){
	 ColorNuevo(ix,iy,t);
      }
}
}
 double LatticeBoltzmann::Colisionador2(int ix, int iy,int t,int j, bool Rojo){
   double GradienteColor[2]={0,0};
     for(int i=1;i<9;i++){
 GradienteColor[0]+=(V[0][i]*(Rho1((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly,t,false,true)-Rho1((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly,t,false,false)));
 GradienteColor[1]+=(V[1][i]*(Rho1((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly,t,false,true)-Rho1((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly,t,false,false)));
    }
     double NormaGradiente=sqrt((GradienteColor[0]*GradienteColor[0])+(GradienteColor[1]*GradienteColor[1]));
      double VjdotGradiente=V[0][j]*GradienteColor[0]+V[1][j]*GradienteColor[1];
      double Forcing;
      if(NormaGradiente>Err){
	Forcing=0.5*NormaGradiente*(w[j]*VjdotGradiente*VjdotGradiente/(NormaGradiente*NormaGradiente)-Col2[j]);
      }
      else Forcing=0;
      if(Rojo) Forcing*=Ar;
      else Forcing*=Ab;
      return Forcing;
  }
 void LatticeBoltzmann::ColorNuevo(int ix,int iy, int t){
   double GradienteColor[2]={0,0};
     for(int i=0;i<9;i++){
 GradienteColor[0]+=(V[0][i]*(Rho1((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly,t,false,true)-Rho1((ix+V[0][i]+Lx)%Ly,(iy+V[1][i]+Ly)%Ly,t,false,false)));
 GradienteColor[1]+=(V[1][i]*(Rho1((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly,t,false,true)-Rho1((ix+V[0][i]+Lx)%Lx,(iy+V[1][i]+Ly)%Ly,t,false,false)));
    }
     double NormaGradiente=sqrt((GradienteColor[0]*GradienteColor[0])+(GradienteColor[1]*GradienteColor[1]));
     double Ene[2];
     if(NormaGradiente>Err){
       Ene[0]=GradienteColor[0]/NormaGradiente;       Ene[1]=GradienteColor[1]/NormaGradiente;
       double Rho=Rho1(ix,iy,t,true,true)+Rho1(ix,iy,t,true,false);
       double Rhor=Rho1(ix,iy,t,false,true);
       double Rhob=Rho1(ix,iy,t,false,false);
       double VidotEne,fi;
       for(int i=0;i<9;i++){
       VidotEne=V[0][i]*Ene[0]+V[1][i]*Ene[1];
       fi=frnew[ix][iy][i]+fbnew[ix][iy][i];
       frnew[ix][iy][i]=Rhor*fi/Rho + Betta*w[i]*Rhor*Rhob*VidotEne/(Rhor+Rhob);
       fbnew[ix][iy][i]=Rhob*fi/Rho - Betta*w[i]*Rhor*Rhob*VidotEne/(Rhor+Rhob);
     }
}
}
void LatticeBoltzmann::Adveccione1(int t){
  int ix, iy, i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<9;i++){
	  fr[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=frnew[ix][iy][i];
          fb[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fbnew[ix][iy][i];
      }
  //Condiciones de frontera sobre las paredes
  //Pared izquierda
  /*
  for(ix=0,iy=1;iy<Ly-1;iy++){
    fr[ix][iy][1]=fr[ix][iy][3];
    fr[ix][iy][8]=fr[ix][iy][6]-0.5*(fr[ix][iy][4]-fr[ix][iy][2]);
    fr[ix][iy][5]=fr[ix][iy][7]+0.5*(fr[ix][iy][4]-fr[ix][iy][2]);
    fb[ix][iy][1]=fb[ix][iy][3];
    fb[ix][iy][8]=fb[ix][iy][6]-0.5*(fb[ix][iy][4]-fb[ix][iy][2]);
    fb[ix][iy][5]=fb[ix][iy][7]+0.5*(fb[ix][iy][4]-fb[ix][iy][2]);
  }
  //Pared derecha
  for(ix=Lx-1,iy=1;iy<Ly-1;iy++) {
    fr[ix][iy][3]=fr[ix][iy][1];
    fr[ix][iy][6]=fr[ix][iy][8]-0.5*(fr[ix][iy][2]-fr[ix][iy][4]);
    fr[ix][iy][7]=fr[ix][iy][5]+0.5*(fr[ix][iy][2]-fr[ix][iy][4]);
     fb[ix][iy][3]=fb[ix][iy][1];
    fb[ix][iy][6]=fb[ix][iy][8]-0.5*(fb[ix][iy][2]-fb[ix][iy][4]);
    fb[ix][iy][7]=fb[ix][iy][5]+0.5*(fb[ix][iy][2]-fb[ix][iy][4]);
  }
  //Pared inferior
  for(ix=1,iy=0;ix<Lx-1;ix++) {
    fr[ix][iy][2]=fr[ix][iy][4];
    fr[ix][iy][5]=fr[ix][iy][7]-0.5*(fr[ix][iy][1]-fr[ix][iy][3]);
    fr[ix][iy][6]=fr[ix][iy][8]+0.5*(fr[ix][iy][1]-fr[ix][iy][3]);
     fb[ix][iy][2]=fb[ix][iy][4];
    fb[ix][iy][5]=fb[ix][iy][7]-0.5*(fb[ix][iy][1]-fb[ix][iy][3]);
    fb[ix][iy][6]=fb[ix][iy][8]+0.5*(fb[ix][iy][1]-fb[ix][iy][3]);
  }
  //Pared superior
  for(ix=1,iy=Ly-1;ix<Lx-1;ix++){
    fr[ix][iy][4]=fr[ix][iy][2];
    fr[ix][iy][7]=fr[ix][iy][5]-0.5*(fr[ix][iy][3]-fr[ix][iy][1]);
    fr[ix][iy][8]=fr[ix][iy][6]+0.5*(fr[ix][iy][3]-fr[ix][iy][1]);
     fb[ix][iy][4]=fb[ix][iy][2];
    fb[ix][iy][7]=fb[ix][iy][5]-0.5*(fb[ix][iy][3]-fb[ix][iy][1]);
    fb[ix][iy][8]=fb[ix][iy][6]+0.5*(fb[ix][iy][3]-fb[ix][iy][1]);
  }
  //Condiciones de frontera en las esquinas
  double sumar,sumab;
  //Inferior izquierda
   fr[0][0][1]= fr[0][0][3];  fb[0][0][1]= fb[0][0][3];
  fr[0][0][2]= fr[0][0][4];  fb[0][0][2]= fb[0][0][4];
  fr[0][0][5]= fr[0][0][7];  fb[0][0][5]= fb[0][0][7];
 for(i=0,sumar=0,sumab=0;i<8;i++) {
    sumar+=fr[0][0][i]; sumab+=fb[0][0][i];
  }
  fr[0][0][6]=0.5*(Rho1(1,0,t,false,true)-sumar+fr[0][0][6]);  fb[0][0][6]= 0.5*(Rho1(1,0,t,false,false)-sumab+fb[0][0][6]); 
  fr[0][0][8]= fr[0][0][6];  fb[0][0][8]= fb[0][0][6];
   //Superior izquierda
   fr[0][Ly-1][1]= fr[0][Ly-1][3];  fb[0][Ly-1][1]= fb[0][Ly-1][3];
  fr[0][Ly-1][4]= fr[0][Ly-1][2];  fb[0][Ly-1][4]= fb[0][Ly-1][2];
  fr[0][Ly-1][8]= fr[0][Ly-1][6];  fb[0][Ly-1][8]= fb[0][Ly-1][6];
 for(i=0,sumar=0,sumab=0;i<9;i++) {
    sumar+=fr[0][Ly-1][i]; sumab+=fb[0][Ly-1][i];
  }
fr[0][Ly-1][5]=0.5*(Rho1(0,Ly-2,t,false,true)-sumar+fr[0][Ly-1][5]+fr[0][Ly-1][7]); fb[0][Ly-2][5]= 0.5*(Rho1(0,Ly-1,t,false,false)-sumab+fb[0][Ly-1][5]+fb[0][Ly-1][7]); 
  fr[0][Ly-1][7]= fr[0][Ly-1][5];  fb[0][Ly-1][7]= fb[0][Ly-1][5];
   //Superior derecha
   fr[Lx-1][Ly-1][3]= fr[Lx-1][Ly-1][1];  fb[Lx-1][Ly-1][3]= fb[Lx-1][Ly-1][1];
  fr[Lx-1][Ly-1][4]= fr[Lx-1][Ly-1][2];  fb[Lx-1][Ly-1][4]= fb[Lx-1][Ly-1][2];
  fr[Lx-1][Ly-1][7]= fr[Lx-1][Ly-1][5];  fb[Lx-1][Ly-1][7]= fb[Lx-1][Ly-1][5];
 for(i=0,sumar=0,sumab=0;i<9;i++) {
    sumar+=fr[Lx-1][Ly-1][i]; sumab+=fb[Lx-1][Ly-1][i];
  }
fr[Lx-1][Ly-1][6]=0.5*(Rho1(Lx-1,Ly-2,t,false,true)-sumar+fr[Lx-1][Ly-1][6]+fr[Lx-1][Ly-1][8]); fb[Lx-1][Ly-1][5]= 0.5*(Rho1(Lx-1,Ly-2,t,false,false)-sumab+fb[Lx-1][Ly-1][6]+fb[Lx-1][Ly-1][8]); 
  fr[Lx-1][Ly-1][8]= fr[Lx-1][Ly-1][6];  fb[Lx-1][Ly-1][8]= fb[Lx-1][Ly-1][6];
   //Inferior derecha
  fr[Lx-1][0][3]= fr[Lx-1][0][1];  fb[Lx-1][0][3]= fb[Lx-1][0][1];
  fr[Lx-1][0][2]= fr[Lx-1][0][4];  fb[Lx-1][0][2]= fb[Lx-1][0][4];
  fr[Lx-1][0][6]= fr[Lx-1][0][8];  fb[Lx-1][0][6]= fb[Lx-1][0][8];
 for(i=0,sumar=0,sumab=0;i<9;i++) {
    sumar+=fr[Lx-1][0][i]; sumab+=fb[Lx-1][0][i];
  }
fr[Lx-1][0][5]=0.5*(Rho1(Lx-2,0,t,false,true)-sumar+fr[Lx-1][0][5]+fr[Lx-1][0][7]);  fb[Lx-1][0][5]= 0.5*(Rho1(Lx-2,0,t,false,false)-sumab+fb[Lx-1][0][5]+fb[Lx-1][0][7]); 
  fr[Lx-1][0][7]= fr[Lx-1][0][5];  fb[Lx-1][0][7]= fb[Lx-1][0][5];
  
  */}
void LatticeBoltzmann::ImprimaseVelocidad(char const * NombreArchivo,int t,bool Rojo){
  ofstream MiArchivo(NombreArchivo);
 
  for(int ix=0;ix<Lx;ix+=4)
    for(int iy=0;iy<Ly;iy+=4)
      MiArchivo<<ix<<" "<<iy<<" "<<Ux1(ix,iy,t,false)<<" "<<Uy1(ix,iy,t,false)<<endl;
  MiArchivo.close();

}
void LatticeBoltzmann::ImprimaseDensidad(char const * NombreArchivo,int t,bool Rojo){
  ofstream MiArchivo(NombreArchivo);
 
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0; iy<Ly;iy++)
      //    MiArchivo<<ix<<" "<<iy<<" "<<Tau(ix,iy,t)<<endl;
      MiArchivo<<ix<<" "<<iy<<" "<<Rho1(ix,iy,t,false,Rojo)<<endl;
    MiArchivo<<endl;
  }
      MiArchivo.close();
 }

//------------------------ Funciones Globales ---------------


int main(void){
  LatticeBoltzmann Fluidos;
  int t,tmax=400;

   Fluidos.ConstruyaLaGeometria();
      Fluidos.Inicie();
      t=0;
      
for(t=0;t<tmax;t++){
   Fluidos.Colisione1(t);
   Fluidos.Adveccione1(t);
      }
  Fluidos.ImprimaseVelocidad("VelocidadRojo.dat",t,true);
  Fluidos.ImprimaseDensidad("DensidadRojo.dat",t,true);
  Fluidos.ImprimaseDensidad("DensidadBlue.dat",t,false);   
 
  return 0;
}
