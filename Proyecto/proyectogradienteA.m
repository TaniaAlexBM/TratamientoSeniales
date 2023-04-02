%PROGRAMA QUE OBTIENE EL TENSOR (T) GADIENTE DE GRAVEDAD A PARTIR DE LA
%COMPONENTE VERTICAL DEL CAMPO GRAVITACIONAL gz
clear all
close all
clc
%SE CARGAN LOS DATOS CORRESPONDIENTES A gz DOMINIO ESPACIAL
gz=load('Gz_prisma2.RES');
%SE FORMA UNA MATRIZ CON EL NUMERO DE DATOS DE "x" e "y"
%DEL TAMAÑO DE gz
[Nx Ny]=size(gz);
%SE DETERMINA DELTAS Y DOMINIO DE "x" e "y"
dx=1/20;
dy=1/20;
x=[0:dx:10-dx];
y=[0:dy:10-dy];
%SE CONTRUYE UNA MALLA ESPACIAL
[xx yy]=meshgrid(x,y);   

%SE OBTIENEN LAS FRECUENCIAS DE NYQUIST Y LA FRECUENCIA FUNDAMENTAL EN
%FUNCION DE "x" e "y"
FNX=1/(2*dx);
F0X=1/(Nx*dx);
FNY=1/(2*dy);    
F0Y=1/(Ny*dy);
%DETRMINAMOS EL DOMINIO DE FRECUENCIAS
X=[-FNX:F0X:FNX-F0X];
Y=[-FNY:F0Y:FNY-F0Y];
%SE CONTRUYE UNA MALLA EN ESPECTRO DE FRECUENCIAS
[XX YY]=meshgrid(X,Y);
%MAGNITUD DEL VECTOR NUMERO DE ONDA
K=sqrt(XX.^2+YY.^2);  
%SE OBTIENE LA TRANFORMADA DE FOURIER DE LOS DATOS gz
GZ=fftshift(fft2(gz));

%BASADOS EN EL ARTICULO, OBTENDREMOSGX Y GY, EN FUNCION DE GZ
GX=((-1i)*XX).*GZ./(K+0.00001);
GY=-1i*YY.*GZ./(K+0.00001);
%REGRESAMOS LAS FUNCION EN FRECUENCIAS AL DOMINO ESPACIAL
gx=real(ifft2(fftshift(GX)));
gy=real(ifft2(fftshift(GY)));
%SE GRAFICAN LAS COMPONENTES DEL GRADIENTE DE GRAVEDAD ?g= (gx,gy,gz)
figure
mesh(xx,yy,gx)
title('\fontsize{16}gx')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15) 
figure
mesh(xx,yy,gy)
title('\fontsize{16}gy')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15)
figure
mesh(xx,yy,gz)
title('\fontsize{16}gx')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15)

%EN BASE AL ARTICULO SE OBTIENE LAS COMPONENTES DEL TENSOR DE GRAVEDAD EN
%EL DOMINIO DEL FOURIER Y POR LA SIMETRIA DEL TENSOR SOLO SE OBTENDRAN
%LAS COMPONENTES GXX, GXY, GXZ, GYY, GYZ, GZZ

GXX= ((-XX.^2)./(K+0.00001)).*GZ;
GXY= -1i*YY.*GX;
GXZ=-1i*XX.*GZ;
%GYX= -1i*YY.*GX;
GYY= ((-YY.^2)./(K+0.00001)).*GZ;
GYZ=-1i*YY.*GZ;
%GZX=-1i*XX.*GZ;
%GZY=-1i*YY.*GZ;
GZZ=K.*GZ;

%CON LA INVESRA DE FOURIER PASAMOS LAS COMPONENTES AL DOMINIO ESPACIAL
gxx=real(ifft2(fftshift(GXX)));
gxy=real(ifft2(fftshift(GXY)));
gxz=real(ifft2(fftshift(GXZ)));
gyy=real(ifft2(fftshift(GYY)));
gyz=real(ifft2(fftshift(GYZ)));
gzz=real(ifft2(fftshift(GZZ)));
%DEFINIMOS EL TENSOR
T=[gxx gxy gxz;gxy gyy gyz;gxz gyz gzz];
%GRAFICAMOS LAS COMPONENTES DEL TENSOR
figure
mesh(xx,yy,gxx)
title('\fontsize{16}gxx')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15)
figure
mesh(xx,yy,gxy)
title('\fontsize{16}gxy     gyx')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15)
figure
mesh(xx,yy,gxz)
title('\fontsize{16}gxz     gzx')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15)
figure
mesh(xx,yy,gyy)
title('\fontsize{16}gyy')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15)
figure
mesh(xx,yy,gyz)
title('\fontsize{16}gyz       gzy')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15)
figure
mesh(xx,yy,gzz)
title('\fontsize{16}gzz')
xlabel('x[km]','FontSize',15),ylabel('y[km]','FontSize',15),zlabel('z[km]','FontSize',15)