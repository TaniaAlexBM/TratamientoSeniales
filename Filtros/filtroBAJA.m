close all
clc

%ANALISIS DE FRECUENCIAS Y FILTRADO FB DE UNA SEÑAL DISCRETA

dt=0.01;

%Eje Temporal
t=[0:dt:20-dt];
f=6*sin(2*pi*t*4)+2*sin(2*pi*t*8);
SFD(t,f)


%Definimos el domimio espacio-tiempo
dt=t(2)-t(1);
%Creamos los elementos necesarios para descomponer la señal en series
%de Fourier
N=length(t);
%Calculamos A0 que es el promedio de la señal
A0=sum(f)/N;
%Definimos el vector de enteros de n que va de 0 a N-1
n=[0:N-1];
%Generamos los vectores que contienen los factores de peso para la
%recontrucción
Ak=[];
Bk=[];
%k va de 0 a N-1, pero k=0 se calcula aparte y es A0
for k=1:N-1
    Ak(k)=(2/N)*sum(f.*cos(2*pi*n*k/N));
    Bk(k)=(2/N)*sum(f.*sin(2*pi*n*k/N));
end
AMP=sqrt(Ak.^2+Bk.^2);
%Partimos la señal para saber en que frecuencias está oscilando la señal
%Reacomodamos el vector de amplitudes
AMPD=AMP(1:N/2);
AMPI=AMP(1+N/2:end);
%Agregamos el armónico cero para completar el espectro de amplitudes
AMP=[AMPI A0 AMPD];
%Para caracterizar el eje de frecuencia con los conceptos de F0 y FN
%F0 Que tan espaciados estan
F0=1/(dt*N);
FN=1/(2*dt);
F=[-FN+F0:F0:FN];
figure
subplot(3,1,1)
stem(F,AMP)
title('Espectro de amplitudes')
xlabel('f','fontsize', 10)
ylabel('Espectro de amplitudes','fontsize', 10)
%Generamos el FILTRO BASABAJAS
%Es en el dominio de las frecuencias
FPB=zeros(1,N);
%Definimos la frecuencia de corte FC
for c1=1:N
    if F(c1)>=-FC && F(c1)<=FC
        FPB(c1)=1;
    end
end
subplot(3,1,2)
stem(F,FPB)
title('Filtro Pasa Bajas')
xlabel('f','fontsize', 10)
ylabel('AMP','fontsize', 10)

%Sólo para ver, filtramos el vector AMP
AMP=AMP.*FPB;
subplot(3,1,3)
stem(F,AMP)
title('Espectro de amplitudes filtrados')
xlabel('f')
ylabel('Espectro de amplitudes')
FPBI=FPB(1:N/2-1);
FPBD=FPB(N/2+1:end);
%Se quita el cero porque Ak y Bk no lo contienen
FPB=[FPBD FPBI];
%Convolucion
Ak=Ak.*FPB;
Bk=Bk.*FPB;
%Realizamos la reconstrucción con Ak y Bk filtrados
armonicos=0;
for k=1:N/2
    armonicos=armonicos+Ak(k)*cos(2*pi*n*k/N)+Bk(k)*sin(2*pi*n*k/N); %Suma de vectores
end


%Sumamos finalmente el armónico cero
fr=A0+armonicos;