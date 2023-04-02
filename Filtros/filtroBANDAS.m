close all
clc

% Análisis en frecuencias y filtro pasa bandas
% Filtrado de las frecuencias en función de un intervalo

dt=0.01;

%Eje Temporal
t=[0:dt:20-dt];
f=6*sin(2*pi*t*4)+2*sin(2*pi*t*8);
SFD(t,f)


%definimos el dominio en el espacio tiempo
dt=t(2)-t(1);

%creamos los elementos necesarios para descomponer la señal en series de Fourier
%calculamos A0, que es el promedio de la señal
N=length(t);
A0=sum(f)/N;  %sumamos los elementos del vector f y los dividimos entre N

%definimos el vector de enteros de n que va de 0 a N-1
n=[0:N-1];
%generamos dos vectores que contienen dos factores de peso para la contrucción 
Ak=[];
Bk=[];

for k=1:N-1;
     Ak(k)=(2/N)*sum(f.*cos(2*pi*n*k/N));
     Bk(k)=(2/N)*sum(f.*sin(2*pi*n*k/N));
end

%Amplitud
AMP=sqrt(Ak.^2+Bk.^2);
 
%reacomodamos el vector de amplitudes
AMPD=AMP(1:N/2);
AMPI=AMP(N/2+1:end);
 
%agregamos el armónico cero para completar el espectro de amplitudes
AMP=[AMPI A0 AMPD];

%creamos el eje de frecuencia con la Frecuencia fundamental
F0=1/(dt*N);
FN=1/(2*dt);
F=[-FN+F0:F0:FN];  %eje de frecuencias avanzando cada F0
%graficamos F contra AMP 
figure
subplot(3,1,1) %espectro de amplitudes original
stem(F,AMP)
title('Espectro de amplitudes')
xlabel('t')
ylabel('Amp')

%generamos el FPB
FPB=zeros(1,N); 

% Frecuencias de corte
FC1 = 3;
FC2 = 7;

for c1=1:N
   if abs(F(c1))>=FC1 && abs(F(c1))<=FC2
       FPB(c1)=1;
    end
end 

subplot(3,1,2)
stem(F,FPB) 
title('Filtro pasa bandas')
xlabel('f')
ylabel('Amp')
%filtro pasa bajas, eje de frecuencias, amplitud del filtro
AMP=AMP.*FPB;

subplot(3,1,3)
stem(F,AMP) 
title('Espectro de amplitudes FILTRADO')
xlabel('f')
ylabel('Amp')

FPBI=FPB(1:(N/2)-1);
FPBD=FPB((N/2)+1:end);
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