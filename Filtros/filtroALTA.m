close all
clc

%ANÁLISIS EN FRECUENCIAS Y FILTRADO PASA ALTAS


dt=0.01;

%Eje Temporal
t=[0:dt:20-dt];
f=6*sin(2*pi*t*4)+2*sin(2*pi*t*8);
SFD(t,f)


%Definimos la discretización en tiempo
%Definimos el dominio espacio-tiempo
dt=t(2)-t(1);
N=length(t);
A0=sum(f)/N;
%Definimos el vector de enteros n que va de 0 a N-1
n=[0:N-1]; %n debe tener el elemento 0
%Generamos los vectores que contienen los factores de peso para la
%reconstrucción

Ak=[];%Almacena Ak
Bk=[];%Almacena Bk

%k va de 0 a N-1 porque el K=0 se calcula aparte, ya se calculó en el A0
for k=1:N-1
    Ak(k)=(2/N)*sum(f.*cos(2*pi*n*k/N));%Sumamos para generar un escalar
    Bk(k)=(2/N)*sum(f.*sin(2*pi*n*k/N));%Sumamos para generar un escalar
end
%Amplitud
AMP=sqrt(Ak.^2+Bk.^2);

%Se parte la señal en dos porque nos interesa interpretar en que frecuencia
%esta oscilando y nos interesa hablar de un espectro de amplityud y de fase

%Reacomodamos el vector de amplitudes
AMPD=AMP(1:N/2);
AMPI=AMP(N/2+1:end);

%Agregamos el armónico cero para completar el espectro de amplitudes
AMP=[AMPI A0 AMPD];

%Creamos el eje de frecuencia con los conceptos Frecuencia fundamental (que
%tan separado esta un armonico de otro) y frecuencia de Nyquist
F0=1/(dt*N);
FN=1/(2*dt);

F=[-FN+F0:F0:FN]; %Dominio de las Frecuencias
figure
subplot(3,1,1)
stem(F,AMP)
title('Espectro de Amplitudes')
xlabel('f')
ylabel('Amp')
grid on
%Generamos el filtro pasa altas, en el dominio de las frecuencias F
FPA=ones(1,N); %En dominio de la frecuencias

%Definimos la frecuencia de corte FC
FC = 7;
%de -Fc a la Fc es igual a 1
%De donde a donde se va a filtrar la señal

for c1=1:N
    if F(c1)>=-FC && F(c1)<=FC %Si F en ese punto es menor o igual a FC y si F 
        %en ese punto es menor o igual a FC entonces es igual a 1
        FPA(c1)=0;
    end
end
subplot(3,1,2)
stem(F,FPA)
title('Filtro pasa altas')
xlabel('f')
ylabel('Amp')
grid on

%Sólo para visualizar filtramos el vector AMP
AMP=AMP.*FPA;
subplot(3,1,3)
stem(F,AMP)
title('Espectro de Amplitud filtrado')
xlabel('f')
ylabel('Amp')
grid on
%Para reconstruir usamos Ak y Bk por lo tanto hay que filtrarlos 
%Reacomodamos el filtro pasa bajas

%La parte izquierda va desde 1 hasta N/2 -1, 
FPBI=FPA(1:N/2-1);
%Nos saltamos el 500 porque ya esta en el A0
FPBD=FPA(N/2+1:end);

FPA=[FPBD FPBI]; %Le quitamos el 0 para que empate con Ak y Bk

%Como estamos en el dominio de la frecuencia el filtrado en el tiempo se
%define como una convolución
Ak=Ak.*FPA;
Bk=Bk.*FPA;

%Realizamos la reconstucion con Ak y Bk filtrados
armonicos=0;
for k=1:N/2
    armonicos=armonicos+Ak(k)*cos(2*pi*n*k/N)+Bk(k)*sin(2*pi*n*k/N); %Suma de vectores
end

%Sumamos finalmente el armónico cero
Fr=A0+armonicos;