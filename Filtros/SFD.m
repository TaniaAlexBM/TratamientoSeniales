function [FASE]=SFD(t,f)

%ANALISIS DE FRECUENCIAS Y FILTRADO FB DE UNA SEÑAL DISCRETA
%dt es constante y es la separación de cada muestra
%Definimos el domimio espacio-tiempo
dt=t(2)-t(1);
figure
subplot(3,1,1)
plot(t,f)
title('Señal original')
xlabel('t')
ylabel('f(t)')
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
%Amplitud
AMP=sqrt(Ak.^2+Bk.^2);
%Se le suma un numero muy pequeño que hace que la división dé un número
%grande y distorsione las gráficas
Ak=Ak+0.01;
FASE=atan(-Bk./Ak);
%Partimos la señal para saber en que frecuencias está oscilando la señal
%Reacomodamos el vector de amplitudes
AMPD=AMP(1:N/2);
AMPI=AMP(1+N/2:end);

FASED=FASE(1:N/2);
FASEI=FASE(N/2+1:end);

%Agregamos el armónico cero para completar los espectros
AMP=[AMPI A0 AMPD];
FASE=[FASEI 0 FASED];
%Para caracterizar el eje de frecuencia con los conceptos de F0 y FN
%F0 Que tan espaciados estan
F0=1/(dt*N);
FN=1/(2*dt);
F=[-FN+F0:F0:FN];
subplot(3,1,2)
stem(F,AMP)
title('Espectro de amplitudes')
xlabel('f')
ylabel('Amp')
subplot(3,1,3)
stem(F,FASE)
title('Espectro de fase')
xlabel('f','fontsize', 10)
ylabel('Fase')