clear all, close all, clc

phi=pi-0.1;       % mÃ¡ximo cambio de fase,
p = 1.1; % tamaño medio del grano de speckle
R = 2; % relación Iref/Iobj

%% Calcula fase teÃ³rica para comparar (igual que en las funciones)
% Simula interferogramas de cada exposiciÃ³n y haces objeto y referencia


[I1,I2,I1s,I2s,Iref,phi_esp]=simulacion_correlacion_full_con_scaling(phi,p,R);


%% 

L=24; % Tamaño de la ventana de análisis
Fo=genero_matriz_de_correccion(I1,I2,L);

%% Calcula fase en base al coeficiente de correlacion para varios coeficientes 
% EvalÃºa usando ventanas sin procesamiento paralelo
% Go Order Statistics Correlation Coefficient
     
[Go]=corr_order_statis_aprox(phi_esp,I1,I2,L);

save coeficientes

%% Calcula la fase ingresando los parametros externos

alfa=0.013; % parámetros de ajuste lineal del factor de corrección
beta=0;

tic
[F_ideal,Goc,Gocs]=correction(phi_esp,alfa,beta,L,Fo,Go);
toc

save fase

%% Procesa usando phase-shifting

% Simula secuencia de interferogramas para phase-shifting
sh=phi;
phia=rand(512)*2*pi-pi;

% Simula errores en phase shift
shift1=[0 sh*(1+0.012) 2*sh*(1-0.012) 3*sh*(1+0.012)];
shift2=shift1;
[I1s1,I2s1]=simulacion_4step_oop(phi_esp,[shift1;shift2],phia,p,R);

fase_carre1=mide_fase_carre(I1s1,I2s1);

n=size(Go,1);
Fase_Carre_recuperada(1:n,1:n)=fase_carre1(L/2+1:L/2+n,L/2+1:L/2+n);

save carre Fase_Carre_recuperada

figure
plot(1:n,squeeze(F_ideal(:,n/2)),'b',1:n,squeeze(Fase_Carre_recuperada(:,n/2)),'k',1:n,squeeze(Gocs(:,n/2)),'r')











