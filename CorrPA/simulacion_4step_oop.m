function [I1,I2]=simulacion_4step_oop(phi_esp,alfa,phia,p,R)
% phi_esp es el corrimiento en fase desde la posicion inicial hasta la
% final y fuera del plano en forma simétrica.
% alfa: matriz con desfasajes para casa paso (4 pasos en 2 filas)
% I1: matriz 3D con datos temporales de primera exposición
% I2: matriz 3D con datos temporales de segunda exposición

% Nn=512;
% Nm=512;
% phia=rand(size(phi_esp))*2*pi-pi;

%%
% p = input('Ingrese tamaño promedio del grano en pixeles (1 x def): ');
% if isempty(p)
%     p=1;
% end
% 
% R = input('Ingrese relacion Iref/Iobj (2 x def): '); %R=2,4,6
% if isempty(R)
%     R=2;
% end

% Reserva memoria para los interferogramas
% I1=zeros([size(phi_esp) 4]);      % primera exposición
% I2=zeros([size(phi_esp) 4]);      % segunda exposición
I1=zeros(512,512,4);        % primera exposición
I2=zeros(512,512,4);        % segunda exposición

%% Calcula filtro H e intensidad haz referencia Iref

Ua=1*exp(1i*(phia));
Uf=fftshift(fft2(Ua));
[f c]=size(Uf); % ejes de frec
[FU FV]=freqspace(size(Uf),'meshgrid');     %ejes de frecuencias espaciales relativas entre -1 y 1
FU=FU*(c/2);                                % los ejes de frec esp quedan entre -128 y 127, paso de frec relat a frec no relativas
FV=FV*(f/2);

fc=512/(2*p);         % frecuencia de corte
H=zeros(512);                             % filtro H(u,v)
H(sqrt(FU.^2+FV.^2)<=fc)=1; 

Uf1 = Uf.*H;
Uobj = ifft2(ifftshift(Uf1));  % campo del obj
Iobj = (abs(Uobj)).^2;
Iref = R*mean2(Iobj);  Uref = sqrt(Iref);  % campo haz de ref. Si p=1 y R=2, Uref=1.2538

%% Calcula secuencia de primera exposición
for paso=1:4
    Ua=1*exp(1i*(phia+alfa(1,paso)));
    Uf=fftshift(fft2(Ua));
    Uf1 = Uf.*H;
    Uobj = ifft2(ifftshift(Uf1));  % campo del obj
    U = Uobj + Uref; %campo obj+ref que llega al CCD
    I1(:,:,paso)=abs(U).^2;               % normaliza secuencia completa al final
end

%% Calcula secuencia de segunda exposición
for paso=1:4
    Ua=1*exp(1i*(phia+phi_esp+alfa(2,paso)));
    Uf=ifftshift(fft2(Ua));
    Uf1=Uf.*H;
    Uobj=ifft2(ifftshift(Uf1));  % campo del obj
    U=Uobj + Uref; %campo obj+ref que llega al CCD
    I2(:,:,paso)=abs(U).^2;               % normaliza secuencia completa al final
end

%% Normaliza secuencia completa
I_maximo=max([max(I1(:)) max(I2(:))]);
I_minimo=min([min(I1(:)) min(I2(:))]);

I1=uint8((I1-I_minimo)/(I_maximo-I_minimo)*255);
I2=uint8((I2-I_minimo)/(I_maximo-I_minimo)*255);


return