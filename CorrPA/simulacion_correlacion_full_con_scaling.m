% FUNCTION [I1,I2,I1s,I2s,Iref,phi_esp]=simulacion_correlacion_full_con_scaling(N,phi,p,R,tipo)
% Simula las imágenes de interferometría de speckle deestado inicial y
% final para distintas situaciones de simulación.
% Entradas:
%   - N: tamaño de las imágenes
%   - phi: valor pico de desplazamiento
%   - p: tamaño medio de grano de speckle en pixeles
%   - R: relación entre intensidad de referencia y objeto: Iref/Iobj
%   - tipo: tipo de señal [1: sombrero gaussiano; 2: tilt; 3: sombrero senoidal; 4: peaks; 5: damero]
%
% Salidas:
%   - I1: interferograma de estado inicial
%   - I2: interferograma de estado final
%   - I1s: imagen de haz de objeto en estado inicial
%   - I2s: imagen de haz de objeto en estado final
%   - Iref: imagen del haz de referencia
%   - phi_esp: distribución de fase de desplazamientos simulada

function [I1,I2,I1s,I2s,Iref,phi_esp]=simulacion_correlacion_full_con_scaling(N,phi,p,R,tipo)
% phi es el corrimiento en fase desde la pos inicial hasta la final y fuera
% del plano en forma simï¿½trica. 
Nn = N;
Nm = N;
phia = rand(N)*2*pi-pi;

% ---------------
% INTERFEROGRAMA DEL OBJ INSTANTE t=0
H = zeros(N); % filtro H(u,v)

phi_esp = zeros(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIPO
switch tipo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Para un sombrero
    case 1
        A=50*N; B=64;
        for m=1:N
            for n=1:N
                x2=(n-Nn/2)^2;
                y2=(m-Nm/2)^2;
                phi_esp(m,n) = B*2*pi*exp((-x2-y2)/A);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Para un tilt
    case 2
        for m=1:N
            for n=1:N
                x2=(n-Nn/2);
                y2=(m-Nm/2);
                phi_esp(m,n) = pi*(x2^2+y2^2+x2*y2*0.5);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Para un sombrero
    case 3
        [X,Y] = meshgrid(1:N);
        X=X/N;
        Y=Y/N;
        R = sqrt((2*X-1).^2 + (2*Y-1).^2) + eps;
        phi_esp = sin(8*R)./R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% peaks
    case 4
        u = linspace(-3.5,3.5,N);
        [x, y] = meshgrid(u);
        
        phi_esp = 5*(1-x.^2).*exp(-x.^2-(y+1).^2)-20*(x/2-x.^3-y.^5).*exp(-x.^2-y.^2)...
            -exp(-(x+1).^2-y.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% damero
    case 5
        sizeCell = ceil(N/8);
        phi_esp = checkerboard(sizeCell);
        phi_esp = phi_esp(1:N,1:N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalización de distribución de fase
phi_esp_min=min(phi_esp(:));
phi_esp_max=max(phi_esp(:));
phi_esp = (phi_esp-phi_esp_min)/phi_esp_max*phi;

% Simulación del campo de speckle
Ua = 1*exp(1i*phia);
Uf = ifftshift(fft2(Ua));
[f, c] = size(Uf); % ejes de frec
[FU, FV] = freqspace(size(Uf),'meshgrid');     %ejes de frecuencias espaciales relativas entre -1 y 1
FU = FU*(c/2); % los ejes de frec esp quedan entre -128 y 127, paso de frec relat a frec no relativas
FV = FV*(f/2);
fc = N/(2*p);
H(sqrt(FU.^2+FV.^2)<=fc) = 1; 
Uf1 = Uf.*H;
Uobj = ifft2(ifftshift(Uf1));  % campo del obj

% Intensidad de objeto
Iobj = (abs(Uobj)).^2;
I1s = Iobj;

% Intensidad de referencia
Iref = R*mean(Iobj(:));  
Uref = sqrt(Iref);  % campo haz de ref. Si p=1 y R=2, Uref=1.2538

% Interferencia
U = Uobj + Uref; %campo obj+ref que llega al CCD
I1 =abs(U).^2; % intensidad incidente en el CCD

% Simulación del desplazamiento
Ua = 1*exp(1i*(phia+phi_esp));
Uf = ifftshift(fft2(Ua));
Uf1 = Uf.*H;
Uobj = ifft2(ifftshift(Uf1));  % campo del obj
U = Uobj + Uref; %campo obj+ref que llega al CCD
I2 = (abs(U)).^2; % intensidad incidente en el CCD
Iref = Iref*ones(N);
I2s = abs(Uobj).^2;

%Normalizacion y discretizacion
M=max(max([I1 I2 I1s I2s Iref]));
m=min(min([I1 I2 I1s I2s Iref]));

I1=uint8(255/(M-m)*(I1-m));
I2=uint8(255/(M-m)*(I2-m));
I1s=uint8(255/(M-m)*(I1s-m));
I2s=uint8(255/(M-m)*(I2s-m));
Iref=uint8(255/(M-m)*(Iref-m));