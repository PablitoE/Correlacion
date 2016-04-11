function F=mide_fase_carre(I1,I2)
% Calcula fase usando el método de Carré
%   I1: primera exposición, 4 pasos (ver simulacion_4step() )
%   I2: segunda exposición, 4 pasos

I1=double(I1);
I2=double(I2);

%% Calcula numerador y denominador de atan() de cada exposición
[sen_f1,cos_f1]=numden_atan(I1);
[sen_f2,cos_f2]=numden_atan(I2);

%% Calcula tan() de la diferencia de fases directamente
sen_fi=sen_f2.*cos_f1-cos_f2.*sen_f1;
cos_fi=cos_f2.*cos_f1+sen_f2.*sen_f1;

% Filtra parte del ruido
sen_fi=conv2(sen_fi,ones(16)/256,'same');
cos_fi=conv2(cos_fi,ones(16)/256,'same');

%% Fase resultante
F=atan2(sen_fi,cos_fi);
% F=F-min(F(:));
%%
return

%%

function [sen_fi,cos_fi]=numden_atan(I)
% Calcula numerador y denominador de atan()

I23=I(:,:,2)-I(:,:,3);
I14=I(:,:,1)-I(:,:,4);

sen_fi=sqrt(abs((I14+I23).*(3*I23-I14)));       % agregado sqrt(abs())
signo_I=sign(I23);              % función sign()
signo_I(signo_I==0)=1;          % corrección x>=0 => sign=1 (ver Rastoghi pg 67)
sen_fi=signo_I.*sen_fi;

cos_fi=I(:,:,2)+I(:,:,3)-I(:,:,1)-I(:,:,4);

return
