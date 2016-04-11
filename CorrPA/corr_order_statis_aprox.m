

function Go=corr_order_statis_aprox(phi_esp,I1,I2,L)
% function [Gp,Gk,Gs,Go]=corr_order_statis_aprox(phi_esp,I1,I2,L)
% Calcula fase en base al coeficiente de correlaci�n
% Evalúa usando ventanas sin procesamiento paralelo

I1=double(I1);
I2=double(I2);

[Nfil Ncol]=size(I1);

tic
for col=fix(L/2)+1:Ncol-fix(L/2)
    parfor fil=fix(L/2)+1:Nfil-fix(L/2)
        
        % Coeficiente de correlacion
        X=I1(fil-fix(L/2):fil+fix(L/2),col-fix(L/2):col+fix(L/2));
        Y=I2(fil-fix(L/2):fil+fix(L/2),col-fix(L/2):col+fix(L/2));
        
        X=X(:);
        Y=Y(:);
        
%         % Coeficiente de Pearson
          %cp=corr(X,Y,'type','Pearson'); % Sensible a la transformaci�n de X e Y
       
%         % Coeficiente de Kendall
%          c=corr(X,Y,'type','Kendall'); % Anda muy bien en las partes estables de la fase       
%          ck=sin(pi*c/2);
%         
%         % Coeficiente de Spearman
%          c=corr(X,Y,'type','Spearman');
%          cs=2*sin(pi*c/6);
% 
        % Order Statistics Correlation Coefficient
         co=OSCC(X,Y);
         

% Paso a matriz:
%              gp(fil,col)=cp;
%              gk(fil,col)=ck;
%              gs(fil,col)=cs;
              go(fil,col)=co;
            
    end
%     col
end
toc

% Gp=gp(fix(L/2)+1:Nfil-fix(L/2),fix(L/2)+1:Ncol-fix(L/2));
% Gk=gk(fix(L/2)+1:Nfil-fix(L/2),fix(L/2)+1:Ncol-fix(L/2));
% Gs=gs(fix(L/2)+1:Nfil-fix(L/2),fix(L/2)+1:Ncol-fix(L/2));
Go=go(fix(L/2)+1:Nfil-fix(L/2),fix(L/2)+1:Ncol-fix(L/2));


% save Go Go
%figure; imagesc(Go); colorbar

