

function Fo=genero_matriz_de_correccion(I1,I2,L)

I1=double(I1);
I2=double(I2);

[Nfil Ncol]=size(I1);

tic
for col=fix(L/2)+1:Ncol-fix(L/2)
    for fil=fix(L/2)+1:Nfil-fix(L/2)
        
        X=I1(fil-fix(L/2):fil+fix(L/2),col-fix(L/2):col+fix(L/2));
        Y=I2(fil-fix(L/2):fil+fix(L/2),col-fix(L/2):col+fix(L/2));
        
        X=X(:); 
        Y=Y(:); 
              
        Correct(fil,col)=sqrt(mean(X.^2)*mean(Y.^2))/mean(X.*Y);

    end
end
toc

Fo=Correct(fix(L/2)+1:Nfil-fix(L/2),fix(L/2)+1:Ncol-fix(L/2));
% save Fo Fo