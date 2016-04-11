

function [F_ideal,Goc,Gocs]=correction(phi_esp,alfa,beta,L,Fo,Go)
%function [F_ideal,Goc,Gocs,Gpc,Gpcs,Gkc,Gkcs,Gsc,Gscs]=correction(phi_esp,alfa,beta,L,Fo,Gp,Gs,Gk,Go)

Nfil=size(phi_esp,1);
n=size(Fo,1);

fos=suavizado_spline(Fo,n,0.00003,ones(n));

% Gpc=Gp.*sqrt(1+alfa*mat2gray(fos)+beta);
% Gkc=Gk.*sqrt(1+alfa*mat2gray(fos)+beta);
% Gsc=Gs.*sqrt(1+alfa*mat2gray(fos)+beta);
Goc=Go.*sqrt(1+alfa*mat2gray(fos)+beta);
%

%  Gp=acos(Gp);
%  Gpc=real(acos(Gpc));
%  Gk=acos(Gk);
%  Gkc=real(acos(Gkc));
%  Gs=acos(Gs);
%  Gsc=real(acos(Gsc));
 Go=acos(Go);
 Goc=real(acos(Goc));

n=size(Go,1);
F_ideal(1:n,1:n)=phi_esp(L/2+1:L/2+n,L/2+1:L/2+n);

% figure
% imagesc(Gpc-F_ideal); colorbar, title('fase con Pearson - fase ideal')
% figure
% imagesc(Gkc-F_ideal); colorbar, title('fase con Kendall - fase ideal')
% figure
% imagesc(Gsc-F_ideal); colorbar, title('fase con Spearman - fase ideal')
% figure
% imagesc(Goc-F_ideal); colorbar, title('fase con OSCC - fase ideal')




N=size(Goc,1);
% Gpcs=suavizado_spline(Gpc,N,0.1,ones(N));
% Gkcs=suavizado_spline(Gkc,N,0.1,ones(N));
% Gscs=suavizado_spline(Gsc,N,0.1,ones(N));
Gocs=suavizado_spline(Goc,N,0.1,ones(N));

% figure
% plot(1:n,F_ideal(1:n,(Nfil-L)/2),'b',1:n,Gp(1:n,(Nfil-L)/2),'k',1:n,Gk(1:n,(Nfil-L)/2),'g',1:n,Gs(1:n,(Nfil-L)/2),'r',1:n,Go(1:n,(Nfil-L)/2),'m'),legend('Fase ideal','Pearson','Kendall', 'Spearman', 'OSCC'), title('Sin el suavisado fc');

% figure
% plot(1:n,F_ideal(1:n,(Nfil-L)/2),'b',1:n,Gpc(1:n,(Nfil-L)/2),'k',1:n,Gp(1:n,(Nfil-L)/2),'g',1:n,Gpcs(1:n,(Nfil-L)/2),'r'),legend('Fase ideal','Pearson con alfa y beta','Pearson', 'Pearson suavisado');
% figure
% plot(1:n,F_ideal(1:n,(Nfil-L)/2),'b',1:n,Gkc(1:n,(Nfil-L)/2),'k',1:n,Gk(1:n,(Nfil-L)/2),'g',1:n,Gkcs(1:n,(Nfil-L)/2),'r'),legend('Fase ideal','Kendall con alfa y beta','Kendall', 'Kendall suavisado');
% figure
% plot(1:n,F_ideal(1:n,(Nfil-L)/2),'b',1:n,Gsc(1:n,(Nfil-L)/2),'k',1:n,Gs(1:n,(Nfil-L)/2),'g',1:n,Gscs(1:n,(Nfil-L)/2),'r'),legend('Fase ideal','Spearman con alfa y beta','Spearman', 'Spearman suavisado');
figure
plot(1:n,F_ideal(1:n,(Nfil-L)/2),'b',1:n,Goc(1:n,(Nfil-L)/2),'k',1:n,Go(1:n,(Nfil-L)/2),'g',1:n,Gocs(1:n,(Nfil-L)/2),'r'),legend('Fase ideal','OSCC con alfa y beta','OSCC', 'OSCC suavisado');


