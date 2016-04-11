
function O=suavizado_spline(G,N,p0,lambda)

x=1:N;
pv(1)=p0;
w(1)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:N

for i=2:N
pv(i)=lambda(k,i)/(1-p0);
w(i)=1;
end

    y=G(k,1:N);
    pp = csaps(x,y,pv,[],w);
    v_fila(k,:) = fnval(pp,x);   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:N

for i=2:N
pv(i)=lambda(i,k)/(1-p0);
w(i)=1;
end

    y=v_fila(1:N,k);
    pp = csaps(x,y,pv,[],w);
    v(k,:) = fnval(pp,x);   
end

O=v';
