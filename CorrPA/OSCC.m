

function c=OSCC(X,Y)

% Implementacion Order Statistics Correlation Coefficient (OSCC):

    [x,orden_x] = sort(X);
    y = sort(Y);
    
    x_inv=fliplr(x');
    
    alfa=x-x_inv';
    num=sum(alfa.*Y(orden_x));
    den=sum(alfa.*y);
    c=num/den;