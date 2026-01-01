function window=rcoswindow(beta,SymLen)
    % t=0:(1+beta)*SymLen;
    % window = zeros(1,(1+beta)*SymLen);
    % for i=1:beta*SymLen
    %     window(i)=0.5+0.5*cos(pi+t(i)*pi/(beta*SymLen));
    % end
    % window(beta*SymLen+1:SymLen)=1;
    % for j=SymLen+1:(1+beta)*SymLen+1
    %     window(j-1)=0.5+0.5*cos((t(j)-SymLen)*pi/(beta*SymLen));
    % end
    t=0:SymLen;
    window = zeros(1,SymLen);
    for i=1:beta*SymLen
        window(i)=0.5+0.5*cos(pi+t(i)*pi/(beta*SymLen));
    end
    window(beta*SymLen+1:SymLen)=1;
end