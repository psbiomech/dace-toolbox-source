function  [f, df] = regpoly3(S)
%REGPOLY3  Third order polynomial regression function
% Call:    f = regpoly2(S)
%          [f, df] = regpoly2(S)
%          
% S : m*n matrix with design sites
% f =  [1 S S(:,1)*S S(:,2)S(:,2:n) ... S(:,n)^3]
% df : Jacobian at the first point (first row in S) 

% westono@ufl.edu  
% Last update September 15, 2009

[m n] = size(S); %m sets of n design variables
nn = (n+1)*(n+2)*(n+3)/factorial(3);  % Number of columns in f  
% Compute  f
f = [ones(m,1) S zeros(m,nn-n-1)]; %[constant linear other_terms]
dumy1 = 0;
% quadratic terms: x1x1 or x1x2
for  i = 1 : n
    for j = i : n
        dumy1 = dumy1 + 1;
        f(:,n+dumy1+1) = S(:,i).*S(:,j);
    end
end
%cubic terms
dumy2 = 0;
for i = 1 : n
    for j = i : n
        for k = j : n
            dumy2 = dumy2 + 1;
            f(:,n+dumy1+dumy2+1)=S(:,i).*S(:,j).*S(:,k);
        end
    end
end

if nargout > 1
    df = [zeros(n,1)  eye(n)  zeros(n,nn-n-1)];
    dumy1 = 0;
    % quadratic terms: x1x1 or x1x2
    for  i = 1 : n
        for j = i : n
            dumy1 = dumy1 + 1;
            for p = 1 : n
                if p == i && p == j
                    df(p,n+dumy1+1) = 2*S(1,i);
                elseif p == i && p ~=j
                    df(p,n+dumy1+1) = S(1,j);
                elseif p == j && p ~=i
                    df(p,n+dumy1+1) = S(1,i);                    
                end
            end
        end
    end
    %cubic terms
    dumy2 = 0;
    for i = 1 : n
        for j = i : n
            for k = j : n
                dumy2 = dumy2 + 1;
                for p = 1 : n
                    if p == i && p == j && p == k
                        df(p,n+dumy1+dumy2+1) = 3*S(1,i)^2;
                    elseif p == i && p == j && p ~= k
                        df(p,n+dumy1+dumy2+1) = 2*S(1,i)*S(1,k);
                    elseif p == i && p ~= j && p == k
                        df(p,n+dumy1+dumy2+1) = 2*S(1,i)*S(1,j);
                    elseif p ~= i && p == j && p == k
                        df(p,n+dumy1+dumy2+1) = 2*S(1,j)*S(1,i);
                    elseif p == i
                        df(p,n+dumy1+dumy2+1) = S(1,j)*S(1,k);                        
                    elseif p == j
                        df(p,n+dumy1+dumy2+1) = S(1,i)*S(1,k);
                    elseif p == k
                        df(p,n+dumy1+dumy2+1) = S(1,i)*S(1,j);                        
                    end
                end
            end
        end
    end
end