function[a, b, mn]= minInd(A)
    [n,m] = size(A);
    u =ones(n,1);
    for i = 1:n
        for j = 1:m
            if(A(i,j) == 0)
                A(i,j) = 16;
            else
                 A(i,j) = A(i,j);
            end
        end
    end
    for i =1:n 
        u(i) = min(A(1:n, i));
    end
    mn = min(u);
    for i =1:n
        for j=1:m
            if A(i,j) == mn
                a = i;
                b = j;
            end
        end
    end
  end
  
  