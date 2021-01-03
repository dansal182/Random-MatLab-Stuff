% Roberto Tello 154706
% Jos� Pliego 157103
% Daniel Salnikov 159817

% This algorithm implements the simplex method for a problem of the type
% max c^T x subject to Ax=b.
% To do this, we first use the phase 1 algorithm to check if the problem
% is feasible, and solve an auxiliary program to obtain an initial basic
% feasible solution and its feasible basis.
% We then use this basic feasible solution and the feasible basis to 
% implement the regular simplex method, using Bland's Rule as the pivoting
% method to avoid cycling.
% 
% In case the feasible set is empty, all the outputs besides status are
% set to NaN (Not a Number):
% If the problem is feasible but unbounded, the optimal value is set to Inf
% and the other outputs excluding status are set to NaN.

function[status, obasis, obfs, oval] = bothPhases(A, b, c)
% maximise c^T x
% subject to Ax = b, x >= 0, b >=0
%
% Input:
% A mxn matrix with m <= n and rank of A is m
% b column vector with m rows
% c column vector with n rows
%
% Output:
% status = -1 if the feasible set is empty
% status = 0 if the feasible set is non-empty but the problem is unbounded
%   (there is no optimal solution)
% status = 1 if the problem is bounded (there is an optimal solution)
% obasis = a vector of size m of indices of an optimal feasible basis for 
%   the problem if the feasible set is non-empty and the problem is 
%   bounded (in terms of a set of indices of column vectors)
% obfs = a vector of size n which is the optimal basic feasible solution 
%   corresponding to this optimal basis if the feasible set is non-empty
%   and the problem is bounded
% oval = the objective value of this optimal basic feasible solution
% (if the feasible set is non-empty and the problem is bounded)


function[nvac, basis, bfs] = phaseOne(A, b, c)
% maximise c^T x
% subject to Ax = b, x >= 0, b >=0
%
% Input:
% A mxn matrix with m <= n and rank of A is m
% b column vector with m rows
% c column vector with n rows
%
% Output:
% nvac = 0 if the feasible set is empty
% nvac = 1 if the feasible set is non-empty
% basis = a vector of size m of indices of column vectors for a 
%   feasible basis for the problem if the feasible set is non-empty
% bfs = a vector of size n of the basic feasible solution corresponding
%   to this basis (if the feasible set is non-empty)
    
    v = b;
    
    % All variables of vector b are set no non negative
    negatives = b < 0;
    for i = 1:length(b)
        if negatives(i) == 1
            b(i) = b(i)*-1;
        end
    end
    
    % The dimensions of matrix A are taken
    [m,n] = size(A);
    B = A;
    
    % Matrix A is set expanded
    A = [A,eye(m)];
    
    % The new vector c for the auxiliary program is created
    c = [zeros(1,n),-1*ones(1,m)];
    
    % The auxiliary linear program is solved
    [~, basis, bfs, val] = phaseTwo(A,b,c',n+1:m+n,[zeros(1,n),b']');
    
    bfs = bfs(1:n);
    
    % Feasibility is checked
    if val==0
        nvac = 1;
    else
        nvac = 0;
    end
end

function[bound, obasis, obfs, oval] = phaseTwo(A, b, c, sbasis, sbfs)
% maximise c^T x
% subject to Ax = b, x >= 0, b >=0
%
% Input:
% A mxn matrix with m <= n and rank of A is m
% b column vector with m rows
% c column vector with n rows
% sbasis a vector of size m of indices of column vectors for a feasible 
%   basis for this problem from which to start the simplex method
% sbfs a vector of size n which is the basic feasible solution 
%   corresponding to this basis
%
% Output:

% bound = 0 if the problem is unbounded (there is no optimal solution)
% bound = 1 if the problem is bounded (there is an optimal solution)
% obasis = a vector of size m of indices of column vectors which gives an
%   optimal feasible basis for the problem if the problem is bounded
% obfs = a vector of size n which is the optimal basic feasible solution
%   corresponding to this optimal basis if the problem is bounded
% oval = the objective value of this optimal basic feasible solution
%   (if the problem is bounded) 

        obasis = sbasis;
        obfs = sbfs;
        oval = dot(c,obfs);
        bound = 1;
    
        % Matrix Ab (containing the basic columns is created)
        Ab =[];
        for i = obasis
            Ab = [Ab, A(:,i)];
        end
        
        % It is verified that the matrix Ab is non singular, and if it is
        % the simplex tableaux is created. 
    
        if det(Ab) ~= 0
            [m,n] = size(A);
            
            % Matrix An is created
            basen = 1:n;
            N = [];
            for i = basen
                if ~ismember(i,obasis)
                    N = [N, i];
                end
            end
            An = [];
            for i = N
                An = [An, A(:,i)];
            end
            
            % Matrix Q is created
            Q = -Ab\An;
            
            % vector p is created
            p = Ab\b;
            
            if any(p<0)
                return
            end
            
            % Auxiliary vector cb and cn are created
            cb = [];
            for i = obasis
                cb = [cb, c(i)];
            end
            cb = cb';
            
            cn = [];
            for i = N
                cn = [cn, c(i)];
            end
            cn = cn';
            
            % Objective value is obtained
            z0 = dot(cb,Ab\b); 
            
            % Coefficient vector r
            r = cn - (cb'*(Ab\An))';
            
            % If all the coefficients are nonpositive the optimum has been
            % reached
            while any(r > 0)
                % The variable that increases the optimum value is 
                %   selected using Bland's rule        
                
                t = r > 0;
                index = find(t,1,'first');              
                enters = N(index);
                
                % m is the number of rows from the original matrix
                % Here the epsilon that determines which variable exits 
                %   the basis is determined
                minimus = Inf;
                Qr = Q(:,index);
                leaves = 0;
                
                for i = 1:m
                   if Qr(i) < 0                    
                       if minimus > -p(i)/Qr(i)
                           minimus = -p(i)/Qr(i);
                           leaves = obasis(i);                          
                       end
                   end
                end
                
                % If no element is elegible to leave the basis then the
                % program is unbounded
                if leaves == 0
                    bound = 0;                 
                    oval = Inf;
                    return
                else
                    bound = 1;
                end
                
                
                % We take from the basis the one that leaves
                leaves = ismember(obasis,leaves);
                obasis(leaves) = [];
                
                % We insert the new element into the basis
                obasis = [obasis,enters];
                obasis = sort(obasis);
                
                Ab =[];
                for i = obasis
                    Ab = [Ab, A(:,i)];
                end                
               
                inverAb = inv(Ab);
            
                % Matrix An is created
                basen = 1:n;
                N = [];
                for i = basen
                    if ~ismember(i,obasis)
                        N = [N, i];
                    end
                end
                N = sort(N);
                
                An = [];
                for i = N
                    An = [An, A(:,i)];
                end
            
                % Matrix Q is created
                Q = -Ab\An;
            
                % vector p is created
                p = Ab\b;
            
                % Auxiliary vector cb and cn are created
                cb = [];
                for i = obasis
                    cb = [cb, c(i)];
                end
                cb = cb';
            
                cn = [];
                for i = N
                    cn = [cn, c(i)];
                end
                cn = cn';
            
            
                % Objective value is obtained
                z0 = dot(cb,Ab\b);
                obfs = zeros(1,n);
                
                j = 0;
                for i = obasis
                    j = j+1;
                    obfs(i) = p(j);   
                end
                
                oval = z0;
            
                % Coefficient vector r
                r = cn - (cb'*(Ab\An))';  
            end            
            
        else
            return          
        end
end

% We implement phase 1 of the algorithm.
[nvac, basis, bfs] = phaseOne(A, b, c);

% We implement phase two if the feasible set is not empty.
if nvac == 1
    [bound, obasis, obfs, oval] = phaseTwo(A, b, c, basis, bfs);
    if bound == 0 % The problem is feasible but unbounded.
        status = 0;
        obasis = NaN;
        obfs = NaN;
        oval = Inf;
        disp("The feasible set is not empty but the problem is unbounded.");
        return;
    else
        status = 1; % The problem is feasible and bounded.
        return;
    end
else
    status = -1;
    obasis = NaN;
    obfs = NaN;
    oval = NaN;
    disp("The feasible set is empty.");
    return;
end
    



end