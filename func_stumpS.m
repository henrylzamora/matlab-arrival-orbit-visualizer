function S = func_stumpS(Z)
        % This function evaluates the Stumpff function S(z).
        if size(Z,1)==1
            Z=Z';
        end
        S=zeros(size(Z));    
        for i=1:size(Z)
            z=Z(i,1);
            if z>0
                S(i,1) = (sqrt(z) - sin(sqrt(z)))./(sqrt(z)).^3;
            elseif z<0
                S(i,1) = (sinh(sqrt(-z)) - sqrt(-z))./(sqrt(-z)).^3;
            else
                S(i,1) = 1/6;
            end
        end
    end