function C = func_stumpC(Z)
    warning('off', 'MATLAB:colon:nonRealScalar');
    warning('off', 'MATLAB:colon:operandsNotReal'); % Suppresses the common numeric warning
    warning('off', 'MATLAB:future');
        % This function evaluates the Stumpff function C(z).
        if size(Z,1)==1
            Z=Z';
        end
        C=zeros(size(Z));
        for i=1:size(Z)
            z=Z(i,1);
            if z>0
                C(i,1) = (1 - cos(sqrt(z)))./z;
            elseif z<0
                C(i,1) = (cosh(sqrt(-z)) - 1)./(-z);
            else
                C(i,1) = 1/2;
            end
        end
    end