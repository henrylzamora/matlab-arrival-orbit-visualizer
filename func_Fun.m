    function F = Fun(z,t,mu,r1,r2,A) % Lamberts Problem
        warning('off', 'MATLAB:colon:nonRealScalar');
        warning('off', 'MATLAB:colon:operandsNotReal'); % Suppresses the common numeric warning
        warning('off', 'MATLAB:future');
        if size(z,1)==1
            z=z';
        end
        if size(A,1)==1
            A=A';
        end
        if size(r1,1)==1
            r1=r1';
        end
        if size(r2,1)==1
            r2=r2';
        end
        AU = 1.49597871E8; % Astronomical Unit
        TU = 5022643;
        y = r1/AU + r2/AU + A./AU.*(z.*func_stumpS(z) - 1)./sqrt(func_stumpC(z));
        F = (y./func_stumpC(z)).^1.5.*func_stumpS(z) + A./AU.*sqrt(y) - sqrt(mu/AU^3*TU^2)*(t/TU);
        F = F';
    end