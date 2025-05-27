function [norm_e,a,i,w,Om,TA,gam,norm_h] = func_RVtoOE(r,v,mu)
        % Calculates Orbital Elements from Position and Inertial Velocity Vectors.
        i1 = [1,0,0]';     i3 = [0,0,1]';
        
        h = cross(r,v);
        p3 = h/norm(h);
        e = cross(v,h)/mu - r/norm(r);
        norm_e=norm(e);
        p1 = e/norm(e);
        N = cross(i3,p3)/norm(cross(i3,p3));
        if N(2)>=0
            Om = acosd(i1'*N);
        else
            Om = 360-acosd(i1'*N);
        end
        i = acosd(i3'*p3);
        if p1(3)>=0
            w = acosd(p1'*N);
        else
            w = 360-acosd(p1'*N);
        end
        vr = r'*v;
        if vr>=0
            TA = acosd(p1'*(r/norm(r)));
        else
            TA = 360-acosd(p1'*(r/norm(r)));
        end
        gam = 90 - acosd((r/norm(r))'*(v/norm(v)));
        a = norm(h)^2/mu/abs(1-norm_e^2);
        norm_h = norm(h);
    end