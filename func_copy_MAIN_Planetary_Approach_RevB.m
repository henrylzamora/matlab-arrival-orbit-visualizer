function [coe,r_vec_SOI_Entry_a_I, RA_Entry_a, DEC_Entry_a]  = func_copy_MAIN_Planetary_Approach_RevB(v_inf_Entry_I,zp,var,Primary,input)
% Choose between:
% Orbital Inclination (i),
% Right Ascension of Ascending Node (RAAN),
% Argument of Periapsis (ARGP).
% var   : Numeric value of angle.
% input : Text ('i','RAAN','ARGP') used to control which loop is executed.
%
%
% When Multiple Solutions are possible, the first choice, "a" is used as
% coe = [e, a, ia, ARGPa, RAANa,TA_SOI_Entry, h];
%
% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% INPUTS


User_Input = input;
fprintf('%s is the selected angle\n',input)
%% ORBIT DETERMINATION
% Planet's Properties
planetparam = func_Primary_Parameters_RevC(Primary);
RE = planetparam.RE;
RP = planetparam.RP;
mu = planetparam.mu;                Color_Primary     = planetparam.Color_Primary;
r_SOI = planetparam.r_SOI;          Color_Primary_Neg = planetparam.Color_Primary_Neg;



v1 = v_inf_Entry_I(1);
v2 = v_inf_Entry_I(2);
v3 = v_inf_Entry_I(3);
vinf = norm(v_inf_Entry_I);   %[km/s] Excess Speed
rp = zp + RE;                 %[km] Periapsis Radius
a = mu/vinf^2;                %[km] Semi-major Axis
e = rp/a+1;                   %[--] Eccentricity
h = mu/vinf*sqrt(e^2-1);      %[km^2/s] Angular Momentum Vector Norm
En = vinf^2/2;                %[km^2/s^2] Orbital Energy
Turn_Angle = 2*asind(1/e);    %[deg] Turn Angle
Aiming_Radius = a*sqrt(e^2-1);%[km] Aiming Radius
TA_inf_Entry = -acosd(-1/e);  %[deg]
TA_inf_Exit  =  acosd(-1/e);  %[deg]

% Inclination Limits given current v_inf_Entry_I
i_min = asind(abs(v3)/vinf);
i_max = 180-asind(abs(v3)/vinf);
fprintf(['Orbital Inclination:\n' ...
         '         %.1f < i < %.1f deg,\n\n'],i_min,i_max)


% ARGP Limits given current v_inf_Entry_I
w(1) = mod((      asind(abs(v3)/vinf) - TA_inf_Entry),360);
w(2) = mod((180 - asind(abs(v3)/vinf) - TA_inf_Entry),360);
w(3) = mod((180 + asind(abs(v3)/vinf) - TA_inf_Entry),360);
w(4) = mod((360 - asind(abs(v3)/vinf) - TA_inf_Entry),360);

fprintf(['Argument of Periapsis:\n' ...
         '      %4.1f <  ARGP < %4.1f deg  or\n' ...
         '      %4.1f <  ARGP < %4.1f deg  or\n\n'],w(1),w(2),w(3),w(4))

%% RAAN
if strcmp(User_Input,'RAAN')
    RAAN = var;
    % RAAN = 0:10:350;
    K = size(RAAN,2);
    for k=1:K
        fprintf('For RAAN = %.1f deg, a unique trajectory exists.\n',RAAN(k))
    
        % Right Ascension of the Ascending Node
        RAANa(k) = RAAN(k);
        RAANb(k) = RAAN(k);
    
        % Line of Nodes
        N_vec_a_I = [cosd(RAANa(k)); sind(RAANa(k)); 0];
        N_vec_b_I = [cosd(RAANb(k)); sind(RAANb(k)); 0];
        
        % Orbital Inclination
        if (v2*cosd(RAANa(k)) - v1*sind(RAANa(k)))/v3>=0
            ia(k) = asind( abs(v3) / sqrt( ( v2*cosd(RAANa(k)) - v1*sind(RAANa(k)) )^2 + v3^2 ) );
        else
            ia(k) = 180 - asind( abs(v3) / sqrt( ( v2*cosd(RAANa(k)) - v1*sind(RAANa(k)) )^2 + v3^2 ) );
        end
        ib(k) = ia(k);

        % p3 versor and momentum vector
        h1a =  h*sind(ia(k))*sind(RAANa(k));
        h2a = -h*sind(ia(k))*cosd(RAANa(k));
        h3a =  h*cosd(ia(k));
        
        h_vec_a_I = [h1a;h2a;h3a];
        h_vec_b_I = h_vec_a_I;

        p3a_I = h_vec_a_I/h;
        p3b_I = h_vec_b_I/h;

    
        % p1 versor and eccentricity vector
        p1a_I = cross(v_inf_Entry_I,h_vec_a_I)/mu/e + v_inf_Entry_I/vinf/e;
        p1b_I = cross(v_inf_Entry_I,h_vec_b_I)/mu/e + v_inf_Entry_I/vinf/e;
    
        % Argument of Periapsis
        if p1a_I(3)>0
            ARGPa(k) =     acosd(N_vec_a_I'*p1a_I);
        else
            ARGPa(k) = 360-acosd(N_vec_a_I'*p1a_I);
        end
        if p1b_I(3)>0
            ARGPb(k) =     acosd(N_vec_b_I'*p1b_I);
        else
            ARGPb(k) = 360-acosd(N_vec_b_I'*p1b_I);
        end  

    end
end
%% INCLINATION
if strcmp(User_Input,'i')
%     if i>i_max % Check if user values are possible and correct if needed.
%     fprintf('Orbital Inclination to high, lowering to max possible \n')
%     i = i_max;
% elseif i <i_min
%     fprintf('Orbital Inclination to low , raising to min possible \n')
%     i = i_min;
% end
i = var;
    % i = i_min:(i_max-i_min)/20:i_max;
    K = size(i,2);
    for k=1:K
        if i(k)<i_min || i(k)>i_max 
            fprintf('For i = %.1f deg, no possible trajectory exists. Choose a valid orbital inclination.\n',i(k))
            return
        else
            fprintf('For i = %.1f deg, two possible trajectories exist.\n',i(k))  
        end
        % Orbital Inclination
        ia(k) = i(k);
        ib(k) = i(k);
    
        % RAAN
        % if v1~=0
        %     A = (v2/v1)^2+1;
        %     B = -2*v2*v3/v1^2*cosd(i(k))/sind(i(k));
        %     C = (v3/v1)^2*(cosd(i(k))/sind(i(k)))^2-1;
        %     discr = B^2-4*A*C;
        %     if abs(discr)<1e-12
        %         discr=0;
        %     end
        %     cosRAANa = -B/(2*A) + 1/(2*A)*sqrt(discr);
        %     sinRAANa = v2/v1*cosRAANa - v3/v1*cosd(ia(k))/sind(ia(k));
        %     cosRAANb = -B/(2*A) - 1/(2*A)*sqrt(discr);
        %     sinRAANb = v2/v1*cosRAANb - v3/v1*cosd(ib(k))/sind(ib(k));
        % end
        % if v2~=0
        %     A = (v1/v2)^2+1;
        %     B = 2*v1*v3/v2^2*cosd(i(k))/sind(i(k));
        %     C = (v3/v2)^2*(cosd(i(k))/sind(i(k)))^2-1;
        %     discr = B^2-4*A*C;
        %     if abs(discr)<1e-12
        %         discr=0;
        %     end
        %     sinRAANa = -B/(2*A) - 1/(2*A)*sqrt(discr);
        %     cosRAANa = v1/v2*sinRAANa + v3/v2*cosd(ia(k))/sind(ia(k));
        %     sinRAANb = -B/(2*A) + 1/(2*A)*sqrt(discr);
        %     cosRAANb = v1/v2*sinRAANb + v3/v2*cosd(ib(k))/sind(ib(k));      
        % end
        % RAANa(k) = atan2d(sinRAANa,cosRAANa);
        % if RAANa(k)<0
        %     RAANa(k)=RAANa(k)+360;
        % end
        % RAANb(k) = atan2d(sinRAANb,cosRAANb);
        % if RAANb(k)<0
        %     RAANb(k)=RAANb(k)+360;
        % end

        h3(k)  = h*cosd(i(k));
        discr = vinf^2*sind(i(k))^2-v3^2;
        if abs(vinf^2*sind(i(k))^2-v3^2) < 1E-12
            discr = 0;
        end
        if v1~=0
            h2a(k) = h*(-v2*v3*cosd(i(k)) + v1*sqrt(discr))/(v1^2+v2^2);
            h1a(k) = -v2/v1*h2a(k) - v3/v1*h3(k);
        
            h2b(k) = h*(-v2*v3*cosd(i(k)) - v1*sqrt(discr))/(v1^2+v2^2);
            h1b(k) = -v2/v1*h2b(k) - v3/v1*h3(k);
         
            h_vec_a_I(:,k) = [h1a(k);h2a(k);h3(k)];
            h_vec_b_I(:,k) = [h1b(k);h2b(k);h3(k)];
        end
        if v2~=0
            h1a(k) = h*(-v1*v3*cosd(i(k)) - v2*sqrt(discr))/(v1^2+v2^2);
            h2a(k) = -v1/v2*h1a(k) - v3/v2*h3(k);
        
            h1b(k) = h*(-v1*v3*cosd(i(k)) + v2*sqrt(discr))/(v1^2+v2^2);
            h2b(k) = -v1/v2*h1b(k) - v3/v2*h3(k);
         
            h_vec_a_I(:,k) = [h1a(k);h2a(k);h3(k)];
            h_vec_b_I(:,k) = [h1b(k);h2b(k);h3(k)];
        end
        
        if v1==0 && v2==0
            i = 90;
            h3 = 0;
        end
        
        % p3-versor (Perifocal RF)
        p3a_I(:,k) = h_vec_a_I(:,k)/h;
        p3b_I(:,k) = h_vec_b_I(:,k)/h;
        
        % Line of Nodes
        N_vec_a_I(:,k) = [-h2a(k)/(h*sind(i(k)));h1a(k)/(h*sind(i(k)));0];
        N_vec_b_I(:,k) = [-h2b(k)/(h*sind(i(k)));h1b(k)/(h*sind(i(k)));0];
        
        % Right Ascension of the Ascending Node (RAAN)
        if N_vec_a_I(2,k)>0
            RAANa(k) = acosd(-h2a(k)/(h*sind(i(k))));
        else
            RAANa(k) = 360-acosd(-h2a(k)/(h*sind(i(k))));
        end
        if N_vec_b_I(2,k)>0
            RAANb(k) = acosd(-h2b(k)/(h*sind(i(k))));
        else
            RAANb(k) = 360-acosd(-h2b(k)/(h*sind(i(k))));
        end

         (v1*sind(RAANa(k)) - v2*cosd(RAANa(k)))*sind(ia(k)) + v3*cosd(ia(k))

        % Eccentricity Vector
        e_vec_a_I(:,k) = [(v2*h3(k) -v3*h2a(k))/mu+v1/vinf
                          (v3*h1a(k)-v1*h3(k) )/mu+v2/vinf
                          (v1*h2a(k)-v2*h1a(k))/mu+v3/vinf];
       
        e_vec_b_I(:,k) = [(v2*h3(k) -v3*h2b(k))/mu+v1/vinf
                          (v3*h1b(k)-v1*h3(k) )/mu+v2/vinf
                          (v1*h2b(k)-v2*h1b(k))/mu+v3/vinf];
        
        % p1-versor (Perifocal RF)
        p1a_I(:,k) = e_vec_a_I(:,k)/e;
        p1b_I(:,k) = e_vec_b_I(:,k)/e;
        
        % p2-versor (Perifocal RF)
        p2a_I(:,k) = cross(p3a_I(:,k),p1a_I(:,k));
        p2b_I(:,k) = cross(p3b_I(:,k),p1b_I(:,k));
        
        % Argument of Periapsis
        if e_vec_a_I(3,k)>0
            ARGPa(k) =     acosd(N_vec_a_I(:,k)'*p1a_I(:,k));
        else
            ARGPa(k) = 360-acosd(N_vec_a_I(:,k)'*p1a_I(:,k));
        end

        if e_vec_b_I(3,k)>0
            ARGPb(k) =     acosd(N_vec_b_I(:,k)'*p1b_I(:,k));
        else
            ARGPb(k) = 360-acosd(N_vec_b_I(:,k)'*p1b_I(:,k));
        end

    end
end
%% ARGP
if strcmp(User_Input,'ARGP')
    % ARGP = [w(1):5:w(2), w(3):5:w(4)];
    ARGP = var;
    K = size(ARGP,2);
    for k=1:K
        if (ARGP(k)>w(2) && ARGP(k)<w(3)) || (ARGP(k)>w(4) && ARGP(k)<w(1))
            fprintf('For ARGP = %.1f deg, no possible trajectory exists. Choose a valid argument of periapsis.\n',ARGP(k))
            return
        else
            fprintf('For ARGP = %.1f deg, two possible trajectories exist.\n',ARGP(k))
        end
        % Argument of Periapsis
        ARGPa(k) = ARGP(k);
        ARGPb(k) = ARGP(k);
    
        % RAAN
        if v1~=0
            A = (v2/v1)^2+1;
            B = 2*(v2/v1)*mu/(h*v1)*((e + cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)));
            C = (mu/(v1*h))^2*((e + cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)))^2 - 1;
            discr = B^2-4*A*C;
            if abs(discr)<1e-12
                discr=0;
            end
            sinRAANa = -B/(2*A)+1/(2*A)*sqrt(discr);
            cosRAANa = -v2/v1*sinRAANa-mu/(v1*h)*((e + cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)));
            % cosRAANa = -v2/v1*sinRAANa-vinf/v1/e*(sqrt(e^2-1)*sind(ARGP(k))-cosd(ARGP(k)));
            sinRAANb = -B/(2*A)-1/(2*A)*sqrt(discr);
            cosRAANb = -v2/v1*sinRAANb-mu/(v1*h)*((e + cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)));
            % cosRAANb = -v2/v1*sinRAANb-vinf/v1/e*(sqrt(e^2-1)*sind(ARGP(k))-cosd(ARGP(k)));
        end
        if v2~=0
            A = (v1/v2)^2+1;
            B = 2*(v1/v2)*mu/(h*v2)*((e + cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)));
            C = (mu/(v2*h))^2*((e + cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)))^2 - 1;
            discr = B^2-4*A*C;
            if abs(discr)<1e-12
                discr=0;
            end
            cosRAANa = -B/(2*A)-1/(2*A)*sqrt(discr);
            sinRAANa = -v1/v2*cosRAANa-mu/(v2*h)*((e + cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)));
            % sinRAANa = -v1/v2*cosRAANa-vinf/v2/e*(sqrt(e^2-1)*sind(ARGP(k))-cosd(ARGP(k)));
            cosRAANb = -B/(2*A)+1/(2*A)*sqrt(discr);
            sinRAANb = -v1/v2*cosRAANb-mu/(v2*h)*((e + cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)));
            % sinRAANb = -v1/v2*cosRAANb-vinf/v2/e*(sqrt(e^2-1)*sind(ARGP(k))-cosd(ARGP(k)));
        end
        RAANa(k) = atan2d(sinRAANa,cosRAANa);
        if RAANa(k)<0
            RAANa(k)=RAANa(k)+360;
        end
        RAANb(k) = atan2d(sinRAANb,cosRAANb);
        if RAANb(k)<0
            RAANb(k)=RAANb(k)+360;
        end
    
        % Orbital Inclination
        sinia(k) =                                     v3/(mu/h*(cosd(ARGP(k))*(e+cosd(TA_inf_Entry)) - sind(ARGP(k))*sind(TA_inf_Entry)));
        cosia(k) = (-v1*sind(RAANa(k))+v2*cosd(RAANa(k)))/(mu/h*(cosd(ARGP(k))*(e+cosd(TA_inf_Entry)) - sind(ARGP(k))*sind(TA_inf_Entry)));
        ia(k)    = atan2d(sinia(k),cosia(k));

        sinib(k) =                                     v3/(mu/h*(cosd(ARGP(k))*(e+cosd(TA_inf_Entry)) - sind(ARGP(k))*sind(TA_inf_Entry)));
        cosib(k) = (-v1*sind(RAANb(k))+v2*cosd(RAANb(k)))/(mu/h*(cosd(ARGP(k))*(e+cosd(TA_inf_Entry)) - sind(ARGP(k))*sind(TA_inf_Entry)));
        ib(k)    = atan2d(sinib(k),cosib(k));


        % % Orbital Inclination
        % if (v2*cosd(RAANa(k)) - v1*sind(RAANa(k)))/v3>=0
        %     ia(k) =       asind( abs(v3) / sqrt( ( v2*cosd(RAANa(k)) - v1*sind(RAANa(k)) )^2 + v3^2 ) );
        % else
        %     ia(k) = 180 - asind( abs(v3) / sqrt( ( v2*cosd(RAANa(k)) - v1*sind(RAANa(k)) )^2 + v3^2 ) );
        % end
        % if (v2*cosd(RAANb(k)) - v1*sind(RAANb(k)))/v3>=0
        %     ib(k) =       asind( abs(v3) / sqrt( ( v2*cosd(RAANb(k)) - v1*sind(RAANb(k)) )^2 + v3^2 ) );
        % else
        %     ib(k) = 180 - asind( abs(v3) / sqrt( ( v2*cosd(RAANb(k)) - v1*sind(RAANb(k)) )^2 + v3^2 ) );
        % end

        v1*cosd(RAANa(k))+v2*sind(RAANa(k)) + mu/h*( (e+cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)) );
        v1*cosd(RAANb(k))+v2*sind(RAANb(k)) + mu/h*( (e+cosd(TA_inf_Entry))*sind(ARGP(k)) + sind(TA_inf_Entry)*cosd(ARGP(k)) );

        (v1*sind(RAANa(k))-v2*cosd(RAANa(k)))*sind(ia(k)) + v3*cosd(ia(k));
        (v1*sind(RAANb(k))-v2*cosd(RAANb(k)))*sind(ib(k)) + v3*cosd(ib(k));
    end
end
% return
for k = 1;
% for k=1:K
    % Coordinate Transformation Matrix between PCI and P-RF
    R_PaI(:,:,k) = R3(ARGPa(k))*R1(ia(k))*R3(RAANa(k));
    R_PbI(:,:,k) = R3(ARGPb(k))*R1(ib(k))*R3(RAANb(k));
    
    %% SOI Entrance
    TA_SOI_Entry = -acosd(1/e*(h^2/mu/r_SOI -1));
    r_vec_SOI_Entry_P = r_SOI*[cosd(TA_SOI_Entry); sind(TA_SOI_Entry); 0]; 
    r_vec_SOI_Entry_a_I(:,k) = R_PaI(:,:,k)'*r_vec_SOI_Entry_P;
    r_vec_SOI_Entry_b_I(:,k) = R_PbI(:,:,k)'*r_vec_SOI_Entry_P;

    DEC_Entry_a(k) =  90 - acosd(dot(r_vec_SOI_Entry_a_I(:,k)/r_SOI,[0;0;1]));
    RA_Entry_a(k)   = atan2d(r_vec_SOI_Entry_a_I(2,k),r_vec_SOI_Entry_a_I(1,k));
    if RA_Entry_a(k)<0
        RA_Entry_a(k)=RA_Entry_a(k)+360;
    end
    DEC_Entry_b(k) =  90 - acosd(dot(r_vec_SOI_Entry_b_I(:,k)/r_SOI,[0;0;1]));
    RA_Entry_b(k)   = atan2d(r_vec_SOI_Entry_b_I(2,k),r_vec_SOI_Entry_b_I(1,k));
    if RA_Entry_b(k)<0
        RA_Entry_b(k)=RA_Entry_b(k)+360;
    end
    % Velocity Vector at SOI Exit
    v_vec_SOI_Entry_P = mu/h*[-sind(TA_SOI_Entry);e+cosd(TA_SOI_Entry);0];
    v_vec_SOI_Entry_a_I(:,k) = R_PaI(:,:,k)'*v_vec_SOI_Entry_P;
    v_vec_SOI_Entry_b_I(:,k) = R_PbI(:,:,k)'*v_vec_SOI_Entry_P;

    % V Inf at Entrance
    v_inf_Entry_P = mu/h*[-sind(TA_inf_Entry);e+cosd(TA_inf_Entry);0];
    v_inf_Entry_a_I(:,k) = R_PaI(:,:,k)'*v_inf_Entry_P;
    v_inf_Entry_b_I(:,k) = R_PbI(:,:,k)'*v_inf_Entry_P;

    %% SOI Exit
    TA_SOI_Exit  =  acosd(1/e*(h^2/mu/r_SOI -1));
    r_vec_SOI_Exit_P = r_SOI*[cosd(TA_SOI_Exit); sind(TA_SOI_Exit); 0];
    r_vec_SOI_Exit_a_I(:,k) = R_PaI(:,:,k)'*r_vec_SOI_Exit_P;
    r_vec_SOI_Exit_b_I(:,k) = R_PbI(:,:,k)'*r_vec_SOI_Exit_P;

    DEC_Exit_a(k) =  90 - acosd(dot(r_vec_SOI_Exit_a_I(:,k)/r_SOI,[0;0;1]));
    RA_Exit_a(k)   = atan2d(r_vec_SOI_Exit_a_I(2,k),r_vec_SOI_Exit_a_I(1,k));
    if RA_Exit_a(k)<0
        RA_Exit_a(k)=RA_Exit_a(k)+360;
    end
    DEC_Exit_b(k) =  90 - acosd(dot(r_vec_SOI_Exit_b_I(:,k)/r_SOI,[0;0;1]));
    RA_Exit_b(k)   = atan2d(r_vec_SOI_Exit_b_I(2,k),r_vec_SOI_Exit_b_I(1,k));
    if RA_Exit_b(k)<0
        RA_Exit_b(k)=RA_Exit_b(k)+360;
    end
    % Velocity Vector at SOI Exit
    v_vec_SOI_Exit_P = mu/h*[-sind(TA_SOI_Exit);e+cosd(TA_SOI_Exit);0];
    v_vec_SOI_Exit_a_I(:,k) = R_PaI(:,:,k)'*v_vec_SOI_Exit_P;
    v_vec_SOI_Exit_b_I(:,k) = R_PbI(:,:,k)'*v_vec_SOI_Exit_P;

    % V Inf at Exit
    v_inf_Exit_P = mu/h*[-sind(TA_inf_Exit);e+cosd(TA_inf_Exit);0];
    v_inf_Exit_a_I(:,k) = R_PaI(:,:,k)'*v_inf_Exit_P;
    v_inf_Exit_b_I(:,k) = R_PbI(:,:,k)'*v_inf_Exit_P;
    coe = [e, a, ia, ARGPa, RAANa,TA_SOI_Entry, h];
    fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    fprintf(' Altitude of Periapsis:    zp = %.1f km\n'   ,zp)
    fprintf('          Excess Speed: v_inf = %.3f km/s\n',vinf)
    fprintf(' Characteristic Energy:    C3 = %.1f km^2/s^2\n',vinf^2)
    fprintf('        Orbital Energy:     E = %.1f km^2/s^2\n',En)
    fprintf('       Semi-Major Axis:     a = %.1f km\n'    ,a)
    fprintf('          Eccentricity:     e = %.3f\n'       ,e)
    fprintf(' Angular Momentum Norm:     h = %.1f km^2/s\n',h)
    fprintf('            Turn Angle:     d = %.1f deg\n' ,Turn_Angle)
    fprintf('         Aiming Radius:     D = %.1f km\n\n',Aiming_Radius)
    fprintf('Right Asc.of Asc. Node:  RAAN = %.1f  or  %.1f deg\n'  ,RAANa(k),RAANb(k))
    fprintf('   Orbital Inclination:     i = %.1f  or  %.1f deg\n'  ,ia(k), ib(k))
    fprintf(' Argument of Periapsis:  ARGP = %.1f  or  %.1f deg\n\n',ARGPa(k),ARGPb(k))
    fprintf('        @ SOI Entrance:    RA = %.2f  or  %.2f deg\n'  ,RA_Entry_a(k), RA_Entry_b(k))
    fprintf('                          DEC = %.2f  or  %.2f deg\n'  ,DEC_Entry_a(k),DEC_Entry_b(k))
    fprintf('        @ SOI     Exit:    RA = %.2f  or  %.2f deg\n'  ,RA_Exit_a(k),  RA_Exit_b(k))
    fprintf('                          DEC = %.2f  or  %.2f deg\n'  ,DEC_Exit_a(k), DEC_Exit_b(k))
    fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')

    TA = TA_SOI_Entry:(TA_SOI_Exit-TA_SOI_Entry)/1000:TA_SOI_Exit;

    for j=1:size(TA,2)
        r(j) = h^2/mu/(1+e*cosd(TA(j)));
        r_P(:,j) = r(j)*[cosd(TA(j)); sind(TA(j)); 0];
        r_a_I(:,j,k) = R_PaI(:,:,k)'*r_P(:,j);
        r_b_I(:,j,k) = R_PbI(:,:,k)'*r_P(:,j);
    end
    r_PER_P  = [h^2/mu/(1+e);0;0];
    r_PER_a_I(:,k) = R_PaI(:,:,k)'*r_PER_P;
    r_PER_b_I(:,k) = R_PbI(:,:,k)'*r_PER_P;
end

Color_Entry = [0,0.7,0];
Color_Exit  = [1,0.4,0];

% AU = 1E6;
% figure(1),hold on
% set(gcf,'units','normalized','position',[0.0,0.0,1.0,0.9])
% 
% subplot(1,2,1),hold on
% tlt1 = title('Sphere of Influence (SOI) & Celestial Equatorial Plane');
% set(tlt1,'Color','w')
% set(gca,'Color','k')
% set(gcf,'Color','k')    
% % set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
% ax = gca;
% ax.XColor = 'w'; ax.YColor = 'w'; ax.ZColor = 'w';
% 
% axis('equal')
% view([-240,16])
% [xSOI,ySOI,zSOI] = ellipsoid(0,0,0,r_SOI/AU,r_SOI/AU,r_SOI/AU,60);
% SOI_Plot = surface(xSOI,ySOI,zSOI);
% SOI_Plot.FaceColor = 'b';
% SOI_Plot.FaceAlpha = 0.15;
% SOI_Plot.EdgeColor = 'none';
% xlim([-r_SOI/AU,r_SOI/AU])
% ylim([-r_SOI/AU,r_SOI/AU])
% zlim([-r_SOI/AU,r_SOI/AU])
% 
% xlabel('X [Mm]')
% ylabel('Y [Mm]')
% zlabel('Z [Mm]')

        % PCI Axes
        % line(r_SOI/AU*[0,1],[0,0],[0,0],'LineWidth',3,'color',[1,0,0.3]);
        % line([0,0],r_SOI/AU*[0,1],[0,0],'LineWidth',3,'color',[1,0,0.3]);
        % line([0,0],[0,0],r_SOI/AU*[0,1],'LineWidth',3,'color',[1,0,0.3]);
        % 
        % text(r_SOI/AU,0,0,'$\mathbf{\hat{I}_1}$','Color',[1,0,0.3],'FontSize',18,'interpreter','latex');
        % text(0,r_SOI/AU,0,'$\mathbf{\hat{I}_2}$','Color',[1,0,0.3],'FontSize',18,'interpreter','latex');
        % text(0,0,r_SOI/AU,'$\mathbf{\hat{I}_3}$','Color',[1,0,0.3],'FontSize',18,'interpreter','latex');
% 
% phi = 0:6:360;
% for j=1:size(phi,2)
%     A_I(:,j) = r_SOI/AU*[cosd(phi(j)),sind(phi(j)),0];
% end
% Eq_Plane = patch(A_I(1,:),A_I(2,:),A_I(3,:),[1,0,1]); 
% Eq_Plane.FaceColor = [0,0.5,1];
% Eq_Plane.EdgeColor = 'none';
% Eq_Plane.FaceAlpha = 0.4; 
% 
% 
% for k=1:K
%     plot3(r_a_I(1,:,k)/AU,r_a_I(2,:,k)/AU,r_a_I(3,:,k)/AU,'color',ColorTraj)
%     plot3(r_b_I(1,:,k)/AU,r_b_I(2,:,k)/AU,r_b_I(3,:,k)/AU,'color',ColorTraj)
% end
% 
% plot3(r_vec_SOI_Entry_a_I(1,:)/AU,r_vec_SOI_Entry_a_I(2,:)/AU,r_vec_SOI_Entry_a_I(3,:)/AU,'marker','s','markeredgecolor',Color_Entry,'markerfacecolor',Color_Entry,'markersize',6,'linestyle','none');
% plot3(r_vec_SOI_Entry_b_I(1,:)/AU,r_vec_SOI_Entry_b_I(2,:)/AU,r_vec_SOI_Entry_b_I(3,:)/AU,'marker','s','markeredgecolor',Color_Entry,'markerfacecolor',Color_Entry,'markersize',6,'linestyle','none');
% 
% plot3(r_vec_SOI_Exit_a_I(1,:)/AU,r_vec_SOI_Exit_a_I(2,:)/AU,r_vec_SOI_Exit_a_I(3,:)/AU,'marker','s','markeredgecolor',Color_Exit,'markerfacecolor',Color_Exit,'markersize',6,'linestyle','none');
% plot3(r_vec_SOI_Exit_b_I(1,:)/AU,r_vec_SOI_Exit_b_I(2,:)/AU,r_vec_SOI_Exit_b_I(3,:)/AU,'marker','s','markeredgecolor',Color_Exit,'markerfacecolor',Color_Exit,'markersize',6,'linestyle','none');
% 
% subplot(1,2,2),hold on
% tlt2 = title('Zoom-In on Planet with Locus of Periapsis');
% set(tlt2,'Color','w')
% set(gca,'Color','k')
% set(gcf,'Color','k')
% ax = gca;
% ax.XColor = 'w'; ax.YColor = 'w'; ax.ZColor = 'w';
% 
% axis('equal')
% [xPRI,yPRI,zPRI] = ellipsoid(0,0,0,RE/RE,RE/RE,RP/RE,60);
% PRI_Plot = surface(xPRI,yPRI,-zPRI);
% PRI_Plot.FaceColor = 'texturemap';
% PRI_Plot.FaceAlpha = 1.0;
% PRI_Plot.EdgeColor = 'none';
% PRI_Plot.CData     = Primary.Figure;
% xlim([-2*rp/RE,2*rp/RE])
% ylim([-2*rp/RE,2*rp/RE])
% zlim([-2*rp/RE,2*rp/RE])
% view([0,0])
% 
% xlabel('X - Planet''s Radii')
% ylabel('Y - Planet''s Radii')
% zlabel('Z - Planet''s Radii')
% 
% phi = 0:6:360;
% for j=1:size(phi,2)
%     A_I(:,j) = 2*rp/RE*[cosd(phi(j)),sind(phi(j)),0];
% end
% Eq_Plane = patch(A_I(1,:),A_I(2,:),A_I(3,:),[1,0,1]); 
% Eq_Plane.FaceColor = [0,0.5,1];
% Eq_Plane.EdgeColor = 'none';
% Eq_Plane.FaceAlpha = 0.4;
% 
% for k=1:K
%     plot3(r_a_I(1,:,k)/RE,r_a_I(2,:,k)/RE,r_a_I(3,:,k)/RE,'color',ColorTraj)
%     plot3(r_b_I(1,:,k)/RE,r_b_I(2,:,k)/RE,r_b_I(3,:,k)/RE,'color',ColorTraj)
% end
% plot3(r_PER_a_I(1,:)/RE,r_PER_a_I(2,:)/RE,r_PER_a_I(3,:)/RE,'marker','s','markeredgecolor',ColorPeri,'markerfacecolor',ColorPeri,'markersize',5,'linestyle','none');
% plot3(r_PER_b_I(1,:)/RE,r_PER_b_I(2,:)/RE,r_PER_b_I(3,:)/RE,'marker','s','markeredgecolor',ColorPeri,'markerfacecolor',ColorPeri,'markersize',5,'linestyle','none');

end


%% Extra Functions

function R1 = R1(a)
    R1 = [1,      0,       0;
          0, cosd(a),sind(a);
          0,-sind(a),cosd(a)];
end
function R2 = R2(a)
    R2 = [cosd(a),0,-sind(a);
                0,1,       0;
          sind(a),0, cosd(a)];
end
function R3 = R3(a)
    R3 = [ cosd(a),sind(a),0;
          -sind(a),cosd(a),0;
                 0,      0,1];
end