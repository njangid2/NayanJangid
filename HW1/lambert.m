function [V1, V2] =lambert(R1, R2, t, string,mu)
    r1 = norm(R1);
    r2 = norm(R2);
    c12 = cross(R1, R2);
    theta = acos(dot(R1,R2)/r1/r2);

    if strcmp(string, 'pro')
        if c12(3) <= 0
            theta = 2*pi - theta;
        end
    elseif strcmp(string,'retro')
        if c12(3) >= 0
            theta = 2*pi - theta;
        end
    else
        string = 'pro'
        fprintf('\n ** Prograde trajectory assumed.\n')
    end
    A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));
    z = -100;
    if z > 0
        S=(sqrt(z) - sin(sqrt(z)))/((sqrt(z))^3);
    elseif z < 0
        S=(sinh(sqrt(-z)) - sqrt(-z))/((sqrt(-z))^3);
    else
        S=1/6;
    end
    if z > 0
        C=(1 - cos(sqrt(z)))/z;
    elseif z < 0
        C=(cosh(sqrt(-z)) - 1)/(-z);
    else
        C=1/2;
    end
    y=(r1 + r2 + A*(z*S - 1)/sqrt(C));
    F=(y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*t;
    while F < 0
        z=z+ 0.1;
            if z > 0
        S=(sqrt(z) - sin(sqrt(z)))/((sqrt(z))^3);
        elseif z < 0
            S=(sinh(sqrt(-z)) - sqrt(-z))/((sqrt(-z))^3);
        else
            S=1/6;
        end
        if z > 0
            C=(1 - cos(sqrt(z)))/z;
        elseif z < 0
            C=(cosh(sqrt(-z)) - 1)/(-z);
        else
            C=1/2;
        end
        y=(r1 + r2 + A*(z*S - 1)/sqrt(C));
        F=(y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*t;
    end
    tol = 1*10^(-8);
    nmax = 1000;
    ratio = 1;
    n = 0;
    while (abs(ratio) > tol) && (n <= nmax)
        n = n + 1;
        if z > 0
            S=(sqrt(z) - sin(sqrt(z)))/((sqrt(z))^3);
        elseif z < 0
            S=(sinh(sqrt(-z)) - sqrt(-z))/((sqrt(-z))^3);
        else
            S=1/6;
        end
        if z > 0
            C=(1 - cos(sqrt(z)))/z;
        elseif z < 0
            C=(cosh(sqrt(-z)) - 1)/(-z);
        else
            C=1/2;
        end
        y = r1 + r2 + A*(z*S - 1)/sqrt(C);
        F=(y/C)^1.5*S + A*sqrt(y) - sqrt(mu)*t;
        if z == 0
            dfz = sqrt(2)/40*yo^1.5 + A/8*(sqrt(y)+ A*sqrt(1/2/y));
        else
            dfz = (y/C)^1.5*(1/2/z*(C - 3*S/2/C)+ 3*S^2/4/C)+ A/8*(3*S/C*sqrt(y)+ A*sqrt(C/y));
        end
        ratio = F/dfz;
        z = z - ratio;
    end
    if n >= nmax
        fprintf('\n\n **Number of iterations exceeds')
        fprintf(' %g \n\n ', nmax)
    end
    if z > 0
        S=(sqrt(z) - sin(sqrt(z)))/((sqrt(z))^3);
    elseif z < 0
        S=(sinh(sqrt(-z)) - sqrt(-z))/((sqrt(-z))^3);
    else
        S=1/6;
    end
    if z > 0
        C=(1 - cos(sqrt(z)))/z;
    elseif z < 0
        C=(cosh(sqrt(-z)) - 1)/(-z);
    else
        C=1/2;
    end
    f = 1 - (r1 + r2 + A*(z*S - 1)/sqrt(C))/r1;
    g = A*sqrt((r1 + r2 + A*(z*S - 1)/sqrt(C))/mu);
    gdot = 1 - (r1 + r2 + A*(z*S - 1)/sqrt(C))/r2;
    V1 = (1/g)*(R2 - f*R1);
    V2 = (1/g)*(gdot*R2 - R1);
return