function [s, j] = nr(f_, df_, x1, x2, xacc, nmax ,u, d, lambda, one_min_z, s, mu)
    x1 = x1 + 1000*eps;
    x2 = x2 - 1000*eps;
    err = 0;
    j = 0;
    fl = f_(x1,u, d, lambda, one_min_z, s, mu); fh = f_(x2,u, d,lambda, one_min_z, s, mu);


    if fl == inf
            fl = -inf;
    end

    if 0 && ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))

        fprintf('Root must be bracketed in rtsafe');
        return;
    end
    if (fl == 0.0)
        s = x1;
        return 
    end
    if (fh == 0.0) 
        s = x2;
        return 
    end
    if (fl < 0.0) 
        xl=x1;
        xh=x2;
    else
        xh=x1;
        xl=x2;
    end

    rts=0.5*(x1+x2); 
    %rts = 0.99 * x1 + 0.01 * x2;
    dxold=abs(x2-x1);
    dx=dxold;

    f = f_(rts,u, d,lambda, one_min_z, s, mu); df = df_(rts,u, d, lambda, one_min_z, s, mu);

    for j=1:nmax

        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (abs(2.0*f) > abs(dxold*df))) 
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            if (xl == rts) 
                s = rts;
                return;
            end
        else
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts = rts - dx;
            if (temp == rts)
                s = rts;
                return;
            end
        end

%         
%         if (abs(dx) < xacc)
%             s = rts;
%             %fprintf("# itr/ = %d\n", j);
%             return;
%         end

        f = f_(rts,u, d,lambda, one_min_z, s, mu); df = df_(rts,u, d, lambda, one_min_z, s, mu);
        
        if (abs(f) < xacc)
            s = rts;
            return
        end
        
        if (f < 0.0)
            xl=rts;
        else
            xh=rts;
        end
    end
    
    %fprintf('\n [!] Algorithm did not converge.\n');
    err = 0;
    s = x1 - 1000*eps;
        
end