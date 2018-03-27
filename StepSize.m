function [alfa, xnew] = StepSize(fun, x, d, y,z, alfa, params)

% <input>
% fun    : pointer to a function
% x      : x.p is the starting point
% d      : search direction
% alfa   : initial step length
% params : contains c1, c2, maxit
%
% <output>
% alfa   : resulting alfa
% xnew   : xnew.p = x(k+1), xnew.f = f(x(k+1)), xnew.g = f'(x(k+1))

% set initial value of alpha_{i-1}, and corresponding phi value
aprev=0.0;
x.g=feval(fun,x.p,1,y,z);
phiprev=x.f; phidprev=x.g'*d;

%store phi and phi' evaluated at alpha=0
phi0 = x.f; phid0=x.g'*d;

% if d is not a descent direction, return with alfa=0 (*** check for this in
% the calling code ***) 
if phid0 >= 0.0 
  alfa = 0.0; 
  xnew = x;
  return;
end

amax=realmax;

%x.f = feval(fun,x.p,1);
%x.g = feval(fun,x.p,2);
xL = [];                % to reuse f,g values in INT
xR = [];                % to reuse f,g values in INT

%fprintf('initial alfa=%8.3e\n', alfa);

for i=1:params.maxit
    % Calculate phi(alfa)
    xnew.p = x.p + alfa*d; 
    xnew.f = feval(fun,xnew.p,1,y,z);
    phi = xnew.f;

    %  Test for a final zoom
    if( phi > phi0 + params.c1*phid0*alfa ) | (i>1 & phi >= phiprev)
      alo = aprev; ahi = alfa;  
      philo=phiprev; phihi=phi; phidlo=phidprev;
      break;   % jump out of loop, to ZOOM section
    end

    % evaluate phi' at the latest alfa
    xnew.g=feval(fun,xnew.p,2,y,z);
    phid=xnew.g'*d;

    % test second strong Wolfe condition
    if ( abs(phid) <= -params.c2*phid0 )
      %fprintf('**accepted\n');
      return;
    end

   % if slope is positive, do a final zoom and exit
   if (phid >= 0.0)
     alo=alfa; ahi=aprev;
     philo=phi; phihi=phiprev; phidlo=phid;
     break;   % jump out of loop, to ZOOM section
   end

   % otherwise double the trial step
   aprev = alfa; phiprev = phi; phidprev = phid;
   alfa = 3*alfa;
   %fprintf('increasing alfa=%8.3e\n', alfa);

end

if i==params.maxit
  return;
end


%% ZOOM section

% set safeguarding parameter
mu = 0.1;

% keep iterating, starting at the last iteration of the earlier loop
for ii=i:params.maxit

% do the quadratic interpolation
        atemp = (phihi - philo - phidlo*(ahi-alo)) / (ahi-alo)^2;
        aint = -phidlo/(2*atemp) + alo;
        
        % safeguarding
        if (alo < ahi) & (aint < alo + mu*(ahi-alo))
          aint = alo + mu*(ahi-alo);
        elseif (alo < ahi) & (aint > ahi - mu*(ahi-alo))
          aint = ahi - mu*(ahi-alo);
        elseif (alo > ahi) & (aint < ahi + mu*(alo-ahi))
          aint = ahi + mu*(alo-ahi);
        elseif (alo > ahi) & (aint > alo - mu*(alo-ahi))
          aint = alo - mu*(alo-ahi);
        end

        % set alfa and evaluate the function
        alfa = aint;
	% fprintf(' trying alfa=%e\n', alfa);
        xnew.p=x.p+alfa*d;
        xnew.f=feval(fun,xnew.p,1,y,z);  phi = xnew.f;

        if (phi > phi0 + params.c1*alfa*phid0) | (phi >= philo)
          ahi = alfa; phihi=phi;
        else
% evaluate the gradient
          xnew.g=feval(fun,xnew.p,2,y,z); 
          phid = xnew.g'*d;
          if (abs(phid) <= -params.c2*phid0)
  % we are done!
            return;
          end
          if (phid*(ahi-alo) >= 0.0)
             ahi = alo; phihi = philo;
          end
          alo = alfa; philo = phi; phidlo = phid;
        end
end


