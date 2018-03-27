% We want to solve and optimization problem with Interior point method
% f(x)
% such that CIneq>=0 and Ceq=0; 
clear
clc
global numf numg numH ProblemNum
global numfCIneq numgCIneq numHCIneq
global numfCeq numgCeq numHCeq

ProblemNum=2
% Initial guess for the function value


switch ProblemNum
    case(0)
        x = struct('p',[2; 3]);
        x_sol=[3.00000000, 0.66666]';
    case(1)
        x = struct('p',[1; 5; 5; 1]);
        x_sol=[1.00000000, 4.74299963, 3.82114998, 1.37940829]';
    case(2)
        n=10;
        x = struct('p',zeros(n,1));
        x_sol=1:1:n;
end

params = struct('maxNewtownIter', 100,'maxtotalIter',10,'c1',0.01,'c2',0.3,...
                'epsilon',1e-4,'taw',0.995,'Final_toler',1.0e-6);

numf=0; numg=0;
[inform,xnew] = IPM_3(@objF,@Ceq,@CIneq,@MeritF,x,params);
norm(xnew.p' -x_sol)
if norm(xnew.p -x_sol)>1e-3
    fprintf('CONVERGENCE FAILURE: %d steps were taken without reaching to solution\n', inform.Newton);
    fprintf('gradient size decreasing below %10.6g.\n', params.Final_toler);
else
    fprintf('Success: %d steps taken\n', inform.Newton);
end
fprintf('  Ending point: '); fprintf('%10.6g ',xnew.p);
fprintf('\n  Ending function value: %10.6g\n', xnew.f);
fprintf('  No. function evaluations: %d, No. gradient evaluations %d\n',...
    numf, numg);
fprintf('  Norm of ending gradient: %10.6g\n\n\n', norm(xnew.g));
% [X,Y]=meshgrid(-5:0.05:5, -5:0.05:5);
% Z=Y.*(5+X);
% surf(X,Y,Z)
% Z=20-(X.^2+Y.^2);
% hold on
% surf(X,Y,Z,'FaceAlpha',0.2)
% shading flat
% Z=X.*Y-5;
% hold on
% surf(X,Y,Z,'FaceAlpha',0.2)
% shading flat
