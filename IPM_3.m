function [inform,x] = IPM_3(fun,Ceq,CIneq,MeritF,x,trparams)

% This function solves an optimization problem for :
% min f(x) s.t. Ceq=0; CIneq>=0;
% This is posed as a linear system of the form (Page 570 Nocedel and Wright)
% [ Hess_xx(L)+A_Ineq^T*(S^-1*Z)*A_Ineq  -A_Eq^T(x)] [Dx]  [Grad(f)-A_Eq^T(x)*y-A_Ineq^T(x)*z;]
% [-A_Eq(x)                                0       ] |Dy|=-[-Ceq(x);                          ]
% with         Z*ds+S*dz=-(S*z-mu*e);
% AND          dz=-(S^-1*Z)*(CIneq-mu*Z^-1*e+A_Ineq*dx)
% Sigma=S^-1*Z
% x is the input variable (n*1)
% where Hess_xx(L) is the Hessian of the Lagrangian function f-CEq^T(x)*y-CIneq^T(x)*z (n*n);
% Grad(f) is the gradient of the object function (n*1).
% A_Eq(x) and A_Ineq(x) are the jacobian of the Eq and Ineq
% constraints w.r.t. x and are respectively n_eq*n and n_Ineq*n; y and z are the lagrange multiplier associated with
% the Ceq and CIneq,so they have the same size as the number of Eq and Ineq constraints (n_eq*1) and (n_Ineq*1).
% Y, Z are basically y and z that are placed on the diagonal of a matrix : (n_eq*n_eq) and (n_Ineq*n_Ineq)
% s is the slack variable to convert the inequality constraints into the
% standard format (n_Ineq*1). S is the diagonal matrix build from s (n_Ineq*n_Ineq).

global numf numg numH;
numf=0; numg=0; numH=0;
status=0;
n=length(x.p)
n_eq=length(feval(Ceq,x.p,1))
n_Ineq=length(feval(CIneq,x.p,1))

mu=100;
s=feval(CIneq,x.p,1)+0.001;
z=mu./s;
y=zeros(n_eq,1);
e=ones(n_Ineq,1);

fprintf('n=%d, n_eq=%d, n_Ineq=%d\n',n,n_eq,n_Ineq);
fprintf('\n\t\t#Newton     F1       F2      alpha_s   alpha_z     mu       merit\n');
fprintf('\t\t-------------------------------------------------------------------------------------\n');
mu=s'*z/n_Ineq;
for Newton=1:trparams.maxtotalIter
    for k=1:trparams.maxNewtownIter
        Hf = feval(fun,x.p,4);
        HCeq=feval(Ceq,x.p,4);
        HCIneq=feval(CIneq,x.p,4);
        HCIneqT_z=zeros(n,n);
        HCeqT_y=zeros(n,n);
        
        for i = 1:n_Ineq
            HCIneqT_z =HCIneqT_z+ z(i) * HCIneq(:,:,i);
        end
        for i = 1:n_eq
            HCeqT_y =HCeqT_y+ y(i) * HCeq(:,:,i);
        end
        H_Lag=Hf-HCeqT_y-HCIneqT_z;
        A_Eq=feval(Ceq,x.p,2);
        A_Ineq=feval(CIneq,x.p,2);
        S=diag(s);
        Z=diag(z);
        
        A_Matrix=...
            [H_Lag+A_Ineq'*(inv(S)*Z)*A_Ineq  -A_Eq'            ;
            -A_Eq                            zeros(n_eq,n_eq)  ;
            ];
        F1=feval(fun,x.p,2)-A_Eq'*y-A_Ineq'*z+...
           A_Ineq'*(inv(S)*Z)*(feval(CIneq,x.p,1)-mu*inv(Z)*e);
        F2=-(feval(Ceq,x.p,1));
        
        b=-[F1;
            F2];
        delta=inv(A_Matrix)*b;
        start_idx=1; end_idx=n;
        dx=delta(start_idx:end_idx);
        start_idx=n+1; end_idx=n+n_eq;
        dy=delta(start_idx:end_idx);
        dz=diag(z./s)*(-(feval(CIneq,x.p,1)-mu*inv(Z)*e)-A_Ineq*dx);
        ds=inv(Z)*(-(S*z-mu*e)-S*dz);

        [alpha_s, alpha_z] = FindAlphas(x.p,s,y,z,dx,ds,dy,dz,fun,Ceq,CIneq,MeritF,trparams.taw);
        x.p=x.p+alpha_s*dx;
        y=y+alpha_z*dy;
        z=z+alpha_z*dz;
        s=s+alpha_s*ds;

        merit=feval(MeritF,x.p,fun,Ceq,CIneq,y,z,1);
        Error=norm(b,2);
        if(Error<trparams.epsilon)
            break;
        end
    end
    mu=mu*0.2;
    fprintf('Iter %d\t\t %d\t %2.2e, %2.2e, %2.2e, %2.2e, %2.2e, %2.2e\n',Newton,k,norm(F1,inf), norm(F2,inf),alpha_s,alpha_z,mu,merit);
end
if(Error<trparams.Final_toler)
    status=1;
end

x.f = feval(fun,x.p,1);
x.g = feval(fun,x.p,2);
x.h = feval(fun,x.p,4);
x = struct('p',x.p,'f',x.f,'g',x.g,'h',x.h);
inform = struct('status',status,'Newton',Newton);
end

