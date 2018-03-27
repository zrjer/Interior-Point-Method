function varargout = Ceq(x,mode)

global numfCeq numgCeq numHCeq ProblemNum
argout = 0;
if bitand(mode,1)
    numfCeq  = numfCeq  + 1;
    argout = argout + 1;
    switch ProblemNum
        case(0)
            varargout(argout) = {x(1)-3};
        case(1)
            varargout(argout) = {x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2-40};
        case(2)
            varargout(argout) = {x'*[ones(length(x),1)]-10};
        case(3)
            varargout(argout) = {x(1)^2 - 2*x(1)*x(2) + 4*x(2)^2};
        case(4)
            varargout(argout) = {x(2)*(5+x(1))};
        case(5)
            varargout(argout) = {x(1)^2-2*x(1)*x(2)+4*x(2)^2};
        case(6)
            varargout(argout) = {x(1)^2 + 2 * x(2)^2};
        case(7)
            varargout(argout) = {(x(1)-5)^2 + (x(2)-5)^2};
        case(8)
            varargout(argout) = {(x(1)-5)^2};
        case(9)
            varargout(argout) = {(x(1)-5)^2};
        case(10)
            varargout(argout) = {x(1)^2-2*x(1)*x(2)+4*x(2)^2};
    end
end
if bitand(mode,2)
    numgCeq  = numgCeq  + 1;
    argout = argout + 1;
    
    switch ProblemNum
        case(0)
            varargout(argout) = {[1 0]};
        case(1)
            g=[2*x(1) 2*x(2) 2*x(3) 2*x(4)];
            varargout(argout) = {g};
        case(2)
            n = length(x);
            g=ones(1,n);
            varargout(argout) = {g};
        case(3)
            g=zeros(2,1);
            g(1) = 2*x(1)-2*x(2);
            g(2) = -2*x(1)+8*x(2);
            varargout(argout) = {g};
        case(4)
            g=zeros(2,1);
            g(1) = x(2);
            g(2) = (5+x(1));
            varargout(argout) = {g};
        case(5)
            g=zeros(2,1);
            g(1) =  2*x(1)-2*x(2);
            g(2) = -2*x(1)+8*x(2);
            varargout(argout) = {g};
        case(6)
            g=zeros(2,1);
            g(1) = 2.0 * x(1);
            g(2) = 4.0 * x(2);
            varargout(argout) = {g};
        case(7)
            g=zeros(2,1);
            g(1) = 2 * (x(1)-5);
            g(2) = 2 * (x(2)-5);
            varargout(argout) = {g};
        case(8)
            g=zeros(2,1);
            g(1) = 2 * (x(1) - 5);
            g(2) = 0;
            varargout(argout) = {g};
        case(9)
            g=zeros(1,1);
            g(1) = 2 * (x(1) - 5);
            varargout(argout) = {g};
        case(10)
            g=zeros(2,1);
            g(1) =  2*x(1)-2*x(2);
            g(2) = -2*x(1)+8*x(2);
            varargout(argout) = {g};
    end
end
if bitand(mode,4)
    numHCeq  = numHCeq  + 1;
    argout = argout + 1;
    % H_eq is % n*n*neq, so we calculate the hessian for each
    % inequality and extend it.
    % H_Ineq = [0 1; 1 0;]; %for x1=2
    switch ProblemNum
        case(0)
            H_eq=[0 0; 0 0];
        case(1)
            H_eq=2*eye(4);
        case(2)
            n = length(x);
            H_eq=zeros(n,n);
            
    end
    varargout(argout) = {H_eq};
end

return;
