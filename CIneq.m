function varargout = CIneq(x,mode)
% x1*x2>5 and x1^2+x2^2<20
global numfCIneq numgCIneq numHCIneq ProblemNum

argout = 0;
if bitand(mode,1)
    numfCIneq = numfCIneq + 1;
    argout = argout + 1;
    switch ProblemNum
        case(0)
            varargout(argout) = {...
                [x(1)*x(2)-2;
                20-(x(1)^2+x(2)^2);
                x(1)+1;]};
        case(1)
            varargout(argout) = {...
                [x(1)*x(2)*x(3)*x(4)-25;
                x(1)-1;
                x(2)-1;
                x(3)-1;
                x(4)-1;
                5-x(1);
                5-x(2);
                5-x(3);
                5-x(4);]};
        case(2)
            varargout(argout) = {x(1)-1000};
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
    numgCIneq = numgCIneq + 1;
    argout = argout + 1;
    
    switch ProblemNum
        case(0)
            varargout(argout) = {...
                [x(2)    x(1);
                -2*x(1) -2*x(2);
                1.0        0]};
        case(1)
            varargout(argout) = {...
                [x(2)*x(3)*x(4) x(1)*x(3)*x(4) x(1)*x(2)*x(4) x(1)*x(2)*x(3);
                1 0 0 0;
                0 1 0 0;
                0 0 1 0;
                0 0 0 1;
                -1 0 0 0;
                0 -1 0 0;
                0 0 -1 0;
                0 0 0 -1;]};
        case(2)
            g=zeros(1,length(x));
            g(1)=1;
            varargout(argout) = {g};

    end
    
    
    
end
if bitand(mode,4)
    numHCIneq = numHCIneq + 1;
    argout = argout + 1;
    % H_Ineq is % n*n*nIneq, so we calculate the hessian for each
    % inequality and extend it.

    
    switch ProblemNum
        case(0)
            H_Ineq = [0 1; 1 0;]; %for x1*x2>5 and 20-(x1^2+x2^2)>0
            H_Ineq(:,:,2) = [-2 0; 0 -2];
            H_Ineq(:,:,3) = [0 0; 0 0];
        case(1)
            H_Ineq=[0           x(3)*x(4)   x(2)*x(4)   x(2)*x(3);
                    x(3)*x(4)   0           x(1)*x(4)   x(1)*x(3);
                    x(2)*x(4)   x(1)*x(4)   0           x(1)*x(2);
                    x(2)*x(3)   x(1)*x(3)   x(1)*x(2)   0;];
            H_Ineq(:,:,2)=zeros(4,4);
            H_Ineq(:,:,3)=zeros(4,4);
            H_Ineq(:,:,4)=zeros(4,4);
            H_Ineq(:,:,5)=zeros(4,4);
            H_Ineq(:,:,6)=zeros(4,4);
            H_Ineq(:,:,7)=zeros(4,4);
            H_Ineq(:,:,8)=zeros(4,4);
            H_Ineq(:,:,9)=zeros(4,4);
        case(2)
            H_Ineq=zeros(length(x));

            
    end
    
    
    
    varargout(argout) = {H_Ineq};
end
return;
