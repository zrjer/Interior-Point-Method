function varargout = objF(x,mode)

global numf numg numH ProblemNum

argout = 0;
if bitand(mode,1)
    numf = numf + 1;
    argout = argout + 1;
    switch ProblemNum
        case(0)
            varargout(argout) = {x(2)*(5+x(1))};
        case(1)
            varargout(argout) = {x(1)*x(4)*(x(1) + x(2) + x(3))  +  x(3)};
        case(2)
            ob = 0;
            n = size(x,2);
            for i = 1:n
                ob = ob + (x(i)-i)^2;
            end
            varargout(argout) = {ob};
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
    numg = numg + 1;
    argout = argout + 1;
    switch ProblemNum
        case(0)
            g=[x(2); (5+x(1))];
        case(1)
            g=zeros(4,1);
            g(1) = x(4)*(2*x(1)+x(2)+x(3));
            g(2) = x(1)*x(4);
            g(3) = x(1)*x(4) + 1;
            g(4) = x(1)*(x(1)+x(2)+x(3));
        case(2)
            n = length(x);
            g=zeros(n,1);
            for i = 1:n
                g(i) = 2*(x(i)-i);
            end
            
        case(3)
            g=zeros(2,1);
            g(1) = 2*x(1)-2*x(2);
            g(2) = -2*x(1)+8*x(2);
        case(4)
            g=zeros(2,1);
            g(1) = x(2);
            g(2) = (5+x(1));
        case(5)
            g=zeros(2,1);
            g(1) =  2*x(1)-2*x(2);
            g(2) = -2*x(1)+8*x(2);
        case(6)
            g=zeros(2,1);
            g(1) = 2.0 * x(1);
            g(2) = 4.0 * x(2);
        case(7)
            g=zeros(2,1);
            g(1) = 2 * (x(1)-5);
            g(2) = 2 * (x(2)-5);
        case(8)
            g=zeros(2,1);
            g(1) = 2 * (x(1) - 5);
            g(2) = 0;
        case(9)
            g=zeros(1,1);
            g(1) = 2 * (x(1) - 5);
        case(10)
            g=zeros(2,1);
            g(1) =  2*x(1)-2*x(2);
            g(2) = -2*x(1)+8*x(2);
    end
    varargout(argout) = {g};
    
end
if bitand(mode,4)
    numH = numH + 1;
    argout = argout + 1;
    
    switch ProblemNum
        case(0)
            Hf=[0 1; 1 0];
        case(1)
            Hf=[2*x(4) x(4) x(4) 2*x(1)+x(2)+x(3);
                x(4)    0   0    x(1);
                x(4)    0   0    x(1);
                2*x(1)+x(2)+x(3) x(1)   x(1)    0;];
        case(2)
            Hf=2*eye(length(x));
    end
    varargout(argout) = {Hf};
    
end
return;
