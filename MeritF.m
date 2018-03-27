function varargout = MeritF(x,Obj_F,C_eq,C_Ineq,y,z,mode)
% x is the input, y and z are the lagrange multiplier ascociated with 
% the equality and inequality constraints.

global numf numg numH

argout = 0;
if bitand(mode,1) 
  numf = numf + 1;
  argout = argout + 1;
  varargout(argout) = {feval(Obj_F,x,1)+y'*abs(feval(C_eq,x,1))+z'*abs(feval(C_Ineq,x,1))};
end
if bitand(mode,2) 
  numg = numg + 1;
  argout = argout + 1;
  varargout(argout) = {feval(Obj_F,x,2)+y'*abs(feval(C_eq,x,2))+z'*abs(feval(C_Ineq,x,2))};
end
if bitand(mode,4) 
  numH = numH + 1;
  argout = argout + 1;
  varargout(argout) = {[0 1; 1 0]};
end
return;
