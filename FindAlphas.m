function [alpha_s, alpha_z] = FindAlphas(x,s,y,z,px,ps,py,pz,objF,Ceq,CIneq,MeritF,taw)

alpha_s=1.0;
alpha_z=1.0;
N=10;
for i=1:N
    if(sum(s+alpha_s*ps-(1-taw)*s>=0)==length(s))
        break;
    else
        alpha_s=alpha_s*0.9;
    end
end

for i=1:N
    if(sum(z+alpha_z*pz-(1-taw)*z>=0)==length(z))
        break;
    else
        alpha_z=alpha_z*0.9;
    end
end


merit_old=feval(MeritF,x,objF,Ceq,CIneq,y,z,1);

for i=1:10
        x_new=x+alpha_s*px;
        s_new=s+alpha_s*ps;
        y_new=y+alpha_z*py;
        z_new=z+alpha_z*pz;
        merit_new=feval(MeritF,x_new,objF,Ceq,CIneq,y_new,z_new,1);
        if(merit_new>merit_old)
            break;
        else
            alpha_s=alpha_s*0.95;
            alpha_z=alpha_z*0.95;
        end
end


end