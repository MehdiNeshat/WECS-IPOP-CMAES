function [Coordinate,flag] = repair_function(arx,Opt)
global array1
global Index
global sumv1

sumv1=inf;

flag      = 1;
Coordinate= arx;
array1    = Coordinate;
%% fminsearchbnd parameters
opts                  = optimset('fminsearch');
opts.Display          =  'off';
opts.MaxFunEvals      = 10000;
fitnessf              = @Eval_NM;
opts.TolFun           = 0.001;
opts.OutputFcn        = @Stop_output_Func;
%%
in=0;
for i=1:Opt.Buoy_Num-1
    
    for j=i+1:Opt.Buoy_Num
        DisX=(Coordinate(1,i)-Coordinate(1,j))^2;
        DisY=(Coordinate(2,i)-Coordinate(2,j))^2;
        
        R   = sqrt(DisX+DisY);
        if R < 50
            in=in+1;
            Index(in)=i;
            break
        end
    end
    
end
Xstart =[Coordinate(1,Index),Coordinate(2,Index)];
LB     = ones(1,2*length(Index))*Opt.VarMin;
UB     = ones(1,2*length(Index))*Opt.VarMax;
[xsol,fval,exitflag,output] = fminsearchbnd(fitnessf,Xstart,LB,UB,opts);
Solution           = [xsol(1:length(Index));xsol(length(Index)+1:end)];
Coordinate(:,Index)= Solution;


for i=1:Opt.Buoy_Num-1
    
    for j=i+1:Opt.Buoy_Num
        DisX=(Coordinate(1,i)-Coordinate(1,j))^2;
        DisY=(Coordinate(2,i)-Coordinate(2,j))^2;
        
        R   = sqrt(DisX+DisY);
        if R < 50
            flag=0;
            break
        end
    end
    if flag==0
        break
    end
end


clear global 
end