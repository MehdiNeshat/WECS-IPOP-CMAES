function [power]=Eval_Power(position)

global siteOpts
global buoy

if size(position,2)==1
    position=position';
end
Opt.Buoy_Num = round(size(position,2)/2);
Opt.VarMin   = 0;
Opt.VarMax   = round( sqrt(Opt.Buoy_Num *20000));

array.number               = Opt.Buoy_Num;
array.radius               = 5* ones(1,array.number);
array.sphereCoordinate(1,:)= position(1:2:end);
array.sphereCoordinate(2,:)= position(2:2:end);
array.sphereCoordinate(3,:)= -8;

[flag, violation] = Control_feasibility(array.sphereCoordinate(1:2,:));

 [Parray_AAP, ~, ~] = arrayBuoyPlacement_v20180601(array, siteOpts, buoy);
if flag==1
   
    
    power=-Parray_AAP;
else
    disp(['sum violation=',num2str(violation)])
    power_temp=Parray_AAP-(violation*5*10^4);
    if power_temp <= 0
        power = 0;
    else
        power = -power_temp;
    end
    
end
end