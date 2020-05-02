function [flag, violation] = Control_feasibility(position)

violation=0;
for i=1:size(position,2)-1
    for j=i+1:size(position,2)
        
        D1= (position(1,i) - position(1,j))^2;
        D2= (position(2,i) - position(2,j))^2;
        Dis= sqrt (D1 + D2);
        if (Dis < 50)
            violation=violation+ (50-Dis);
        end
        
        
    end
end
if violation == 0
    flag=1;
else
    flag=0;
end

end