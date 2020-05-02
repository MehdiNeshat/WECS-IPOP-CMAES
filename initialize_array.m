function [position]=initialize_array(Opt)

while(1)
    
   position       = rand(Opt.nVar,1)*Opt.UBounds;
   coordinate(1,:)= position(1:2:end);
   coordinate(2,:)= position(2:2:end);
   
   [flag, violation] = Control_feasibility(coordinate);
   
   if flag==1
       break
   end
end

end