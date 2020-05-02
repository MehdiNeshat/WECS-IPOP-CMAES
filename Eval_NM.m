function [SumV]=Eval_NM(infeasible)
global array1
global Index
global sumv1
% close all
% scatter(array1(1,:),array1(2,:))
% axis([0 566 0 566])
N        = size(array1,2);
M        = size(Index,2);
Infe(1,:)= infeasible(1:M);
Infe(2,:)= infeasible(M+1:end);
SumV     = 0;

for i=1:M
    
    for j=1:N
        if j==Index(i)
            continue
        end
        X= (Infe(1,i)-array1(1,j))^2;
        Y= (Infe(2,i)-array1(2,j))^2;
        R=  sqrt(X+Y);
        if R <50
           SumV=SumV+(50-R) ;
        end
    end
end

if SumV < sumv1
   sumv1          = SumV;
   array1(:,Index)= Infe;
%    hold on
%    scatter(array1(1,:),array1(2,:),'filled')
end

end
