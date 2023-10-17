function [b] = ViP_fons(Q_inf,n,t,alpha)
for i = 1:size(n,1)
    for j = 1:size(n,1)
        b(i) = -Q_inf*n(i);


    end
end


end