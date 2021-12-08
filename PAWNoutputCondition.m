function [idx] = PAWNoutputCondition(y,param)
% function for a more sofisticated condition for the computation of the
% PAWN sensitivity indices

for ii = 1:length(y)
    
    if y(ii) >= param(1) && y(ii) <= param(2)
        idx(ii) = true(1);
        
    else
        idx(ii) = false(1);
    end
    
end



end