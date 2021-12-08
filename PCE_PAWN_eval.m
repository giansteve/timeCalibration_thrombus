function [struct_MetaModel] = PCE_PAWN_eval(k_NonConst,Ns_const,output,int_num)
% compute the PCE with constant RV

struct_MetaModel = struct;
pt_evaluation_intervals = 0:int_num;

% do it for each non_const RVs
for m = 1:max(pt_evaluation_intervals)
    % Get input sample
    [struct_MetaModel(m).ED,~] = InputConst_PAWN(k_NonConst,Ns_const,pt_evaluation_intervals(m),max(pt_evaluation_intervals));
    % uq evaluation Model
    try % check if PCE is in time
        for tt = 1:size(output.time,1)
            struct_MetaModel(m).eval(:,tt) = uq_evalModel(output.time(tt),struct_MetaModel(m).ED);
        end
    catch
        struct_MetaModel(m).eval = uq_evalModel(output,struct_MetaModel(m).ED);
    end
    % eval = mean(eval_time);
    
end

end % function


