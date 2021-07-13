% input data curation, replace NA with 0,  ignore column with all NA
function T=tb_fix(T)
for i=size(T,2):-1:1
    if iscell(T{1,i})
        if isempty(T{1,i}{1})
            T(:,end)=[];
        end
    end    
end
for i=1:size(T,1)
    for j=1:size(T,2)
        if ~iscell(T{i,j})
            if isnan(T{i,j})
                T{i,j}=0;
            end
        end                
    end
end