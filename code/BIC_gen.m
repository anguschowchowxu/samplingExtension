function ret = BIC_gen(p_post,p,patterns,nVar)
% N: maybe the # of elements of p
% D: maybe the # of constraint category tuples
    D = 0;
    for i=1:size(patterns,1)
        D = D+sum(find(patterns(i,:)==1),2)+1;
    end
    N = 1;  % prod(nVar)
    ret = -2*max_entropy(p_post,p,patterns,nVar,D) + N*log(D);
end

function ret = max_entropy(p_post,p,patterns,nVar,D)
    nPattern = length(p_post);
    post = cell(nPattern,1);

    for i=1:nPattern
        if sum(p_post{i}(:)) == 0
            post{i} = 0;
        else
            post{i} = p_post{i}./sum(p_post{i}(:));
        end
    end
    
    ret = 0;
    for i=1:nPattern
        pattern = patterns(i,:);
        [p_cond, p, perm, patternLen] = p_cond_gen(p,pattern);
        
        patternVars = prod(nVar(find(pattern==1)));
        nVarNew1 = [nVar(find(pattern==1)),nVar(find(pattern==0))];
        nVarNew2 = [patternVars,nVar(find(pattern==0))];
        p = reshape(p,nVarNew2);

        for j=1:patternVars
            if p_cond(j) ~= 0
                ret = ret + post{i}(j)* log(p_cond(j));
            end
        end
        p = reshape(p,nVarNew1);
        p = ipermute(p,perm);
    end
    ret = ret*D;
end

function [p_cond, p, perm, patternLen] = p_cond_gen(p,pattern)
    perm = [find(pattern==1) find(pattern==0)];
    p = permute(p,perm);
    p_cond = p;
    patternLen = length(find(pattern==1));
    for i=patternLen+1:size(pattern,2)  % i-th dim
        p_cond = sum(p_cond,i);
    end
end