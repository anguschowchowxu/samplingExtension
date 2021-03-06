function ret = KL_gen(p_post,p,patterns, nVar, varSubs)
    nPattern = length(p_post);
    post = cell(nPattern,1);

    for i=1:nPattern
        post{i} = p_post{i}';
    end
    ret = 0;
    for i=1:nPattern
        pattern = patterns(i,:); 
        [p_cond, p, perm, patternLen] = p_cond_gen(p,pattern);
        ret = ret + KL(post{i},p_cond,pattern,nVar,varSubs{i});
        p = ipermute(p,perm);
    end
end

function ret = KL(p_post,p,pattern,nVar,varSub)
    inds = [];
    for j=1:size(varSub,1)
        dim = size(varSub,2);
        if dim == 2
            % default 2D
            ind = sub2ind(nVar(pattern==1), varSub(j,1),varSub(j,2));
        elseif dim == 1
            ind = varSub(j,1);
        end
        inds = [inds,ind];
    end
    p = p(inds');
    ret = p.*log(p./p_post) + (1-p).*log((1-p)./(1-p_post));
    ret(isinf(ret)) = 0;
    ret(isnan(ret)) = 0;
    ret = sum(ret(:));
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