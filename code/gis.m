function [p,out_p_cond, kl] = gis(p_post, patterns, nVar, varSub, varargin)
% p_post:   the posterior probability of constraint, a (1,nPatterns) cell, 
%           each cell with size(1,nVar)
% patterns: with 1 or 0 to express the constraint pattern, with size
%           (nPatterns, nFeature)
% nVar:     # of  variables in each constraint, with size(1,nPatterns)
% varSub:   the sub of each p_post, a (1,nPatterns) cell, each cell with 
%           size(nVar, ndim) - varSub{i}(j,:) --> p_post{i}(j)
% varargin: output patterns like argument patterns

nPattern = size(patterns,1);  % # of patterns
patternList = cell(1,nPattern);  % as {1} {2} {4,5} {4,5,6}
for i=1:nPattern
   patternList{i} =  find(patterns(i,:)==1);
end

% computed given
post = cell(nPattern,1);
for i=1:nPattern
    post{i} = p_post{i}';
end

disp(['possible value counts for income, race and occupation are:',...
    num2str(nVar(1)), ',',num2str(nVar(2)),',',num2str(nVar(3))]) 

% initialize
p = ones(nVar);
p = p/sum(p(:));
disp(['KL divergence at ',num2str(0), ' : ',...
            num2str(KL_gen(p_post,p,patterns,nVar,varSub))])

iterMAX = 5;
kl = zeros(1,iterMAX+1);
kl(1) = KL_gen(p_post,p,patterns,nVar,varSub);

tic

for iter=1:iterMAX
    p_old = p;
    for i=1:nPattern
        pattern = patterns(i,:);  % pattern [3 5]
        [p_cond, p, perm, patternLen] = p_cond_gen(p,pattern);  
        % perm=[3 5 1 2 4 6]
        patternVars = prod(nVar(pattern==1));
        nVarNew1 = [nVar(pattern==1),nVar(pattern==0)];
        nVarNew2 = [patternVars,nVar(pattern==0)];
        p = reshape(p,nVarNew2);

        for j=1:size(varSub{i},1)
            if size(varSub{i},2) == 2
                % default 2D
                ind = sub2ind(nVar(pattern==1), varSub{i}(j,1),varSub{i}(j,2));
            elseif size(varSub{i},2) == 1
                ind = varSub{i}(j,1);
            end
            if p_cond(ind) == 0
                p(ind,:) = 0;
            else
                p(ind,:) = p(ind,:)*post{i}(j)/p_cond(ind);  % *(1-p_cond(j))/(1-post{i}(j))
                % u0 = u0*(1-post{i}(j))/(1-p_cond(j))
            end
        end
        p = reshape(p,nVarNew1);
        p = ipermute(p,perm);
    end
%     p=u0*p;  % [PROBLEM] an optional algo with offseted with factor
%             (1-p_cond(j))/(1-post{i}(j))
    res = abs(p-p_old);
    residual = sum(res(:));
    kl(iter+1) = KL_gen(p_post, p, patterns, nVar, varSub);
    if mod(iter,5) == 0
        disp(['KL divergence at ',num2str(iter), ' : ',...
            num2str(KL_gen(p_post, p, patterns,nVar, varSub))])
%         disp(['BIC:',num2str(BIC_gen(p_post,p,patterns,nVar))])
%         disp(['residual at ',num2str(iter), ' : ',...
%             num2str(residual)])
    end
end
toc

% disp(['BIC:',num2str(BIC_gen(p_post,p,patterns,nVar))])
% figure
% p_cond = p_cond_gen(p,patterns(1,:));
% subplot(3,1,1);bar([post{1}';p_cond']','grouped');
% legend('p\_post','p\_MEM','Location','NE');
% title({['marginal distribution of '];['with KL-divergence ',...
% num2str(KL(p_post{i},p_cond))]});
% p_cond = p_cond_gen(p,patterns(2,:));
% subplot(3,1,2);bar([post{2}';p_cond']');
% p_cond = p_cond_gen(p,patterns(3,:));
% subplot(3,1,3);bar([post{3}';p_cond']');
% 
% pattern = zeros(1,7);
% pattern([1 5]) = 1;
% figure
% p_cond = p_cond_gen(p,pattern);
% bar3(p_cond)

figure(10)
semilogy(abs(kl(1:end)))
xlabel('iteration');
xlim([1,iterMAX+1]);
title('Kullback-Leibler (KL) divergence');

hold on

if nargin == 6
    out_p_cond = cell(1,size(varargin{1},1));
    patterns = varargin{1};
    for i=1:size(patterns,1)
        pattern = varargin{1}(i,:);
        [p_cond, p, perm, patternLen] =  p_cond_gen(p,pattern);
        varSub = varargin{2};
        inds = [];
        for j=1:size(varSub{i},1)
            dim = size(varSub{i},2);
            if dim == 2
                % default 2D
                ind = sub2ind(nVar(pattern==1), varSub{i}(j,1),varSub{i}(j,2));
            elseif dim == 1
                ind = varSub{i}(j,1);
            end
            inds = [inds,ind];
        end
        out_p_cond{i} = p_cond(inds);
        p = ipermute(p,perm);
    end
end
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

function ret = KL(p_post,p)
    ret = p.*log(p./p_post) + (1-p).*log((1-p)./(1-p_post));
    ret(isinf(ret)) = 0;
    ret(isnan(ret)) = 0;
    ret = sum(ret(:));
end

