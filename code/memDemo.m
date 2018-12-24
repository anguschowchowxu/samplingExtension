function [p,p_pop,category,patterns,varSub,p_post] = memDemo(varargin)
    global str prefix columns cols

    if nargin == 0
        idx_sample = 3;
        idx_population = 2;
    else
        idx_sample = varargin{1};
        idx_population = idx_sample-1;
    end
       
    params = struct;

    nFeature = length(columns);
    p_post = cell(1, nFeature);
    p_pop = cell(1, nFeature);
    nVar = zeros(1, nFeature);
    category = cell(1, nFeature);
    for i=1:nFeature
        if str=='1'
            df = readtable(['../data/sample/population_sample_',columns{i},'.csv']);
        elseif str=='2'
            df = readtable(['../data/population_sample_',columns{i},'.csv']);
        elseif str=='0'
            path = strjoin([prefix,'_',columns{i},'.csv'],'_');
            df = readtable(path);
        end
        if sum(df.(idx_sample),1) ~= 0
            p_post{i} = df.(idx_sample)'/sum(df.(idx_sample),1);
        else
            p_post{i} = zeros(size(df.(idx_sample)))';
        end
        p_pop{i} = df.(idx_population)'/sum(df.(idx_population),1);
        nVar(i) = size(df,1);
        category{i} = df.Var1';
    end
    patterns = eye(nFeature);
    
    varSub = cell(1,nFeature);  % 1d [[1;2;3]] 2d - [[1;1;4] [1;2;5]]
    for i=1:nFeature
        varSub{i} = (1:nVar(i))';
    end
    
    dim1 = idivide(nFeature, int32(2), 'ceil');
    outpatterns = zeros(dim1,nFeature);
    for i=1:dim1
        if 2*i < nFeature
            outpatterns(i,2*i-1:2*i) = 1;
        else
            outpatterns(i,2*i-2:2*i-1) = 1;
        end
    end
   
%     tmp = cell(1,11);i=1;
    
    if nargin == 0
        for MIN = 1
            outVarSub = gen_outVarSub(outpatterns,nVar,MIN);
            tic
%             [p,pcond] = mem(p_post, patterns, nVar, varSub);
            params.outpatterns = outpatterns;
            params.outVarSub = outVarSub;
            params.graph = 1;
            params.MaxIter = 1;
            [p,pcond] = mem(p_post, patterns, nVar, varSub,...
                            params);
%             mem(pcond, outpatterns, nVar, outVarSub);
            toc

    %         tmp{i} = num2str(MIN);
    %         i = i+1;
    %         pause()
        end

    %     ax = gca();
    %     legend(ax,tmp);



        for i=1:nFeature
            figure(i)
            p_cond = p_cond_gen(p,patterns(i,:));
            bar([p_pop{i};p_cond']');
            ax = gca();
            legend('p\_pop','p\_MEM','Location','NE');
            title({['marginal distribution of ',cols{i}];['with KL-divergence ',...
                num2str(KL(p_pop{i}',p_cond,patterns(i,:),nVar,varSub{i}))]});
            ax.XTick = 1:nVar(i);
            ax.XTickLabel = category{i};
            ax.XTickLabelRotation = 90;
%             saveas(gcf,['../report/image/partial_marginal_homtype_maxKL_',num2str(i),'.jpg']);
        end
    else
        tic
        p = mem(p_post, patterns, nVar, varSub,params);
        toc
    end
end

function outVarSub = gen_outVarSub(outpatterns,nVar,MIN)
    outVarSub = cell(1,size(outpatterns,1));
    for i=1:size(outpatterns,1)
        pattern = outpatterns(i,:);
        tmp = nVar;
        tmp(pattern==0) = 0;
        [maxLen,ind] = max(tmp);
        dim = sum(pattern==1);
        if dim == 1
            outVarSub{i} = (1:nVar(ind))';
        elseif dim == 2
            nVarPattern = prod(nVar(pattern==1));
            sample = randperm(nVarPattern,nVar(ind));  % nVar(ind),min([nVar(ind),MIN])
            [I1,I2] = ind2sub(nVar(pattern==1),sample);
            outVarSub{i} = [I1',I2'];
        elseif dim == 3
            nVarPattern = prod(nVar(pattern==1));
            sample = randperm(nVarPattern,nVar(ind));
            [I1,I2,I3] = ind2sub(nVar(pattern==1),sample);
            outVarSub{i} = [I1',I2',I3'];
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