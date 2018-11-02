function gisDemo()

    columns = {'Pax_ID', 'P1_Age', 'P2_Gender', 'P8_Income', 'P5_EconActivity'};
    cols = {'ID', 'Age', 'Gender', 'Income', 'EconActivity'};
    nFeature = length(columns);
    nVar = zeros(1, nFeature);
    p_post = cell(1, nFeature);
    p_pop = cell(1, nFeature);
    category = cell(1, nFeature);
    for i=1:nFeature
        df = readtable(['../data/sample/population_sample_',columns{i},'.csv']);
        p_post{i} = df.sample'/sum(df.sample,1);
        p_pop{i} = df.population'/sum(df.population,1);
        nVar(i) = size(df,1);
        category{i} = df.Var1';
    end
    patterns = eye(nFeature);
    
    varSub = cell(1,nFeature);  % 1d [[1;2;3]] 2d - [[1;1;4] [1;2;5]]
    for i=1:nFeature
        varSub{i} = (1:nVar(i))';
    end
    
    outpatterns = eye(5);
    tmp = eye(4);
    tmp = [zeros(4,1) tmp];
    tmp = [tmp;zeros(1,5)];
    outpatterns = outpatterns + tmp;

    [p,pcond] = gis(p_post, patterns, nVar, varSub, outpatterns);
    
%     pause()
%     gis(pcond, outpatterns, nVar, varSub);
    
    for i=1:nFeature
        figure(i)
        p_cond = p_cond_gen(p,patterns(i,:));
        bar([p_post{i};p_cond']');
        ax = gca();
        legend('p\_post','p\_MEM','Location','NE');
        title({['marginal distribution of ',cols{i}];['with KL-divergence ',...
            num2str(KL(p_post{i},p_cond))]});
        KL(p_post{i},p_cond);
        xticklabels(ax,category{i});
        ax.XTick=1:nVar(i);
        ax.XTickLabel=category{i};
        ax.XTickLabelRotation = 90;
%         saveas(gcf,['../report/image/marginal_2_',num2str(i),'.jpg']);
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