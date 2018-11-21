function homeDemo()
    clc; clear all; clf;

    df = readtable(['../data/sample/population_sample_hometype.csv']);
    figure(6)
%     yyaxis left
    bar([df.(2),df.(3)]);
    xlabel('Hometype of different members')
    ylabel('Counts')
    ax = gca();
    legend('hometype\_pop','hometype\_sample','Location','NW');
    title(['sum of hometypes with different members counts']);
    ax.XTick=1:length(df.(1)');
    ax.XTickLabel=df.(1)';
    ax.XTickLabelRotation = 90;
%     saveas(gcf,['../report/image/hometype_distribution.jpg']);
    proportion = sum(df.(3),1)/sum(df.(2),1)
    p_p = df.(2)'/sum(df.(2),1);
    p_s = df.(3)'/sum(df.(3),1);
    nVar = [11,17,2,13,13];
    p = zeros(nVar);
    hometype = df.(1);
    sum(p_s)

    columns = {'Pax_ID', 'P1_Age', 'P2_Gender', 'P8_Income', 'P5_EconActivity'};
    cols = {'ID', 'Age', 'Gender', 'Income', 'EconActivity'};
    nFeature = length(columns);
    [tmp,p_pop,nVar,category,patterns,varSub,p_post] = gisDemo();
    pause();
    

    kl_div = zeros(1,length(hometype));
    for i=1:length(hometype)
        p_tmp = gisDemo(3+2*i);
        p = p + p_s(i).*abs(p_tmp);  
        for j=1:nFeature
            p_cond = p_cond_gen(p_tmp,patterns(j,:));
            kl_div(i) = kl_div(i) + ...
                KL(p_pop{j}',p_cond,patterns(j,:),nVar,varSub{j});
        end
    end
    p = p/sum(sum(sum(sum(sum(p,1),2),3),4),5);

    figure(7)
%     yyaxis left
    plot(abs(kl_div));
    xlabel('Hometype of different members');
    ylabel('KL-divergence');
    title('KL-divergence of different hometype in distribution');
    ax = gca();
    ax.XTick=1:length(df.(1)');
    ax.XTickLabel=df.(1)';
    ax.XTickLabelRotation = 90;
%     saveas(gcf,['../report/image/hometype_distribution_KL.jpg']);
    
    for i=1:nFeature
        figure(i)
        p_cond = p_cond_gen(p,patterns(i,:));
        bar([p_pop{i};p_post{i};p_cond']');
        ax = gca();
        legend('p\_pop','p\_post','p\_MEM','Location','NE');
        title({['marginal distribution of ',cols{i}];['with KL-divergence ',...
            num2str(KL(p_pop{i}',p_cond,patterns(i,:),nVar,varSub{i}))]});
        ax.XTick=1:nVar(i);
        ax.XTickLabel=category{i};
        ax.XTickLabelRotation = 90;
        saveas(gcf,['../report/image/partial_marginal_hometype_',num2str(i),'.jpg']);
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
            ind = sub2ind(nVar(pattern==1), varSub(j,1),varSub(j,2))
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