function [p,prefix,columns] = homeMEM(varargin)
    clc; clf;
    
    if nargin==0
        prompt = ['Which demo do you want? 1/2 [1]: \n1 - for 5 dimension population data; \n2 - for 8 dimesion population data;\n'];
        global str
        str = input(prompt,'s');
        if isempty(str)
            str = '1';
        end
    else
        str = '0';
        path = varargin{1};
    end
        
    global prefix
    if str=='1'
        df = readtable('../data/sample/population_sample_hometype.csv');
        df2 = readtable('../data/sample/population_sample.csv');
        tmp = strsplit('../data/sample/population_sample.csv','.csv');
        prefix = tmp(1);
    elseif str=='2'
        df = readtable('../data/population_sample_hometype.csv');
        df2 = readtable('../data/population_sample.csv');
        tmp = strsplit('../data/population_sample.csv','.csv');
        prefix = tmp(1);
    elseif str=='0'
        tmp = strsplit(path,'.csv'); 
        prefix = tmp(1);
        path2 = strjoin([prefix,'hometype.csv'],'_');
        df2 = readtable(path);
        df = readtable(path2);
    end
    
    figure(11)
%     yyaxis left
    bar([df.(2),df.(3)]);
    xlabel('Hometype of different members')
    ylabel('Counts')
    ax = gca();
    legend('hometype\_pop','hometype\_sample','Location','NW');
    title(['sum of hometypes with different members counts']); 
    ax.XTick=1:length(df.(1)');
    xlim([0,length(df.(1)')+2]);
    ax.XTickLabel=df.(1)';
    hometype = df.(1)';
    ax.XTickLabelRotation = 90;
%     saveas(gcf,['../report/image/hometype_distribution_2.jpg']);


    proportion = sum(df.(3),1)/sum(df.(2),1);
    p_p = df.(2)'/sum(df.(2),1);
    p_s = df.(3)'/sum(df.(3),1);
    disp(['sampled ',num2str(proportion*100), '% population data ',...
        'with 5% homes'])

    global columns cols
    columns = df2.Properties.VariableNames(2:end);
    cols = cell(1,length(columns));
    for i=1:length(columns)
        tmp = strsplit(columns{i},'_');
        if strcmp(tmp{1},'Area')
            cols{i} = 'Area Name';
        else
            cols{i} = tmp{2};
        end
    end
    
    nFeature = length(columns);
    nVar = zeros(1, nFeature);
    for i=1:nFeature
        if str=='1'
            df = readtable(['../data/sample/population_sample_',columns{i},'.csv']);
        elseif str=='2'
            df = readtable(['../data/population_sample_',columns{i},'.csv']);
        elseif str=='0'
            tmp = strjoin([prefix,'_',columns{i},'.csv'],'');
            df = readtable(tmp);
        end
        nVar(i) = size(df,1);
    end
    
    p = zeros(nVar);
    sum(p_s);
    
    [tmp,p_pop,category,patterns,varSub,p_post] = memDemo();
    disp(['MEM result without considering hometype, [ENTER] to continue'])
    pause();
    
    prompt = ['How much size want to extend?[default 1000]:'];
    exsize = input(prompt);
    if isempty(exsize)
        exsize = 1000;
    end
    ret_df = table;
    for i=1:nFeature+1
       ret_df.(i) = 1;
    end
    ret_df(1,:) = [];
    kl_div = zeros(1,length(hometype));
    for i=1:length(hometype)
        df_tmp = table;
        p_tmp = memDemo(3+2*i);
        p = p + p_s(i).*abs(p_tmp); 
        size_tmp = floor(exsize*p_s(i));
        df_tmp.(1) = repmat(hometype(i),size_tmp,1);
        for j=1:nFeature
            p_cond = p_cond_gen(p_tmp,patterns(j,:));
            df_tmp.(j+1) = randsample(category{j},size_tmp,...
                true,abs(p_cond)')';
            kl_div(i) = kl_div(i) + ...
                KL(p_pop{j}',p_cond,patterns(j,:),nVar,varSub{j});
        end
        ret_df = [ret_df;df_tmp];
    end
    p = p/sum(p(:));
    
    size(ret_df)
    columns_tmp(2:nFeature+1) = columns;
    columns_tmp{1} = 'hometype'
    ret_df.Properties.VariableNames = columns_tmp;
    strjoin([prefix,'extension.csv'],'_')
    writetable(ret_df, strjoin([prefix,'extension.csv'],'_'));
    
    disp('MEM result without considering hometype')
    
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
        legend('p\_pop','p\_post','p\_MEM\_hometype','Location','NE');
        title({['marginal distribution of ',cols{i}];['with KL-divergence ',...
            num2str(KL(p_pop{i}',p_cond,patterns(i,:),nVar,varSub{i}))]});
        ax.XTick=1:nVar(i);
        ax.XTickLabel=category{i};
        ax.XTickLabelRotation = 90;
%         saveas(gcf,['../report/image/partial_marginal_hometype_',num2str(i),'.jpg']);
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