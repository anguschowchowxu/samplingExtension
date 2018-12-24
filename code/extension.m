function extension()
    prompt = ['Please input your location of sample csv file:\n'];
    path = input(prompt,'s');
    [p,prefix,columns] = homeMEM(path);
    
    nFeature = length(columns)+1;
    p_post = cell(1, nFeature);
    category = cell(1, nFeature);
    patterns = eye(nFeature);
    
    path2 = strjoin([prefix,'hometype.csv'],'_');
    df = readtable(path2);
    category{1} = df.(1)';
    p_post{1} = df.(3)'/sum(df.(3),1);
    

    for i=2:nFeature
        path = strjoin([prefix,'_',columns{i-1},'.csv'],'');
        df = readtable(path);
        p_post{i} = p_cond_gen(p,patterns(i-1,:))';
        category{i} = df.Var1';
    end
    
    tmp = table;
    for i=1:nFeature
        tmp.(i) = randsample(category{i},100,true,abs(p_post{i}))';
    end
    columns(2:nFeature) = columns;
    columns{1} = 'hometype';
    tmp.Properties.VariableNames = columns;
    strjoin([prefix,'extension.csv'],'_')
    writetable(tmp, strjoin([prefix,'extension.csv'],'_'));
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