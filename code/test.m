nVar = [2 10 5 5 9 5 8];
A = rand(nVar); 
pattern = [0 1 1 0 0 0 0];
varSub = cell(1,7);
% varSub{1} = [[1 1 4];[2 3 5]];
% A_sub = cell(1,3);

tmp = nVar;
tmp(pattern==0) = 0;
[maxLen,ind] = max(tmp);
dim = sum(pattern==1)
if dim == 1
    outVarSub{1} = (1:nVar(ind))';
elseif dim == 2
    nVarPattern = prod(nVar(pattern==1));
    sample = randperm(nVarPattern,nVar(ind));
    [I1,I2] = ind2sub(nVar(pattern==1),sample);
    outVarSub{1} = [I1',I2'];
elseif dim == 3
    nVarPattern = prod(nVar(pattern==1));
    sample = randperm(nVarPattern,nVar(ind));
    [I1,I2,I3] = ind2sub(nVar(pattern==1),sample);
    outVarSub{1} = [I1',I2',I3'];
end
 
if 0
    nVar = [2 10 5 5 9 5 8];
    A = rand(nVar); 
    pattern = [0 1 1 0 0 0 0];
    varSub = cell(1,7);
    varSub{1} = [[1 1 4];[2 3 5]];
    A_sub = cell(1,3);

    for i=1:size(varSub{1},2)
        ind = sub2ind(nVar(pattern==1), varSub{1}(1,i),varSub{1}(2,i));
        perm = [find(pattern==1) find(pattern==0)];
        A = permute(A,perm);
        patternVars = prod(nVar(find(pattern==1)));
        nVarNew1 = [nVar(find(pattern==1)),nVar(find(pattern==0))];
        nVarNew2 = [patternVars,nVar(find(pattern==0))];

        A = reshape(A,nVarNew2);
        A_sub{i} = reshape(A(ind,:),nVar(find(pattern==0)));
        A = reshape(A,nVarNew1);
        A = ipermute(A,perm);
    end
    
    squeeze(A(:,1,2,:,:,:,:)) == A_sub{1};
end
    


% subs = zeros(numel(varSub{1}),size(varSub{1},1));
% ix = 1;
% for i=1:length(size(varSub{1}))
%     tmp = [];
%     for j=varSub{1}(i,:)
%         varSub{1}(i,:)
%         tmp = [tmp j];
%     end
%     size(tmp);
%     subs(i,:) = tmp;
% end
% subs

if 0
    nVar = [5 9 5 10 5];
    p = rand(nVar);
    pattern=[0 0 1 0 1];
    perm = [find(pattern==1) find(pattern==0)];
    nVarNew = [nVar(find(pattern==1)),nVar(find(pattern==0))];

    p_perm = permute(p,perm);
    p_reshape = reshape(p_perm,nVarNew);
    p_ans = ipermute(p_reshape,perm);

    test_ans = p_ans == p;
    sum(test_ans(:))
end