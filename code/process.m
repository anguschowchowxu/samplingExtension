columns = {'Pax_ID', 'P1_Age', 'P2_Gender', 'P8_Income', 'P5_EconActivity'};
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

