#-*- encoding: utf-8 -*-

import pandas as pd
import os

PROPORTION = 0.1
RANDOM_SEED = 1
df = pd.read_csv('../data/sample/population_sample.csv')
df_sample = df.sample(frac=PROPORTION, replace=False,  random_state=RANDOM_SEED, axis=0)

columns = df.columns[1:]

def get_statistic(df, df_sample, columns):
	
	for col in columns:
		tmp = df[col].value_counts()
		tmp_sample = df_sample[col].value_counts()
		tmp_df = pd.DataFrame(tmp.values,columns=['population'],\
							index=tmp.index.tolist())
		tmp_df['sample'] = 0.0
		for i in tmp_sample.index:
			tmp_df.loc[i, 'sample'] = tmp_sample[i]
		tmp_df.to_csv('../data/sample/population_sample_{}.csv'.format(col))

if __name__ == '__main__':
	get_statistic(df, df_sample, columns)


	
