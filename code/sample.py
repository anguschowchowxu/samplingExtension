#-*- encoding: utf-8 -*-

import pandas as pd
import numpy as np
import os
import sys

PROPORTION = 0.05
RANDOM_SEED = 1

def get_path():
	if len(sys.argv) == 1:
		global demo
		demo = input("which demo want to choose (1/2): [1]\n")
		if demo == '2':
			PATH = '../data/population_sample.csv'
		else:
			PATH = '../data/sample/population_sample.csv'
	elif len(sys.argv) > 1:
		PATH = sys.argv[1]
	PREFIX = PATH.replace('.csv','')
	return PATH, PREFIX

def gen_df(PATH, PREFIX):
	df = pd.read_csv(PATH, index_col=0)
	idx = np.random.choice(df.index.unique(), int(PROPORTION*len(df)), \
							replace=False)
	df_sample = df.loc[idx,:]
	print('{} entries of data'.format(len(df)))
	print('{} entries of sampled data'.format(len(df_sample)))
	columns = df.columns

	global df_hometype
	df_hometype = df.index.value_counts()
	if len(sys.argv) == 1:
		if demo == '2':
			df_hometype[df_hometype>7] = '8 or above'
		else:
			df_hometype[df_hometype>20] = '21 or above'
	list_hometypes = df_hometype.unique()

	tmp_df = pd.DataFrame(df_hometype.value_counts().values,columns=['population'],\
						index=df_hometype.value_counts().index.tolist())
	tmp_df['sample'] = 0
	tmp = df_sample.index.value_counts()
	if len(sys.argv) == 1:
		if demo == '2':
			tmp[tmp>20] = '8 or above'
		else:
			tmp[tmp>7] = '21 or above'
	df_hometype_sample = tmp.value_counts()
	for i in df_hometype_sample.index:
		tmp_df.loc[i, 'sample'] = df_hometype_sample[i]
		tmp_df = tmp_df.loc[list_hometypes,:]
	# tmp_df.to_csv('../data/sample/population_sample_hometype.csv')
	tmp_df.to_csv(PREFIX+'_hometype.csv')
	print('home will be seperated in the following types:\n')
	print(list_hometypes)
	return df, df_sample, columns, list_hometypes

def get_statistic(df, df_sample, columns, list_hometypes):
	
	for col in columns:
		tmp = df[col].value_counts()
		tmp_sample = df_sample[col].value_counts()
		tmp_df = pd.DataFrame(tmp.values,columns=['population'],\
							index=tmp.index.tolist())
		tmp_df['sample'] = 0.0

		for i in tmp_sample.index:
			tmp_df.loc[i, 'sample'] = tmp_sample[i]

		for hometype in list_hometypes:

			tmp_df['population_{}'.format(hometype)] = 0.0
			tmp_df['sample_{}'.format(hometype)] = 0.0
			df_home, df_home_sample = get_df_by_home(df, df_sample,\
								col, hometype)

			for i in df_home.index:
				tmp_df.loc[i,'population_{}'.format(hometype)] = df_home[i]
				if i in df_home_sample.index:
					tmp_df.loc[i,'sample_{}'.format(hometype)] = df_home_sample[i]

		# tmp_df.to_csv('../data/sample/population_sample_{}.csv'.format(col))
		tmp_df.to_csv(PREFIX+'_{}.csv'.format(col))

def get_df_by_home(df, df_sample, col, hometype):

	idx = df_hometype[df_hometype==hometype].index.tolist()
	df_home = df.loc[idx,col]
	idx = set(idx) & set(df_sample.index) 
	df_home_sample = df_sample.loc[list(idx),col]

	return df_home.value_counts(), df_home_sample.value_counts()

if __name__ == '__main__':
	PROPORTION = input('How much PROPORTION will to sample? (default 5% in all homes):\n'
		+'integer format in percentage, e.g., 6 for 6%\n')
	if not PROPORTION:
		PROPORTION = 0.05
	else:
		PROPORTION = float(PROPORTION)/100.0
	print('you mean {}%'.format(100*PROPORTION))
	PATH, PREFIX = get_path()
	df, df_sample, columns, list_hometypes,  = gen_df(PATH, PREFIX)
	get_statistic(df, df_sample, columns, list_hometypes)


	
