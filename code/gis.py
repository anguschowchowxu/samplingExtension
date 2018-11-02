#-*- encoding: utf-8 -*-

import pandas as pd
import numpy as np

population = pd.read_csv('../data/population.csv', sep='\t', index_col=0)
race = pd.read_csv('../data/race.csv', sep='\t', index_col=0)

population['Total Population'] = population['Total Population']\
					.apply(lambda x: float(x.replace(',','')), 0)\
					.values.astype('float64')
values = race.applymap(lambda x: 0 if x=='N' else float(x.replace(',','')))\
					.values.astype('float64')
race = pd.DataFrame(values, index=race.index.tolist(), columns=race.columns)
race.index = race.index.set_names('Geography')

counts = len(population)
assert (population.index == race.index).sum() == counts, \
				"The Geography info are not same"

test = race.iloc[:,1:3].div(race.iloc[:,1:3].sum().values, axis=1).values
# 0th col for p(x|'Hispanic or Latino'),
# 1th col for p(x|'White alone, Not Hispanic or Latino')
p = np.ones([counts,counts], dtype=np.float64)
p /= np.sum(p)

res = 10
iterN = 0
tol = 1e-5
basic = False
advance = True
u0 = 1

def factor(p_post, p_pri):
	return p_post/p_pri * (1-p_pri)/(1-p_post)

def factor2(p_post, p_pri):
	return (1-p_post)/(1-p_pri)

while res>tol:
	p_old = np.copy(p)

	for i in range(counts):
		for j in range(counts):

			# basic algo
			if basic:
				p[i,j] = p_old[i,j] * test[i,0] / np.sum(p_old[i,:]) * \
							test[j,1] / np.sum(p_old[:,j])

			# [PROBLEM: don't know why] algo in wu at el.
			if advance:
				p[i,j] = p_old[i,j] * factor(test[i,0],np.sum(p_old[i,:])) \
									* factor(test[j,1],np.sum(p_old[:,j])) \

				u0 = u0 * factor2(test[i,0], np.sum(p_old[i,:])) \
						* factor2(test[j,1], np.sum(p_old[:,j]))

			if np.isnan(p[i,j]):  # !!
				p[i,j] = 0.0

	if advance:
		p = np.multiply(p, u0)
	res = np.sum(np.abs(np.subtract(p,p_old)))
	iterN += 1
	print('iter {}: residual is {}'.format(iterN, res))

