#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd

""" 
https://www.nature.com/scitable/definition/hardy-weinberg-equilibrium-122/ 
"""


def hwe(
		obs_hets: int | float, obs_hom1: int | float, obs_hom2: int | float
) -> float:
	""" Python interpretation hwe - https://github.com/jeremymcrae/snphwe

	:param obs_hets: - количество наблюдаемых гетерозигот
	:param obs_hom1: - количество наблюдаемыех гомозигот1
	:param obs_hom2: - количество наблюдаемыех гомозигот2
	:return: - сдесь возвращается p-value
	"""

	obs_hets = round(obs_hets)
	obs_hom1 = round(obs_hom1)
	obs_hom2 = round(obs_hom2)

	if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
		raise ValueError("snphwe: negative allele count")

	obs_homr = min(obs_hom1, obs_hom2)
	obs_homc = max(obs_hom1, obs_hom2)

	rare = 2 * obs_homr + obs_hets
	genotypes = obs_hets + obs_homc + obs_homr

	if genotypes == 0:
		raise ValueError("snphwe: zero genotypes")

	probs = np.zeros(round(rare) + 1)

	# get distribution midpoint, but ensure midpoint and rare alleles have
	# same parity
	mid = int(rare * (2 * genotypes - rare) / (2 * genotypes))
	if mid % 2 != rare % 2:
		mid += 1

	probs[mid] = 1.0
	_sum = probs[mid]

	curr_homr = (rare - mid) / 2
	curr_homc = genotypes - mid - curr_homr
	curr_hets = mid
	while curr_hets > 1:
		probs[curr_hets - 2] = (
			probs[curr_hets] * curr_hets * (curr_hets - 1.0)
			/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
		)
		_sum += probs[curr_hets - 2]

		# fewer heterozygotes -> add one rare, one common homozygote
		curr_homr += 1
		curr_homc += 1
		curr_hets -= 2

	# calculate probabilities from midpoint up
	curr_homr = (rare - mid) / 2
	curr_homc = genotypes - mid - curr_homr

	curr_hets = mid
	while curr_hets <= rare - 2:
		probs[curr_hets + 2] = \
			(probs[curr_hets] * 4.0 * curr_homr * curr_homc
			 / ((curr_hets + 2.0) * (curr_hets + 1.0)))
		_sum += probs[curr_hets + 2]

		# add 2 heterozygotes -> subtract one rare, one common homozygote
		curr_homr -= 1
		curr_homc -= 1
		curr_hets += 2

	# p-value calculation for p_hwe
	target = probs[obs_hets]
	p_hwe = 0.0

	for p in probs:
		if p <= target:
			p_hwe += p / _sum

	return min(1.0, p_hwe)


def hwe_test(
		seq_snp: pd.Series, freq: float, crit_chi2: float = 3.841
) -> bool:
	""" The Hardy-Weinberg equilibrium is a principle stating that the genetic
	variation in a population will remain constant from one generation to the
	next in the absence of disturbing factors

	:param seq_snp: -
	:param freq: -
	:param crit_chi2: -  The critical value for a test ("either / or":
		observed and expected values are either one way or the other),
		therefore with degrees of freedom = df = 1 is 3.84 at p = 0.05
	:return: - Возвращается решение исключить или оставить проверяемый снп
	"""

	_seq = seq_snp.replace(5, np.nan)

	if _seq.nunique() == 1:
		return True

	n_genotypes = _seq.count()

	observed = {
		0: (_seq == 0).sum(),
		1: (_seq == 1).sum(),
		2: (_seq == 2).sum()
	}

	expected = {
		0: ((1 - freq) ** 2) * n_genotypes,
		1: (2 * ((1 - freq) * freq)) * n_genotypes,
		2: (freq ** 2) * n_genotypes
	}

	chi = sum([
		((obs - exp) ** 2) / exp
		for obs, exp in zip(observed.values(), expected.values())
	])

	if chi > crit_chi2:
		return False
	else:
		return True