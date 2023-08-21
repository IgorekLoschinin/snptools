#!/usr/bin/env python
# coding: utf-8

import numpy as np


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
	:return: -
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


def hwe_test(self):
	allele_freq = self.calculate_allele_freq() if self.alleles_freq.empty else self.alleles_freq
	replace_missing_markers = self.replace(5, np.nan)
	critical_value_chi = 3.841

	def calc_hwe(row):
		n_unique_values = row.nunique()
		if n_unique_values == 1:
			return True
		if n_unique_values == 0:
			return np.nan

		n_genotypes = row.count()

		obs = {
			0: (row == 0).sum(),
			1: (row == 1).sum(),
			2: (row == 2).sum()
		}
		b_allele_freq = allele_freq[row.name]

		est = {
			0: ((1 - b_allele_freq) ** 2) * n_genotypes,
			1: (2 * ((1 - b_allele_freq) * b_allele_freq)) * n_genotypes,
			2: (b_allele_freq ** 2) * n_genotypes
		}

		chi = sum([
			((o - e) ** 2) / e
			for o, e in zip(obs.values(), est.values())
		])

		if chi > critical_value_chi:
			return False
		else:
			return True

	return replace_missing_markers.apply(calc_hwe, axis=1)
