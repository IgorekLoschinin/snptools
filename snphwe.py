import numpy as np


def snphwe(obs_hets: int, obs_hom1: int, obs_hom2: int) -> float:
	""" https://www.nature.com/scitable/definition/hardy-weinberg
	-equilibrium-122/ """
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
		probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (
				curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (
				curr_homc + 1.0))
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
			probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
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


import unittest


class TestHWE(unittest.TestCase):
	def test_snphwe(self):
		''' check snphwe gives expected p-values
		'''
		self.assertAlmostEqual(snphwe(500, 10, 5000), 0.65157189991)
		self.assertAlmostEqual(snphwe(1000, 20, 5000), 1.26598491e-5)

	def test_snphwe_odd_inputs(self):
		''' check snphwe with odd inputs
		'''
		# should raise errors with odd inputs
		with self.assertRaises(ValueError):
			snphwe(0, 0, 0)
		with self.assertRaises(ValueError):
			snphwe(-5, 10, 1000)

	def test_snphwe_large_input(self):
		''' check snphwe doesn't give errors with large sample sizes
		'''
		self.assertEqual(snphwe(200000, 200000, 200000), 0.0)

	def test_snphwe_uncertain_genotypes(self):
		''' check uncertain genotypes give correct p-values

		Imputed genotypes need HWE checks too. I checked this against the method
		from:
		Shriner D. (2011) Genet Epidemiol. 35:632‐637. doi:10.1002/gepi.20612
		Approximate and exact tests of Hardy-Weinberg equilibrium using uncertain genotypes.

		Instead of running the full exact test from that, I simply round the
		genotype totals before calling the underlying c++ function. This should
		still be reasonably accurate and fast.
		'''
		self.assertAlmostEqual(snphwe(4989.99999, 494999.999, 9.9999999),
							   0.570223198305)