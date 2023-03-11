import pandas as pd
import numpy as np


class Snp(pd.DataFrame):
    alleles_freq = pd.Series(dtype=object)

    def __init__(self, filename, precalculate_freq=False):

        super().__init__(pd.read_csv(filename, index_col=0).astype(np.int8))
        if precalculate_freq:
            self.alleles_freq = self.calculate_allele_freq()

    def get_markers(self, markers):
        markers_intersection = self.index.intersection(markers)
        return self.loc[markers_intersection]

    def one_parent_verification(self, calf_number, parent_number):
        calf_parent_table = self[[calf_number, parent_number]].replace(5, np.nan).dropna()
        only_homozygous = calf_parent_table.loc[calf_parent_table.isin([0, 2]).all(axis=1)]
        same_genotype_bool = only_homozygous.apply(lambda x: x[0] == x[1], axis=1)
        return len(same_genotype_bool) - same_genotype_bool.sum()

    def animal_call_rate(self):
        return self.apply(lambda x: 1 - ((x == 5).sum() / len(x)))

    def marker_call_rate(self):
        return self.apply(lambda x: 1 - ((x == 5).sum() / len(x)), axis=1)

    def calculate_allele_freq(self):
        drop_missing_markers = self.replace(5, np.nan)

        def _calc_freq_for_row(row):
            num_of_ref_alleles = row.sum()
            if num_of_ref_alleles == 0:
                return num_of_ref_alleles
            return num_of_ref_alleles / (row.count() * 2)

        return drop_missing_markers.apply(_calc_freq_for_row, axis=1)

    def calculate_maf(self):
        allele_freq = self.calculate_allele_freq() if self.alleles_freq.empty else self.alleles_freq

        def _replace_freq(val):
            if val > 0.5:
                return round(1 - val, 3)
            return round(val, 3)

        return allele_freq.apply(lambda x: _replace_freq(x))

    def linkage_disequilibrium(self):
        df = self.replace(5, np.nan).T.astype(np.float16)
        return df.corr().applymap(lambda x: round(x * x, 3))

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
            obs = {0: (row == 0).sum(), 1: (row == 1).sum(), 2: (row == 2).sum()}
            b_allele_freq = allele_freq[row.name]
            est = {0: ((1 - b_allele_freq) ** 2) * n_genotypes,
                   1: (2 * ((1 - b_allele_freq) * b_allele_freq)) * n_genotypes,
                   2: (b_allele_freq ** 2) * n_genotypes}
            chi = sum([((o - e) ** 2) / e for o, e in zip(obs.values(), est.values())])

            if chi > critical_value_chi:
                return False
            else:
                return True

        return replace_missing_markers.apply(calc_hwe, axis=1)


if __name__ == '__main__':
    path_to_bovine_csv = 'all_50K.csv'
    snp_data = Snp(path_to_bovine_csv, precalculate_freq=False)
    print(snp_data.hwe_test().to_csv('hwe.csv'))
