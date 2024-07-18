"""
Unit tests for ld.py.  This should be run via pytest.
"""

import os
import sys
main_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
test_directory = os.path.abspath(os.path.join(main_directory, 'test'))
data_directory = os.path.abspath(os.path.join(test_directory, 'data'))
sys.path.append(main_directory)

import numpy as np
import pandas as pd
import pytest

from mama.pop_info import PopInfo

from mama.util.bed import BED_SUFFIX, write_bed_file
from mama.util.bim import (BIM_SUFFIX, BIM_COL_TYPES, BIM_CHR_COL, BIM_RSID_COL, BIM_CM_COL,
                           BIM_BP_COL, BIM_A1_COL, BIM_A2_COL)
from mama.util.fam import (FAM_SUFFIX, FAM_FID_COL, FAM_IID_COL, FAM_IIDF_COL, FAM_IIDM_COL,
                           FAM_SEX_COL, FAM_PHEN_COL)
from mama.util.sumstats import COMPLEMENT


def bedbimfam(bedbimfam_prefix, M, N, G, rsids, bps, a1s, a2s):
    assert G.shape == (M, N)
    assert len(rsids) == len(bps) == len(a1s) == len(a2s) == M
    bed_filename = f"{bedbimfam_prefix}{BED_SUFFIX}"

    bim_df = pd.DataFrame(
        {
            BIM_CHR_COL : [1] * M,
            BIM_RSID_COL : rsids,
            BIM_CM_COL : [0] * M,
            BIM_BP_COL : bps,
            BIM_A1_COL : a1s,
            BIM_A2_COL : a2s
        }
    )
    bim_filename = f"{bedbimfam_prefix}{BIM_SUFFIX}"

    fam_df = pd.DataFrame(
        {
            FAM_FID_COL : list(range(10, 10*N+10, 10)),
            FAM_IID_COL : list(range(1, N+1)),
            FAM_IIDF_COL : [0] * N,
            FAM_IIDM_COL : [0] * N,
            FAM_SEX_COL : [0] * N,
            FAM_PHEN_COL : [1] * N
        }
    )
    fam_filename = f"{bedbimfam_prefix}{FAM_SUFFIX}"

    write_bed_file(bed_filename, G)
    bim_df.to_csv(bim_filename, sep="\t", header=False, index=False)
    fam_df.to_csv(fam_filename, sep="\t", header=False, index=False)

    return bedbimfam_prefix

def is_swap(pop1_a1, pop1_a2, pop2_a1, pop2_a2):
    in_agreement1 = (pop1_a1 == pop2_a1 and pop1_a2 == pop2_a2)
    in_agreement2 = (pop1_a1 == COMPLEMENT[pop2_a1] and pop1_a2 == COMPLEMENT[pop2_a2])

    return not (in_agreement1 or in_agreement2)


LD_SCORE_DATA_DIR = os.path.join(data_directory, "ld_scores")


POP1_G = np.array(
    [[0,1,2,0],
     [1,2,0,1],
     [2,0,1,2],
     [0,1,2,0],
     [1,1,2,2],
     [0,0,1,1],
     [0,1,0,1]], dtype=float)
POP1_COV = np.cov(POP1_G)
POP1_M, POP1_N = POP1_G.shape
POP1_RSIDs = ['RS01', 'RS03', 'RS04', 'RS05', 'RS06', 'RS07', 'RS09']
POP1_BPs = [10, 30, 40, 50, 60, 70, 90]
POP1_A1s = ['A'] * POP1_M
POP1_A2s = ['G'] * POP1_M

POP2_G = np.array(
    [[2,2,1,0],
     [2,2,0,1],
     [0,0,1,2],
     [0,1,2,1],
     [1,1,0,2],
     [0,0,0,1],
     [0,2,1,1]], dtype=float)
POP2_COV = np.cov(POP2_G)
POP2_M, POP2_N = POP2_G.shape
POP2_RSIDs = ['RS01', 'RS02', 'RS04', 'RS05', 'RS07', 'RS08', 'RS09']
POP2_BPs = [10, 20, 40, 50, 70, 80, 90]
POP2_A1s = ['T', 'A', 'G', 'A', 'C', 'A', 'A']
POP2_A2s = ['C', 'G', 'A', 'G', 'T', 'G', 'G']


POP1_POP2_SNP_INTERSECTION = set(POP1_RSIDs).intersection(set(POP2_RSIDs))
POP1_POP2_INTERSECTION_LIST = sorted(list(POP1_POP2_SNP_INTERSECTION))
POP1_SNPS_VALID_FOR_POP2 = [rsid in POP1_POP2_SNP_INTERSECTION for rsid in POP1_RSIDs]
POP2_SNPS_VALID_FOR_POP1 = [rsid in POP1_POP2_SNP_INTERSECTION for rsid in POP2_RSIDs]
POP1_G_FILTERED_FOR_POP2 = POP1_G[POP1_SNPS_VALID_FOR_POP2, :]
POP2_G_FILTERED_FOR_POP1 = POP2_G[POP2_SNPS_VALID_FOR_POP1, :]

POP1_POP2_SIGN_DIFFS = [is_swap(POP1_A1s[POP1_RSIDs.index(rsid)],
                                     POP1_A2s[POP1_RSIDs.index(rsid)],
                                     POP2_A1s[POP2_RSIDs.index(rsid)],
                                     POP2_A2s[POP2_RSIDs.index(rsid)])
                        for rsid in POP1_POP2_INTERSECTION_LIST]


POP1_G_FILTERED_FOR_POP2_COV = np.cov(POP1_G_FILTERED_FOR_POP2)
POP2_G_FILTERED_FOR_POP1_COV = np.cov(POP2_G_FILTERED_FOR_POP1)

SWAP_FACTOR = np.where(POP1_POP2_SIGN_DIFFS, -1.0, 1.0)[:, np.newaxis]
POP1_POP2_COV_PRODUCT_UNSWAPPED = POP1_G_FILTERED_FOR_POP2_COV * POP2_G_FILTERED_FOR_POP1_COV
POP1_POP2_COV_PRODUCT = SWAP_FACTOR.T * POP1_POP2_COV_PRODUCT_UNSWAPPED * SWAP_FACTOR
POP1_POP2_COV_PRODUCT_DIAG = np.diag(POP1_POP2_COV_PRODUCT)


@pytest.fixture(scope='module')
def pop1_bedbimfam():
    return bedbimfam(
        bedbimfam_prefix=os.path.join(LD_SCORE_DATA_DIR, "POP1"),
        M=POP1_M,
        N=POP1_N,
        G=POP1_G,
        rsids=POP1_RSIDs,
        bps=POP1_BPs,
        a1s=POP1_A1s,
        a2s=POP1_A2s
    )

@pytest.fixture(scope='module')
def pop2_bedbimfam():
    return bedbimfam(
        bedbimfam_prefix=os.path.join(LD_SCORE_DATA_DIR, "POP2"),
        M=POP2_M,
        N=POP2_N,
        G=POP2_G,
        rsids=POP2_RSIDs,
        bps=POP2_BPs,
        a1s=POP2_A1s,
        a2s=POP2_A2s
    )


###########################################



class TestCalculateLdScores:


    #########
    # This is just a test that it returns the same as manually-calculated for these inputs
    def test__twopop__max_extents__return_expected(self, pop1_bedbimfam, pop2_bedbimfam, tmp_path):
        expected_ld_scores = POP1_POP2_COV_PRODUCT.sum(axis=0) * np.reciprocal(POP1_POP2_COV_PRODUCT_DIAG)

        pop1_r = os.path.join(tmp_path, "pop1_r.npy")
        pop2_r = os.path.join(tmp_path, "pop2_r.npy")

        pop1 = PopInfo('POP1', bedbimfam_prefix=pop1_bedbimfam, dist_col=BIM_BP_COL,
                               win_size=10**6, standardize=False, r_band_filename=pop1_r)
        pop2 = PopInfo('POP2', bedbimfam_prefix=pop2_bedbimfam, dist_col=BIM_BP_COL,
                               win_size=10**6, standardize=False, r_band_filename=pop2_r)

        pop1.calc_cross_pop_indices(pop2)
        actual_ld_scores = pop1.calc_ldscores(pop2).to_numpy()
        assert np.allclose(expected_ld_scores, actual_ld_scores)


    # This is just a test that it returns the same as manually-calculated for these inputs
    def test__twopop__min_extents__return_expected(self, pop1_bedbimfam, pop2_bedbimfam, tmp_path):
        pass
