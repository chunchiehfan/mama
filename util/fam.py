# TODO(jonbjala) Comment code / include function descriptions

# FAM file suffix
FAM_SUFFIX = ".fam"

# FAM columns
FAM_FID_COL = 'FID'
FAM_IID_COL = 'IID'
FAM_IIDF_COL = 'IIDF'
FAM_IIDM_COL = 'IIDM'
FAM_SEX_COL = 'SEX'
FAM_PHEN_COL = 'PHEN'
FAM_COL_TYPES = {
    FAM_FID_COL : str,
    FAM_IID_COL : str,
    FAM_IIDF_COL : str,
    FAM_IIDM_COL : str,
    FAM_SEX_COL : int,
    FAM_PHEN_COL : str
    }
FAM_COLUMNS = tuple(FAM_COL_TYPES.keys())


def get_sample_size_from_fam_file(fam_filename: str):
    with open(fam_filename) as f:
        N = sum(1 for line in f)

    return N
