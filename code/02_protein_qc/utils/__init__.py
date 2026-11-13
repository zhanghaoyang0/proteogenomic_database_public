"""
Utils package for data processing
"""
from .process_data import (
    norm_id,
    qc_olink,
    format_proteins_for_plink,
    format_covariates_for_plink,
    map_covar_to_genotype_iid,
)

__all__ = [
    'norm_id',
    'qc_olink',
    'format_proteins_for_plink',
    'format_covariates_for_plink',
    'map_covar_to_genotype_iid',
]

