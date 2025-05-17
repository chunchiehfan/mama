#!/usr/bin/env bash

# Merge the EAS and EUR genotype files into a single PLINK fileset.

plink \
  --bfile chr22_mind02_geno02_maf01_EAS \
  --bmerge chr22_mind02_geno02_maf01_EUR \
  --make-bed \
  --out chr22_mind02_geno02_maf01_EAS_EUR

