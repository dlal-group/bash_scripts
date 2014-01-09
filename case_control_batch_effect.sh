  --noweb
  --bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed
  --bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY
  --fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam
  --allow-no-sex
  --keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples.keeplist
  --pheno /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples_pheno.centre
  --maf 0.01
  --logistic
  --covar /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_JOINT_samples_pheno.pop
  --adjust
  --out chr4.txt

#analysis on cluster as phenotype and center as covariate - TWINS
mkdir -p logistic_middle_vs_low
chr=8

plink --noweb \
--bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed \
--bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY \
--fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam \
--allow-no-sex \
--maf 0.01 \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_TWINSUK_samples.keeplist \
--pheno /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/twins/PRUNED/MDS/UNFILTERED_1/CHR${chr}/TWINS_${chr}_cluster_middle_vs_low.pheno \
--covar /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_TWINSUK_samples_pheno.centre \
--logistic \
--adjust \
--out logistic_middle_vs_low/chr${chr}_cluster_center.txt

#analysis on cluster as phenotype and center as covariate - ALSPAC
mkdir -p logistic_middle_vs_low
chr=8

plink --noweb \
--bed /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bed \
--bim /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.bim.COPY \
--fam /lustre/scratch113/projects/uk10k/users/jh21/references_panel/uk10k/plink/chr${chr}.pruned.fam \
--allow-no-sex \
--maf 0.01 \
--keep /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_ALSPAC_samples.keeplist \
--pheno /lustre/scratch113/projects/uk10k/cohorts/REL-2012-06-02/v2/batch_effect/alspac/PRUNED/MDS/UNFILTERED_1/CHR${chr}/ALSPAC_${chr}_cluster_middle_vs_low.pheno \
--covar /nfs/users/nfs_m/mc14/Work/SANGER/UK10K/BATCHEFF/samples/3621_ALSPAC_samples_pheno.centre \
--logistic \
--adjust \
--out logistic_middle_vs_low/chr${chr}_cluster_center.txt

