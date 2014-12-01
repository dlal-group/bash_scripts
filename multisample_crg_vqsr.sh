#!/bin/bash

#modified script to concat together all the called regions and perform the vqsr filtering
OUTF=$1

#variation type
VARTYPE=$2

#check if we have all chr called:
chr_num=`ls $OUTF/*.multisampleinitial.allregions.${VARTYPE}.done | wc -l`

if [ ${chr_num} -eq "23" ]
then
  touch $OUTF/1.call.done
fi


if [ -f $OUTF/1.call.done ]
then

  #REF=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19/ucsc.hg19.fasta <- this file generate errors during the contig header check: mismatch of contig names
  REF=/nfs/users/GD/resource/human/hg19/hg19.fasta
  DBSNP=/users/GD/resource/human/hg19/databases/dbSNP/dbsnp_138.hg19.vcf
  #GATK=/users/GD/tools/GATK/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar
  GATK=/users/GD/tools/GATK/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar
  #PROBE=/users/GD/resource/human/probesets/nimblegene/v3/Target_Regions/SeqCap_EZ_Exome_v3_capture.bed_plus_50
  #PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/nimblegen_plus50.bed
  #PROBE=/nfs/users/xe/ggirotto/multisample/test_multisample_chr17_WES_WGS/nimblegen_plus50_chr${chr}_r${reg}.bed
  GATKRS=/users/GD/resource/human/hg19/databases/GATK_resources/bundle/2.8/hg19
  CPU=8

  #this has to be general, so we need to concat all chr together to have a full file with all the data to filter
  for c in {1..22} X
  do
  ls $OUTF/${c}.multisampleinitial.allregions.${VARTYPE}.vcf
  done > $OUTF/all_vcf.${VARTYPE}.list

  ##concat all files together
  #if [ ! -s $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf ]
  #then
    #echo "Already created concat file for ${VARTYPE}.."
  #else
    echo "Create concat file for ${VARTYPE}.."
    bcftools2 concat -f $OUTF/all_vcf.${VARTYPE}.list -O v -o $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf
    echo "..DONE!"
  #fi
if [ ! -s $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.tranches ]
then
  ## Variant Recalibration : we need to use different criteria for SNPs and INDELs
  case ${VARTYPE} in
    SNP)
java -jar $GATK \
-T VariantRecalibrator -R $REF -input $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf \
-recalFile $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.recal \
-tranchesFile $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.tranches \
-U LENIENT_VCF_PROCESSING --maxGaussians 6 \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATKRS/hapmap_3.3.hg19.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 $GATKRS/1000G_omni2.5.hg19.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=6.0 /nfs/users/xe/ggirotto/annotations/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz \
-resource:1000g,known=false,training=true,truth=false,prior=10.0 $GATKRS/1000G_phase1.snps.high_confidence.hg19.vcf \
-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff \
-mode ${VARTYPE} \
--target_titv 3.0 \
--TStranche 70.0 --TStranche 71.0 --TStranche 72.0 --TStranche 73.0 --TStranche 74.0 --TStranche 75.0 --TStranche 76.0 --TStranche 77.0 --TStranche 78.0 --TStranche 79.0 --TStranche 80.0 --TStranche 81.0 --TStranche 82.0 --TStranche 83.0 --TStranche 84.0 --TStranche 85.0 --TStranche 86.0 --TStranche 87.0 --TStranche 88.0 --TStranche 89.0 --TStranche 90.0 --TStranche 91.0 --TStranche 92.0 --TStranche 93.0 --TStranche 94.0 --TStranche 95.0 --TStranche 96.0 --TStranche 96.2 --TStranche 96.4 --TStranche 96.6 --TStranche 96.8 --TStranche 97.0 --TStranche 97.2 --TStranche 97.4 --TStranche 97.6 --TStranche 97.8 --TStranche 98.0 --TStranche 98.2 --TStranche 98.4 --TStranche 98.6 --TStranche 98.8 --TStranche 99.0 --TStranche 99.2 --TStranche 99.4 --TStranche 99.6 --TStranche 99.8 --TStranche 100.0 \
--rscript_file $OUTF/All.multisampleinitial.allregions.${VARTYPE}.r
# -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATKRS/dbsnp_138.hg19.excluding_sites_after_129.vcf \
    ;;
    INDEL)
java -jar $GATK \
-T VariantRecalibrator -R $REF -input $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf \
-recalFile $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.recal \
--tranchesFile $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.tranches \
-U LENIENT_VCF_PROCESSING --maxGaussians 6 \
-resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $GATKRS/Mills_and_1000G_gold_standard.indels.hg19.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /nfs/users/xe/ggirotto/annotations/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz \
-an QD -an FS -an ReadPosRankSum -an MQRankSum \
-mode ${VARTYPE} --target_titv 3.0 \
--TStranche 10 --TStranche 15 --TStranche 20 --TStranche 70 --TStranche 75 --TStranche 80 --TStranche 85 --TStranche 90 --TStranche 95 --TStranche 98.0 --TStranche 99.0 --TStranche 99.9 --TStranche 100 \
--rscript_file $OUTF/All.multisampleinitial.allregions.${VARTYPE}.r
# -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKRS/dbsnp_138.hg19.excluding_sites_after_129.vcf \
    ;;
  esac
fi
  if [ -s $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.tranches ]
  then
    echo "Recalibration of ${VARTYPE} done"
    ## Apply Recalibration
    echo "Apply Recalibration on ${VARTYPE} data.."

    java -jar $GATK -T ApplyRecalibration -R $REF -input $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf -tranchesFile $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.tranches -recalFile $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.recal -o $OUTF/All.multisampleinitial.allregions.${VARTYPE}.recalibrated.filtered.vcf -ts_filter_level 99.0 -mode ${VARTYPE}
    ## grep PASS snps from the recalibration
    egrep 'PASS|^#' $OUTF/All.multisampleinitial.allregions.${VARTYPE}.recalibrated.filtered.vcf > $OUTF/All.multisampleinitial.allregions.${VARTYPE}.recalibrated.filtered.clean.vcf
    
    if [ -s $OUTF/All.multisampleinitial.allregions.${VARTYPE}.recalibrated.filtered.clean.vcf ]
    then
      touch $OUTF/2.vqsr.apply.${VARTYPE}.done

      ## variant Evaluation
      java -jar $GATK -T VariantEval -R $REF --dbsnp $DBSNP -o $OUTF/report.All.multisampleinitial.allregions.${VARTYPE}.recalibrated.filtered.clean.vcf --eval $OUTF/All.multisampleinitial.allregions.${VARTYPE}.recalibrated.filtered.clean.vcf -l INFO --downsampling_type NONE

    	if [ -s $OUTF/report.All.multisampleinitial.allregions.${VARTYPE}.recalibrated.filtered.clean.vcf ]
    	then
    	touch $OUTF/3.report.${VARTYPE}.done
    	fi
    fi
  else
    echo "Recalibration of ${VARTYPE} not performed!!!!Why??"
  fi
fi

