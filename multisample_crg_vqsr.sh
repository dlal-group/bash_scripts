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
  done > $OUTF/all_vcf.list

  ##concat all files together
  vcf-concat -f $OUTF/all_vcf.list > $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf
  echo "Created concat file for ${VARTYPE}.."

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
      -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATKRS/dbsnp_138.hg19.excluding_sites_after_129.vcf \
      -resource:1000g,known=false,training=true,truth=false,prior=10.0 $GATKRS/1000G_phase1.snps.high_confidence.hg19.vcf \
      -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff \
      -mode ${VARTYPE} \
      --target_titv 3.0 \
      -tranche 70.0 -tranche 71.0 -tranche 72.0 -tranche 73.0 -tranche 74.0 -tranche 75.0 -tranche 76.0 -tranche 77.0 -tranche 78.0 -tranche 79.0 -tranche 80.0 -tranche 81.0 -tranche 82.0 -tranche 83.0 -tranche 84.0 -tranche 85.0 -tranche 86.0 -tranche 87.0 -tranche 88.0 -tranche 89.0 -tranche 90.0 -tranche 91.0 -tranche 92.0 -tranche 93.0 -tranche 94.0 -tranche 95.0 -tranche 96.0 -tranche 96.2 -tranche 96.4 -tranche 96.6 -tranche 96.8 -tranche 97.0 -tranche 97.2 -tranche 97.4 -tranche 97.6 -tranche 97.8 -tranche 98.0 -tranche 98.2 -tranche 98.4 -tranche 98.6 -tranche 98.8 -tranche 99.0 -tranche 99.2 -tranche 99.4 -tranche 99.6 -tranche 99.8 -tranche 100.0
    ;;
    INDEL)
      java -jar $GATK \
      -T VariantRecalibrator -R $REF -input $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf \
      -recalFile $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.recal \
      -tranchesFile $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.tranches \
      -U LENIENT_VCF_PROCESSING --maxGaussians 6 \
      -resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $GATKRS/Mills_and_1000G_gold_standard.indels.hg19.vcf \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATKRS/dbsnp_138.hg19.excluding_sites_after_129.vcf \
      -an DP -an FS -an ReadPosRankSum -an MQRankSum \
      -mode ${VARTYPE} --target_titv 3.0 \
      -tranche 10 -tranche 15 -tranche 20 -tranche 25 -tranche 30 -tranche 35 -tranche 40 -tranche 45 -tranche 50 -tranche 55 -tranche 60 -tranche 65 -tranche 70 -tranche 75 -tranche 80 -tranche 85 -tranche 90 -tranche 95 -tranche 98.0 -tranche 99.0 --tranche 99.9 -tranche 100
    ;;
  esac

  if [ ! -s $OUTF/All.multisampleinitial.allregions.${VARTYPE}.vcf.tranches ]
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

