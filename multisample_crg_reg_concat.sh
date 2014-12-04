#!/bin/bash

#concat step: concat all regions together
filepath=$1
TYPE=$2

ls ${filepath}.*.multisampleinitial.allregions.${TYPE}.vcf > ${filepath}.all_reg.list

echo "Create concat file for ${filepath} of ${TYPE}.."
bcftools2 concat -f ${filepath}.all_reg.list -O v -o ${filepath}.multisampleinitial.allregions.${TYPE}.vcf
echo "..DONE!"