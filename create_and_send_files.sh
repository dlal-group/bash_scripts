#!/bin/bash 
#ottimizzato per velocita' e spazio
#un solo ciclo che fa' tutto
#a patto che i file ci siano....
#copyleft 2008 "SistemistaMannaro©"  for "Tonioloids©"  

#DICHIARAZIONE VARS=PATHS
#FILE_DAT="TEST2.dat"
#FILE_PED="TEST2.ped"
#FILE_TRAIT="traits.txt"

##FUNZIONE BUILD SKELETON
function build_template() {                                
cat << EOF        
#PBS -S /bin/bash 
#PBS -V
#PBS -q toniolo
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -o OUT
PBS_O_INITDIR=\$PBS_O_WORKDIR
cd \$PBS_O_WORKDIR
/usr/local/cluster/toniolo/bin/ghost -d $FILE_DAT -p $FILE_PED --poly --bivariate --balance --trait $FIRST,$SECOND --covariate sex,age --norm --prefix $FIRST_$SECOND > $FIRST_$SECOND.ghout.txt
EOF
}

###CONTROLLO PARAMETRI INGRESSO SE OK COSTRUISCI E LANCIA
if [ -e $FILE_TRAIT ] && [ -e $FILE_DAT ] && [ -e $FILE_PED ] 
then
while read line 
do
        FIRST=${line% *}
        SECOND=${line#* }
        #togli i 3 commenti per creare e lanciare
        build_template > $FIRST_$SECOND.torun
        #qsub $FIRST_$SECOND.torun
        #rm $FIRST_$SECOND.torun &> /dev/null
done # < $FILE_TRAIT
else
  echo "Controlla i parametri, hai fornito:
          FILE_DAT=$FILE_DAT
        FILE_PED=$FILE_PED
        FILE_TRAIT=$FILE_TRAIT"
        exit
fi
#LA FINE (registered trademark)
