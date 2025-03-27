#!/bin/bash

#SBATCH --partition batch
#SBATCH -w compute06
#SBATCH -n 8
##SBATCH --gres=gpu:v100:1
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --job-name=phylo-chikv
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=J.Juma@cgiar.org

# activate conda environment
hostname=$(hostname -s)
ostype=$(uname)

# local
if [ "${hostname}" == "jjuma" ] && [ "${ostype}" == 'Linux' ]; then
    echo -e "running pipeline locally: ${hostname}"
    source ${HOME}/miniconda3/etc/profile.d/conda.sh
    conda activate ${HOME}/miniconda3/envs/rvfv-phylo

elif [ "${hostname}" == "jjuma-2" ] && [ "${ostype}" == 'Darwin' ]; then
    echo -e "running pipeline locally: ${hostname}"
    source ${HOME}/anaconda3/etc/profile.d/conda.sh
    conda activate ${HOME}/anaconda3/envs/rvfv-phylo

elif [ "${hostname}" == "ilrike-23195" ] && [ "${ostype}" == 'Darwin' ]; then
    echo -e "running pipeline locally: ${hostname}"
    source ${HOME}/anaconda3/etc/profile.d/conda.sh
    conda activate ${HOME}/anaconda3/envs/rvfv-phylo
else
    echo -e "running pipeline on cluster: ${hostname}"
    # activate conda environment
    source ${HOME}/miniconda3/etc/profile.d/conda.sh
    conda activate ${HOME}/miniconda3/envs/rvfv-phylo
fi



# # M segment glycoprotein
# BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/M/complete/global"
# refs=("DQ380208" "DQ380213" "DQ380193")
# for ref in ${refs[*]}; do
#     nextflow run main.nf \
#         --subworkflow rvfvcirculatingstrains \
#         --fasta ${BASEDIR}/merged-sequences/RVFV-M.fasta \
#         --metadata ${BASEDIR}/merged-sequences/RVFV-M.csv \
#         --segment M \
#         --prefix riftM \
#         --start 20 \
#         --end 3611 \
#         --remove-duplicates false \
#         --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
#         --vaccine_reference ${ref} \
#         --outliers ./assets/RVFV-M-outliers-circulating.txt \
#         --recombinants ./assets/rvfv-M-potential-recombinants.txt \
#         --outdir ./output-dir/riftM \
#         -work-dir ./work/riftM \
#         -resume
# done

# # python bin/split_fasta.py \
# #    --fasta output-dir/riftM/netglyc/riftM.human.fasta \
# #    --num 85 --outDir output-dir/riftM/netglyc/

# NETGLYC_DIR="output-dir/riftM/netglyc"
# for f in ${NETGLYC_DIR}/*.netglyc.txt;
# do
# fname="$(basename $f)"
# prefix="${fname%.netglyc.txt}.parsed.netglyc"

# python bin/parseNetnglycOutput.py \
#     --input $f \
#     --prefix $prefix \
#     --outdir $NETGLYC_DIR
# done


# # S segment non-structural gene
# BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/S/complete/global"
# refs=("DQ380154" "DQ380182" "DQ380157")
# for ref in ${refs[*]}; do
#     nextflow run main.nf \
#         --subworkflow rvfvcirculatingstrains \
#         --fasta ${BASEDIR}/merged-sequences/RVFV-S.fasta \
#         --metadata ${BASEDIR}/merged-sequences/RVFV-S.csv \
#         --segment S-NSS \
#         --prefix riftS-NSS \
#         --start 34 \
#         --end 829 \
#         --remove-duplicates false \
#         --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
#         --vaccine_reference ${ref} \
#         --outliers ./assets/RVFV-S-NSS-outliers-circulating.txt \
#         --recombinants ./assets/rvfv-S-potential-recombinants.txt \
#         --outdir ./output-dir/riftS-NSS \
#         -work-dir ./work/riftS-NSS \
#         -resume
# done

# NETGLYC_DIR="output-dir/riftS-NSS/netglyc"
# for f in ${NETGLYC_DIR}/*.netglyc.txt;
# do
# fname="$(basename $f)"
# prefix="${fname%.netglyc.txt}.parsed.netglyc"

# python bin/parseNetnglycOutput.py \
#     --input $f \
#     --prefix $prefix \
#     --outdir $NETGLYC_DIR
# done

# # S segment nucleocapsid gene
# BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/S/complete/global"
# refs=("DQ380154" "DQ380182" "DQ380157")
# for ref in ${refs[*]}; do
#     nextflow run main.nf \
#         --subworkflow rvfvcirculatingstrains \
#         --fasta ${BASEDIR}/merged-sequences/RVFV-S.fasta \
#         --metadata ${BASEDIR}/merged-sequences/RVFV-S.csv \
#         --segment S-NP \
#         --prefix riftS-NP \
#         --start 925 \
#         --end 1660 \
#         --remove-duplicates false \
#         --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
#         --vaccine_reference ${ref} \
#         --outliers ./assets/RVFV-S-NSS-outliers-circulating.txt \
#         --recombinants ./assets/rvfv-S-potential-recombinants.txt \
#         --outdir ./output-dir/riftS-NP \
#         -work-dir ./work/riftS-NP \
#         -resume
# done

# NETGLYC_DIR="output-dir/riftS-NP/netglyc"
# for f in ${NETGLYC_DIR}/*.netglyc.txt;
# do
# fname="$(basename $f)"
# prefix="${fname%.netglyc.txt}.parsed.netglyc"

# python bin/parseNetnglycOutput.py \
#     --input $f \
#     --prefix $prefix \
#     --outdir $NETGLYC_DIR
# done


# # L segment RdRp
# BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/L/complete/global"
# refs=("DQ375404" "DQ375417" "DQ375430")

# for ref in ${refs[*]}; do
#     nextflow run main.nf \
#         --subworkflow rvfvcirculatingstrains \
#         --fasta ${BASEDIR}/merged-sequences/rvfv-L.fasta \
#         --metadata ${BASEDIR}/merged-sequences/rvfv-L.csv \
#         --segment L \
#         --prefix riftL \
#         --start 18 \
#         --end 6294 \
#         --remove-duplicates false \
#         --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
#         --vaccine_reference ${ref} \
#         --outliers ./assets/RVFV-L-outliers-circulating.txt \
#         --recombinants ./assets/rvfv-L-potential-recombinants.txt \
#         --outdir ./output-dir/riftL \
#         -work-dir ./work/riftL \
#         -resume
# done


# # python bin/split_fasta.py \
# #     --fasta output-dir/riftL/netglyc/riftL.human.fasta \
# #     --num 85 --outDir output-dir/riftL/netglyc/

# NETGLYC_DIR="output-dir/riftL/netglyc"
# for f in ${NETGLYC_DIR}/*.netglyc.txt;
# do
# fname="$(basename $f)"
# prefix="${fname%.netglyc.txt}.parsed.netglyc"

# python bin/parseNetnglycOutput.py \
#     --input $f \
#     --prefix $prefix \
#     --outdir $NETGLYC_DIR
# done

################################################################################
###############  Mutational profiles with reference genomes ####################
################################################################################

# L segment RdRp
BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/L/complete/global"

nextflow run main.nf \
    --subworkflow rvfvmutationalprofiling \
    --fasta ${BASEDIR}/merged-sequences/rvfv-L.fasta \
    --metadata ${BASEDIR}/merged-sequences/rvfv-L.csv \
    --segment L \
    --prefix riftL \
    --start 18 \
    --end 6294 \
    --remove-duplicates false \
    --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
    --vaccine_reference NC_014397 \
    --outliers ./assets/RVFV-L-outliers-circulating.txt \
    --recombinants ./assets/rvfv-L-potential-recombinants.txt \
    --outdir ./output-dir/mutational-profiling/riftL \
    -work-dir ./work/mutational-profiling/riftL \
    -resume

# both strands
python bin/getSequenceContexts.py \
    --positive-strand ./output-dir/mutational-profiling/riftL/strain-types/strain_type.NC_014397.sorted.alignment.fasta \
    --negative-strand ./output-dir/mutational-profiling/riftL/strain-types-rev/strain_type.NC_014397.sorted.alignment.fasta \
    --positive-strand-snps ./output-dir/mutational-profiling/riftL/strain-types/riftL.strain_type.NC_014397.mutations.per.strain.singleton.csv \
    --negative-strand-snps ./output-dir/mutational-profiling/riftL/strain-types-rev/riftL.strain_type.NC_014397.mutations.per.strain.singleton.csv \
    --num 2 \
    --column host \
    --prefix riftL \
    --outdir ./output-dir/mutational-profiling/riftL/sequence-contexts



# M segment RdRp
BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/M/complete/global"

nextflow run main.nf \
    --subworkflow rvfvmutationalprofiling \
    --fasta ${BASEDIR}/merged-sequences/rvfv-M.fasta \
    --metadata ${BASEDIR}/merged-sequences/rvfv-M.csv \
    --segment M \
    --prefix riftM \
    --start 20 \
    --end 3611 \
    --remove-duplicates false \
    --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
    --vaccine_reference NC_014396 \
    --outliers ./assets/RVFV-M-outliers-circulating.txt \
    --recombinants ./assets/rvfv-M-potential-recombinants.txt \
    --outdir ./output-dir/mutational-profiling/riftM \
    -work-dir ./work/mutational-profiling/riftM \
    -resume

# both strands
python bin/getSequenceContexts.py \
    --positive-strand ./output-dir/mutational-profiling/riftM/strain-types/strain_type.NC_014396.sorted.alignment.fasta \
    --negative-strand ./output-dir/mutational-profiling/riftM/strain-types-rev/strain_type.NC_014396.sorted.alignment.fasta \
    --positive-strand-snps ./output-dir/mutational-profiling/riftM/strain-types/riftM.strain_type.NC_014396.mutations.per.strain.singleton.csv \
    --negative-strand-snps ./output-dir/mutational-profiling/riftM/strain-types-rev/riftM.strain_type.NC_014396.mutations.per.strain.singleton.csv \
    --num 2 \
    --column host \
    --prefix riftM \
    --outdir ./output-dir/mutational-profiling/riftM/sequence-contexts

# NSS 

BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/S/complete/global"

nextflow run main.nf \
    --subworkflow rvfvmutationalprofiling \
    --fasta ${BASEDIR}/merged-sequences/rvfv-S.fasta \
    --metadata ${BASEDIR}/merged-sequences/rvfv-S.csv \
    --segment S-NSS \
    --prefix riftS-NSS \
    --start 34 \
    --end 829 \
    --remove-duplicates false \
    --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
    --vaccine_reference NC_014395 \
    --outliers ./assets/RVFV-S-NSS-outliers-circulating.txt  \
    --recombinants ./assets/rvfv-S-potential-recombinants.txt \
    --outdir ./output-dir/mutational-profiling/riftS-NSS \
    -work-dir ./work/mutational-profiling/riftS-NSS \
    -resume

# both strands
python bin/getSequenceContexts.py \
    --positive-strand ./output-dir/mutational-profiling/riftS-NSS/strain-types/strain_type.NC_014395.sorted.alignment.fasta \
    --negative-strand ./output-dir/mutational-profiling/riftS-NSS/strain-types-rev/strain_type.NC_014395.sorted.alignment.fasta \
    --positive-strand-snps ./output-dir/mutational-profiling/riftS-NSS/strain-types/riftS-NSS.strain_type.NC_014395.mutations.per.strain.singleton.csv \
    --negative-strand-snps ./output-dir/mutational-profiling/riftS-NSS/strain-types-rev/riftS-NSS.strain_type.NC_014395.mutations.per.strain.singleton.csv \
    --num 2 \
    --column host \
    --prefix riftS-NSS \
    --outdir ./output-dir/mutational-profiling/riftS-NSS/sequence-contexts



# NP

BASEDIR="${HOME}/projects/vaccine-and-circulating-strains-analysis/segments/S/complete/global"

nextflow run main.nf \
    --subworkflow rvfvmutationalprofiling \
    --fasta ${BASEDIR}/merged-sequences/rvfv-S.fasta \
    --metadata ${BASEDIR}/merged-sequences/rvfv-S.csv \
    --segment S-NP \
    --prefix riftS-NP \
    --start 925 \
    --end 1660 \
    --remove-duplicates false \
    --lineages ${BASEDIR}/assignment/output-dir/report/lineages.csv \
    --vaccine_reference NC_014395 \
    --outliers ./assets/RVFV-S-NSS-outliers-circulating.txt  \
    --recombinants ./assets/rvfv-S-potential-recombinants.txt \
    --outdir ./output-dir/mutational-profiling/riftS-NP \
    -work-dir ./work/mutational-profiling/riftS-NP \
    -resume

# both strands
python bin/getSequenceContexts.py \
    --positive-strand ./output-dir/mutational-profiling/riftS-NP/strain-types-rev/strain_type.NC_014395.sorted.alignment.fasta \
    --negative-strand ./output-dir/mutational-profiling/riftS-NP/strain-types/strain_type.NC_014395.sorted.alignment.fasta \
    --positive-strand-snps ./output-dir/mutational-profiling/riftS-NP/strain-types-rev/riftS-NP.strain_type.NC_014395.mutations.per.strain.singleton.csv \
    --negative-strand-snps ./output-dir/mutational-profiling/riftS-NP/strain-types/riftS-NP.strain_type.NC_014395.mutations.per.strain.singleton.csv \
    --num 2 \
    --column host \
    --prefix riftS-NP \
    --outdir ./output-dir/mutational-profiling/riftS-NP/sequence-contexts

# BASEDIR="${HOME}/projects/RVFV/continuous"

# nextflow run main.nf \
# --subworkflow rvfvphylocontinuous \
# --fasta ${BASEDIR}/segments/M/complete/global/merged-sequences/M-global.fasta \
# --metadata ${BASEDIR}/segments/M/complete/global/merged-sequences/M-global.csv \
# --segment M \
# --prefix riftM \
# --filter_column location \
# --start 20 \
# --end 3611 \
# --lineages ${BASEDIR}/segments/M/complete/global/assignment/output-dir/report/lineages.csv \
# --outliers ${BASEDIR}/segments/M/complete/RVFV-M-outliers-circulating.txt \
# --outdir ./output-dir/continuous/riftM \
# -work-dir ./work-dir/continuous/riftM \
# -resume


