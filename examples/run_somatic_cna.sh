#! /bin/bash
# End-to-end somatic centromere CNA on the HG008 matched tumor/normal pair (GIAB).
# Companion to docs/somatic_cna_tutorial.md.
#
# Prerequisites:
#   pip install ascairn[s3,plot]       # boto3 (S3 BAM access) + matplotlib (PDF)

set -eu

# Step 0a: reference panel (ascairn_resource).
if [ ! -d ascairn_resource ]; then
    git clone https://github.com/friend1ws/ascairn_resource.git
fi
RESOURCE_DIR=ascairn_resource/resource/panel/ascairn_paper_2025
COMMON_DIR=ascairn_resource/resource/common
# The reference the BAMs are aligned to. Passing it lets samtools skip the slow
# MD5-based reference lookup when reading the S3 BAM (~5x faster extraction).
REF=GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

S3=s3://giab/data_somatic/HG008/Liss_lab/Element_AVITI_20240626
NORMAL_BAM=${S3}/HG008-N-D_Element-StdInsert_72x_GRCh38-GIABv3.bam
TUMOR_BAM=${S3}/HG008-T_Element-StdInsert_102x_GRCh38-GIABv3.bam

OUT=output
mkdir -p ${OUT}


# Step 0b: fetch the GRCh38-GIABv3 reference the HG008 BAMs are aligned to.
REF_URL=https://42basepairs.com/download/web/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
if [ ! -f ${REF} ]; then
    curl -L ${REF_URL} | gunzip > ${REF}
    samtools faidx ${REF}
fi


# Step 1: depth on the normal (drives coverage / sex for cen_type).
ascairn check_depth ${NORMAL_BAM} \
    -o ${OUT}/HG008-N.depth.txt \
    --baseline_region ${COMMON_DIR}/chr22_long_arm_hg38.bed \
    --x_region        ${COMMON_DIR}/chrX_short_arm_hg38.bed \
    -r ${REF} -t 4


# Step 2: rare k-mer counts for both normal and tumor.
for LABEL_BAM in "N:${NORMAL_BAM}" "T:${TUMOR_BAM}"; do
    LABEL=${LABEL_BAM%%:*}
    BAM=${LABEL_BAM#*:}
    ascairn kmer_count ${BAM} \
        -o ${OUT}/HG008-${LABEL}.kmer_count.txt \
        --kmer_file  ${RESOURCE_DIR}/rare_kmer_list.fa \
        --cen_region ${COMMON_DIR}/cen_region_curated_margin_hg38.bed \
        -r ${REF} -t 4
done


# Step 3+4: per chromosome, type the normal then call somatic CNA.
SEX=$(grep '^Sex:' ${OUT}/HG008-N.depth.txt | cut -d ' ' -f 2)

CHR_LIST="$(seq 1 22) X"
if [ "${SEX}" = "male" ]; then
    CHR_LIST="${CHR_LIST} Y"
fi

for CHR in ${CHR_LIST}; do

    cen_type_cmd=(
        ascairn cen_type
        ${OUT}/HG008-N.kmer_count.txt
        -o ${OUT}/HG008-N.chr${CHR}
        --kmer_info  ${RESOURCE_DIR}/kmer_info/chr${CHR}.kmer_info.txt.gz
        --hap_info   ${RESOURCE_DIR}/hap_info/chr${CHR}.hap_info.txt
        --depth_file ${OUT}/HG008-N.depth.txt
    )
    if { [ "${CHR}" = "X" ] && [ "${SEX}" = "male" ]; } || [ "${CHR}" = "Y" ]; then
        cen_type_cmd+=(--single_hap)
        "${cen_type_cmd[@]}"
        # single-haplotype chromosomes are not plotted by somatic_cna; skip.
        continue
    fi
    "${cen_type_cmd[@]}"

    # Normal marker_prob defines the proxy hap pair; normal/tumor counts are paired per marker.
    ascairn somatic_cna \
        --marker_prob  ${OUT}/HG008-N.chr${CHR}.haplotype.marker_prob.txt \
        --normal_count ${OUT}/HG008-N.kmer_count.txt \
        --tumor_count  ${OUT}/HG008-T.kmer_count.txt \
        --kmer_info    ${RESOURCE_DIR}/kmer_info/chr${CHR}.kmer_info.txt.gz \
        --hap_info     ${RESOURCE_DIR}/hap_info/chr${CHR}.hap_info.txt \
        -o ${OUT}/HG008.chr${CHR}
done
