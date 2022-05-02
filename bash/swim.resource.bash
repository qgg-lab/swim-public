# ==============================
# = prepare resources for swim =
# ==============================

# 1. prepare file for illumina chips, this file does not contain duplicate in the first column
# but may in the second and third
# there are also duplicates between Illumina and Affy arrays with the same Affx ID, remove these as well
# ============================================================

cat PorcineSNPillumina_final.info PorcineSNPAffymetrix_final.info KPS_PorcineBreeding_52312.info | awk '{print $3"\t"$1"\t"$2"\t"$4"\t"$5}' | sort -k1,1 | /opt/software/bedtools2/bin/bedtools groupby -g 1 -c 2,3,4,5 -o distinct,distinct,distinct,distinct | awk '!($2 ~ /,/ || $3 ~ /,/ || $4 ~ /,/ || $5 ~ /,/)' > /var/data/swim/resource/chip.info
