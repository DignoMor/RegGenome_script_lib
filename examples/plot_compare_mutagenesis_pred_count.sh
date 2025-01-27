
INPATH=sample_data/sample_exogeneous_sequence
OPATH=wdir/plot_compare_mutagenesis_pred_count

mkdir -p ${OPATH}

scripts/exogeneous_tool.py compare_mutagenesis \
    --fasta ${INPATH}/sample.exogeneous.chr2.fa \
    --region_file_path ${INPATH}/sample.exogeneous.chr2.predicted_peaks.bed \
    --region_file_type bed3 \
    --regex_ref_seqs "(.*)_ref" \
    --regex_mut_seqs "(.*)_[0-9]*:[ATCG]2[ATCG]" \
    --pl_track_npy ${INPATH}/sample.exogeneous.chr2.pred_profiles_pl.npy \
    --mn_track_npy ${INPATH}/sample.exogeneous.chr2.pred_profiles_mn.npy \
    --opath ${OPATH} \
    --total_count_plot_path ${OPATH}/total_count.png \