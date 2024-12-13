
INPATH=sample_data/sample_exogeneous_sequence
OPATH=wdir/plot_exogeneous_metaplot

mkdir -p ${OPATH}

scripts/exogeneous_tool.py metaplot \
    --inpath ${INPATH}/sample.exogeneous.chr2.fa \
    --region_path ${INPATH}/sample.exogeneous.chr2.predicted_peaks.bed \
    --outpath ${OPATH}/sample.png \
    --signal_tracks_pl ${INPATH}/sample.exogeneous.chr2.pred_profiles_pl.npy \
    --signal_tracks_mn ${INPATH}/sample.exogeneous.chr2.pred_profiles_mn.npy
