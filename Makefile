# Top-level shortcuts
.PHONY: t_qc t_norm t_filter p_qc p_norm p_filter surv uni penal

t_qc:
\tRscript 01_transcriptomics/pipeline/01_qc.R
t_norm:
\tRscript 01_transcriptomics/pipeline/02_norm_batch.R
t_filter:
\tRscript 01_transcriptomics/pipeline/03_feature_filter.R

p_qc:
\tRscript 02_proteomics/pipeline/01_qc.R
p_norm:
\tRscript 02_proteomics/pipeline/02_norm_batch.R
p_filter:
\tRscript 02_proteomics/pipeline/03_feature_filter.R

uni:
\tRscript 03_survival/scripts/01_univariate_cox.R
penal:
\tRscript 03_survival/scripts/02_penalized_cox.R
