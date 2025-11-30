import yaml
from pathlib import Path

CONFIG = yaml.safe_load(open("analysis/00_config/cancers.yaml"))
PATHS  = yaml.safe_load(open("analysis/00_config/paths.yaml"))
CANCERS = CONFIG["cancers"]
DATATYPES = ["gene", "iso", "iso_frac"]

def t_path(kind, cancer):
    return PATHS["transcriptomics"][kind].format(cancer=cancer)

rule all:
    input:
        expand("analysis/03_survival/01_univariate/out/{cancer}.transcriptomics.{dtype}.univariate.csv",
               cancer=CANCERS, dtype=DATATYPES)

# ---------- QC ----------
rule qc_transcriptomics:
    input:
        expr=lambda wc: t_path(wc.dtype, wc.cancer),
        clinical=lambda wc: t_path("clinical", wc.cancer),
        thresholds="analysis/00_config/qc_thresholds.yaml"
    output:
        out_csv="analysis/01_transcriptomics/01_qc/out/TCGA_{cancer}_{dtype}.qc_passed.csv"
    params:
        script="analysis/01_transcriptomics/scripts/01_qc.R"
    shell:
        "Rscript {params.script} {input.expr} {input.clinical} {input.thresholds} {output.out_csv}"

# ---------- Normalize + batch ----------
rule norm_transcriptomics:
    input:
        qc_csv="analysis/01_transcriptomics/01_qc/out/TCGA_{cancer}_{dtype}.qc_passed.csv",
        clinical=lambda wc: t_path("clinical", wc.cancer),
        batches="analysis/00_config/batches.yaml",
        covars="analysis/00_config/covariates.yaml"
    output:
        out_csv="analysis/01_transcriptomics/02_norm_batch/out/TCGA_{cancer}_{dtype}.normalized.csv"
    params:
        script="analysis/01_transcriptomics/scripts/02_norm_batch.R"
    shell:
        "Rscript {params.script} {input.qc_csv} {input.clinical} {input.batches} {input.covars} {output.out_csv}"

# ---------- (Simple) pass-through filter to survival input ----------
# If you later add a real filter script, swap the shell to call 03_filter.R
rule prepare_for_survival:
    input:
        norm_csv="analysis/01_transcriptomics/02_norm_batch/out/TCGA_{cancer}_{dtype}.normalized.csv",
        clinical=lambda wc: t_path("clinical", wc.cancer)
    output:
        expr="analysis/01_transcriptomics/03_feature_filter/out/TCGA_{cancer}_{dtype}.expr_normalized.csv",
        manifest="analysis/01_transcriptomics/03_feature_filter/out/TCGA_{cancer}_{dtype}.sample_manifest.csv"
    run:
        import pandas as pd
        df = pd.read_csv(input.norm_csv)
        df.to_csv(output.expr, index=False)
        clin = pd.read_csv(input.clinical)
        clin.to_csv(output.manifest, index=False)

# ---------- Univariate Cox ----------
rule univariate_cox:
    input:
        expr="analysis/01_transcriptomics/03_feature_filter/out/TCGA_{cancer}_{dtype}.expr_normalized.csv",
        clinical=lambda wc: t_path("clinical", wc.cancer),
        covars="analysis/00_config/covariates.yaml"
    output:
        out_csv="analysis/03_survival/01_univariate/out/{cancer}.transcriptomics.{dtype}.univariate.csv"
    params:
        script="analysis/03_survival/scripts/01_univariate_cox.R"
    shell:
        "Rscript {params.script} {input.expr} {input.clinical} {input.covars} {output.out_csv}"
