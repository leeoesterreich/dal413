#!/bin/sh
# properties = {"type": "single", "rule": "GC_bias", "local": false, "input": ["results/GC_counts/sample_name_1.GC_counts.txt"], "output": ["results/GC_bias/sample_name_1.GC_bias.txt", "results/GC_plots/sample_name_1.GC_bias.summary.pdf"], "wildcards": {"out_dir": "results", "samples": "sample_name_1"}, "params": {"sample_name": "sample_name_1", "mappable_name": "k100_minus_exclusion_lists.mappable_regions.hg38", "genome_GC_frequency": "/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/Ref/genome_GC_frequency/", "size_range": "15 500", "griffin_GC_bias_script": "/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/scripts//griffin_GC_bias.py"}, "log": [], "threads": 1, "resources": {}, "jobid": 2, "cluster": {}}
 cd /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction && \
/ihome/alee/rak373/.conda/envs/griffin_demo/bin/python \
-m snakemake results/GC_bias/sample_name_1.GC_bias.txt --snakefile /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/griffin_GC_and_mappability_correction.snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/.snakemake/tmp.f3m98s6t results/GC_counts/sample_name_1.GC_counts.txt --latency-wait 60 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/config/cluster_slurm.yaml  --allowed-rules GC_bias --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/.snakemake/tmp.f3m98s6t/2.jobfinished || (touch /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/.snakemake/tmp.f3m98s6t/2.jobfailed; exit 1)

