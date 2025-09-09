#!/bin/sh
# properties = {"type": "single", "rule": "GC_counts", "local": false, "input": ["/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/5_mkdup/TP19-M480_FOL6151A5_S12rg.mkdp.bam"], "output": ["results/GC_counts/sample_name_1.GC_counts.txt"], "wildcards": {"out_dir": "results", "samples": "sample_name_1"}, "params": {"sample_name": "sample_name_1", "mappable_regions_path": "/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/Ref/k100_minus_exclusion_lists.mappable_regions.hg38.bed", "ref_seq": "/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/Ref/hg38.fa", "chrom_sizes": "/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/Ref/hg38.standard.chrom.sizes", "map_q": 20, "size_range": "15 500", "CPU": 8, "griffin_GC_counts_script": "/ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/scripts//griffin_GC_counts.py"}, "log": [], "threads": 1, "resources": {}, "jobid": 1, "cluster": {}}
 cd /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction && \
/ihome/alee/rak373/.conda/envs/griffin_demo/bin/python \
-m snakemake results/GC_counts/sample_name_1.GC_counts.txt --snakefile /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/griffin_GC_and_mappability_correction.snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/.snakemake/tmp.10wmzxa_ /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/5_mkdup/TP19-M480_FOL6151A5_S12rg.mkdp.bam --latency-wait 60 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/config/cluster_slurm.yaml  --allowed-rules GC_counts --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/.snakemake/tmp.10wmzxa_/1.jobfinished || (touch /ix1/alee/LO_LAB/Personal/Daisong/Test_sample/Griffin/Griffin_updated/griffin_GC_and_mappability_correction/.snakemake/tmp.10wmzxa_/1.jobfailed; exit 1)

