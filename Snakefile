configfile: "config.yaml"

rule all:
    input:
        os.path.join(config["output_dir_RELACS"],"Consensus_peaks/Consensus_H3K4me3/Consensus_H3K4me3.bed"),
        os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K4me3/Counts/H3K4me3_counts.mat.gz"),
        os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K4me1/Counts/H3K4me1_counts.mat.gz"),
        os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K27ac/Counts/H3K27ac_counts.mat.gz")



rule create_consensus_peaks:
    input:
        H3K4me3=expand(os.path.join(config["snakepipe_base_dir"],
        "MACS2/hippoN_{treatment}_H3K4me3_{replicate}.filtered.BAMPE_peaks.narrowPeak"),
        treatment=config["treatment"], replicate=config["replicate"]),
        H3K4me1=expand(os.path.join(config["snakepipe_base_dir"],
        "MACS2/hippoN_{treatment}_H3K4me1_{replicate}.filtered.BAMPE_peaks.broadPeak"),
        treatment=config["treatment"], replicate=config["replicate"]),
        H3K27ac=expand(os.path.join(config["snakepipe_base_dir"],
        "MACS2/hippoN_{treatment}_H3K27ac_{replicate}.filtered.BAMPE_peaks.narrowPeak"),
        treatment=config["treatment"], replicate=config["replicate"])
    output:
        H3K4me3=os.path.join(config["output_dir_RELACS"],"Consensus_peaks/Consensus_H3K4me3/Consensus_H3K4me3.bed"),
        H3K4me1=os.path.join(config["output_dir_RELACS"],"Consensus_peaks/Consensus_H3K4me1/Consensus_H3K4me1.bed"),
        H3K27ac=os.path.join(config["output_dir_RELACS"],"Consensus_peaks/Consensus_H3K27ac/Consensus_H3K27ac.bed")
    params:
        merge_distance_H3K4me3 = 500,
        merge_distance_H3K4me1 = 500,
        merge_distance_H3K27ac = 500
    shell:
        "cat {input.H3K4me3} | bedtools sort | bedtools merge -d {params.merge_distance_H3K4me3} > {output.H3K4me3} && \
        cat {input.H3K4me1} | bedtools sort | bedtools merge -d {params.merge_distance_H3K4me1} > {output.H3K4me1} && \
        cat {input.H3K27ac} | bedtools sort | bedtools merge -d {params.merge_distance_H3K27ac} > {output.H3K27ac}"


rule count_reads_multiBamSummary_H3K4me3:
    input:
        bams=expand(os.path.join(config["snakepipe_base_dir"],
        "filtered_bam/hippoN_{treatment}_H3K4me3_{replicate}.filtered.bam"),
        treatment=config["treatment"], replicate=config["replicate"]),
        bed_peaks=os.path.join(config["output_dir_RELACS"],"Consensus_peaks/Consensus_H3K4me3/Consensus_H3K4me3.bed"),
        blacklist = config["blacklist_ChIP-Seq"]
    output:
        matrix_binary=os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K4me3/Counts/H3K4me3_counts.mat.gz"),
        matrix_count=os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K4me3/Counts/H3K4me3_counts.counts")
    params:
        p=20,
        mapq=3
    shell:
        """ multiBamSummary BED-file -b {input.bams} \
            -o {output.matrix_binary} --BED {input.bed_peaks} -bl {input.blacklist} -p {params.p} \
            --outRawCounts {output.matrix_count} -e \
            --minMappingQuality {params.mapq} """


rule count_reads_multiBamSummary_H3K4me1:
    input:
        bams=expand(os.path.join(config["snakepipe_base_dir"],
        "filtered_bam/hippoN_{treatment}_H3K4me1_{replicate}.filtered.bam"),
        treatment=config["treatment"], replicate=config["replicate"]),
        bed_peaks=os.path.join(config["output_dir_RELACS"],"Consensus_peaks/Consensus_H3K4me1/Consensus_H3K4me1.bed"),
        blacklist = config["blacklist_ChIP-Seq"]
    output:
        matrix_binary=os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K4me1/Counts/H3K4me1_counts.mat.gz"),
        matrix_count=os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K4me1/Counts/H3K4me1_counts.counts")
    params:
        p=20,
        mapq=3
    shell:
        """ multiBamSummary BED-file -b {input.bams} \
            -o {output.matrix_binary} --BED {input.bed_peaks} -bl {input.blacklist} -p {params.p} \
            --outRawCounts {output.matrix_count} -e \
            --minMappingQuality {params.mapq} """


rule count_reads_multiBamSummary_H3K27ac:
    input:
        bams=expand(os.path.join(config["snakepipe_base_dir"],
        "filtered_bam/hippoN_{treatment}_H3K27ac_{replicate}.filtered.bam"),
        treatment=config["treatment"], replicate=config["replicate"]),
        bed_peaks=os.path.join(config["output_dir_RELACS"],"Consensus_peaks/Consensus_H3K27ac/Consensus_H3K27ac.bed"),
        blacklist = config["blacklist_ChIP-Seq"]
    output:
        matrix_binary=os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K27ac/Counts/H3K27ac_counts.mat.gz"),
        matrix_count=os.path.join(config["output_dir_RELACS"],"Differential_analysis/Differential_H3K27ac/Counts/H3K27ac_counts.counts")
    params:
        p=20,
        mapq=3
    shell:
        """ multiBamSummary BED-file -b {input.bams} \
            -o {output.matrix_binary} --BED {input.bed_peaks} -bl {input.blacklist} -p {params.p} \
            --outRawCounts {output.matrix_count} -e \
            --minMappingQuality {params.mapq} """



# rule differential_peak_analysis_H3K4me3:
#     input:
#         bams=expand(os.path.join(config["snakepipe_base_dir"],
#         "filtered_bam/hippoN_{treatment}_H3K4me3_{replicate}.filtered.bam"),
#         treatment=config["treatment"], replicate=config["replicate"]),
#         bed_peaks="Consensus_peaks/Consensus_H3K4me3/Consensus_H3K4me3.bed"
#     output:
#         directory("Differential_H3K4me3/")
#     shell:
#         """ bash ./scripts/pipe_DiffPeaks.sh H3K4me3 {input.bed_peaks} {output} {input.bams} """
