# Snakemake workflow to generate a diploid genome based on VCFs
# Juan Caballero
# (C) 2024

import yaml

# some subroutines needed
def dump_config_to_yaml(config):
    output_file = "run_config.yaml"
    with open(output_file, 'w') as configfile:
        yaml.dump(config, configfile)

# setting configurations
configfile: "config.yaml"

# before running the process
onstart:
    print("\n==== Variant calling pipeline starts ====")
    print("Configuration:")
    print(config)
    print("=" * 80)
    print()
    dump_config_to_yaml(config)

# main workflow
rule all:
    input:
        expand("{sample}_diploid_genome.fasta", sample=config["sample_id"])

rule merge_mom_phased:
    input:
        hom_vcf = config["hom_vcf"],
        par_vcf = config["mom_phased_vcf"]
    output:
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["mom_label"]), 
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["mom_label"])
    params:
        par = "-O z --allow-overlaps"
    log:
        "logs/merge_mom_phased.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            concat  \
            {params.par} \
            -o {output.vcf} \
            {input.hom_vcf} \
            {input.par_vcf} 

        bcftools index {output.vcf}
        """
        
rule merge_dad_phased:
    input:
        hom_vcf = config["hom_vcf"],
        par_vcf = config["dad_phased_vcf"]
    output:
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["dad_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["dad_label"])
    params:
        par = "-O z --allow-overlaps"
    log:
        "logs/merge_dad_phased.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            concat  \
            {params.par} \
            -o {output.vcf} \
            {input.hom_vcf} \
            {input.par_vcf} 

        bcftools index {output.vcf}
        """

rule consensus_mom:
    input:
        ref = config["reference_genome"],
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["mom_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["mom_label"])
    output:
        fas = expand("{sample}_{label}_phased_genome.fasta", sample=config["sample_id"], label=config["mom_label"])
    params:
        id = config["sample_id"],
        lab = config["mom_label"],
        par = "-H A",
        chr = config["phased_chrom"]
    log:
        "logs/consensus_mom_phased.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2",
        "samtools/1.12"
    shell:
        """
        samtools \
            faidx \
            {input.ref} \
            -r {params.chr} \
        | \
        bcftools \
            consensus \
            {params.par} \
            --sample {params.id} \
            {input.vcf} \
        | \
        perl -pe 's/^>(.*)/>$1_{params.lab}/' > {output.fas}
        """

rule consensus_dad:
    input:
        ref = config["reference_genome"],
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["dad_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["dad_label"])
    output:
        fas = expand("{sample}_{label}_phased_genome.fasta", sample=config["sample_id"], label=config["dad_label"])
    params:
        id = config["sample_id"],
        lab = config["dad_label"],
        par = "-H A",
        chr = config["phased_chrom"]
    log:
        "logs/consensus_dad_phased.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2",
        "samtools/1.12"
    shell:
        """
        samtools \
            faidx \
            {input.ref} \
            -r {params.chr} \
        | \
        bcftools \
            consensus \
            {params.par} \
            --sample {params.id} \
            {input.vcf} \
        | \
        perl -pe 's/^>(.*)/>$1_{params.lab}/' > {output.fas}
        """

rule consensus_unphased_X:
    input:
        ref = config["reference_genome"],
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["mom_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["mom_label"])
    output:
        fas = expand("{sample}_unphased_genome_X.fasta", sample=config["sample_id"])
    params:
        id  = config["sample_id"],
        lab = config["mom_label"],
        par = "-H A",
        chr = config["unphased_chrom_X"]
    log:
        "logs/consensus_unphased_X.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2",
        "samtools/1.12"
    shell:
        """
        samtools \
            faidx \
            {input.ref} \
            -r {params.chr} \
        | \
        bcftools \
            consensus \
            {params.par} \
            --sample {params.id} \
            {input.vcf} \
        | \
        perl -pe 's/^>(.*)/>$1_{params.lab}/' > {output.fas}
        """

rule consensus_unphased_Y:
    input:
        ref = config["reference_genome"],
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["dad_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["dad_label"])
    output:
        fas = expand("{sample}_unphased_genome_Y.fasta", sample=config["sample_id"])
    params:
        id  = config["sample_id"],
        lab = config["dad_label"],
        par = "-H A",
        chr = config["unphased_chrom_Y"]
    log:
        "logs/consensus_unphased_Y.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2",
        "samtools/1.12"
    shell:
        """
        samtools \
            faidx \
            {input.ref} \
            -r {params.chr} \
        | \
        bcftools \
            consensus \
            {params.par} \
            --sample {params.id} \
            {input.vcf} \
        | \
        perl -pe 's/^>(.*)/>$1_{params.lab}/' > {output.fas}
        """

rule combine_fasta:
    input:
        mom = expand("{sample}_{label}_phased_genome.fasta", sample=config["sample_id"], label=config["mom_label"]),
        dad = expand("{sample}_{label}_phased_genome.fasta", sample=config["sample_id"], label=config["dad_label"]),
        upx = expand("{sample}_unphased_genome_X.fasta", sample=config["sample_id"]),
        upy = expand("{sample}_unphased_genome_Y.fasta", sample=config["sample_id"])
    output:
        fas = expand("{sample}_diploid_genome.fasta", sample=config["sample_id"])
    threads: 1
    log:
        "logs/combine_fasta.log"
    shell:
        """
        cat \
            {input.mom} \
            {input.dad} \
            {input.upx} \
            {input.upy} > {output.fas}
        """

onsuccess:
    print("\n==== Workflow finished successfully! ====\n")
