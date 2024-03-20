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
        par_vcf = config["mom_phased_vcf"]
    output:
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["dad_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["dad_label"])
    params:
        par = "-O z  --allow-overlaps"
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
        id  = config["sample_id"],
        ref = config["reference_genome"],
        lab = config["mom_label"],
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["mom_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["mom_label"])
    output:
        fas = expand("{sample}_{label}_phased_genome.fasta", sample=config["sample_id"], label=config["mom_label"])
    params:
        par = "-H A"
        chr = config["phased_chroms"]
    log:
        "logs/consensus_mom_phased.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
        "samtools/1.12"
    shell:
        """
        samtools \
            faidx \
            {input.ref} \
            {params.chr} \
        | \
        bcftools \
            consensus \
            {params.par} \
            --sample {input.id} \
            {input.vcf} \
        | \
        sed 's/^>\(.*\)$/\>\1_{input.lab}/' \
        > {output.fas}
        """

rule consensus_dad:
    input:
        id  = config["sample_id"],
        ref = config["reference_genome"],
        lab = config["dad_label"],
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["dad_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["dad_label"])
    output:
        fas = expand("{sample}_{label}_phased_genome.fasta", sample=config["sample_id"], label=config["dad_label"])
    params:
        par = "-H A"
        chr = config["phased_chroms"]
    log:
        "logs/consensus_dad_phased.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
        "samtools/1.12"
    shell:
        """
        samtools \
            faidx \
            {input.ref} \
            {params.chr} \
        | \
        bcftools \
            consensus \
            {params.par} \
            --sample {input.id} \
            {input.vcf} \
        | \
        sed 's/^>\(.*\)$/\>\1_{input.lab}/' \
        > {output.fas}
        """

rule consensus_unphased:
    input:
        id  = config["sample_id"],
        ref = config["reference_genome"],
        lab = config["unphased_label"],
        vcf = expand("{sample}_{label}_phased_chroms.vcf.gz", sample=config["sample_id"], label=config["unphased_label"]),
        idx = expand("{sample}_{label}_phased_chroms.vcf.gz.csi", sample=config["sample_id"], label=config["unphased_label"])
    output:
        fas = expand("{sample}_unphased_genome.fasta", sample=config["sample_id"])
    params:
        par = "-H A"
        chr = config["unphased_chroms"]
    log:
        "logs/consensus_unphased.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
        "samtools/1.12"
    shell:
        """
        samtools \
            faidx \
            {input.ref} \
            {params.chr} \
        | \
        bcftools \
            consensus \
            {params.par} \
            --sample {input.id} \
            {input.vcf} \
        > {output.fas}
        """

rule combine_fasta:
    input:
        id  = config["sample_id"],
        mom = expand("{sample}_{label}_phased_genome.fasta", sample=config["sample_id"], label=config["mom_label"]),
        dad = expand("{sample}_{label}_phased_genome.fasta", sample=config["sample_id"], label=config["dad_label"]),
        unp = expand("{sample}_unphased_genome.fasta", sample=config["sample_id"])
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
            {input.unp} \
            > {output.fas}
        """

onsuccess:
    print("\n==== Workflow finished successfully! ====\n")
