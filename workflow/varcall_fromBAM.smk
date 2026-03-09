# Snakemake workflow for variant calling
# Juan Caballero
# (C) 2024

import yaml

#wildcard_constraints:
#    sample=r".+\d+"

threads: 40

# some subroutines needed
def dump_config_to_yaml(config):
    output_file = "run_config.yaml"
    with open(output_file, 'w') as configfile:
        yaml.dump(config, configfile)

# setting configurations
configfile: "config.yaml"
bam_dir = str(config['bam_dir'])

# before running the process
onstart:
    print("\n==== Variant calling pipeline starts ====")
    print("Configuration:")
    print(config)
    print("=" * 80)
    print()
#    dump_config_to_yaml(config)

# main workflow
rule all:
    input:
        expand("03_calls/{sample}_filt.bcf.csi", sample=config["samples"]) 


rule bcftools_mpileup:
    input:
        ref = config["genome"],
        bam = "{sample}.bam",
        bai = "{sample}.bam.bai"
    output:
        bcf_pil = "{sample}_pileup.bcf"
    params:
        par = config["bcf_mpileup_param"]
    log:
        "logs/bcftoos_mpileup/{sample}.log"
    threads: 10
    #conda:
    #    "envs/environment_varcall.yaml"
    envmodules:
        "bcftools/1.21"
    shell:
        """
        bcftools \
            mpileup \
            -f {input.ref} \
            {params.par} \
            --threads {threads} \
            -o {output.bcf_pil} \
            {input.bam} 
        """


rule bcftools_call:
    input:
        bcf_pal = "{sample}_pileup.bcf"
    output:
        bcf_raw = "{sample}_raw.bcf"
    params:
        par = config["bcf_call_param"]
    log:
        "logs/bcftoos_call/{sample}.log"
    threads: 10
    #conda:
    #    "envs/environment_varcall.yaml"
    envmodules:
        "bcftools/1.21"
    shell:
        """
        bcftools \
            call \
            {params.par} \
            --threads {threads} \
            -o {output} \
            {input} 
        """


rule bcftools_norm:
    input:
        ref = config["genome"],
        bcf_raw = "{sample}_raw.bcf"
    output:
        bcf_nor = "{sample}_norm.bcf"
    params:
        par = config["bcf_norm_param"]
    log:
        "logs/bcftoos_norm/{sample}.log"
    threads: 10
    #conda:
    #    "envs/environment_varcall.yaml"
    envmodules:
        "bcftools/1.21"
    shell:
        """
        bcftools \
            norm \
            -f {input.ref} \
            {params.par} \
            --threads {threads} \
            -o {output.bcf_nor} \
            {input.bcf_raw} 
        """


rule bcftools_filter:
    input:
        bcf_nor = "{sample}_norm.bcf"
    output:
        bcf_fil = "{sample}_filt.bcf"
    params:
        par = config["bcf_filter_param"]
    log:
        "logs/bcftoos_filter/{sample}.log"
    threads: 10
    #conda:
    #    "envs/environment_varcall.yaml"
    envmodules:
        "bcftools/1.21"
    shell:
        """
        bcftools \
            filter \
            {params.par} \
            --threads {threads} \
            -o {output.bcf_fil} \
            {input.bcf_nor} 
        """


rule bcftools_index:
    input:
        bcf_fil = "{sample}_filt.bcf"
    output:
        bcf_idx = "{sample}_filt.bcf.csi"
    log:
        "logs/bcftoos_index/{sample}.log"
    threads: 1
    #conda:
    #    "envs/environment_varcall.yaml"
    envmodules:
        "bcftools/1.21"
    shell:
        """
        bcftools \
            index \
            --threads {threads} \
            {input.bcf_fil} 
        """


onsuccess:
    print("\n==== Workflow finished successfully! ====\n")
