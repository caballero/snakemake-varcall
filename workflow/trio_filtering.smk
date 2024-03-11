# Snakemake workflow for filtering variants in trios
# Juan Caballero
# (C) 2024

import os
import yaml

wildcard_constraints:
    sample=r".+\d+"

# some subroutines needed
def dump_config_to_yaml(config, output_dir):
    output_file = os.path.join(output_dir, "run_config.yaml")
    with open(output_file, 'w') as configfile:
        yaml.dump(config, configfile)

# setting configurations
configfile: "config.yaml"
outdir = workflow.overwrite_workdir #get outdir from -d option of snakemake

# before running the process
onstart:
    print("\n==== Variant calling pipeline starts ====")
    print("Configuration:")
    print(config)
    print("=" * 80)
    print()
    dump_config_to_yaml(config, '.')

# main workflow
rule all:
    input:
        config["child_id"] + "_heterozygous_" + config["mother"] + "_phased.vcf.gz",
        config["child_id"] + "_heterozygous_" + config["father"] + "_phased.vcf.gz"

rule filter_het:
    input:
        config["child_bcf"]
    output:
        config["child_id"] + "_heterozygous.vcf.gz"
    params:
        par = "-i 'GT=\"0/1\"' -O z"
    log:
        "logs/filter_het.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            filter \
            {params.par} \
            --threads {threads} \
            -o {output} \
            {input}

        bcftools index {output}
        """

rule filter_hom:
    input:
        config["child_bcf"]
    output:
        vcf_hom = config["child_id"] + "_homozygous.vcf.gz"
    params:
        par = "-i 'GT=\"1/1\"' -O z"
    log:
        "logs/filter_hom.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            filter \
            {params.par} \
            --threads {threads} \
            -o {output} \
            {input}

        bcftools index {output}
        """

rule intersect_child_hom_mother:
    input:
        mother = config["mother_bcf"],
        child  = config["child_id"] + "_homozygous.vcf.gz"
    output:
        config["child_id"] + "_homo-vs-" + config["mother"]
    params:
        par = "-O z"
    log:
        "logs/isec_child_hom_mother.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            isec \
            {input.child} \
            {input.mother} \
            {params.par} \
            --threads {threads} \
            -p {output}
        """

rule intersect_child_het_mother:
    input:
        mother = config["mother_bcf"],
        child  = config["child_id"] + "_heterozygous.vcf.gz"
    output:
        config["child_id"] + "_hetero-vs-" + config["mother"]
    params:
        par = "-O z"
    log:
        "logs/isec_child_het_mother.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            isec \
            {input.child} \
            {input.mother} \
            {params.par} \
            --threads {threads} \
            -p {output}
        """

rule intersect_child_hom_father:
    input:
        father = config["father_bcf"],
        child  = config["child_id"] + "_homozygous.vcf.gz"
    output:
        config["child_id"] + "_homo-vs-" + config["father"]
    params:
        par = "-O z"
    log:
        "logs/isec_child_hom_father.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            isec \
            {input.child} \
            {input.father} \
            {params.par} \
            --threads {threads} \
            -p {output}
        """

rule intersect_child_het_father:
    input:
        mother = config["father_bcf"],
        child  = config["child_id"] + "_heterozygous.vcf.gz"
    output:
        config["child_id"] + "_hetero-vs-" + config["father"]
    params:
        par = "-O z"
    log:
        "logs/isec_child_het_father.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            isec \
            {input.child} \
            {input.father} \
            {params.par} \
            --threads {threads} \
            -p {output}
        """

rule intersect_child_het_mother_not_father:
    input:
        child_and_mother = config["child_id"] + "_hetero-vs-" + config["mother"] + "/0002.vcf.gz",
        father_not_child = config["child_id"] + "_hetero-vs-" + config["father"] + "/0000.vcf.gz"
    output:
        isec_dir = config["child_id"] + "_hetero-and-" + config["mother"] + "-and-not-" + config['father'],
        vcf_out  = config["child_id"] + "_heterozygous_" + config["mother"] + "_phased.vcf.gz"
    params:
        par = "-O z"
    log:
        "logs/isec_child_het_and_mother_not_father.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            isec \
            {input.child_and_mother} \
            {input.father_not_child} \
            {params.par} \
            --threads {threads} \
            -p {output.isec_dir}

        bcftools \
            view \
            {output.isec_dir}/0002.vcf.gz \
            -O z \
            -o {output.vcf_out}

        """

rule intersect_child_het_father_not_mother:
    input:
        child_and_father = config["child_id"] + "_hetero-vs-" + config["father"] + "/0002.vcf.gz",
        mother_not_child = config["child_id"] + "_hetero-vs-" + config["mother"] + "/0000.vcf.gz"
    output:
        isec_dir = config["child_id"] + "_hetero-and-" + config["father"] + "-and-not-" + config['mother'],
        vcf_out  = config["child_id"] + "_heterozygous_" + config["father"] + "_phased.vcf.gz"
    params:
        par = "-O z"
    log:
        "logs/isec_child_het_and_father_not_mother.log"
    threads: 1
    envmodules:
        "bcftools/1.10.2"
    shell:
        """
        bcftools \
            isec \
            {input.child_and_father} \
            {input.mother_not_child} \
            {params.par} \
            --threads {threads} \
            -p {output.isec_dir}

        bcftools \
            view \
            {output.isec_dir}/0002.vcf.gz \
            -O z \
            -o {output.vcf_out}

        """

onsuccess:
    print("\n==== Workflow finished successfully! ====\n")
