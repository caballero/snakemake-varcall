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
        expand("{child}_heterozygous_{mother}_phased.vcf.gz", child=config["child_id"], mother=config["mother_id"]),
        expand("{child}_heterozygous_{father}_phased.vcf.gz", child=config["child_id"], mother=config["father_id"])
 
rule filter_het:
    input:
        child = config["child_bcf"]
    output:
        vcf = expand("{child}_heterozygous.vcf.gz", child=config["child_id"])
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
            -o {output.vcf} \
            {input.child}

        bcftools index {output}
        """

rule filter_hom:
    input:
        child = config["child_bcf"]
    output:
        vcf = expand("{child}_homozygous.vcf.gz", child=config["child_id"])
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
            -o {output.vcf} \
            {input.child}

        bcftools index {output}
        """

rule intersect_child_hom_mother:
    input:
        mother = config["mother_bcf"],
        child  = expand("{child}_homozygous.vcf.gz", child=config["child_id"])
    output:
        dirout = directory(expand("{child}_homo-vs-{mother}", child=config["child_id"], mothr=config["mother_id"]))
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
            -p {output.dirout}
        """

rule intersect_child_het_mother:
    input:
        mother = config["mother_bcf"],
        child  = expand("{child}_heterozygous.vcf.gz", child=config["child_id"])
    output:
        dirout = directory(expand("{child}_hetero-vs-{mother}", child=config["child_id"], mother=config["mother_id"]))
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
            -p {output.dirout}
        """

rule intersect_child_hom_father:
    input:
        father = config["father_bcf"],
        child  = expand("{child}_homozygous.vcf.gz", child=config["child_id"])
    output:
        dirout = directory(expand("{child}_homo-vs-{father}", child=config["child_id"], father=config["father_id"]))
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
            -p {output.dirout}
        """

rule intersect_child_het_father:
    input:
        mother = config["father_bcf"],
        child  = expand("{child}_heterozygous.vcf.gz", child=config["child_id"])
    output:
        dirout = directory(expand("{child}_hetero-vs-{father}", child=config["child_id"], father=config["father_id"]))
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
            -p {output.dirout}
        """

rule intersect_child_het_mother_not_father:
    input:
        child_and_mother = expand("{child}_hetero-vs-{mother}/0002.vcf.gz", child=config["child_id"], mother=config["mother_id"]),
        father_not_child = expand("{child}_hetero-vs-{father}/0000.vcf.gz", child=config["child_id"], father=config["father_id"])
    output:
        isec_dir = directory(expand("{child}_hetero-and-{mother}-and-not-{father}", child=config["child_id"], mother=config["mother_id"], father=config['father_id'])),
        vcf_out  = expand("{child}_heterozygous_{mother}_phased.vcf.gz", child=config["child_id"], mother=config["mother_id"])
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
        child_and_father = expand("{child}_hetero-vs-{father}/0002.vcf.gz", child=config["child_id"], father=config["father_id"]),
        mother_not_child = expand("{child}_hetero-vs-{mother}/0000.vcf.gz", child=config["child_id"], mother=config["mother_id"])
    output:
        isec_dir = directory(expand("{child}_hetero-and-{father}-and-not-{mother}", child=config["child_id"], mother=config["mother_id"], father=config["father_id"])),
        vcf_out  = expand("{child}_heterozygous_{father}_phased.vcf.gz", child=config["child_id"], father=config["father_id"])
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
