rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species="drosophila_melanogaster",
        build="BDGP6.46",
        release="110",
    log:
        "logs/vep/cache.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v2.6.0/bio/vep/cache"