from Utils import *
configfile: "config.yaml"


rule all:
    input:
        expand("Sequences/{ref}.fa", ref=config["ref_sequences"]),
        
rule download_ref:
    params:
        ref = lambda wildcards: config["ref_sequences"][wildcards.ref]
    output:
        "Sequences/{ref}.fa"
    run:
        download_seqs(expand_range(params['ref']), output[0], email)