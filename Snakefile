from Utils import *
configfile: "config.yaml"


rule all:
    input:
        "blast_db/ref.nhr",
        "blast_db/ref.nin",
        "blast_db/ref.nsq"
        
rule download_ref:
    params:
        ref = lambda wildcards: config["ref_sequences"][wildcards.ref]
    output:
        "Sequences/{ref}.fa"
    run:
        download_seqs(expand_range(params['ref']), output[0], email)
        
rule cat_seqs:
    input:
        expand("Sequences/{ref}.fa", ref=config["ref_sequences"])
    output:
        "Sequences/all.fa"
    shell:
        "cat {input} > {output}"

rule make_blast_db:
    input:
        "Sequences/all.fa"
    output:
        "blast_db/ref.nhr",
        "blast_db/ref.nin",
        "blast_db/ref.nsq"
    shell:
        "makeblastdb -dbtype nucl -in {input} -out blast_db/ref"
