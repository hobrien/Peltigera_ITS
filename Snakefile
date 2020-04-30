from Utils import *
import os
configfile: "config.yaml"

email = os.environ.get('EMAIL')

rule all:
    input:
        "blast_db/ref.nhr",
        "blast_db/ref.nin",
        "blast_db/ref.nsq"
        
rule download_ref:
    params:
        ref = lambda wildcards: config["ref_sequences"][wildcards.ref]
    output:
        "RefSequences/{ref}.fa"
    run:
        download_seqs(expand_range(params['ref']), output[0], email)
        
rule cat_seqs:
    input:
        expand("RefSequences/{ref}.fa", ref=config["ref_sequences"])
    output:
        "RefSequences/all.fa"
    shell:
        "cat {input} > {output}"

rule make_blast_db:
    input:
        "RefSequences/all.fa"
    output:
        "blast_db/ref.nhr",
        "blast_db/ref.nin",
        "blast_db/ref.nsq"
    shell:
        "makeblastdb -dbtype nucl -in {input} -out blast_db/ref"
