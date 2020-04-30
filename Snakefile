from Utils import *
import os
configfile: "config.yaml"

email = os.environ.get('EMAIL')

queries = glob_wildcards("QuerySequences/{query}.fa")

rule all:
    input:
        expand("blast/{query}.bl", query=queries[0])

rule blast_queries:
    input:
        "QuerySequences/{query}.fa",
        "blast_db/ref.nhr",
        "blast_db/ref.nin",
        "blast_db/ref.nsq"
    output:
        "blast/{query}.bl"
    params:
        format="'6 qseqid qlen stitle slen pident length evalue bitscore'"
    shell:
        "blastn -outfmt {params.format} -query {input[0]} -db blast_db/ref -out {output}"

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
        "cat {input} | sed 's/18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 26S ribosomal RNA gene, partial sequence//' > {output}"

rule make_blast_db:
    input:
        "RefSequences/all.fa"
    output:
        "blast_db/ref.nhr",
        "blast_db/ref.nin",
        "blast_db/ref.nsq"
    shell:
        "makeblastdb -dbtype nucl -in {input} -out blast_db/ref"
