from Utils import *
import os
configfile: "config.yaml"

email = os.environ.get('EMAIL')

queries = glob_wildcards("QuerySequences/{query}.fa")

rule all:
    input:
        expand("blast/{query}.bl", query=queries[0]),
        expand("SpeciesComplexes/{complex}/{complex}_ref.fa", complex=config['complexes']),

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
        "cat {input} | perl -pe 's/((internal)|(18S)|(5.8S)).*//' > {output}"

rule make_blast_db:
    input:
        "RefSequences/all.fa"
    output:
        "blast_db/ref.nhr",
        "blast_db/ref.nin",
        "blast_db/ref.nsq"
    shell:
        "makeblastdb -dbtype nucl -in {input} -out blast_db/ref"
        
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

rule separate_seqs:
    input:
        rules.cat_seqs.output
    output:
        expand("SpeciesComplexes/{complex}/{complex}_ref.fa", complex=config['complexes'])
    run:
        complexes = invert_list_dict(config['complexes'])
        separated_seqs = separate_seqs(input[0], complexes)
        for sp_complex in separated_seqs:
            outfilename = os.path.join("SpeciesComplexes", sp_complex, "%s_ref.fa" % sp_complex)
            SeqIO.write(separated_seqs[sp_complex], outfilename, "fasta")
