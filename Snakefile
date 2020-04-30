from Utils import *
import os
configfile: "config.yaml"

email = os.environ.get('EMAIL')

queries = glob_wildcards("QuerySequences/{query}.fa")

rule all:
    input:
        expand("SpeciesComplexes/{complex}/{complex}.phy", complex=config['complexes'])

rule download_ref:
    params:
        ref = lambda wildcards: config["ref_sequences"][wildcards.ref]
    output:
        "RefSequences/{ref}.fa"
    run:
        download_seqs(expand_range(params['ref']), output[0], email)
        
rule cat_ref_seqs:
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
        format = "'6 qseqid qlen stitle qstart qend sstart send slen pident length evalue bitscore'",
        eval = '1e-150'
    shell:
        "blastn -outfmt {params.format} -query {input[0]} -db blast_db/ref -evalue {params.eval} -out {output}"

rule separate_ref_seqs:
    input:
        rules.cat_ref_seqs.output
    output:
        expand("SpeciesComplexes/{complex}/{complex}_ref.fa", complex=config['complexes'])
    run:
        complexes = invert_list_dict(config['complexes'])
        separated_seqs = separate_seqs(input, complexes)
        for sp_complex in separated_seqs:
            outfilename = os.path.join("SpeciesComplexes", sp_complex, "%s_ref.fa" % sp_complex)
            SeqIO.write(separated_seqs[sp_complex], outfilename, "fasta")
            
rule separate_query_seqs:
    input:
        blast = expand("blast/{query}.bl", query=queries[0]),
        seqs = expand("QuerySequences/{query}.fa", query=queries[0])
    output:
        expand("SpeciesComplexes/{complex}/{complex}_queries.fa", complex=config['complexes'])
    run:
        complexes = make_lookup_table(input['blast'], invert_list_dict(config['complexes']))
        separated_seqs = separate_seqs(input['seqs'], complexes, 'description')
        for sp_complex in config['complexes']:
            outfilename = os.path.join("SpeciesComplexes", sp_complex, "%s_queries.fa" % sp_complex)
            if sp_complex in separated_seqs:
                SeqIO.write(separated_seqs[sp_complex], outfilename, "fasta")
            else:
                open(outfilename, 'a').close()

rule cat_seqs:
    input:
        "SpeciesComplexes/{complex}/{complex}_ref.fa", 
        "SpeciesComplexes/{complex}/{complex}_queries.fa"
    output:
        "SpeciesComplexes/{complex}/{complex}_all.fa"
    shell:
        "cat {input} > {output}"
        
rule align_seqs:
    input:
        rules.cat_seqs.output
    output:
        "SpeciesComplexes/{complex}/{complex}_aln.fa"
    shell:
        "mafft {input} > {output}"

rule trim_aln:
    input:
        rules.align_seqs.output
    output:
        "SpeciesComplexes/{complex}/{complex}.phy"
    shell:
        "trimal -in {input} -out {output} -phylip -gappyout"
