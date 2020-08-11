import os
import re
import pandas as pd
import warnings
from Bio import SeqIO
from Bio import Entrez
from collections import defaultdict

"""takes a range of alpha-numeric strings separated by `delim` and creates
a `out_delim`-separated list by expanding the range of the numeric portion
eg; "AJ1:AJ4,BG5" -> ["AJ1 AJ2 AJ3 AJ4 BG5"]
numeric poriton must be at the end of the string, and the number of
non-numeric characters must be specified (`num_letters`)
"""
def expand_range(range_strings, num_letters=2, range_delim=':', list_delim=',', out_delim=' '):
    out_list = []
    for range_string in range_strings.split(','):
        if ':' in range_string:
            letters=range_string[:num_letters]
            digits_length=len(range_string.split(range_delim)[0])-num_letters
            (start, end) = [int(string[num_letters:]) for string in range_string.split(range_delim)]
            out_list = out_list + [letters + str(num).rjust(digits_length, '0') for num in range(start, end+1)]
        else:
            out_list.append(range_string)
    return out_delim.join(out_list)
        
"""eg download_seqs(expand_range("FJ708820 FJ708821 FJ708822 FJ708823"), "my_example.fa", email)
"""
def download_seqs(accessions, filename, email, db='nucleotide', filetype='fasta'):
    Entrez.email=email
    handle = Entrez.efetch(db=db, id=accessions, rettype=filetype, retmode="text")
    return SeqIO.write(SeqIO.parse(handle, filetype), filename, filetype)

"""takes dict with lists as values and returns a dict where each list member is a key
and the values are the original dic keys; eg {'A':[1,2,3], 'B':[4,5]} -> {1:'A',2:'A',3:'A',4:'B',5:'B'}
"""
def invert_list_dict(list_dict):
    inverted_dict = {}
    for key in list_dict:
        for value in list_dict[key]:
            assert value not in inverted_dict
            inverted_dict[value] = key
    return inverted_dict

"""Parse blast outputs, extract species of top hit, compare to species complex list
and return a dict with sequence descriptions as keys and species complexes as values

returns a two-item tuple. first item includes all seqs. second includes only sequences to be reverse-complemented
"""
def make_lookup_table(blast_file_names, complexes_dict):
    col_names=["qseqid", "qlen", "stitle", "qstart", "qend", "sstart", "send", "slen", "pident", "length", "evalue", "bitscore"]
    df = pd.concat([pd.read_table(file, names=col_names) for file in blast_file_names])
    df['sspecies'] = df['stitle'].str.extract('(Peltigera \w+)').str.slice_replace(1, 9, '.')
    df = df[(df['sspecies'] != 'Peltigera sp') & (df['pident'] > 95)]
    df = df.sort_values('bitscore').groupby(['qseqid']).last()
    for name in list(df[df['sstart'] > df['send']].index.values):
        warnings.warn("%s in reverse strand" % name)
    df_rc = df[df['sstart'] > df['send']]
    df_rc = df_rc.join(pd.DataFrame.from_dict(complexes_dict, orient='index'), 'sspecies', how='inner')
    df_rc = df_rc.rename(columns={0: 'complex'})
    df = df.join(pd.DataFrame.from_dict(complexes_dict, orient='index'), 'sspecies', how='inner')
    df = df.rename(columns={0: 'complex'})
    return (df['complex'].to_dict(), df_rc['complex'].to_dict())
      
"""separate sequences into species complexes accoring to dict of species-complex
- I think the easiest thing to do here is to just hold all sequences in memory
I don't think there's ever going to be enough of them that this approach will
be a problem. If it is, I can always refactor

complexes is a two-item tuple. first item includes all seqs. second includes only sequences to be reverse-complemented
"""
def separate_seqs(sequence_files, complexes, lookup, excluded, manual_refs={}):
    separated_seqs = defaultdict(list)
    excluded = set(excluded)
    for sequence_file in sequence_files:
        record_iterator = SeqIO.parse(sequence_file, "fasta")
        for seq_record in record_iterator:
            if seq_record.id in excluded:
                continue
            seq_record.description = seq_record.description.replace("aff. ", "")
            if seq_record.id.split('.')[0] in manual_refs:
                key = manual_refs[seq_record.id.split('.')[0]]
            else:
                if lookup == 'species':
                    id = ' '.join(['P.', seq_record.description.split(' ')[2]])
                elif lookup == 'description':
                    id = seq_record.description
                else:
                    raise Exception("lookup not recognised. should be one of `species`, `description`")
                try:
                    key = complexes[0][id]
                except KeyError:
                    warnings.warn("No complex for %s (%s)" % (id, seq_record.id))
                if id in complexes[1]:
                    seq_record.seq = seq_record.seq.reverse_complement()
            seq_record.description = re.sub(r'\([^)]*\)', '', seq_record.description)
            seq_record.id = seq_record.description.replace(' ', '_')
            seq_record.description = ''
            separated_seqs[key].append(seq_record)
    return separated_seqs


if __name__ == '__main__':
    email = os.environ.get('EMAIL')
    download_seqs(expand_range("KC139749:KC139751,KC139753:KC139761,KJ413197,KJ413220:KJ413226,KJ413246:KJ413250,KJ616384:KJ616388,KT861422:KT861425,KT948755:KT948756,KU954062:KU954063"), "RefSequences/Lichenol16.fa", email)
    complexes = {'Pcan': ['P. canina', 'P. koponenii', 'P. praetextata', 'P. islandica', 'P. evansiana', 'P. fuscopraetextata'],
             'Pleu': ['P. leucophlebia'], 
             'Pcin': ['P. cinnamomea', 'P. neocanina'], 
             'Paphth': ['P. aphthosa', 'P. britannica', 'P. malacea', 'P. chionophila'],
             'Ppono': ['P. ponojensis', 'P. monticola']
    }
    complexes = invert_list_dict(complexes)
    lookup = make_lookup_table(["blast/ITS_Aug17.bl", "blast/ITS_Aug18.bl"], complexes)
    for key in lookup:
        print(key, lookup[key])
        
    print(separate_seqs(["RefSequences/Lichenol16.fa"], complexes, lookup='species'))
