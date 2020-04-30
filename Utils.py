import os
import pandas as pd
import warnings
from Bio import SeqIO
from Bio import Entrez
from collections import defaultdict

"""takes a range of alpha-numeric strings separated by `delim` and creates
a `out_delim`-separated list by expanding the range of the numeric portion
eg; "AJ1:AJ4" -> ["AJ1 AJ2 AJ3 AJ4"]
numeric poriton must be at the end of the string, and the number of
non-numeric characters must be specified (`num_letters`)
"""
def expand_range(range_string, num_letters=2, delim=':', out_delim=' '):
    letters=range_string[:num_letters]
    digits_length=len(range_string.split(delim)[0])-num_letters
    (start, end) = [int(string[num_letters:]) for string in range_string.split(delim)]
    return out_delim.join(letters + str(num).rjust(digits_length, '0') for num in range(start, end+1))
        
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
"""
def make_lookup_table(blast_file_names, complexes_dict):
    col_names=["qseqid", "qlen", "stitle", "qstart", "qend", "sstart", "send", "slen", "pident", "length", "evalue", "bitscore"]
    df = pd.concat([pd.read_table(file, names=col_names) for file in blast_file_names])
    df['sspecies'] = df['stitle'].str.extract('(Peltigera \w+)').str.slice_replace(1, 9, '.')
    df = df[(df['sspecies'] != 'Peltigera sp') & (df['pident'] > 95)]
    df = df.sort_values('bitscore').groupby(['qseqid']).last()
    for name in list(df[df['sstart'] > df['send']].index.values):
        warnings.warn("%s in reverse strand. please reverse-complement" % name)
    df = df[df['sstart'] < df['send']]
    df = df.join(pd.DataFrame.from_dict(complexes_dict, orient='index'), 'sspecies', how='inner')
    df = df.rename(columns={0: 'complex'})
    return df['complex'].to_dict()  
      
"""separate sequences into species complexes accoring to dict of species-complex dict
- I think the easiest thing to do here is to just hold all sequences in memory
I don't think there's ever going to be enough of them that this approach will
be a problem. If it is, I can always refactor
"""
def separate_seqs(sequence_files, complexes, lookup='species'):
    separated_seqs = defaultdict(list)
    for sequence_file in sequence_files:
        record_iterator = SeqIO.parse(sequence_file, "fasta")
        for seq_record in record_iterator:
            if lookup == 'species':
                key = ' '.join(['P.', seq_record.description.split(' ')[2]])
            elif lookup == 'description':
                key = seq_record.description
            else:
                raise Exception("lookup not recognised. should be one of `species`, `description`")
            try:
                separated_seqs[complexes[key]].append(seq_record)
            except KeyError:
                warnings.warn("No complex for %s" % key)
    return separated_seqs


if __name__ == '__main__':
    email = os.environ.get('EMAIL')
    download_seqs(expand_range("MF067362:MF067365"), "my_example.fa", email)
    complexes = {'Pcan': ['P. canina', 'P. koponenii', 'P. praetextata', 'P. islandica', 'P. evansiana', 'P. fuscopraetextata'],
             'Pleu': ['P. leucophlebia'], 
             'Pcin': ['P. cinnamomea', 'P. neocanina'], 
             'Paphth': ['P. aphthosa', 'P. britannica', 'P. malacea', 'P. chionophila'],
             'Ppono': ['P. ponojensis', 'P. monticola']
    }
    lookup = make_lookup_table(["blast/ITS_Aug17.bl", "blast/ITS_Aug18.bl"], invert_list_dict(complexes))
    for key in lookup:
        print(key, lookup[key])
