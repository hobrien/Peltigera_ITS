import os
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

"""separate sequences into species complexes accoring to dict of species-complex dict
- I think the easiest thing to do here is to just hold all sequences in memory
I don't think there's ever going to be enough of them that this approach will
be a problem. If it is, I can always refactor
"""
def separate_seqs(sequence_file, complexes):
    separated_seqs = defaultdict(list)
    record_iterator = SeqIO.parse(sequence_file, "fasta")
    for seq_record in record_iterator:
        species = ' '.join(['P.', seq_record.description.split(' ')[2]])
        try:
            separated_seqs[complexes[species]].append(seq_record)
        except KeyError:
            print("No complex for %s" % species)
    return separated_seqs
    

#seq_dict = SeqIO.index("seqs.fasta", "fasta")

if __name__ == '__main__':
    email = os.environ.get('EMAIL')
    download_seqs(expand_range("MF067362:MF067365"), "my_example.fa", email)