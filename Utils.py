from Bio import SeqIO
from Bio import Entrez
email = "heath.obrien@gmail.com"  # Always tell NCBI who you are

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

#seq_dict = SeqIO.index("seqs.fasta", "fasta")

if __name__ == '__main__':
    download_seqs(expand_range("MF067362:MF067365"), "my_example.fa", email)