Peltigera ITS species identification pipeline

- This workflow will place unidentified sequences in "species complexes" using a top blast hit approach, then align them to all other members of that species complex and build a tree.
- Sequences are automatically downloaded from NCBI and used to build a blast DB
- The following species complexes are currently supported:
    - Section Peltidea (Paphth)
        - P. aphthosa, P. britannica, P. malacea, P. chionophila
    - Section Chloropeltigera (Pleu)
        - P. leucophlebia Clade I, P. leucophlebia Clade II, P. leucophlebia Clade III, 
    - Section Peltigera group C (Pcan)
        - P. canina, P. koponenii, P. praetextata, P. islandica, P. evansiana, P. "fuscopraetextata"
    - Section Peltigera group D (Pcin)
        - P. cinnamomea, P. "neocanina"

- `EMAIL` environmental variable needs to be set for this to work. The easiest way to do this is to add "export EMAIL=..." to `~/.bash_profile`