# PythonII_Portfolio

## Sequence Objects 1-4

Portfolio of python code for BISC 450C (Python 2)
```python



from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
# does not count over lapping occurrences
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# total number bases in sequence
len(my_seq)
```




    32




```python
my_seq.count("G")
```




    9




```python
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    46.875




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
gc_fraction(my_seq)
```




    0.46875




```python
my_seq[4:12]
```




    Seq('GATGGGCC')




```python
# Slices with start and stop (Stride)
my_seq[0::3]
```




    Seq('GCTGTAGTAAG')




```python
my_seq[1::3]
```




    Seq('AGGCATGCATC')




```python
my_seq[2:3]
```




    Seq('T')




```python
# Use negatives to start sequence from opposite ends (prints backwards)
my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')




```python
str(my_seq)
```




    'GATCGATGGGCCTATATAGGATCGAAAATCGC'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
# Create name for fasta string
print(fasta_format_string)
```

    >Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" * 10)
```


```python
# joins 10 Ns between each sequence
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
# shortcut to change sequence to uppercase
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# shortcut to change sequence to lowercase
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# false because of capitalization
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# complementary base pairs to go with sequence
my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
# complementary base pairs with sequence in reverse (same as -1 function)
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
# V can mean A C or G
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
# Turns Ts into Us for rna
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
# using specific table to generate sequence
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
# stoping at stop codon from line 58
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table =2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGCCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# prints stop codons from mito table
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
# prints start codons from mito table
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
# double equal signs used to double check sequence
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# None because not defined in sequence
seq[1000:1020]
```




    Seq(None, length=20)




```python
# Defined in Line 96 so able to produce Sequence
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# Not putting number after : results from entered number through end of sequence
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length =10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGTGCCCGA")
```


```python
# For randomly mutating genes
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGTGCCCGA')




```python
# after importing mutable seq, able to mutate sequence as desired, this changes 5th letter to C, seen in Line 109
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGTGCCCGA')




```python
# In mutable seq, removes first desired letter from sequence, seen in Line 111
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGTGCCCGA')




```python
# Prints in reverse, Line 113
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGAAAGTCGCCGGGTAATGCACCG')




```python
# Changes back to sequence that is no longer mutable "protected"
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```




    'AVMGRWKGGRAAG*'




```python

```

## Sequence Annotations 1-4

```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
# Storing sequences with anotations
simple_seq_r.description = "Made up sequence for teh VDB Computational Biology Class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for teh VDB Computational Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
# Now holds annotations listed from Lines 10 & 11
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for teh VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
# prints annotations
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
# per letter annotations (scores for sequences)
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
# in line 22 record saved as file, printed with annotations
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# Can pull only id
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
# Can pull out only description
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# Lines 28-30 will print as empty based on information not being saved in file imported
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.GB.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# Not printed empty, because information was stored in new file imported (Same with lines 42 & 43)
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
# Source where DNA is from
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
# prints project number
record.dbxrefs
```




    ['Project:58037']




```python
#Shows how many features are in table
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
# Stating start position for gene feature is somewhere after 5
start_pos = SeqFeature.AfterPosition(5)
```


```python
# Stating end position is somewhere between 8 & 9
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
# Lines 59-66 printing annotations created
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
from Bio.SeqRecord import SeqRecord
```


```python
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD" "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK" "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM" "SSAC"), id = "gi|14150838|gb|AAK54648.1|AF376133_1", description = "chalcone synthase [Cucumis sativus]")
```


```python
print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC
    



```python
# Gives information, but not full sequence
print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')



```python
from Bio import SeqIO
```


```python
# Read in record, now that we have downloaded record
record = SeqIO.read("NC_005816.GB.txt", "genbank")
```


```python
# Prints sequence information
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# total number of base pairs
len(record)
```




    9609




```python
# total number of features in record
len(record.features)
```




    41




```python
# Only prints feature number of number in brackets
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# Selecting location from Line 82 and subdividing
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
# prints number of features in selected area
len(sub_record.features)
```




    2




```python
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
# information from Line 81 stays same, but location changes due to selecting only sub record
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
# Can add annotations to record
sub_record.annotations["topology"] = "linear"
```


```python
# topology- plasma circular DNA, but since looking at "chunk" can say its linear
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_record.id
```




    'NC_005816.1'




```python
sub_record.name
```




    'NC_005816'




```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# Can edit sub record description taken from Line, 99
sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
# printing record between points
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...



```python
record = SeqIO.read("NC_005816.GB.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
# pulls up key words from annotations
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# Circular genome, endless, start and stop codon can be any, so can shift sequence
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# Length should be the same, as we only shifted where we started (endless circle)
len(shifted)
```




    9609




```python
# Lost feature after shift
len(shifted.features)
```




    40




```python
# now only has molecule type
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
# now empty. annotations not transferred after shift
shifted.dbxrefs
```




    []




```python
# Chanbes annotations to as they were before shift
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
# Now back to how they were before
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# printing first value of string (s) and next four as intergers (i)
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
rc = record.reverse_complement(id="Testing")
```


```python
# reverse complement (backwards)
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# Same as Line 130, but using reverse complement shorthand
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



```python

```

## Sequence IO 1-3

```python
from Bio import SeqIO
```


```python
# printint multiple items form fasta document
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
# Pulls out all identifers in sequence form ls_orhcid.gbk
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
# Choosing specific record
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
# Prints out next record in order, do not have to specify number
second_record = next(record_iterator)
```


```python
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
# printing multiple
print("The last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records[0]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740



```python
from Bio import SeqIO
```


```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
# printing all and only annotations
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
# annotations with numerical values
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
# tells organism
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
# also tells organism
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
# Creating empty vector
all_species = []
```


```python

for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
# prints all 94 species found in ls_orchid.gbk
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species= [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
# Renaming id
first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
first_record.description = first_record.id + " " + "Mutations induced randomly"
```


```python
# Printing through the first 200 bases adding new description
print(first_record.format("fasta"[:200]))
```

    >new_id Mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord(
Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
       "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
       "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
       "SSAC",), 
    id = "gi|14150838|gb|AAK54658.1| AF376133_1", 
    description = "chalcone synthase [Cucumis sativus]")
```


```python
rec2 = SeqRecord(
Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
       "DMVVVEIPKLGKEAAVKAIKEWGQ",),
    id = "gi|13919613|gb|AAK33142.1|",
    description = "chalcone synthase [Fragaria vesca subsp. bracteata]")
```


```python
rec3 = SeqRecord(
Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
       "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
       "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
       "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
       "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
       "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
       "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
   ), 
    id="gi|13925890|gb|AAK49457.1|",
    description = "chalcone synthase [Nicotiana tabacum]",)
```


```python
my_records = [rec1, rec2, rec3]
```


```python
from Bio import SeqIO
SeqIO.write(my_records, "my_example.faa", "fasta")
```




    3




```python
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
count = SeqIO.write(records, "my_example.fasta", "fasta")
print("Convertd %i recoreds" % count)
```

    Convertd 94 recoreds



```python
# Prints 94 reverse complements 
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement)
```

    Z78533.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')>
    Z78532.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')>
    Z78531.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')>
    Z78530.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')>
    Z78529.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')>
    Z78527.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')>
    Z78526.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')>
    Z78525.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')>
    Z78524.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')>
    Z78523.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')>
    Z78522.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')>
    Z78521.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')>
    Z78520.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')>
    Z78519.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')>
    Z78518.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')>
    Z78517.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')>
    Z78516.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')>
    Z78515.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')>
    Z78514.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')>
    Z78513.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')>
    Z78512.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')>
    Z78511.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')>
    Z78510.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')>
    Z78509.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')>
    Z78508.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')>
    Z78507.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')>
    Z78506.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')>
    Z78505.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')>
    Z78504.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')>
    Z78503.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')>
    Z78502.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')>
    Z78501.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')>
    Z78500.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')>
    Z78499.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')>
    Z78498.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')>
    Z78497.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78496.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')>
    Z78495.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')>
    Z78494.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')>
    Z78493.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')>
    Z78492.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')>
    Z78491.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')>
    Z78490.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78489.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')>
    Z78488.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')>
    Z78487.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')>
    Z78486.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')>
    Z78485.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')>
    Z78484.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')>
    Z78483.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78482.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')>
    Z78481.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78480.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78479.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')>
    Z78478.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')>
    Z78477.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78476.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78475.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')>
    Z78474.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')>
    Z78473.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78472.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78471.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78470.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')>
    Z78469.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')>
    Z78468.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')>
    Z78467.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')>
    Z78466.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78465.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78464.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')>
    Z78463.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')>
    Z78462.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')>
    Z78461.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')>
    Z78460.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')>
    Z78459.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')>
    Z78458.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')>
    Z78457.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')>
    Z78456.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78455.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')>
    Z78454.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')>
    Z78453.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78452.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')>
    Z78451.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')>
    Z78450.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')>
    Z78449.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')>
    Z78448.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78447.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')>
    Z78446.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')>
    Z78445.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')>
    Z78444.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')>
    Z78443.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')>
    Z78442.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')>
    Z78441.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')>
    Z78440.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')>
    Z78439.1
    <bound method _SeqAbstractBaseClass.reverse_complement of Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')>



```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
len(records)
```




    94




```python
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
# From Line 81, only 18 are less than 700
len(records)
```




    18




```python
records = (
rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
)

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18




```python

```

## Multiple Sequence Alignment and Pairwise Alignment

```python
pip install biopython
```

    Requirement already satisfied: biopython in /home/student/anaconda3/lib/python3.7/site-packages (1.81)
    Requirement already satisfied: numpy in /home/student/anaconda3/lib/python3.7/site-packages (from biopython) (1.17.2)
    Note: you may need to restart the kernel to use updated packages.



```python
from Bio import AlignIO
```


```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Alignment length %i" % alignment.get_alignment_length())
```

    Alignment length 52



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
# To look at all annotations
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
# to write an align file
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
align1 = MultipleSeqAlignment(
    [
    SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
    SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"), 
    SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
])

align2 = MultipleSeqAlignment(
    [
    SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
    SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
    SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
])

align3 = MultipleSeqAlignment(
    [
    SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
    SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
    SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
    
])
```


```python
my_alignments = [align1, align2, align3]
```


```python
my_alignments
```




    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7ff012231910>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7ff012231210>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7ff012231c90>]




```python
print(my_alignments)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7ff012231910>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7ff012231210>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7ff012231c90>]



```python
# Creating file my_example.phy
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
# Printing contents of file
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
# Before continuing with Line 33, must set x_align first
last_align = alignments[-1]
```


```python
# Printing last multiple sequence line (align3)
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
first_align = alignments[0]
```


```python
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
from Bio import AlignIO
```


```python
count = AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
# How to see how many alignments have been converted
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
alignments = AlignIO.parse("PF05371_seed.sth", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
# Another method of how to see how many alignments have been converted
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
# Creating file
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
AlignIO.convert("PF05371_seed.sth", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
# Manipulating identifiers/name of files
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" %i
```


```python
# numbered elements in Line 47
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA seq0
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA seq1
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA seq2
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA seq3
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA seq4
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA seq5
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA seq6



```python
alignment = AlignIO.read("PF05371_seed.sth", "stockholm")
```


```python
print("Number of ros: %i" % len(alignment))
```

    Number of ros: 7



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
# Printing specific rows
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
# Selecting from column and row
print(alignment[2,6])
```

    T



```python
# Another way to select from column and row
print(alignment[2].seq[6])
```

    T



```python
# Printing specific column
print(alignment[:, 6])
```

    TTT---T



```python
# Selecting specific alignments (alginments 3-6 and first 6 letters)
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
# All alignments through first six letters
print(alignment[:,:6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
# Opposite of 63. Useful in selecting chunks of alignments to remove or use
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
# : after number selects everything AFTER that number
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
# Spliced
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited.sort()
```


```python
# Sorted based on name of id
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49



```python
pip install biopython
```

    Collecting biopython
    [?25l  Downloading https://files.pythonhosted.org/packages/ad/a4/237edd5f5e5b68d9543c79bcd695ef881e6317fbd0eae1b1e53e694f9d54/biopython-1.81.tar.gz (19.3MB)
    [K     |████████████████████████████████| 19.3MB 20.7MB/s eta 0:00:01
    [?25hRequirement already satisfied: numpy in /home/student/anaconda3/lib/python3.7/site-packages (from biopython) (1.17.2)
    Building wheels for collected packages: biopython
      Building wheel for biopython (setup.py) ... [?25ldone
    [?25h  Created wheel for biopython: filename=biopython-1.81-cp37-cp37m-linux_x86_64.whl size=3114517 sha256=a95347e5c22734cdcb357e0f6a7fee59b946f80627cce11fc549f2ff42b8d721
      Stored in directory: /home/student/.cache/pip/wheels/a8/44/b6/119e779b15a4e3d14d4974defdc045d5500eb08ab6de835e08
    Successfully built biopython
    Installing collected packages: biopython
    Successfully installed biopython-1.81
    Note: you may need to restart the kernel to use updated packages.



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Align import MultipleSeqAlignment
```


```python
alignment = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTCCTA"), id="seq1"),
        SeqRecord(Seq("AAT-CTA"), id="seq2"),
        SeqRecord(Seq("CCTACT-"), id="seq3"),
        SeqRecord(Seq("TCTCCTC"), id="seq4"),
    ]
)
```


```python
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
# Taking all pairs of rows in alignment and counting number of times two letters are aligned together (See Line 9)
substitutions = alignment.substitutions
```


```python
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
# Adding a fourth substitution, but NO "Gs" in alignment so reports as 0
m = substitutions.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
# Rearranging alignments
m = substitutions.select("ACTG")
```


```python
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
import Bio.Align.Applications
```


```python
# Directory for algorithms for aligning sequences
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
from Bio import AlignIO
```


```python
align = AlignIO.read("opuntia.aln", "clustal")
```


```python
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
from Bio import Phylo
```


```python
tree = Phylo.read("opuntia.dnd", "newick")
```


```python
# Creating tree diagram using text, ASCII related to statistics and end products
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    



```python

```

```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
aligner = Align.PairwiseAligner(match_score = 1.0)
```


```python
# Creating target/ what we are aligning to
target = "GAACT"
```


```python
query = "GAT"
```


```python
score = aligner.score(target,query)
```


```python
score
```




    3.0




```python
alignments = aligner.align(target, query)
```


```python
# Showing posibilities for alignments(Match = 1pt, Mismatch = -1pt, Gap = -.5pt,etc)
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
aligner.mode = "local"
```


```python
target = "AGAACTC"
```


```python
query = "GAACT"
```


```python
score = aligner.score(target, query)
```


```python
# Score of 5, ignores first two sequences
score
```




    5.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
# Prints all scores
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
# What we are using
aligner.algorithm
```




    'Smith-Waterman'




```python
# Setting a significance, two scores considered equal for purposes of alignment if difference between is less than epsilon
aligner.epsilon
```




    1e-06




```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
# Setting back to original target
target = "GAACT"
```


```python
query = "GAT"
```


```python
alignments = aligner.align(target, query)
```


```python
# Looking at positions of aligner
alignment = alignments[0]
```


```python
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7fbb787d1610>




```python
# Looking at alignment score for object
alignment.score
```




    3.0




```python
# Looking at what your target is
alignment.target
```




    'GAACT'




```python
# Looking at what your query is 
alignment.query
```




    'GAT'




```python
# Printing so you can see alignment
print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    



```python
# Pulling coordinates of alignment in array form, can see what they align to (based on Line 34)
alignment.coordinates
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
len(alignment)
```




    2




```python
# Shape of alignment (2 rows, 5 columns)
alignment.shape
```




    (2, 5)




```python
# Can set mode 
aligner.mode = "local"
```


```python
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
local_alignment = local_alignments[0]
```


```python
# Because it is local, does not include everything, drops off two end letters from Line 39
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
# Same as Line 37 (2 rows, 4 columns)
local_alignment.shape
```




    (2, 4)




```python
aligner.mode = "global"
```


```python
# Changing mismatch score, see results in Line 47
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
```


```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
len(alignments)
```




    2




```python
# Printing first alignment, introduced gap, would rather gap than mismatch
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
# Printing second alignment stored, same as Line 50, but G & C reversed
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
# Sorting based on number so can print as seen in Line 54
local_alignment.sort()
```


```python
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    



```python
from Bio import Align
```


```python
from Bio.Seq import reverse_complement
```


```python
target = "AAACCC"
```


```python
query = "AACC"
```


```python
aligner = Align.PairwiseAligner(mismatch_score = -1, internal_gap_score = -1)
```


```python
# Aligning to positive strand results in 4
aligner.score(target,query)
```




    4.0




```python
# Aligning to reverse complement strand results in 0
aligner.score(target,reverse_complement(query))
```




    0.0




```python
# Aligning to reverse complement of negative strand results in 4 (basically double negative, cancels out, same as positive strand)
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
# Aligning to negative strand also results in 0 (is reverse complement)
aligner.score(target, query, strand = "-")
```




    0.0




```python
alignments = aligner.align(target,query)
```


```python
len(alignments)
```




    1




```python
# Printing first alignment 
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# Saving as different file format, "bed" type of alignment format used, can also use PSL, SAM, etc.
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
# Results in same as normal alignment due to negative of reverse complement (double negative), See results in Line 71
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, query, strand = "-" )
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
# Reverse of Line 76 (negative score, zero overlap)
print(alignments[1])
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
aligner.left_gap_score = -0.5
```


```python
aligner.right_gap_score = -0.2
```


```python
aligner.score(target, query)
```




    3.3




```python
alignments = aligner.align(target,query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments)
```

    <Bio.Align.PairwiseAlignments object at 0x7fbb78742510>



```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
# Shows printing/mapping DNA to reverse complement gives worse alignment vs aligning to standard string
aligner.score(target, reverse_complement(query), strand = "+")
```




    -2.0




```python
# Better results than Line 90, same as Line 88
aligner.score(target, query, strand = "+")
```




    3.3




```python

```


## Challenge 1

```python
from Bio.Blast import NCBIWWW
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "324")
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("APC_324.fasta", format = "fasta")
```


```python
print(record)
```

    ID: NC_000005.10:112707498-112846239
    Name: NC_000005.10:112707498-112846239
    Description: NC_000005.10:112707498-112846239 Homo sapiens chromosome 5, GRCh38.p14 Primary Assembly
    Number of features: 0
    Seq('GCATTGTAGTCTTCCCACCTCCCACAAGATGGCGGAGGGCAAGTAGCAAGGGGG...GCA')



```python
with open("APC_324.fast", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520B), microRNA
    length: 61
    e value: 4.91307e-23
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b (MIR520B), microRNA
    length: 60
    e value: 1.71483e-22
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133242|ref|NR_032573.1| Macaca mulatta microRNA mir-519a (MIR519A), microRNA
    length: 85
    e value: 2.54503e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| |||||||||||||||||||||||||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTGGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171267|ref|NR_035851.1| Pan troglodytes microRNA mir-519b (MIR519B), microRNA
    length: 80
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205302|ref|NR_030191.1| Homo sapiens microRNA 519b (MIR519B), microRNA
    length: 81
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171259|ref|NR_035850.1| Pan troglodytes microRNA mir-519a (MIR519A), microRNA
    length: 86
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205451|ref|NR_030222.1| Homo sapiens microRNA 519a-2 (MIR519A2), microRNA
    length: 87
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171276|ref|NR_035852.1| Pan troglodytes microRNA mir-519c (MIR519C), microRNA
    length: 86
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205290|ref|NR_030188.1| Homo sapiens microRNA 519c (MIR519C), microRNA
    length: 87
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171354|ref|NR_035860.1| Pan troglodytes microRNA mir-520f (MIR520F), microRNA
    length: 86
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205281|ref|NR_030186.1| Homo sapiens microRNA 520f (MIR520F), microRNA
    length: 87
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 4.60152e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||| ||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAGAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTCTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171394|ref|NR_035865.1| Pan troglodytes microRNA mir-522 (MIR522), microRNA
    length: 86
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205429|ref|NR_030218.1| Homo sapiens microRNA 519a-1 (MIR519A1), microRNA
    length: 85
    e value: 1.60609e-16
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| |||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205423|ref|NR_030217.1| Homo sapiens microRNA 522 (MIR522), microRNA
    length: 87
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171401|ref|NR_035866.1| Pan troglodytes microRNA mir-523 (MIR523), microRNA
    length: 79
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133247|ref|NR_032574.1| Macaca mulatta microRNA mir-519b (MIR519B), microRNA
    length: 81
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  |||||||||||||||||||||||||||||||||| |||| ||| |||||||||...
    CCCTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGTGCATCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205309|ref|NR_030193.1| Homo sapiens microRNA 523 (MIR523), microRNA
    length: 87
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 1.95662e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| ||||||||||||||||  |||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAATAAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| ||   | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTTATTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171437|ref|NR_035870.1| Pan troglodytes microRNA mir-526a-2 (MIR526A-2), microRNA
    length: 68
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||||||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||| |||||||||| | |||||| |||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGGAAGAAAAGAATGCGCTTCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 9.49283e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||| |||| || |   | | || ||||||||||||| | |||||||...
    CCCTCTAAAGGGAAGCGCATTCTTTTCTTCCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171211|ref|NR_035845.1| Pan troglodytes microRNA mir-518b (MIR518B), microRNA
    length: 82
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| ||||||||||| |||||||||||||||||| || |||||||||...
    CCCTCTAGAAGGAAGCACTTTCTGTTGTTTGAAAGAAAAGAAAGTGCATCATTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 3.31332e-06
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |  || || |||||| || | | ||| || ||||||||||||||| |||||||...
    CCTCTAAAATGATGCACTTTCTTTTCTTTCAAACAACAGAAAGTGCTTCCTTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205321|ref|NR_030196.1| Homo sapiens microRNA 518b (MIR518B), microRNA
    length: 83
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171342|ref|NR_035859.1| Pan troglodytes microRNA mir-520e (MIR520E), microRNA
    length: 86
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205265|ref|NR_030183.1| Homo sapiens microRNA 520e (MIR520E), microRNA
    length: 87
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133223|ref|NR_032569.1| Macaca mulatta microRNA mir-518c (MIR518C), microRNA
    length: 101
    e value: 1.01355e-12
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  ||||||||||||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133232|ref|NR_032571.1| Macaca mulatta microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||| |||||||  |||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGATTTCTGTGATCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171300|ref|NR_035855.1| Pan troglodytes microRNA mir-520a (MIR520A), microRNA
    length: 84
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171161|ref|NR_035839.1| Pan troglodytes microRNA mir-516b-1 (MIR516B-1), microRNA
    length: 89
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||| |||||||||    |||||| |||||| |||||||||||||||||||...
    CCCTCCAAAGGGAAGCACTTTCTGTTTGTTGTCTGAGAGAAAACAAAGTGCTTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | |||||| |||| | || ||| |   || ||| ||||||||||||| ||| |||||...
    CTCTAAAAGGAAGCACTTTGTTTTCTCTCAGACAACAAACAGAAAGTGCTTCCCTTTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | |||||||| |||||||||||||||||||||  |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAGCAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||| ||||||||   || | | | | || ||||||||||||| | | |||||...
    CCTCTAAAGGGGAGCGCTTTGCTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132725|ref|NR_032718.1| Macaca mulatta microRNA mir-516 (MIR516), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205396|ref|NR_030212.1| Homo sapiens microRNA 516b-1 (MIR516B1), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205294|ref|NR_030189.1| Homo sapiens microRNA 520a (MIR520A), microRNA
    length: 85
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    CTCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 0.00599063
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGAG...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| |||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171454|ref|NR_035872.1| Pan troglodytes microRNA mir-527 (MIR527), microRNA
    length: 86
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171333|ref|NR_035858.1| Pan troglodytes microRNA mir-520d (MIR520D), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    CTCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171291|ref|NR_035854.1| Pan troglodytes microRNA mir-519e (MIR519E), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||||||||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTTGTCTGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|269847477|ref|NR_031696.1| Homo sapiens microRNA 1283-2 (MIR1283-2), microRNA
    length: 87
    e value: 4.30974e-11
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||  || |||||||||||||||||||||||||||||||  |||||| ||| |||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAATCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205435|ref|NR_030219.1| Homo sapiens microRNA 527 (MIR527), microRNA
    length: 85
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|270133264|ref|NR_032577.1| Macaca mulatta microRNA mir-520a (MIR520A), microRNA
    length: 85
    e value: 1.50425e-10
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||   ||||||||||||||| ||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTATTTTCTGTTGTCTGAAGGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGA...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| ||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGA...
    ****ALIGNMENT****
    sequence: gi|262205358|ref|NR_030204.1| Homo sapiens microRNA 520d (MIR520D), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    TCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| ||||||||||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTTTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171169|ref|NR_035840.1| Pan troglodytes microRNA mir-516b-2 (MIR516B-2), microRNA
    length: 89
    e value: 5.25034e-10
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| |||||||||||||||||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTGAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTCAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171381|ref|NR_035863.1| Pan troglodytes microRNA mir-521-1 (MIR521-1), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133201|ref|NR_032565.1| Macaca mulatta microRNA mir-517 (MIR517), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||| |||| ||||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTGGTCTAAAAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205417|ref|NR_030216.1| Homo sapiens microRNA 521-1 (MIR521-1), microRNA
    length: 87
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133214|ref|NR_032567.1| Macaca mulatta microRNA mir-518a (MIR518A), microRNA
    length: 87
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTC-CTTT...
    ||||  | |||||||| ||||||||||||| || |||||||||||||||| ||||...
    CCCTACAAAGGGAAGCCCTTTCTGTTGTCTAAACGAAAAGAAAGTGCTTCTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205407|ref|NR_030214.1| Homo sapiens microRNA 517c (MIR517C), microRNA
    length: 95
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||||||||  |||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTTGTCT--AAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205377|ref|NR_030208.1| Homo sapiens microRNA 526a-2 (MIR526A2), microRNA
    length: 65
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| ||||||||||    |||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTG----AAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133311|ref|NR_032588.1| Macaca mulatta microRNA mir-523b (MIR523B), microRNA
    length: 89
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCT--GAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| || |||||||||||||||| ||   |||||| || ||| |||||| |||||||...
    CCCTCTAGAGCGAAGCGCTTTCTGTTGGCTAGAAAAGAATAGGAAGCGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133297|ref|NR_032585.1| Macaca mulatta microRNA mir-521 (MIR521), microRNA
    length: 87
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| ||| |||...
    CCCTCCAAAGGGAAGTACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 2.2325e-08
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132721|ref|NR_032717.1| Macaca mulatta microRNA mir-524 (MIR524), microRNA
    length: 85
    e value: 2.71974e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||  | ||| |||| |||||| |||||| ||||||||||| |||||||| |||...
    CCCTACAAAGGCAAGCACTTTCTCTTGTCTAAAAGAAAAGAAGGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205276|ref|NR_030185.1| Homo sapiens microRNA 519e (MIR519E), microRNA
    length: 84
    e value: 9.49283e-07
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||   |||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTT---TGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171363|ref|NR_035861.1| Pan troglodytes microRNA mir-520g (MIR520G), microRNA
    length: 89
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| | || |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-ACGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133196|ref|NR_032564.1| Macaca mulatta microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133191|ref|NR_032563.1| Macaca mulatta microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171387|ref|NR_035864.1| Pan troglodytes microRNA mir-521-2 (MIR521-2), microRNA
    length: 86
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205353|ref|NR_030203.1| Homo sapiens microRNA 521-2 (MIR521-2), microRNA
    length: 87
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133237|ref|NR_032572.1| Macaca mulatta microRNA mir-518f (MIR518F), microRNA
    length: 87
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||  ||||||||||||||| |  | | | | || | ||||||||||| | |||||||...
    CCTCTGAAGGGAAGCGCTTTCTTTCCTTTCACACAAGATAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171283|ref|NR_035853.1| Pan troglodytes microRNA mir-519d (MIR519D), microRNA
    length: 87
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTCAAACAAAGTGCCTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133288|ref|NR_032583.1| Macaca mulatta microRNA mir-520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||  |||| |||||||||    ||||||   |||| ||||||||||| || ||||...
    CCCTCTAGAGA-AAGCACTTTCTGTTTGTTGTCTGAGGAAAAACAAAGTGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205371|ref|NR_030207.1| Homo sapiens microRNA 516b-2 (MIR516B2), microRNA
    length: 85
    e value: 0.00171634
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||     |||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAG-----AAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205349|ref|NR_030202.1| Homo sapiens microRNA 519d (MIR519D), microRNA
    length: 88
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTTAAACAAAGTGCCTCCCTTTAGAG...



```python
# the e-value when compared to chimpanzee(pan troglodytes) is 1.71483e-22
```


```python

```


## Blast 

```python
from Bio.Blast import NCBIWWW
```


```python
# Shorthand for searching for information in NCBI 
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("m_cold.fasta", format = "fasta")
```


```python
print(record)
```

    ID: gi|8332116|gb|BE037100.1|BE037100
    Name: gi|8332116|gb|BE037100.1|BE037100
    Description: gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
    Number of features: 0
    Seq('CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...TTC')



```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
# Use result handle results by reading in and closing to NCBI search
with open("m_cold.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
# Saving as XML file
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
# Looking specific P value or E value at significance of 0.04 or higher
E_VALUE_THRESH = 0.04
```


```python
# Writing loop through hsp if hsp.expect is less than e value thresh print all recorded information
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520B), microRNA
    length: 61
    e value 4.91307e-23
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b (MIR520B), microRNA
    length: 60
    e value 1.71483e-22
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133242|ref|NR_032573.1| Macaca mulatta microRNA mir-519a (MIR519A), microRNA
    length: 85
    e value 2.54503e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| |||||||||||||||||||||||||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTGGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171267|ref|NR_035851.1| Pan troglodytes microRNA mir-519b (MIR519B), microRNA
    length: 80
    e value 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205302|ref|NR_030191.1| Homo sapiens microRNA 519b (MIR519B), microRNA
    length: 81
    e value 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171259|ref|NR_035850.1| Pan troglodytes microRNA mir-519a (MIR519A), microRNA
    length: 86
    e value 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205451|ref|NR_030222.1| Homo sapiens microRNA 519a-2 (MIR519A2), microRNA
    length: 87
    e value 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171276|ref|NR_035852.1| Pan troglodytes microRNA mir-519c (MIR519C), microRNA
    length: 86
    e value 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205290|ref|NR_030188.1| Homo sapiens microRNA 519c (MIR519C), microRNA
    length: 87
    e value 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171354|ref|NR_035860.1| Pan troglodytes microRNA mir-520f (MIR520F), microRNA
    length: 86
    e value 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205281|ref|NR_030186.1| Homo sapiens microRNA 520f (MIR520F), microRNA
    length: 87
    e value 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value 4.60152e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||| ||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAGAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTCTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171394|ref|NR_035865.1| Pan troglodytes microRNA mir-522 (MIR522), microRNA
    length: 86
    e value 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205429|ref|NR_030218.1| Homo sapiens microRNA 519a-1 (MIR519A1), microRNA
    length: 85
    e value 1.60609e-16
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| |||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205423|ref|NR_030217.1| Homo sapiens microRNA 522 (MIR522), microRNA
    length: 87
    e value 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171401|ref|NR_035866.1| Pan troglodytes microRNA mir-523 (MIR523), microRNA
    length: 79
    e value 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133247|ref|NR_032574.1| Macaca mulatta microRNA mir-519b (MIR519B), microRNA
    length: 81
    e value 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  |||||||||||||||||||||||||||||||||| |||| ||| |||||||||...
    CCCTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGTGCATCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205309|ref|NR_030193.1| Homo sapiens microRNA 523 (MIR523), microRNA
    length: 87
    e value 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value 1.95662e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| ||||||||||||||||  |||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAATAAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| ||   | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTTATTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171437|ref|NR_035870.1| Pan troglodytes microRNA mir-526a-2 (MIR526A-2), microRNA
    length: 68
    e value 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||||||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||| |||||||||| | |||||| |||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGGAAGAAAAGAATGCGCTTCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value 9.49283e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||| |||| || |   | | || ||||||||||||| | |||||||...
    CCCTCTAAAGGGAAGCGCATTCTTTTCTTCCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171211|ref|NR_035845.1| Pan troglodytes microRNA mir-518b (MIR518B), microRNA
    length: 82
    e value 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| ||||||||||| |||||||||||||||||| || |||||||||...
    CCCTCTAGAAGGAAGCACTTTCTGTTGTTTGAAAGAAAAGAAAGTGCATCATTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value 3.31332e-06
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |  || || |||||| || | | ||| || ||||||||||||||| |||||||...
    CCTCTAAAATGATGCACTTTCTTTTCTTTCAAACAACAGAAAGTGCTTCCTTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205321|ref|NR_030196.1| Homo sapiens microRNA 518b (MIR518B), microRNA
    length: 83
    e value 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171342|ref|NR_035859.1| Pan troglodytes microRNA mir-520e (MIR520E), microRNA
    length: 86
    e value 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205265|ref|NR_030183.1| Homo sapiens microRNA 520e (MIR520E), microRNA
    length: 87
    e value 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133223|ref|NR_032569.1| Macaca mulatta microRNA mir-518c (MIR518C), microRNA
    length: 101
    e value 1.01355e-12
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  ||||||||||||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133232|ref|NR_032571.1| Macaca mulatta microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||| |||||||  |||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGATTTCTGTGATCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171300|ref|NR_035855.1| Pan troglodytes microRNA mir-520a (MIR520A), microRNA
    length: 84
    e value 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171161|ref|NR_035839.1| Pan troglodytes microRNA mir-516b-1 (MIR516B-1), microRNA
    length: 89
    e value 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||| |||||||||    |||||| |||||| |||||||||||||||||||...
    CCCTCCAAAGGGAAGCACTTTCTGTTTGTTGTCTGAGAGAAAACAAAGTGCTTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | |||||| |||| | || ||| |   || ||| ||||||||||||| ||| |||||...
    CTCTAAAAGGAAGCACTTTGTTTTCTCTCAGACAACAAACAGAAAGTGCTTCCCTTTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | |||||||| |||||||||||||||||||||  |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAGCAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||| ||||||||   || | | | | || ||||||||||||| | | |||||...
    CCTCTAAAGGGGAGCGCTTTGCTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132725|ref|NR_032718.1| Macaca mulatta microRNA mir-516 (MIR516), microRNA
    length: 90
    e value 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205396|ref|NR_030212.1| Homo sapiens microRNA 516b-1 (MIR516B1), microRNA
    length: 90
    e value 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205294|ref|NR_030189.1| Homo sapiens microRNA 520a (MIR520A), microRNA
    length: 85
    e value 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    CTCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value 0.00599063
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGAG...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| |||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171454|ref|NR_035872.1| Pan troglodytes microRNA mir-527 (MIR527), microRNA
    length: 86
    e value 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171333|ref|NR_035858.1| Pan troglodytes microRNA mir-520d (MIR520D), microRNA
    length: 86
    e value 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    CTCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171291|ref|NR_035854.1| Pan troglodytes microRNA mir-519e (MIR519E), microRNA
    length: 86
    e value 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||||||||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTTGTCTGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|269847477|ref|NR_031696.1| Homo sapiens microRNA 1283-2 (MIR1283-2), microRNA
    length: 87
    e value 4.30974e-11
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||  || |||||||||||||||||||||||||||||||  |||||| ||| |||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAATCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205435|ref|NR_030219.1| Homo sapiens microRNA 527 (MIR527), microRNA
    length: 85
    e value 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|270133264|ref|NR_032577.1| Macaca mulatta microRNA mir-520a (MIR520A), microRNA
    length: 85
    e value 1.50425e-10
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||   ||||||||||||||| ||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTATTTTCTGTTGTCTGAAGGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGA...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| ||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGA...
    ****ALIGNMENT****
    sequence: gi|262205358|ref|NR_030204.1| Homo sapiens microRNA 520d (MIR520D), microRNA
    length: 87
    e value 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    TCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| ||||||||||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTTTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171169|ref|NR_035840.1| Pan troglodytes microRNA mir-516b-2 (MIR516B-2), microRNA
    length: 89
    e value 5.25034e-10
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| |||||||||||||||||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTGAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTCAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171381|ref|NR_035863.1| Pan troglodytes microRNA mir-521-1 (MIR521-1), microRNA
    length: 86
    e value 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133201|ref|NR_032565.1| Macaca mulatta microRNA mir-517 (MIR517), microRNA
    length: 86
    e value 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||| |||| ||||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTGGTCTAAAAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205417|ref|NR_030216.1| Homo sapiens microRNA 521-1 (MIR521-1), microRNA
    length: 87
    e value 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133214|ref|NR_032567.1| Macaca mulatta microRNA mir-518a (MIR518A), microRNA
    length: 87
    e value 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTC-CTTT...
    ||||  | |||||||| ||||||||||||| || |||||||||||||||| ||||...
    CCCTACAAAGGGAAGCCCTTTCTGTTGTCTAAACGAAAAGAAAGTGCTTCTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205407|ref|NR_030214.1| Homo sapiens microRNA 517c (MIR517C), microRNA
    length: 95
    e value 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||||||||  |||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTTGTCT--AAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205377|ref|NR_030208.1| Homo sapiens microRNA 526a-2 (MIR526A2), microRNA
    length: 65
    e value 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| ||||||||||    |||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTG----AAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133311|ref|NR_032588.1| Macaca mulatta microRNA mir-523b (MIR523B), microRNA
    length: 89
    e value 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCT--GAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| || |||||||||||||||| ||   |||||| || ||| |||||| |||||||...
    CCCTCTAGAGCGAAGCGCTTTCTGTTGGCTAGAAAAGAATAGGAAGCGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133297|ref|NR_032585.1| Macaca mulatta microRNA mir-521 (MIR521), microRNA
    length: 87
    e value 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| ||| |||...
    CCCTCCAAAGGGAAGTACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value 2.2325e-08
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132721|ref|NR_032717.1| Macaca mulatta microRNA mir-524 (MIR524), microRNA
    length: 85
    e value 2.71974e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||  | ||| |||| |||||| |||||| ||||||||||| |||||||| |||...
    CCCTACAAAGGCAAGCACTTTCTCTTGTCTAAAAGAAAAGAAGGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205276|ref|NR_030185.1| Homo sapiens microRNA 519e (MIR519E), microRNA
    length: 84
    e value 9.49283e-07
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||   |||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTT---TGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171363|ref|NR_035861.1| Pan troglodytes microRNA mir-520g (MIR520G), microRNA
    length: 89
    e value 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| | || |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-ACGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133196|ref|NR_032564.1| Macaca mulatta microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133191|ref|NR_032563.1| Macaca mulatta microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171387|ref|NR_035864.1| Pan troglodytes microRNA mir-521-2 (MIR521-2), microRNA
    length: 86
    e value 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205353|ref|NR_030203.1| Homo sapiens microRNA 521-2 (MIR521-2), microRNA
    length: 87
    e value 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133237|ref|NR_032572.1| Macaca mulatta microRNA mir-518f (MIR518F), microRNA
    length: 87
    e value 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||  ||||||||||||||| |  | | | | || | ||||||||||| | |||||||...
    CCTCTGAAGGGAAGCGCTTTCTTTCCTTTCACACAAGATAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171283|ref|NR_035853.1| Pan troglodytes microRNA mir-519d (MIR519D), microRNA
    length: 87
    e value 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTCAAACAAAGTGCCTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133288|ref|NR_032583.1| Macaca mulatta microRNA mir-520g (MIR520G), microRNA
    length: 90
    e value 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||  |||| |||||||||    ||||||   |||| ||||||||||| || ||||...
    CCCTCTAGAGA-AAGCACTTTCTGTTTGTTGTCTGAGGAAAAACAAAGTGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205371|ref|NR_030207.1| Homo sapiens microRNA 516b-2 (MIR516B2), microRNA
    length: 85
    e value 0.00171634
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||     |||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAG-----AAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205349|ref|NR_030202.1| Homo sapiens microRNA 519d (MIR519D), microRNA
    length: 88
    e value 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTTAAACAAAGTGCCTCCCTTTAGAG...



## Open CV Basics Pt. 1

```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
# imported jpeg of random mushrooms as example
img = cv2.imread("Mushrooms.jpeg")
```


```python
type(img)
```




    numpy.ndarray




```python
# Looks as if correct/has ran, but in next line you can see it is wrong
img_wrong = cv2.imread('wrong/path/doesnot/abcdegh.jpg')
```


```python
# Checking information, see Line 7
type(img_wrong)
```




    NoneType




```python
# Shows image, but matplotlib program shows RGB image
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fdf8ec89850>




![png](output_6_1.png)



```python
# Converting image to original color from Blue Green Red to Red Green Blue
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7fdf8ebba890>




![png](output_8_1.png)



```python
# How to see size of image in pixels
img_gray = cv2.imread("Mushrooms.jpeg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (168, 299)




```python
# Does not directly print as grayscale. Same problem as Line 7, must also be converted as in Line 10
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7fdf8db37150>




![png](output_10_1.png)



```python
# Converting and printing image in true gray scale
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fdf8daa8890>




![png](output_11_1.png)



```python
# Resizing image
fix_img.shape
```




    (168, 299, 3)




```python
# Setting new size of image
new_img = cv2.resize(fix_img, (1000,400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7fdf8ca13cd0>




![png](output_13_1.png)



```python
# Seeing new size of image in pixels after resizing 
new_img.shape
```




    (400, 1000, 3)




```python
# Resizing scale of image, using 0.5, cuts size in half
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7fdf8c9ff1d0>




![png](output_16_1.png)



```python
# Shows new size of image in pixels after reducing w & h by half
new_img.shape
```




    (84, 150, 3)




```python
# Flipping image vertically 
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7fdf8c971690>




![png](output_18_1.png)



```python
# Flips image vertically and reverses it
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7fdf8c8daa90>




![png](output_19_1.png)



```python
type(fix_img)
```




    numpy.ndarray




```python
# Saving new image to files under new name, does NOT save color corrections
cv2.imwrite("Mushrooms_fixed_image.jpeg", flip_img)
```




    True




```python

```


## Open CV Basics Pt. 2 

```python

```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread("Mushrooms.jpeg")
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fbf07c8cad0>




![png](output_4_1.png)



```python
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
# Swaps colors back to original 
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fbf063cb550>




![png](output_6_1.png)



```python
# Changes hue saturation value to a "psychedelic" color scheme seen in Line 11
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7fbf0633edd0>




![png](output_8_1.png)



```python
# Changes hue light saturation. Going with Line 10, ways to go between different color values Seen in Line 13
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7fbf06226510>




![png](output_10_1.png)



```python
# Will print both images in blue hues seen in Line 19
img1 = cv2.imread("DNC.jpeg")
img2 = cv2.imread("Mushrooms.jpeg")
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fbf06215ad0>




![png](output_12_1.png)



```python
# Reverting both images back to original format/color scheme seen in Line 22&23
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fbf0619e990>




![png](output_14_1.png)



```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7fbf06109650>




![png](output_15_1.png)



```python
# Resizing both images
img1 = cv2.resize(img1, (1200,1200))
img2 = cv2.resize(img2, (1200,1200))
```


```python
alpha = 0.5
beta = 0.5
```


```python
# Making images equally transparent and will print together as seen in Line 27
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
# ability to do so important for microscopy as in fluorescence image overlay
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7fbf04030b50>




![png](output_19_1.png)



```python
# Adjusting alpha/beta as such will print DNC img brighter than mushroom img
alpha = 0.8
beta = 0.2

blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7fbef6fb71d0>




![png](output_20_1.png)



```python
# Reversing alpha/beta values from Line 28 will produce Mushrooms brigher than DNC giving watermark effect
alpha = 0.2
beta = 0.8

blended2 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended2)
```




    <matplotlib.image.AxesImage at 0x7fbef6f18150>




![png](output_21_1.png)



```python
img1 = cv2.imread('DNC.jpeg')
img2 = cv2.imread('Mushrooms.jpeg')

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (100,100))
```


```python
# Lines 90-91 converting to original color & resizing DNC img to be smaller and printed on top corner of Mushrooms img
large_img = img2
small_img = img1

x_offset = 0
y_offset = 0

x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]

large_img[y_offset:y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```




    <matplotlib.image.AxesImage at 0x7fbf05df5a90>




![png](output_23_1.png)



```python

```

## Open CV Basics Pt. 3

```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread('rainbow.jpg')
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fc268cb3e50>




![png](output_2_1.png)



```python
# Cancelling out background colors
img = cv2.imread('rainbow.jpg', 0)
```


```python
# Reduced color of image, setting to gray scale
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fc2683e0490>




![png](output_4_1.png)



```python
# Thresholding img. 255 total pixels. Thresholding to 1/2 (127). Everything over 127 set to one color anything below set to another 
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
```


```python
ret1
```




    127.0




```python
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc268348ad0>




![png](output_7_1.png)



```python
# Inversing limits- Similar to Lines 7 & 9, but instead of two colors printed thresh_trunc applies an adaptive threshold to an array to transfer into binary image
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc2682aed10>




![png](output_8_1.png)



```python
#image shown seems almost 3D
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc26821ca90>




![png](output_9_1.png)



```python
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc268201510>




![png](output_10_1.png)



```python
def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
show_pic(img_r)
```


![png](output_12_0.png)



```python
# Making short hand and creating same effect as Line 7 (gets rid of gray only prints black and white image)
ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_13_0.png)



```python
# raising threshold, increased the contrast
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_14_0.png)



```python
# Putting max value of adaptive threshold (outlines boxes, removes filling)
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
```


```python
show_pic(th2)
```


![png](output_16_0.png)



```python
#Layering images on top of each other
blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th2, beta = 0.4, gamma = 0)

show_pic(blended)
```


![png](output_17_0.png)



```python
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th3, beta = 0.4, gamma = 0)

show_pic(blended)
```


![png](output_18_0.png)



```python

```

## Corner Detection

```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
flat_chess = cv2.imread('Checkers.jpeg')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fb7c1494f10>




![png](output_1_1.png)



```python
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fb7c0ca4290>




![png](output_2_1.png)



```python
real_chess = cv2.imread("Chess.jpeg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fb7c0c20350>




![png](output_4_1.png)



```python
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fb7c039d450>




![png](output_5_1.png)



```python
# Using Harris Corner Detection (playing with values to see how detection changes)
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2.dilate(dst, None)
```


```python
#Detecting corners and making red
flat_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fb7c032f650>




![png](output_7_1.png)



```python
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize =2, ksize=3, k=0.04)
dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fb7c029ab90>




![png](output_8_1.png)



```python
#Shi-Tomasi Corner Detection 

corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y),3,(255,0,0), -1)
    
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fb7c0292b10>




![png](output_10_1.png)



```python
# Method worked well with other image in Line 16, looked same as image in Line 10, but had more "conservative" approach here by overlapping and adding green dots with the red using edge detection
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fb7c01faad0>




![png](output_11_1.png)



```python

```


## Edge Detection

```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread("Mushrooms.jpeg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fcc688e5910>




![png](output_3_1.png)



```python
# Setting edgres using threshold (127 = median of threshold)
edges = cv2.Canny(image = img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fcc6889f810>




![png](output_4_1.png)



```python
# Printing/finding median color value
med_value = np.median(img)
med_value
```




    169.0




```python
# Adjusting lower/upper values of threshold  
lower = int(max(0, 0.7*med_value))
upper = int(max(255, 1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fcc61ff96d0>




![png](output_6_1.png)



```python
# Adding to upper threshold to adjust image
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper +100)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fcc61f66690>




![png](output_7_1.png)



```python
# Using blurring effect to adjust image
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fcc61eccc90>




![png](output_8_1.png)



```python
# Same as Line 12 but increasing ksize
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fcc61eb8fd0>




![png](output_9_1.png)



```python
# Same as Line 12, but inreasing upper threshold. All blurring effect images produced not as clear as sucessful as before
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 50)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fcc61e31450>




![png](output_10_1.png)



```python
# Trying Line 12 again but increasing upper threshold by 100, still not as successful as Line 7-11
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 100)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fcc61da07d0>




![png](output_11_1.png)



```python

```


## Feature Match 

```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
def display(img, cmap = 'gray'):
    fig = plt.figure(figsize = (12,10))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
fruity_pebbles = cv2.imread("fruity_pebbles", 0)
display(fruity_pebbles)
```


![png](output_2_0.png)



```python
cereals = cv2.imread('All_Cereals.jpeg', 0)
display(cereals)
```


![png](output_3_0.png)



```python
orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(fruity_pebbles, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches= bf.match(des1,des2)
```


```python
matches = sorted(matches, key = lambda x:x.distance)
```


```python
fruity_pebbles_matches = cv2.drawMatches(fruity_pebbles, kp1, cereals, kp2, matches[:25], None, flags = 2)
```


```python
display(fruity_pebbles_matches)
```


![png](output_8_0.png)



```python
sift = cv2.SIFT_create()
```


```python
kp1, des1 = sift.detectAndCompute(fruity_pebbles, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
# knn = known nearest neighbor; finds best matches from set
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1, des2, k=2)
```


```python
good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
print('length of tatal matches:', len(matches))
print('length of good matches:', len(good))
```

    length of tatal matches: 476
    length of good matches: 1



```python
sift_matches = cv2.drawMatchesKnn(fruity_pebbles, kp1, cereals, kp2, good, None, flags = 2)
display(sift_matches)
```


![png](output_14_0.png)



```python
# flann based matching; fast librarying for nearst neighbors; faster than sift
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(fruity_pebbles, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm=flann_index_KDtree, trees = 5)
search_params = dict(checks=50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)

good = []

for match1, match2 in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
flann_matches = cv2.drawMatchesKnn(fruity_pebbles, kp1, cereals, kp2, good, None, flags = 0)
display(flann_matches)
```


![png](output_18_0.png)



```python
sift = cv2.SIFT_create()

kp1, des1 = sift.detectAndCompute(fruity_pebbles, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
flann_index_KDtree = 0
index_params = dict(algorithm= flann_index_KDtree, trees = 5)
search_param = dict(checks = 50)
```


```python
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k = 2)
```


```python
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
for i, (match1, match2) in enumerate(matches):
    if match1.distance <0.75*match2.distance:
        matchesMask[i] = [1,0]
        
draw_params = dict(matchColor = (0,255,0), 
                  singlePointColor = (255,0,0),
                  matchesMask = matchesMask,
                  flags = 0)
```


```python
#Showing all features identifyied and making single line color
flann_matches = cv2.drawMatchesKnn(fruity_pebbles, kp1, cereals, kp2, matches, None, **draw_params)
display(flann_matches)
```


![png](output_24_0.png)



```python

```


## Object Ditection

```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
```


```python
%matplotlib inline 
```


```python
full = cv2.imread('Training_Rose')
```


```python
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7f92467db710>




![png](output_6_1.png)



```python
test = cv2.imread('Testing_Rose')
```


```python
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7f92467500d0>




![png](output_9_1.png)



```python
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (161, 171, 3)
    Training image shape: (183, 275, 3)



```python
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
```


```python
for m in methods:
    
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
        
    else:
        top_left = max_loc
        
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255,0,0),10)
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title('Detection of template')
    
    plt.suptitle(m)
    
    plt.show()
    print('\n')
    print('\n')
```


![png](output_12_0.png)


    
    
    
    



![png](output_12_2.png)


    
    
    
    



![png](output_12_4.png)


    
    
    
    



![png](output_12_6.png)


    
    
    
    



![png](output_12_8.png)


    
    
    
    



![png](output_12_10.png)


    
    
    
    



```python

```
