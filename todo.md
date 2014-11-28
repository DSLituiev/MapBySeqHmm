
**prior knowledge**:

1. Blast pieces around SNPs for to check conservation

- extract pieces with `bedtools`

   $ cat test.fa
    >chr1
    AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

    $ cat test.bed
    chr1 5 10

    $ bedtools getfasta -fi test.fa -bed test.bed -fo test.fa.out

    $ cat test.fa.out
    >chr1:5-10
    AAACC

OR with `samtools`

    samtools faidx <ref.fasta> [region1 [...]]
    

- run `blast` locally with `-remote` databases

    blastp -query my1000seqs.fasta -db nr -remote

    blastn -remote -db nr -query stuff_to_search.fasta -entrez_query "Homo Sapiens[Organism]" -evalue 1e-20 -num_alignments 10 > stuff_to_search.blastn

[BLAST output formatting](http://www.ncbi.nlm.nih.gov/books/NBK1763/#CmdLineAppsManual.Quick_start)
