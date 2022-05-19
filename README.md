# Evolutionary-rate
How to calculate dN/dS using PAML software

## download and install software

1. get PAML software

       wget http://abacus.gene.ucl.ac.uk/software/
       
   NOTE: the macOS version does not work (paml4.8a.macosx.tgz), so use the unix version for mac: paml4.9j.tgz    
      
2. install paml software:
  
       tar -xf paml4.9j.tgz
       cd paml4.9j
       rm bin/*.exe
       cd src
       make -f Makefile
       ls -lF
       rm *.o
       mv baseml basemlg codeml pamp evolver yn00 chi2 ../bin
       cd ..
       ls -lF bin
       bin/baseml
       bin/codeml
       bin/evolver

3. get pal2nal software. 

   You need this to map nucleotide/cds sequences onto the protein alignment. This greatly reduces misalignment of DNA sequences.
   
   http://www.bork.embl.de/pal2nal/#Download
   
   tar file: pal2nal.v14.tar.gz
   
4. untar pal2nal

       tar -xzf pal2nal.v14.tar.gz
       
5. change directory to where perl script is then try running to see if works

        cd pal2nal.v14
        perl pal2nal.pl 

   Should show help and how to run
   
## build cds tree to run analysis



## dN/dS analysis pairwise analysis

3. copy control file (yn00.ctl) to working directory (should be in paml directory)

   yn00.ctl:

       seqfile = TMP.NT.FA.PHY * sequence data file name
       outfile = TMP.OUT        * main result file
       verbose = 0  * 1: detailed output (list sequences), 0: concise output

       icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

       weighting = 0  * weighting pathways between codons (0/1)?
       commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
       *       ndata = 1


       * Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
       * 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
       * 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
       * 10: blepharisma nu.
       * These codes correspond to transl_table 1 to 11 of GENEBANK.

4. get dna sequences for genes

  i. wget or download from website cds sequence
  
  ii. can simplify names using FastaManager.py (if needed)
  
        for i in *.fa; do python ~/Desktop/post_doc/scripts/FastaManager.py -f simplify -fasta $i; done
  
  iii. check gene names match list- not all cds gene names match protein gene names! If this is the case, use gff file to convert. if cds sequence is not available, use gff and genome fasta:
  
   download genome, gff from NCBI
   
       i.     wget 
       ii.	gunzip *.gz
       
   convert gff to coords
   
       python FastaManager.py -f gff_cds_to_coord -gff <gff file>
       
   get sequence from coords
   
       python FastaManager.py -f get_stretch4 -fasta <genomic fasta file> -coords <coord file from previous step>
       
   if the gff is not available, you may have to use BLAST and find recipricol best match of transcript/protein file
  
  iv.combine all modified fastas
  
        cat *.mod.fa > all_new_cds.fa
        
  v. get gene sequences
  
        python ~/Desktop/post_doc/scripts/FastaManager.py -f getseq2 -fasta all_new_cds_fasta.fa -name genelist.txt

5. get gene pair list; modified to exclude self-pairs

  i. get all by all gene pairs from a gene list
  
        python get_all_genepair_from_list.py <list of genes>
        
          
