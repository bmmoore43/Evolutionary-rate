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
*NOTE: for FastaManager.py and gff_manager.py please see Utilities repository:
https://github.com/bmmoore43/Utilities

1. Make CDS tree by first aligning cds sequences to protein sequence
2. obtain cds sequences:
  i. wget or download from website cds sequence
  
  ii. can simplify names using FastaManager.py (if needed)
  
        for i in *.fa; do python ~/Desktop/post_doc/scripts/FastaManager.py -f simplify -fasta $i; done
        
3. Check that cds names match protein names
       
   i. make file with list of protein names
   ii. use gff file to find protein name and conversion to transcript/cds name
   
       python gff_manager.py -f filter -i <gff_file> -genes <protein_names.txt>
       
   Result: gff_filt file that contains only lines with gene protein name
   
   IF NEEDED:
   *NOTE: for installing/ using BLAST, see https://github.com/bmmoore43/Building-gene-trees
   
   
   iii. use BLAST to match protein file to transcript file:
   
   make BLAST database using protein file:
   
       ncbi-blast-2.10.0+/bin/makeblastdb -in <protein_fasta> -dbtype prot
       
   do BLASTX to blast nucleotide vs protein database
   
       ncbi-blast-2.10.0+/bin/blastx -num_threads <# of computer nodes> -db <protein_fasta> -query <nucle_fasta> -out <output_results.txt> -outfmt 6
       
   iv. parse BLAST results to get recipricol best match
   
       python parse_blastp_files_get_bestmatches.py <directory with BLAST output files ending in .out or _results.txt>
       
4. get all cds sequences

   i. with new gene list, get sequences of all cds genes. First concatenate all cds fastas together. Then run FastaManager.py on the combined fasta to get genes based on your gene list.
   
       cat *.fa > all.fa
       
       python FastaManager.py -f getseq2 -name <gene_list> -fasta all.fa
       
5. convert protein fasta file gene names to cds names so they match

   i. rename your protein fasta where the gene_list file format is: [new_name][old_name]
   
       python FastaManager.py -f rename -fasta <fasta_file> -name <gene_list>
       
6. Align protein file using mafft:
   *NOTE for installing and using MAFFT see: https://github.com/bmmoore43/Building-gene-trees

       mafft --auto --anysymbol <protein_fasta> > <output_name>
       
## use PAL2NAL to align cds sequences to protein alignment
1. In folder with your protein alignment and cds fasta, run pal2nal. Be sure to reference where the pal2nal program is in order to run.

       perl ~/pal2nal.v14/pal2nal <prot_alignment> <cds_fasta> -output fasta > cds_fasta.aln

   -output: output format- fasta is best for making trees afterward
   > cds_fasta.aln: concatenates the output into a file. This can be named anything, but typically I just add .aln to the fasta file name.
   NOTE: if certain genes won't align, check their alignment to the protein gene using tBLASTn. Some may cds genes have small percent ID because of large gaps/ incomplete transcript. Remove these from analysis. If match it very high (~99% ID) it is likely the cds gene is just reversed. Calculate reverse complement and rerun.

## make cds tree 
1. follow raxml-ng instructions to make the tree, see: https://github.com/bmmoore43/Building-gene-trees

## dN/dS models using PAML
These models compare a clade in the tree to all other leaves (genes) in the tree. The models determine the probablity of positive selection on an amino acid residue for that particular clade.

1. designate tree node where you want model to be. You can do this manually (open tree file and edit) or you can use EasyCodeML to label the tree.
2. if labeling manually: label with $1 or #1. Label needs to be after bootstrap value on tree, but before branch length. Example:

       ((Pan_Hal_v3.1_Pahal.1G105200.1:0.047243,(Sevir.1G019300.1:0.003058,Seita.1G019400.1:0.000001)99:0.053396)99:0.038588,(Zm00008a021549_T01:0.071475,Sobic.004G108200.1:0.050152)99:0.051612)91 #1:0.085421,CALSI_Maker00034354:0.178851)98:0.125517);
   
3. if using EasyCodeML:
   download:
       https://github.com/BioEasy/EasyCodeML
   run:
       java -jar EasyCodeML.jar
       
   when GUI pops up, load tree and alignment file, then press label. When finished export tree.

4. make CODEML.ctl file. NOTE: there are many options here and this is just an example. For explanations of options see: http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf. Need to update control file with your alignment file name, tree file name, and output file name.
       
       seqfile = FNSII_clade_monocot_cds.fa.p2n.aln * sequence data filename
       treefile = FNSII_clade_monocot_cds.fa.p2n.aln.raxml.supportFBP_model1     * tree structure file name- must be labeled!
       outfile = FNSII_clade_monocot_cds.fa.p2n_model1_dnds.txt           * main result file name

       noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
       verbose = 1  * 0: concise; 1: detailed, 2: too much
       runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

       seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
       CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

       * ndata = 10
       clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
       aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

       model = 2
       * models for codons:
         * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
       * models for AAs or codon-translated AAs:
         * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
         * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

       NSsites = 2  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0
                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

        fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
            kappa = 2  * initial 	or fixed kappa
        fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
            omega = 1 * initial or fixed omega, for codons or codon-based AAs

        fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
            alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
        Malpha = 0  * different alphas for genes
            ncatG = 8  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
        RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

        Small_Diff = .5e-6
            cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
        fix_blength = 2  * 0: ignore, -1: random, 1: initial, 2: fixed
            method = 0  * Optimization method 0: simultaneous; 1: one branch a time

        * Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
        * 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
        * 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
        * 10: blepharisma nu.
        * These codes correspond to transl_table 1 to 11 of GENEBANK.


5. Run Codeml. NOTE: codeml.ctl file needs to be in the same directory as the alignmnet and the tree files. You then have to run codeml from where the program is on your computer.

       /Users/Beth/Desktop/Github/gene_rate/paml4.9j/bin/codeml

## OLD analysis
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
        
          
