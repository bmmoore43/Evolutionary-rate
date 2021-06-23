# Evolutionary-rate
How to calculate dN/dS using PAML software

## dN/dS analysis- Note this is not updated

1. get PAML software

       wget http://abacus.gene.ucl.ac.uk/software/paml4.8a.macosx.tgz
      
2.get clustal software

       wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64

3.to run python AlnUtility.py, need gene DNA sequence and gene pair list

  d.get all by all gene pairs
  
        python get_all_genepair_from_list.py<list of ADH homologs>

4. install paml software:
  
        tar xf paml4.9j.tgz
        cd paml4.9j
        rm bin/*.exe
        cd src
        make -f Makefile
        ls -lF
        rm *.o
        mv baseml basemlg codeml pamp evolver yn00 chi2 ../bin
        cd ../
        ls -lF bin
        bin/baseml
        bin/codeml
        bin/evolver
        cd ../
        rm paml4.9j.tgz
        
5. install clustalw

        mkdir clustal1.82
        cd clustal1.82/
        wget http://www.clustal.org/download/1.X/ftp-igbmc.u-strasbg.fr/pub/ClustalX/clustalx1.82.linux.tar.gz
        tar zxvf clustalx1.82.linux.tar.gz

6. get dna sequences for genes

  i. wget or download from website cds sequence
  
  ii. simplify names using FastaManager.py
  
        for i in *.fa; do python ~/Desktop/post_doc/scripts/FastaManager.py -f simplify -fasta $i; done
  
  iii.check gene names match list
  
  iv.combine all modified fastas
  
        cat *.mod.fa > all_new_cds.fa
        
  v. get gene sequences
  
        python ~/Desktop/post_doc/scripts/FastaManager.py -f getseq2 -fasta all_new_cds_fasta.fa -name joiv_ADH_orthogroup_genelist_simple.txt

7. get gene pairs with simplified gene names; modify to exclude self-pairs

        python get_all_genepair_from_list.py joiv_ADH_orthogroup_genelist_simple.txt
        
8. copy control file into working directory

        control file: seqfile = TMP.NT.FA.PHY * sequence data file nameoutfile = TMP.OUT        * main result fileverbose = 0  * 1: detailed output (list sequences), 0: concise outputicode = 0  * 0:universal code; 1:mammalian mt; 2-10:see belowweighting = 0  * weighting pathways between codons (0/1)?commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? *       ndata = 1* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., * 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., * 10: blepharisma nu.i.

9. load old python module- need Python/2.7.9- use anaconda environments
  
  
10. need to copy the following into directory with AlnUtility.py

         cp /mnt/home/john3784/Github/Utilities/FastaManager.py /mnt/home/john3784/Github/SM-gene_prediction_Slycopersicum/
         cp /mnt/home/john3784/Github/Utilities/Translation.py /mnt/home/john3784/Github/SM-gene_prediction_Slycopersicum/
         cp /mnt/home/john3784/Github/Utilities/FileUtility.py /mnt/home/john3784/Github/SM-gene_prediction_Slycopersicum/
         cp /mnt/home/john3784/Github/Utilities/ParseBlast.py /mnt/home/john3784/Github/SM-gene_prediction_Slycopersicum/
         cp /mnt/home/john3784/Github/Utilities/BlastUtility.py /mnt/home/john3784/Github/SM-gene_prediction_Slycopersicum/
         cp /mnt/home/john3784/Github/Utilities/SingleLinkage.py /mnt/home/john3784/Github/SM-gene_prediction_Slycopersicum/
         cp /mnt/home/shius/codes/TreeUtility.py* ~/Github/Utilities/
         cp /mnt/home/shius/codes/DrawBlocks.py* ~/Github/Utilities/
         cp /mnt/home/shius/codes/SVGWriter.py* ~/Github/Utilities/
         cp /mnt/home/shius/codes/kaksToolsHacked.py* ~/Github/Utilities/
         cp /mnt/home/shius/codes/DatabaseOp.py* ~/Github/Utilities/
         
  11. need to copy rst file into paml/bin directory
  
          cp ~/paml4.9j/rst1 ~/paml4.9j/bin/
          
  12. run kaks script
    
          python~/Github/SM-gene_prediction_Slycopersicum/AlnUtility.py -f rate_pair -pep <protein fasta> -cds <cds fasta> -p paml -paml ~/paml4.9j/ -clustal ~/clustal/clustal1.82/clustalx1.82.linux/ -pairs <gene pair file>
          
