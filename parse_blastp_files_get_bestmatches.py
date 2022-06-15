import os, sys

start_dir = sys.argv[1] #folder with BLAST file outputs

#function to add BLAST data to dictionary
def add_BLdata_to_dict(inp, D):

    for line in inp:
        if line.startswith("#"):
            pass
        else:
            L = line.strip().split("\t")
            #print (L)
            gene = L[0]
            #gene1 = Atgene.split(".")[0]
            #gene1a = gene1.split("_")[1]
            Hmgene = L[1]
            if Hmgene != "chromosome":
                percentsim = L[2]
                evalue = L[10]
                if float(evalue) <= 1e-6:
                    if gene not in D.keys():
                         D[gene] = [Hmgene, percentsim]
                    else:
                         first_PS = D[gene][1]
                         if float(percentsim) > float(first_PS):
                              D[gene] = [Hmgene, percentsim]
                         else:
                              pass


            
#loop through directory for each file to add input and each filename
title_list = []
for file in os.listdir(start_dir):
    if file.endswith(".out") or file.endswith("_results.txt"): #get BLAST file
        name = file.strip().split(".txt")[0]
        print (name)
        title_list.append(name)
        inp = open(start_dir + "/" + file, 'r')
        D = {}
        add_BLdata_to_dict(inp, D)
        sum_matrix = open(str(name)+"_bestmatch.txt", "w")
        sum_matrix.write('gene\tBLAST_bestmatch\tPercent_Sim\n')
        print (D)
        for gene in D:
            datalist= D[gene]
            datastr= "\t".join(datalist)
            sum_matrix.write('%s\t%s\n' %(gene, datastr))
        sum_matrix.close()
        inp.close()
