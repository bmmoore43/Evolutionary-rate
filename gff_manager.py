"""
Input file, -i

Functions, -f:
	sort		Sort a GFF file
			output: [gff_file].sort
	merge		Merge features in a sorted GFF file
			output: [gff_file].merge
	sort_merge	Sort and then merge a GFF file
			output: [gff_file].sort.merge
	merge_depth	Count the number of features clustered in each line of a merged
		GFF file
			output: [gff_file].depth
	filter	Remove lines where gene and/or motif not in list provided.
		req: 	-i (input gff file)
		opt:	-genes and/or -str
			output: [gff_file].filt
	length	Return the total length and length of features covered in the
		GFF file
			output: print to screen
	lengths	Print the ID and length of features in a GFF file. Include a
		-str flag to only include the lengths of that feature
			opt: -str,-com
			output: print to screen
	min_len	Keep only features longer than a minimum length (-int). Include a 
		-type flag to only keep features of that type.
			req: -int
			opt: -type,-com
			output: [gff_file].[int]_min_len
	len_prcntl	Print the median values and 95th and 99th percentiles of 
		feature lengths in a GFF file. Feature types will be calculated
		separately. Provide an integer value (-int) to calculate a 
		different percentile.
			opt: -int
			output: print to screen
			SciPy module required (module load SciPy)
	compare_lens	Compare the distributions of lengths in two gff files (-i, 
		-gff2) using a Mann Whitney U test. Include a -str flag to only
		test the lengths of a particular feature type.
			req: -gff2; opt: -str
			output: print to screen
			SciPy module required (module load SciPy)
	prefix_id	Adds a prefix to the id in the "ID=" section of description
			req: -str
			output: print to screen
	loc_id	Output location-based identifiers for each gff feature.
		Identifiers take the form of "seq|start-end"
			output: [gff_file].loc_id
	loc_id_keep	Keep feaures based on a list of loc_ids (-str; seq|start-end).
			req: -str
			opt: -out_str
			output: [gff_file].loc_id_keep OR [gff_file].[out_str].loc_id
	loc_id2gff	Convert a 1col file with location IDs into a gff file. Include 
		a -source or -type flag to change the source and type
		fields in the resulting gff file.
			opt: -source,-type
			output: [loc_id_file].gff
	coords2gff	Convert a 2col, 3col, or 4col coords file into a GFF file. 
		Include a -source or -type flag to change the source and type
		fields in the resulting gff file. The script will auto-detect 
		2, 3, or 4-col format
		 2col format: col1=seq/chr col2=start,end
		 3col format: col1=seq/chr col2=start col3=end
		 4col format: col1=seq/chr col2=start col3=end col4=region_name
			opt: -source,-type
			output: [coords_file].gff
	mapped2gff	Convert a mapping file into a GFF file using a coordinate key. 
								req: -i -key
								opt: -source,-type
								output: [coords_file].gff
	mapped2gff_internalkey	Convert a mapping file into a GFF when the first column
								in the mapping file isn't a gene that needs to be converted into coords, but 
								is the coords themselves. 
								req: -i
								opt: -source,-type
								output: [coords_file].gff
	mask		Mask input gff file with a secondary gff file
			req: -gff2
			output: print to screen
	overlap	Identify GFF features that overlap with a secondary GFF file 
		(-gff2); Feature type and location from secondary feature file
		will be added to description in primary file in an "overlap=" 
		section; Include a character string (-type) to only find 
		overlaps with that feature type in the secondary file
			req: -gff2; opt: -str,-type
			output: [infile].[str_]ovrlp
	overlap+	Describe the feature, length, start and end points, and
		proportion of total of overlap regions between a primary (-i)
		and secondary (-gff2) gff file.
			req: -gff2
			output: print to screen

Additional flags:
	-str		character string
	-int		integer
	-gff2		secondary gff file
	-com		comment character (default = #)
	-source	string source input
	-type		string type input
	-out_sffx	output file suffix

	"""
def print_help():
	print(__doc__)

def set_defaults_and_parse_args(argv_l):
	function = input_file = char_str = gff2 = src = typ = intgr = key = genes = o_str = None
	cm_chr = "#"
	import sys
	try:
		for i in range(0,len(sys.argv)):
			if sys.argv[i] == "-i":
				input_file = sys.argv[i+1]
			elif sys.argv[i] == "-f":
				function = sys.argv[i+1]
			elif sys.argv[i] == "-str":
				char_str = sys.argv[i+1]
			elif sys.argv[i] == "-gff2":
				gff2 = sys.argv[i+1]
			elif sys.argv[i] == "-source":
				src = sys.argv[i+1]
			elif sys.argv[i] == "-type":
				typ = sys.argv[i+1]
			elif sys.argv[i] == "-int":
				intgr = sys.argv[i+1]
			elif sys.argv[i] == "-key":
				key = sys.argv[i+1]	
			elif sys.argv[i] == "-genes":
				genes = sys.argv[i+1]	
			elif sys.argv[i] == "-out_str":
				o_str = sys.argv[i+1]
		return input_file,function,char_str,gff2,cm_chr,src,typ,intgr,key,genes,o_str
	except:
		print_help()
		print("Error reading arguments!")
		sys.exit()

def gff2dict(gff_fl,cmm="#",strng=None):
	d = {}
	inp = open(gff_fl)
	for line in inp:
		if not line.startswith(cmm):
			lineLst = line.strip().split("\t")
			lineLst[3] = int(lineLst[3])
			lineLst[4] = int(lineLst[4])
			seq = lineLst[0]
			feat = lineLst[2]
			if strng == None or feat == strng:
				if seq not in d:
					d[seq] = [lineLst]
				else:
					d[seq].append(lineLst)
	inp.close()
	return d

def make_loc_id(gff_lnLst):
	sq = gff_lnLst[0]
	st = gff_lnLst[3]
	en = gff_lnLst[4]
	lid = "%s|%s-%s"%(sq,st,en)
	return lid

def pull_info_from_description(desc_ln,prb):
	desc_lst = desc_ln.split(";")
	info = "NA"
	for item in desc_lst:
		if prb in item:
			info = item.replace(prb,"")
			break
	return info

def pull_lengths(gff_fl,strng):
	len_l = []
	inp = open(gff_fl)
	for line in inp:
		lineLst = line.strip().split("\t")
		feat = lineLst[2]
		if strng == None or feat == strng:
			start = int(lineLst[3])
			end = int(lineLst[4])
			reg_len = end-start+1
			len_l.append(reg_len)
	inp.close()
	return len_l

def seq_max(gff,d):
	inp = open(gff)
	for line in inp:
		lineLst = line.split("\t")
		seq = lineLst[0]
		end = lineLst[4]
		if seq not in d or int(end) > d[seq]:
			d[seq] = int(end)
	inp.close()

def roundup(x,by):
	round_down = int(str(float(x)/float(by)).split(".")[0])*int(by)
	round_up = round_down + int(by)
	return round_up

def rounddown(x,by):
	round_down = int(str(float(x)/float(by)).split(".")[0])*int(by)
	return round_down

def indices_between_values(low,high,span):
	ind_l = []
	min_vals = range(low,high,span)
	for min_val in min_vals:
		start_ind = min_val
		end_ind = start_ind+span-1
		index="%s-%s"%(start_ind,end_ind)
		ind_l.append(index)
	return ind_l

def index_seqs(seq_max_d,span):
	d = {}
	for seq in seq_max_d:
		d[seq]={}
		max = seq_max_d[seq]
		max_ceil = roundup(max,span)
		if max == max_ceil:
			index_list = indices_between_values(0,max_ceil+1,span)
		else:
			index_list = indices_between_values(0,max_ceil,span)
		for index in index_list:
			d[seq][index] = set()
	return d

def coord_index_overlap(gff_line,span):
	lineLst = gff_line.strip().split("\t")
	seq = lineLst[0]
	start = lineLst[3]
	end = lineLst[4]
	
	ind1_start = rounddown(start,span)
	ind1_end = ind1_start+span-1
	
	ind2_start = rounddown(end,span)
	ind2_end = ind2_start+span-1
	
	index1="%s-%s"%(ind1_start,ind1_end)
	index2="%s-%s"%(ind2_start,ind2_end)
	
	id_l = []
	id_l.append(index1)
	if index1 != index2:
		id_l.append(index2)
		if ind1_end+1 != ind2_start:
			mid_index_list = indices_between_values(ind1_end+1,ind2_start-1,span)
			for index in mid_index_list:
				id_l.append(index)
	return id_l

def add_gff_to_index(g,ind_d,span):
	inp = open(g)
	for line in inp:
		lineLst = line.split("\t")
		lineTup = tuple(lineLst)
		seq = lineLst[0]
		all_indices = coord_index_overlap(line,span)
		for index in all_indices:
			ind_d[seq][index].add(lineTup)
	inp.close()

def same_keys_list(g1,g2):
	l = []
	for key in g1.keys():
		if key in g2.keys():
			l.append(key)
	return l

def test_overlap(st1,en1,st2,en2):
	s1 = int(st1)
	e1 = int(en1)
	s2 = int(st2)
	e2 = int(en2)
	overlap = False
	if s2 >= s1 and s2 <= e1:
		overlap = True
	elif s1 >= s2 and s1 <= e2:
		overlap = True
	return overlap

def generate_overlap_dict(in_gff,gff2_d):
	ovrlp_d = {}
	inp = open(in_gff)
	for line in inp:
		lineLst = line.strip().split("\t")
		lineLst_tup = tuple(lineLst)
		ovrlp_d[lineLst_tup] = []
		seq = lineLst[0]
		start = lineLst[3]
		end = lineLst[4]
		if seq in gff2_d:
			comparison_lines = gff2_d[seq]
			for comp_line in comparison_lines:
				start2 = comp_line[3]
				end2 = comp_line[4]
				ovrlp = test_overlap(start,end,start2,end2)
				if ovrlp == True:
					ovrlp_d[lineLst_tup].append(comp_line)
	inp.close()
	return ovrlp_d

def overlap_dict_by_indexing(gf1,gf2,span):
	seq_max_dict = {}
	seq_max(gf1,seq_max_dict)
	seq_max(gf2,seq_max_dict)
	# print seq_max_dict
	seq_index_d = index_seqs(seq_max_dict,span)
	add_gff_to_index(gf2,seq_index_d,span)
	
	ovrlp_d = {}
	inp = open(gf1)
	for line in inp:
		lineLst = line.strip().split("\t")
		lineTup = tuple(lineLst)
		if lineTup in ovrlp_d:
			pass
		else:
			ovrlp_d[lineTup] = []
			
			seq = lineLst[0]
			start = lineLst[3]
			end = lineLst[4]
			
			all_indices = coord_index_overlap(line,span)
			comparison_lines_set = set()
			for index in all_indices:
				try:
					for comp_line in seq_index_d[seq][index]:
						comparison_lines_set.add(comp_line)
				except:
					print(lineLst)
					print(seq,index)
					print(all_indices)
			
			for comp_line in comparison_lines_set:
				start2 = comp_line[3]
				end2 = comp_line[4]
				ovrlp = test_overlap(start,end,start2,end2)
				if ovrlp == True:
					ovrlp_d[lineTup].append(comp_line)
	inp.close()
	return ovrlp_d

def new_coords(gff_line,new_st,new_en):
	new_s = int(new_st)
	new_e = int(new_en)
	new_line = ""
	if new_en-new_st<0:
		pass
	else:
		new_line = gff_line[:]
		new_line[3] = str(new_s)
		new_line[4] = str(new_e)
	return new_line

def adjust_gff_limits(gff1_lineLst,gff2_lineLst):
	# print "HERE:",gff1_lineLst
	# print "HERE:",gff2_lineLst
	s1 = int(gff1_lineLst[3])
	e1 = int(gff1_lineLst[4])
	s2 = int(gff2_lineLst[3])
	e2 = int(gff2_lineLst[4])
	#check if primary is encompassed by secondary
	if s2 <= s1 and e2 >= e1:
		return ["",""]
	#check if secondary is encompassed by primary
	elif s1 <= s2 and e1 >= e2:
		part1 = new_coords(gff1_lineLst,s1,s2-1)
		part2 = new_coords(gff1_lineLst,e2+1,e1)
		return [part1,part2]
	#check if primary starts first, keep front end of primary
	elif s1 < s2:
		new_lineLst = new_coords(gff1_lineLst,s1,s2-1)
		return [new_lineLst,""]
	#check if secondary starts first, keep back end of primary
	elif s2 < s1:
		new_lineLst = new_coords(gff1_lineLst,e2+1,e1)
		return ["",new_lineLst]
	else:
		print("WHAT ELSE IS THERE?",gff1_lineLst,gff2_lineLst)

def region_of_overlap(gff1_line,gff2_line):
	seq0 = gff1_line[0]
	seq1 = gff2_line[0]
	if seq0 != seq1:
		print("Sequences are not identical!:",seq0,seq1)
	
	start0 = int(gff1_line[3])
	end0 = int(gff1_line[4])
	total_length = end0-start0+1
	frag_nm = "%s|%s-%s"%(seq0,start0,end0)
	type = gff1_line[2]
	
	start1 = int(gff2_line[3])
	end1 = int(gff2_line[4])
	feature = gff2_line[2]
	
	if start1 <= start0 and end1 >= end0: #check if txfrag is encompassed by feature
		region_start = start0
		region_end = end0
	elif start0 <= start1 and end0 >= end1:	#check if feature is encompassed by txfrag
		region_start = start1
		region_end = end1
	elif start0 < start1: #check if txfrag starts first, keep back end of frag
		region_start = start1
		region_end = end0
	elif start1 < start0: #check if feature starts first, keep front end of frag, should be last option
		region_start = start0
		region_end = end1
	else:
		print("What's left?")
		print(gff1_line)
		print(gff2_line)
	region_length = region_end-region_start+1
	percent_of_total = float(region_length)/float(total_length)*100
	# oof.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(frag_nm,region_start,region_end,feature,region_length,total_length,percent_of_total))
	out_str = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(frag_nm,type,region_start,region_end,feature,region_length,total_length,percent_of_total)
	print(out_str)

def function_sortGFF(fl,cmmnt,input_str):
	out = open(fl+".sort","w")
	if cmmnt == "#" and input_str == None:
		dict = gff2dict(fl)
	elif cmmnt != "#" and input_str == None:
		dict = gff2dict(fl,cmmnt,None)
	elif cmmnt == "#" and input_str != None:
		dict = gff2dict(fl,"#",input_str)
	else:
		dict = gff2dict(fl,cmmnt,input_str)
	
	for key in dict:
		lines_list = dict[key]
		lines_list.sort(key=lambda k: (k[3]))
		for line in lines_list:
			out_list = []
			for item in line:
				out_list.append(str(item))
			out_line = "\t".join(out_list)
			out.write(out_line+"\n")
			# print out_line
	out.close()

def write_mergeGFF(l,seq,oof):
	str_list = []
	for item in l:
		str_list.append(str(item))
	out_str = seq+"\t"+"\t".join(str_list)
	oof.write(out_str+"\n")
	# print out_str

def merge_lists(ent_l,ky,out_f):
	out_ls = [ent_l[0]]
	for i in range(1,len(ent_l)):
		start = int(ent_l[i][2])
		end = int(ent_l[i][3])
		start_prev = int(out_ls[-1][2])
		end_prev = int(out_ls[-1][3])
		curr_ID = ent_l[i][-1]
		curr_ID
		if start <= end_prev:
			if end > end_prev:
				out_ls[-1][3] = end
			out_ls[-1][-1] = out_ls[-1][-1]+"--merged--"+curr_ID
			feat = ent_l[i][1]
			prev_feat = ent_l[i-1][1]
			if feat == prev_feat:
				out_ls[-1][1] = "%s"%feat
			elif feat != prev_feat:
				flst = [feat,prev_feat]
				flst.sort()
				mrg_nm = "%s_%s"%(flst[0],flst[1])
				out_ls[-1][1] = mrg_nm
		else:
			
			out_ls.append(ent_l[i])
	for list in out_ls:
		write_mergeGFF(list,ky,out_f)

def function_mergeGFF(fl):
	fl_dict = {}
	inp = open(fl)
	for line in inp:
		if not line.startswith("#"):
			lineLst = line.strip().split("\t")
			seq = lineLst[0]
			start = lineLst[3]
			end = lineLst[4]
			IDinfo = lineLst[-1]
			add_to_info="seq_name=%s;region_start=%s;region_end=%s"%(seq,start,end)
			if IDinfo.endswith(";"):
				newIDinfo = IDinfo+add_to_info
			else:
				newIDinfo = IDinfo+";"+add_to_info
			lineLst[-1] = newIDinfo
			if seq not in fl_dict:
				fl_dict[seq] = [lineLst[1:]]
			else:
				fl_dict[seq].append(lineLst[1:])
			
	inp.close()
	
	out = open(fl+".merge","w")
	for key in fl_dict:
		entry_lists = fl_dict[key]
		if len(entry_lists) == 1:
			write_mergeGFF(entry_lists[0],key,out)
		else:
			merge_lists(entry_lists,key,out)
	out.close()

# def featLenDict_and_fullLen(gf):
	# inp = open(gf)
	# all_len = 0
	# d_len = {}
	# for line in inp:
		# lineLst = line.split("\t")
		# feat = lineLst[2]
		# start = int(lineLst[3])
		# end = int(lineLst[4])
		# line_len = end-start+1
		# if line_len < 0:
			# print "NEGATIVE LENGTH:",line.strip()
		# else:
			# all_len += line_len
			# if feat not in d_len:
				# d_len[feat] = line_len
			# else:
				# d_len[feat] += line_len
	# inp.close()
	# return d_len,all_len

def function_filter(gff,chr_str,genes):
	import os.path
	if chr_str == None:
		chr_str = 'pass'
	if genes == None:
		genes = 'pass'
	
	# Read in genes / strings to keep
	if os.path.isfile(genes):
		with open(genes) as gene_file:
			gene_list = gene_file.read().strip().splitlines()
	else:
		gene_list = [genes]
	print('Number of genes to keep: %i' % len(gene_list))

	if os.path.isfile(chr_str):
		with open(chr_str) as str_file:
			str_list = str_file.read().strip().splitlines()
	else:
		str_list = [chr_str]
	print('Number of genes to keep: %i' % len(str_list))

	count_all = 0
	count_keep = 0
	out = open(gff+'_filt', 'w')

	with open(gff) as gff_file:
	  for gff_line in gff_file:
	    count_all += 1
	    gff_split = gff_line.strip().split('\t')

	    gene = gff_split[1]
	    string = gff_split[2]
	    
	    if string in str_list or 'pass' in str_list:
	      if gene in gene_list or 'pass' in gene_list:
	        count_keep +=1
	        out.write(gff_line)

	    if count_all%10000 ==0:
	      print('Finished: %i' % count_all)

	print('Kept %i out of %i lines from the gff file' % (count_keep, count_all))



def function_length(gff):
	# d_ln,ln = featLenDict_and_fullLen(gff)
	inp = open(gff)
	ln = 0
	d_ln = {}
	for line in inp:
		lineLst = line.split("\t")
		feat = lineLst[2]
		start = int(lineLst[3])
		end = int(lineLst[4])
		line_len = end-start+1
		if line_len < 0:
			print("NEGATIVE LENGTH:",line.strip())
		else:
			ln += line_len
			if feat not in d_ln:
				d_ln[feat] = line_len
			else:
				d_ln[feat] += line_len
	inp.close()
	
	fl_nm = gff.split("/")[-1]
	run_ln = 0
	print("#file\tfeature\tlen\t%ofTot")
	for feat in d_ln:
		ft_len = d_ln[feat]
		percn = round(float(ft_len)/float(ln)*100,2)
		print("%s\t%s\t%s\t%s"%(fl_nm,feat,ft_len,percn))
		run_ln += ft_len
	print("%s\ttotal\t%s\t%s"%(fl_nm,ln,run_ln))

def function_prefixID(gff,str):
	inp = open(gff)
	for line in inp:
		lineLst = line.strip().split("\t")
		rpl_desc = lineLst[-1].replace("ID=","ID="+str)
		lineLst[-1] = rpl_desc
		print("\t".join(lineLst))

def function_maskGFF(gff,gff2,feat_type,out_str):
	if out_str == None:
		out = open(gff+".masked","w")
	else:
		out = open("%s.%s_masked"%(gff,out_str),"w")
	index_span = 10000
	overlap_d = overlap_dict_by_indexing(gff,gff2,index_span)
	for key in overlap_d:
		if len(overlap_d[key]) == 0:
			# print "\t".join(key)
			out.write("\t".join(key)+"\n")	
		elif len(overlap_d[key]) == 1:
			adjusted_lines = adjust_gff_limits(list(key),overlap_d[key][0])
			for line in adjusted_lines:
				if line != "":
					# print "\t".join(line)
					out.write("\t".join(line)+"\n")
		elif len(overlap_d[key]) >= 2:
			overlap_lines = []
			for lnTup in overlap_d[key]:
				lnLst = list(lnTup)
				lnLst[3] = int(lnLst[3])
				lnLst[4] = int(lnLst[4])
				overlap_lines.append(lnLst)
			overlap_lines.sort(key=lambda n: (n[3], -n[4]))
			curr_prim = list(key)
			for ovrlp_line in overlap_lines:
				if curr_prim == "":
					break
				else:
					adjusted_lines = adjust_gff_limits(curr_prim,ovrlp_line)
				if adjusted_lines[0] != "":
					out_l = []
					for item in adjusted_lines[0]:
						out_l.append(str(item))
					# print "\t".join(out_l)
					out.write("\t".join(out_l)+"\n")
				curr_prim = adjusted_lines[1]
			if adjusted_lines[1] != "":
				out_l = []
				for item in adjusted_lines[1]:
					out_l.append(str(item))
				# print "\t".join(out_l)
				out.write("\t".join(out_l)+"\n")
	out.close()

def function_overlap(gff,gff2,out_str,feat_type):
	index_span = 10000
	overlap_d = overlap_dict_by_indexing(gff,gff2,index_span)
	inp = open(gff)
	if out_str == None:
		out = open("%s.ovrlp"%(gff),"w")
	else:
		out = open("%s.%s_ovrlp"%(gff,out_str),"w")
	for line in inp:
		lineLst = line.strip().split("\t")
		lineTup = tuple(lineLst)
	# for key in overlap_d:
		# lineTup = key
		overlap_l = []
		for overlap_line in overlap_d[lineTup]:
			ovrlp_ft = overlap_line[2]
			ovrlp_st = overlap_line[3]
			ovrlp_en = overlap_line[4]
			ovrlap_desc = "%s:%s-%s"%(ovrlp_ft,ovrlp_st,ovrlp_en)
			overlap_l.append(ovrlap_desc)
		overlap_str = ",".join(overlap_l)
		if overlap_str == "":
			overlap_str = "NA"
		line_remade = "\t".join(lineTup)
		# print "%s;overlap=%s"%(line_remade,overlap_str)
		out.write("%s;overlap=%s\n"%(line_remade,overlap_str))
	out.close()
	inp.close()

def function_overlapPlus(gff,gff2):
	index_span = 10000
	overlap_d = overlap_dict_by_indexing(gff,gff2,index_span)
	print("#id\tregion_type\tovrlp_start\tovrlp_end\tovrlp_feat\tovrlp_len\ttotal_len\tper_of_total")
	inp = open(gff)
	for line in inp:
		lineLst = line.strip().split("\t")
		lineTup = tuple(lineLst)
		for overlap_line in overlap_d[lineTup]:
			region_of_overlap(lineTup,overlap_line)

def calc_med_len_ave(lst):
	md = numpy.median(lst)
	n = len(lst)
	av = round(sum(lst)/float(n),1)
	return md,n,av

def function_compareLens(gff,gff2,feat_type):
	gff1_len_l = pull_lengths(gff,feat_type)
	gff2_len_l = pull_lengths(gff2,feat_type)
	
	global numpy
	import numpy
	from scipy import stats
	m1,n1,a1 = calc_med_len_ave(gff1_len_l)
	m2,n2,a2 = calc_med_len_ave(gff2_len_l)
	pVal = stats.mannwhitneyu(gff1_len_l,gff2_len_l)[1]
	
	print("stat\tgff1\tgff2")
	print("med\t%s\t%s"%(m1,m2))
	print("n\t%s\t%s"%(n1,n2))
	print("ave\t%s\t%s"%(a1,a2))
	print("MWU p-value:",pVal)

def function_coords2gff(coords,source,type):
	if source == None:
		source = "coords"
	if type == None:
		type = "coords"
	import os
	full_pth_crds = os.path.abspath(coords)
	gff_line = "%s\t%s\t%s\t%s\t%s\t.\t.\t.\tID=%s;coords_file=%s"
	inp = open(coords)
	out = open(coords+".gff","w")
	for line in inp:
		lineLst = line.strip().split("\t")
		if len(lineLst) == 2:
			seq,coords = lineLst
			start,end = coords.split(",")
			id = "%s|%s-%s"%(seq,start,end)
		elif len(lineLst) == 3:
			seq,start,end = lineLst
			id = "%s|%s-%s"%(seq,start,end)
		elif len(lineLst) == 4:
			seq,start,end,id = lineLst
		out_line = gff_line%(seq,source,type,start,end,id,full_pth_crds)
		out.write(out_line+"\n") 
	out.close()
	inp.close()

def function_mapped2gff(coords, key, source, type):
	if source == None:
		source = "coords"
	if type == None:
		type = "coords"
	import os
	full_pth_crds = os.path.abspath(coords)
	
	#Make a dictionary with all coordinates from key
	coord_dic = {}
	for l in open(key, 'r'):
			line = l.strip().split("\t")
			AT = line[3]
			chromo = l[3:4]
			promoter_coords = str(chromo) + "|" + str(line[1]) + "|" + str(line[2])
			coord_dic[AT] = promoter_coords

	gff_line = "%s\t%s\t%s\t%s\t%s\t.\t%s\t.\tID=%s;coords_file=%s"
	inp = open(coords)
	out = open(coords+".gff","w")
	count = 0
	for line in inp:
		if not line[0].isdigit():
				seq,k_start,k_end,direction,kmer,score,pval = line.strip().split("\t")
				try:
					coords = coord_dic[seq]
					x = coords.strip().split("|")
					chr = "chr" + str(x[0])
					if direction == "1":
							d = "+"
					elif direction == "-1":
							d = "-"
					else:
							d = "."
					start = int(x[1]) + int(k_start)
					end = int(x[1]) + int(k_end)
					loc_id = "%s|%s-%s"%(seq,start,end)
					out_line = gff_line%(chr,seq,kmer,start,end,d,loc_id,full_pth_crds)
					out.write(out_line+"\n") 
				except:
								count = count + 1
	print("Genes that were not found in conversion file: " + str(count))
	out.close()
	inp.close()

def function_mapped2gff_internalkey(coords, source, type):
	if source == None:
		source = "coords"
	if type == None:
		type = "coords"
	import os
	full_pth_crds = os.path.abspath(coords)
	

	gff_line = "%s\t%s\t%s\t%s\t%s\t.\t%s\t.\tID=%s;coords_file=%s"
	inp = open(coords)
	out = open(coords+".gff","w")
	count = 0
	for line in inp:
		if not line[0].isdigit():
				seq,k_start,k_end,direction,kmer,score,pval = line.strip().split("\t")
				chro = seq.split("|")[0]
				start_stop = seq.split("|")[1]
				start = start_stop.split("-")[0]

				if direction == "1":
						d = "+"
				elif direction == "-1":
						d = "-"
				else:
						d = "."
				kmer_start = int(start) + int(k_start)
				kmer_end = int(start) + int(k_end)
				loc_id = "%s|%s-%s"%(seq,kmer_start,kmer_end)
				out_line = gff_line%(chro,seq,kmer,kmer_start,kmer_end,d,loc_id,full_pth_crds)
				out.write(out_line+"\n") 

	out.close()
	inp.close()


def function_mapped2gff2(coords, key, source, type):
		import numpy as np
		if source == None:
			source = "coords"
		if type == None:
			type = "coords"
		import os
		full_pth_crds = os.path.abspath(coords)
		
		#Make a dictionary with all coordinates from key
		coord_dic = {}
		for l in open(key, 'r'):
					line = l.strip().split("\t")
					AT = line[3]
					chromo = l[3:4]
					promoter_coords = str(chromo) + "|" + str(line[1]) + "|" + str(line[2])
					coord_dic[AT] = promoter_coords
		gff_line = "%s\t%s\t%s\t%s\t%s\t.\t%s\t.\tID=%s;coords_file=%s"
		inp = open(coords)
		out = open(coords+".gff","w")
		count = 0
		for line in inp:
			if not line[0].isdigit() and not line.startswith('Seq'):
					seq,kmer,motif_seq,hit_seq,k_start,hit_score,threshold = line.strip().split("\t")
					seq = seq[:-3]
					dir = np.sign(int(k_start))
					if dir == 1:
							k_end = int(k_start) + len(hit_seq)
					if dir == -1:
							k_start = 1000 + int(k_start)
							k_end = k_start + len(hit_seq)
					
					try:
									coords = coord_dic[seq]
									x = coords.strip().split("|")
									chr = "chr" + str(x[0])
									if dir == 1:
											d = "+"
									elif dir == -1:
											d = "-"
									else:
											d = "."
									start = int(x[1]) + int(k_start)
									end = int(x[1]) + int(k_end)
									loc_id = "%s|%s-%s"%(seq,start,end)
									out_line = gff_line%(chr,seq,kmer,start,end,d,loc_id,full_pth_crds)
									out.write(out_line+"\n") 
					except:
									count = count + 1
		print("Genes that were not found in conversion file: " + str(count))
		out.close()
		inp.close()


def function_locId2gff(loc_ids,source,type):
	if source == None:
		source = "loc_id"
	if type == None:
		type = "loc_id"
	import os
	full_pth_crds = os.path.abspath(loc_ids)
	gff_line = "%s\t%s\t%s\t%s\t%s\t.\t.\t.\tID=%s;coords_file=%s"
	inp = open(loc_ids)
	out = open(loc_ids+".gff","w")
	for line in inp:
		loc_id = line.strip()
		seq,start_end = loc_id.split("|")
		start,end = start_end.split("-")
		out_line = gff_line%(seq,source,type,start,end,loc_id,full_pth_crds)
		out.write(out_line+"\n") 
	out.close()
	inp.close()

def function_lengths(gff,str,com):
	inp = open(gff)
	for line in inp:
		if not line.startswith(com):
			lineLst = line.strip().split("\t")
			seq = lineLst[0]
			type = lineLst[2]
			start = lineLst[3]
			end = lineLst[4]
			id = pull_info_from_description(lineLst[-1],"ID=")
			# loc_id = "%s|%s-%s"%(seq,start,end)
			loc_id = make_loc_id(lineLst)
			reg_len = int(end)-int(start)+1
			if str == None or type == str:
				print("%s\t%s\t%s\t%s"%(loc_id,id,type,reg_len))
	inp.close()

def function_lenPercentiles(gff,intgr,com):
	# ln_d,allLen = featLenDict_and_fullLen(gff)
	ft_lens_d = {}
	inp = open(gff)
	for line in inp:
		if not line.startswith(com):
			lineLst = line.strip().split("\t")
			ft = lineLst[2]
			start = lineLst[3]
			end = lineLst[4]
			reg_len = int(end)-int(start)+1
			if ft not in ft_lens_d:
				ft_lens_d[ft] = [reg_len]
			else:
				ft_lens_d[ft].append(reg_len)
	inp.close()
	
	if intgr != None:
		print("feat\tmed\t75per\t95per\t99per\t%sper\tmax"%(intgr))
	else:
		print("feat\tmed\t75per\t95per\t99per\tmax")
	import numpy
	for key in ft_lens_d:
		lens_l = ft_lens_d[key]
		fl_l = []
		for len in lens_l:
			fl_l.append(float(len))
		med = numpy.percentile(fl_l,50)
		per75 = numpy.percentile(fl_l,75)
		per95 = numpy.percentile(fl_l,95)
		per99 = numpy.percentile(fl_l,99)
		max_val = max(fl_l)
		if intgr != None:
			perX = numpy.percentile(fl_l,intgr)
			print("%s\t%s\t%s\t%s\t%s\t%s\t%s"%(key,med,per75,per95,per99,perX,max_val))
		else:
			print("%s\t%s\t%s\t%s\t%s\t%s"%(key,med,per75,per95,per99,max_val))
			
	# inp = open(gff)
	# ft_d = {}
	# for line in inp:
		
	# inp.close()

def function_minLen(gff,min_len,type,cmmnt):
	int_len = int(min_len)
	inp = open(gff)
	out = open("%s.%s_min_len"%(gff,min_len),"w")
	for line in inp:
		if not line.startswith(cmmnt):
			lnL = line.split("\t")
			if type == None or lnL[2] == type:
				start = int(lnL[3])
				end = int(lnL[4])
				reg_len = end-start+1
				if reg_len >= int_len:
					out.write(line)
	out.close()
	inp.close()

def function_mergeDepth(gff,com):
	inp = open(gff)
	out = open(gff+".depth","w")
	max_depth = 0
	for line in inp:
		if not line.startswith(com):
			lineLst = line.strip().split("\t")
			loc_id = make_loc_id(lineLst)
			desc = lineLst[-1]
			merge_cnt = desc.count("--merged--")
			depth = merge_cnt+1
			out.write("%s\t%s\n"%(loc_id,depth))
			if depth > max_depth:
				max_depth = depth
	out.close()
	inp.close()
	
	print("Maximum merge depth:",max_depth)

def function_locID(gff,com):
	inp = open(gff)
	out = open(gff+".loc_id","w")
	for line in inp:
		if not line.startswith(com):
			lineLst = line.split("\t")
			loc_id = make_loc_id(lineLst)
			out.write(loc_id+"\n")
	out.close()
	inp.close()

def function_locIDkeep(gff,com,str,out_str):
	id_set = set()
	id_lst_inp = open(str)
	for line in id_lst_inp:
		id = line.strip()
		id_set.add(id)
	id_lst_inp.close()
	
	inp = open(gff)
	if out_str == None:
		out = open(gff+".loc_id_keep","w")
	else:
		out = open("%s.%s"%(gff,out_str),"w")
	for line in inp:
		if line.startswith(com):
			out.write(line)
		else:
			lineLst = line.strip().split("\t")
			loc_id = make_loc_id(lineLst)
			if loc_id in id_set:
				out.write(line)	
	inp.close()
	out.close()

def main():
	import sys
	if len(sys.argv) == 1 or "-h" in sys.argv:
		print_help()
		sys.exit()
	
	infile,func,chr_str,gff2_file,com_char,inp_src,inp_typ,inp_int,key,genes,out_str = set_defaults_and_parse_args(sys.argv)
	
	if func == None or infile == None:
		print_help()
		print("Function (-f) and input file (-i) required!")
		sys.exit()
	elif func == "sort":
		function_sortGFF(infile,com_char,chr_str)
	elif func == "merge":
		function_mergeGFF(infile)
	elif func == "sort_merge":
		print("Sorting GFF file")
		function_sortGFF(infile,com_char,chr_str)
		print("Merging GFF file")
		function_mergeGFF(infile+".sort")
		print("Done!")
	elif func == "filter":
		print("Filtering GFF file")
		function_filter(infile,chr_str,genes)
		print("Done!")
	elif func == "length":
		function_length(infile)
	elif func == "prefix_id":
		if chr_str == None:
			print_help()
			print("Character string (-str) required to run add_prefix_to_id function")
		else:
			function_prefixID(infile,chr_str)
	elif func == "mask":
		if gff2_file == None:
			print_help()
			print("Secondary gff file (-gff2) required to run mask function")
			sys.exit()
		else:
			function_maskGFF(infile,gff2_file,inp_typ,chr_str)
	elif func == "overlap":
		if gff2_file == None:
			print_help()
			print("Secondary gff file (-gff2) required to run overlap function")
			sys.exit()
		else:
			function_overlap(infile,gff2_file,chr_str,inp_typ)
	elif func == "overlap+":
		if gff2_file == None:
			print_help()
			print("Secondary gff file (-gff2) required to run overlap+ function")
			sys.exit()
		else:
			function_overlapPlus(infile,gff2_file)
	elif func == "compare_lens":
		if gff2_file == None:
			print_help()
			print("Secondary gff file (-gff2) required to run compare_lens function")
			sys.exit()
		else:
			function_compareLens(infile,gff2_file,chr_str)
	elif func == "coords2gff":
		function_coords2gff(infile,inp_src,inp_typ)
	elif func == "mapped2gff":
		function_mapped2gff(infile,key, inp_src,inp_typ)
	elif func == "mapped2gff2":
		function_mapped2gff2(infile,key, inp_src,inp_typ)
	elif func == "mapped2gff_internalkey":
		function_mapped2gff_internalkey(infile, inp_src,inp_typ)
	elif func == "lengths":
		function_lengths(infile,chr_str,com_char)
	elif func == "len_prcntl":
		function_lenPercentiles(infile,inp_int,com_char)
	elif func == "min_len":
		if inp_int == None:
			print_help()
			print("Minimum length (-int) required to run min_len function")
			sys.exit()
		else:
			function_minLen(infile,inp_int,inp_typ,com_char)
	elif func == "merge_depth":
		function_mergeDepth(infile,com_char)
	elif func == "loc_id":
		function_locID(infile,com_char)
	elif func == "loc_id_keep":
		function_locIDkeep(infile,com_char,chr_str,out_str)
	elif func == "loc_id2gff":
		function_locId2gff(infile,inp_src,inp_typ)
	else:
		print_help()
		print("Unrecognized function:",func)

if __name__ == "__main__":
	main()
