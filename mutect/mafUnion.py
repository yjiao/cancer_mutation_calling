from glob import glob
from collections import defaultdict
from os import system

# format of MAF files:
#  build chr     start       end ref_allele tum_allele1 tum_allele2       tumor_barcode       normal_barcode  tumor_f init_t_lod t_lod_fstar t_alt_count t_ref_count judgement
#    37   1  12907458  12907458          T           T           C 107T-Tumor-SM-A2WOJ 107T-Normal-SM-A3PYV 0.020548 -68.106301    7.121002           6         286      KEEP
#    37   1  17086183  17086183          T           T           G 107T-Tumor-SM-A2WOJ 107T-Normal-SM-A3PYV 0.025806 -31.913741    6.654500           4         151      KEEP
#    37   1  36229877  36229877          G           G           T 107T-Tumor-SM-A2WOJ 107T-Normal-SM-A3PYV 0.153846   1.899680    6.658286           6          33      KEEP
#    37   1  39879328  39879328          G           G           A 107T-Tumor-SM-A2WOJ 107T-Normal-SM-A3PYV 0.060606 -13.818255    6.464200           6          93      KEEP
#    37   1  44063417  44063417          A           A           C 107T-Tumor-SM-A2WOJ 107T-Normal-SM-A3PYV 0.183333   2.754435    8.579684          11          49      KEEP
#    37   1 142621112 142621112          T           T           A 107T-Tumor-SM-A2WOJ 107T-Normal-SM-A3PYV 0.018450 -64.159155    6.601377           5         266      KEEP


def merge(samps):
    union = set()
    fhList = [open(filename) for filename in samps]
    flag = True;
    cnt = 0

    for fh in fhList: # discard header
	line = fh.readline()

    for fh in fhList:
	for line in fh.readlines():
	    line = line.strip().split('\t')
	    if line[4] != line[5]:
		raise NameError("ERROR: expected ref_allele == tum_allele1:" + line)
	    if line[2] != line[3]:
		raise NameError("ERROR: expected start == end for SNVs:" + line)

	    # a mutation is defined by chr#, start, end, ref_allele, and tum_allele2
	    mutation = line[1:3]
	    mutation.append(line[4])
	    mutation.append(line[6])

	    # convert fields to int for sorting later. Chr X and weird contigs do not convert, so leave as char
	    try:
		mutation[0:2] = map(int, mutation[0:2])
	    except:
		mutation[1:2] = map(int, mutation[1:2])

	    # convert to tuple in order to add to set
	    mutation = tuple(mutation)
	    union.add(mutation)

    return union

sampleGroups = defaultdict(list)

mafFiles = glob('*.maf')
for file in mafFiles:
    pat = file.split('-')[0]
    sampleGroups[pat].append(file)

for pat, samps in sampleGroups.iteritems():
    print pat, "----------------------------"
    union = merge(samps)
    union = list(union)
    union.sort(key=lambda x: (x[0], x[1], x[2]))
    union = map(list, union)
    file = 'union/' + pat + '_mutect_combined.maf'
    with open(file, 'w') as f:
	for u in union:
	    line = ' '.join(map(str, u)) + '\n'
	    f.write(line)