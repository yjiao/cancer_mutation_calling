# group samples by patient, and get the union of mutations from this patient
from glob import glob
from collections import defaultdict
from os import system
from os import path
from os import makedirs

samplespath="/seq/hacohenlab1/yjiao/prepost/bcbio/config/samples/allsamples/allsamples_20160902.csv"
sampleGroups = defaultdict(list)

prefix = '/seq/hacohenlab1/yjiao/prepost/bcbio/vardictRuns/'
suffix = '/final/*/passed.vcf'

destination_dir = '/seq/hacohenlab1/yjiao/prepost/voting/vardict'

if not path.exists(destination_dir):
    makedirs(destination_dir)

with open(samplespath) as f:
    lines = f.readlines()

for line in lines:
    line.strip()
    line = line.split(',')
    if line[3] == 'tumor':
	sampleGroups[line[1]].append(prefix + line[0] + suffix)


for key, val in sampleGroups.iteritems():
    cmd = 'vcfcombine ' + ' '.join(val) + ' > ' + destination_dir + '/' + key + '_vardict_combined.vcf&'
#    qsub = 'qsub -b y -q short -l h_vmem=1g -N vcfcombine_' + key + ' -cwd ' + cmd
#    print qsub
#    system(qsub)
    system(cmd)
    print cmd