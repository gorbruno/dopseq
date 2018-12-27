import pandas as pd
import glob
import sys
import os

'''
Sanity checks of dopseq output
Usage: python test_validate.py {dopseq_dir}
'''


assert pd.__version__.startswith('0.23')

WD = sys.argv[1]
STAT = os.path.join(WD, 'results/stats.xlsx')

sd = pd.read_excel(STAT)

def run_test(test, msg):
    '''If test passes, print message formatted with list of indices'''
    if test.any():
        print(msg.format(', '.join(sd[test].index)))
        
trim_test = (sd.total_reads < sd.trimmed_reads)
run_test(trim_test,
         msg='ERROR: less total reads than trimmed reads for {}')
trim_bp_test = (sd.total_bp < sd.trimmed_bp)
run_test(trim_bp_test,
         msg='ERROR: total bp < trimmed bp for {}')
mapped_test = (sd.trimmed_reads < sd.mapped_reads)
run_test(mapped_test,
         'ERROR: less trimmed reads than mapped reads for {}')
duplicated_test = (sd.mapped_reads < sd.duplicated_reads)
run_test(duplicated_test,
         msg='ERROR: less mapped reads than duplicated reads for {}')
filter_test = ((sd.mapped_reads - sd.duplicated_reads) < sd.mapped_reads_after_filter)
run_test(filter_test,
         msg='ERROR: too many reads remaining after duplicate removal and filtering for {}')
filter_bp_test = (sd.trimmed_bp < sd.mapped_bp_after_filter)
run_test(filter_bp_test,
         msg='ERROR: less trimmed bp than filtered bp for {}')
divergence_test = (sd.error_rate > 10)
run_test(divergence_test,
         msg='ERROR: sequence divergence over 10% for {}')
qual_test = ((sd.average_quality > 40) | (sd.average_quality < 20))
run_test(qual_test,
         msg='ERROR: average mapq not within (20,40) for {}')

print('Checking logs')
for fn in glob.iglob(os.path.join(WD, 'results/logs/**/*.log'), recursive=True):
    with open(fn) as f:
        for l in f:
            lowl = l.lower()
            for s in ('error', 'warning'):
                if s in lowl:
                    print(fn + ':')
                    sys.stdout.write(l)