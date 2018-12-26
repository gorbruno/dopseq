# supress warnings for rpy2-pandas interface deprecation
# and 0-division during log calculation
import warnings
warnings.simplefilter(action='ignore')

import pybedtools
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from scipy.stats import binom_test


def read_fai(gneome_fai): 
    """
    chromosome sizes from fasta index
    chromosome : (0, chromosome_length)
    """
    chrom_lens = dict()
    with open(genome_fai) as f:
        for l in f:
            ll = l.split('\t')
            try:
                chrom_lens[ll[0]] = (0, int(ll[1]))
            except:
                pass
    return chrom_lens

def bam_to_pos_and_dist(in_bam, out_pos):
    # positions and distances between positions
    in_bam = pybedtools.BedTool(in_bam)
    # positions - merged overlapping reads
    pos_bed = in_bam.bam_to_bed().sort().merge(c=1, o='count').saveas(out_pos)
    pos = pos_bed.to_dataframe()
    pos['chrom'] = pos['chrom'].astype(str)
    # distance - complement of positions, i.e. end-to-start distances between positions
    dist = pos_bed.complement(g=chrom_lens).to_dataframe()
    dist['chrom'] = dist['chrom'].astype(str)
    # logscale end-to-start distances
    dist['log.dist'] = np.log10(dist['end'] - dist['start'])

    return (pos, dist)

def segment_genome(dist, sample):
    """
    Segmentation of genome based on logscale distances between read positions
    with DNAcopy circular binary segmentation algorithm.
    Outputs predicted regions with different mean distances between positions (pd_mean)
    """
    # R code for segmentation
    dnacopy = importr('DNAcopy')
    cna = robjects.r['CNA'](robjects.FloatVector(dist['log.dist']),
                            robjects.StrVector(dist['chrom']), 
                            robjects.IntVector(dist['end']), # end of dist = start of pos
                            data_type="logratio", sample=sample)
    cna = robjects.r['smooth.CNA'](cna)
    segm = robjects.r['segment'](cna, verbose=0)

    # convert to pandas dataframe
    pandas2ri.activate()
    regions = robjects.pandas2ri.ri2py(segm[1])
    
    # normalize names
    regions = regions.rename(columns={
        'ID':'sample',
        'loc.start':'reg_start', 
        'loc.end':'reg_end',
        'num.mark':'reg_pos',
        'seg.mean':'lg_pd_mean'})
    # recover mean distances between positions
    regions['pd_mean'] = np.power(10, regions['lg_pd_mean'])
    
    return regions

def shift_regions(regions, pos, chrom_lens):
    """
    Segmentation results correction and annotation
    
    Priority:
    - Start of first region in chromosome - 0
    - End of last region in chromosome - chromosome length
    - Highest pd regions start and end at positions
    - Lower pd regions include distance to higher pd region position.
    
    Implemented as iterrows, speedup possible
    """

    def annotate_region(reg, pos):
        """
        For a single region,
        recalculate number of positions as initial clustering was done 
        on distance complements (+1 position per chromosome).
        Get coverage and position size statistics for positions within a region.
        """
        def fill_empty(reg):
            """Fill data for region without read positions"""
            reg['first_pos_start'] = np.nan
            reg['last_pos_end'] = np.nan
            reg['reg_pos'] = 0
            reg['reg_reads'] = 0
            reg['pos_cov_mean'] = 0
            reg['pos_len_mean'] = 0
            reg['pos_len_sum'] = 0
            # deleting seg_mean and pd_mean, as values produced by DNAcopy are not meaningful
            reg['pd_mean'] = np.nan
            reg['lg_pd_mean'] = np.nan

            return reg

        if reg['reg_pos'] > 1:
            reg_pos = pos[(pos['chrom'] == reg['chrom']) & 
                            (pos['start'] >= reg['reg_start']) & 
                            (pos['start'] <= reg['reg_end'])]
            try:
                reg['first_pos_start'] = reg_pos['start'].iloc[0]
                reg['last_pos_end'] = reg_pos['end'].iloc[-1]
                reg['reg_pos'] = reg_pos.shape[0] # nrows
                reg['reg_reads'] = reg_pos['name'].sum()
                reg['pos_cov_mean'] = reg_pos['name'].mean() # coverage as 'name' column in pos df
                reg['pos_len_mean'] = (reg_pos['end']-reg_pos['start']).mean()
                reg['pos_len_sum'] = (reg_pos['end']-reg_pos['start']).sum()
            except: # no positions in chromosome 
                print('empty')
                fill_empty(reg)
        else:
            fill_empty(reg)
        return reg

    regions = regions.sort_values(['chrom', 'reg_start'])
    shift_regs = []
    # slow iterrows
    for (i, r) in regions.iterrows():
        
        cr = r
        # compare to previous region
        try: 
            pr = regions.iloc[i - 1]
            # first region in chromosome
            if cr.chrom != pr.chrom:  
                cr['reg_start'] = 0
            elif pr.chrom == cr.chrom:
                # low-pd to high-pd
                # shift to start of previous position
                # coincides with previous region end
                if pr.pd_mean > cr.pd_mean:
                    cr['reg_start'] = pr.reg_end
                # high-pd to low-pd
                # shift to end of previous position + 1
                elif pr.pd_mean < cr.pd_mean:
                    prev_pos = pos[(pos['chrom'] == pr.chrom) & (pos['start'] == pr.reg_end)]
                    cr['reg_start'] = prev_pos['end'].iloc[0] + 1
                else:
                    raise ValueError('Consequent regions have same pd_mean:'
                                     '\n{}\n{}'.format(pr, cr))
        # first region in dataset
        except:
            cr['reg_start'] = 0
        # compare to next region
        try:
            nr = regions.iloc[i + 1]
            # last region in chromosome
            if cr.chrom != nr.chrom:
                cr['reg_end'] = chrom_lens[cr.chrom][1]
            elif cr.chrom == nr.chrom: 
                # low-pd to high pd
                # shift 1 bp left (position is correct)
                if cr.pd_mean > nr.pd_mean:
                    cr['reg_end'] = cr.reg_end - 1
                # high-pd to low-pd
                # shift to end of current position (position is correct)
                elif cr.pd_mean < nr.pd_mean: 
                    curr_pos = pos[(pos['chrom'] == cr.chrom) & (pos['start'] == cr.reg_end)]
                    cr['reg_end'] = curr_pos['end'].iloc[0]
                else:
                    raise ValueError('Consequent regions have same pd_mean:'
                                     '\n{}\n{}'.format(cr, nr))
            else:
                pass
        # last region in dataset
        except:
            cr['reg_end'] = chrom_lens[cr.chrom][1]
            
        cr = annotate_region(cr, pos)
        
        shift_regs.append(cr)

    return pd.DataFrame(shift_regs)

def regions_stats(regions, chrom_lens):
    """
    Calculate per-region statistics, 
    log-ratio of enrichment with positions (>0 - enrichment, <0 - depletion),
    p-value for enrichment with positions
    """
    # general stats
    total_pos = regions['reg_pos'].sum()
    genome_len = sum([c[1] for c in chrom_lens.values()])
    regions['reg_len'] = regions['reg_end'] - regions['reg_start']
    regions['chrom_len'] = pd.DataFrame(regions['chrom'].map(chrom_lens).values.tolist())[1]
    regions['log_ratio'] = np.log10( (regions['reg_pos'] / total_pos) / (regions['reg_len'] / genome_len) )
    # enrichment p-value
    regions['p_value'] = regions.apply(lambda row: binom_test(row['reg_pos'],  total_pos,
                                               row['reg_len'] / genome_len, 
                                               alternative='greater'), axis=1)
    return regions

if __name__ == "__main__":
    in_bam = snakemake.input[0]
    genome_fai = str(snakemake.input.genome_fai) # for some reason read as snakemake.io.Namedlist
    out_pos = snakemake.output.pos
    out_reg = snakemake.output.reg
    sample = snakemake.params.sample
    chrom_lens = read_fai(genome_fai)
    (pos, dist) = bam_to_pos_and_dist(in_bam, out_pos)
    regions = segment_genome(dist, sample)
    regions = shift_regions(regions, pos, chrom_lens)
    regions = regions_stats(regions, chrom_lens)
    regions.to_csv(out_reg, sep="\t", index=False)