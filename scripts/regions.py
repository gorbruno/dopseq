# supress warnings for rpy2-pandas interface deprecation
# and 0-division during log calculation
# import warnings
# warnings.simplefilter(action='ignore')

import pybedtools
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
try:
    from scipy.stats import binom_test
except ImportError:
    from scipy.stats import binomtest as binom_test
from collections import OrderedDict


def read_fai(genome_fai):
    """
    chromosome sizes from fasta index
    chromosome : (0, chromosome_length)
    """
    chrom_lens = OrderedDict()
    with open(genome_fai) as f:
        for l in f:
            ll = l.split('\t')
            try:
                chrom_lens[ll[0]] = (0, int(ll[1]))
            except:
                pass
    return chrom_lens

def bam_to_pos_and_dist(in_bam, out_pos, genome_fai, reads_threshold):
    """
    Convert BAM into positions and
    distances between positions
    """

    in_bam = pybedtools.BedTool(in_bam)
    # non-empty BAM
    if in_bam.count() > 0:
        # positions - merged overlapping reads
        pos_bed = in_bam.bam_to_bed().sort(faidx=genome_fai).merge(c=1, o='count').saveas(out_pos)
        pos = pos_bed.to_dataframe()
        pos['chrom'] = pos['chrom'].astype(str)
        # distance - complement of positions, i.e. end-to-start distances between positions
        dist_raw = pos_bed.complement(g=genome_fai).to_dataframe()
        dist_raw['chrom'] = dist_raw['chrom'].astype(str)
        # keep regions with mapped reads
        valid_chroms = pos["chrom"].unique()
        dist_valid = dist_raw[dist_raw["chrom"].isin(valid_chroms)]
        # filter regions based on count of mapped reads
        chrom_count = pos["chrom"].value_counts()
        chroms_to_keep = chrom_count[chrom_count > reads_threshold].index
        if len(chroms_to_keep) == 0:
            print(f'WARNING: No regions with reads above {reads_threshold} threshold, keeping distance and position tables unfiltered')
            dist_filtered = dist_valid
            pos_filtered = pos
        else:
            dist_filtered = dist_valid[dist_valid["chrom"].isin(chroms_to_keep)]
            pos_filtered = pos[pos["chrom"].isin(chroms_to_keep)]
        # logscale end-to-start distances
        dist_filtered['log.dist'] = np.log10(dist_filtered['end'] - dist_filtered['start'])

        return (pos_filtered, dist_filtered)
    # empty BAM
    else:
        # touch positions
        open(out_pos, 'w').close()

        return (None, None)

def segment_genome(dist, sample, do_plot_reg, out_plot, chrom_list, chrom_lens, ncols=8, hipc=2, wipc=2):
    """
    Segmentation of genome based on logscale distances between read positions
    with DNAcopy circular binary segmentation algorithm.
    Outputs predicted regions with different mean distances between positions (pd_mean)
    """
    # R code for segmentation
    print('- importing DNAcopy')
    dnacopy = importr('DNAcopy')
    print('- reading distance data')
    cna = dnacopy.CNA(robjects.FloatVector(dist['log.dist']),
                            robjects.StrVector(dist['chrom']),
                            robjects.IntVector(dist['end']), # end of dist = start of pos
                            data_type="logratio", sample=sample)
    print('- smoothing distance data')
    cna = dnacopy.smooth_CNA(cna)
    print('- segmenting')
    segm = dnacopy.segment(cna, verbose=0)

    # plotting
    if do_plot_reg:
        print('- re-segmenting for plotting')
        # get chromosome list for plotting
        if len(chrom_list) == 0:
            # use faidx order for all chromosomes
            # helps if nrows limit reached
            plot_chrom = list(chrom_lens.keys())
        else:
            plot_chrom = chrom_list
        len_chrom = len(plot_chrom)
        len_dist = len(dist["chrom"].unique())
        nchrom = min(len_chrom, len_dist)
        # set up figure dimensions
        ncols = ncols
        nrows = int(nchrom / ncols) + (nchrom % ncols > 0)
        # need to limit nrows in plot
        if nrows > 200:
            nrows = 200
            plot_nchrom = ncols * nrows
            print('WARNING: too many sequences in reference, '
                  ' only first {} will be displayed on plot'.format(plot_nchrom))
            plot_chrom = plot_chrom[:plot_nchrom]
        # re-segmenting
        # TODO check chromosome set consistency
        plot_dist = dist[dist['chrom'].isin(plot_chrom)]
        plot_cna = dnacopy.CNA(robjects.FloatVector(plot_dist['log.dist']),
                        robjects.StrVector(plot_dist['chrom']),
                        robjects.IntVector(plot_dist['end']), # end of dist = start of pos
                        data_type="logratio", sample=sample)
        plot_cna = dnacopy.smooth_CNA(plot_cna)
        plot_segm = dnacopy.segment(plot_cna, verbose=0)
        print('- plotting')
        dims = robjects.IntVector([nrows, ncols])
        height = nrows * hipc
        width = ncols * wipc
        # set ylimits to meaningful for PD - less than 0 not expected
        ymax = dist['log.dist'].max()
        ylim = robjects.FloatVector([0, ymax])
        grdevices = importr('grDevices')
        grdevices.pdf(file=out_plot, width=width, height=height)
        dnacopy.plot_DNAcopy(plot_segm,
                            plot_type='s',
                            ylim=ylim,
                            sbyc_layout=dims,
                            xmaploc=True)
        grdevices.dev_off()

    # convert to pandas dataframe
    print('- converting results to pandas dataframe')
    try:
        pandas2ri.activate()
        regions = robjects.pandas2ri.ri2py(segm[1])
    except AttributeError:
        with (robjects.default_converter + pandas2ri.converter).context():
            regions = robjects.conversion.get_conversion().rpy2py(segm[1])
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
        on intervals between positions (+1 position per chromosome).
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

        # 1 dist for empty
        if reg['reg_pos'] > 1:
            reg_pos = pos[
                (pos['chrom'] == reg['chrom']) &
                (pos['start'] >= reg['reg_start']) &
                (pos['start'] <= reg['reg_end'])
                ]
            # 
            # assert len(reg_pos) == int(reg['reg_pos']) - 1, \
            #     f'{len(reg_pos) } positions, {reg.reg_pos} intervals for {reg.chrom}'
            # try:
            reg['first_pos_start'] = reg_pos['start'].iloc[0]
            reg['last_pos_end'] = reg_pos['end'].iloc[-1]
            reg['reg_pos'] = reg_pos.shape[0] # recalc from intervals to pos
            reg['reg_reads'] = reg_pos['name'].sum() # coverage as 'name' column in pos df
            reg['pos_cov_mean'] = reg_pos['name'].mean() 
            reg['pos_len_mean'] = (reg_pos['end'] - reg_pos['start']).mean()
            reg['pos_len_sum'] = (reg_pos['end'] - reg_pos['start']).sum()
            # except: # no positions in chromosome
            #     print('empty')
            #     fill_empty(reg)
        else:
            fill_empty(reg)
        return reg

    regions = regions.sort_values(['chrom', 'reg_start'])
    regions = regions.reset_index(drop=True)
    shift_regs = []
    pr = None
    for (i, r) in regions.iterrows():
        
        cr = r
        # first region in dataset
        if pr is None:
            cr['reg_start'] = 0
        else:
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
            
        # compare to next region
        if i < len(regions) - 1:
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
        else:
            cr['reg_end'] = chrom_lens[cr.chrom][1]

        cr = annotate_region(cr, pos)

        shift_regs.append(cr)
        pr = cr

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
    regions['p_value'] = regions.apply(
        lambda row: binom_test(
            row['reg_pos'],
            total_pos,
            row['reg_len'] / genome_len,
            alternative='greater'
            ), axis=1)
    return regions

if __name__ == "__main__":
    in_bam = snakemake.input[0]
    genome_fai = str(snakemake.input.genome_fai) # for some reason read as snakemake.io.Namedlist
    chrom_lens = read_fai(genome_fai)
    out_pos = snakemake.output.pos
    out_reg = snakemake.output.reg
    do_plot_reg = snakemake.params.do_plot_reg
    out_plot = snakemake.params.plot
    sample = snakemake.params.sample
    print('Converting BAM to positions and distances between positions')
    (pos, dist) = bam_to_pos_and_dist(in_bam, out_pos, genome_fai, snakemake.params.reads_threshold)
    # empty BAM - touch remaining outputs
    if pos is None:
        print('No reads - creating empty output files')
        open(out_reg, 'w').close()
        if do_plot_reg:
            open(out_plot, 'w').close()
    else:
        print('Reading chromosome list')
        chrom_list_file = str(snakemake.params.chrom_list)
        print('{}'.format(chrom_list_file))
        chrom_list = []
        # TODO better type casting
        if (chrom_list_file != 'nan') and (len(chrom_list_file) > 0):
            with open(chrom_list_file) as f:
                for l in f:
                    l = l.strip()
                    if len(l) > 0:
                        chrom_list.append(l.split()[0])
        print('Segmenting the genome with DNAcopy')
        regions = segment_genome(
            dist, 
            sample,
            do_plot_reg,
            out_plot,
            chrom_list,
            chrom_lens,
            snakemake.params.plot_ncols,
            snakemake.params.plot_chrom_height,
            snakemake.params.plot_chrom_width
            )
        print('Correcting and annotating the segments')
        regions = shift_regions(regions, pos, chrom_lens)
        print('Collecting statistics')
        regions = regions_stats(regions, chrom_lens)
        print('Writing output')
        regions.to_csv(out_reg, sep="\t", index=False)
        print('Done!')
