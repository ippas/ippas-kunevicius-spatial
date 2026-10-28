#!/usr/bin/env python3
"""Downstream abundance comparison for cell2location mouse-brain 59 types.

DOES NOT filter cell types using disease-group p-values or arbitrary abundance cutoffs.
- Test all available types, with sample/section as observation.
- Full-section BH across all types; per-cluster BH across ALL cluster x cell type tests.
- Choose the 10 most abundant cell types using group-blinded, across-sample mean
  ONLY for secondary figures/tables; FDRs remain those from the full family.

WARNING: sections from the same animal must not be treated as independent
biological replicates for confirmatory claims. Section-level tests are exploratory.
"""
from __future__ import annotations
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(os.environ.get('PROJECT_ROOT', '/home/mateusz/mnt/ippas-kunevicius-spatial'))
HERE = Path(__file__).resolve().parent
OUT = ROOT / 'results/maternalFMT_n16samples/cellTypeDeconvolution_cell2location/mouseBrain59'
META = HERE / 'metadata_autismFMT_for_RCTD.tsv'
SECTION = OUT / 'wholeSection'
CLUSTER = OUT / 'perCluster'
QC = OUT / 'qc'
TOP_N = 10
MIN_GROUP_N = 3  # for a usable test with sample-level variance; not a cell-type filter


def bh(df, pcol, name, scope='all'):
    df[name] = np.nan
    if scope == 'all':
        ok = np.isfinite(df[pcol].to_numpy(dtype=float))
        if ok.any():
            df.loc[ok, name] = multipletests(df.loc[ok, pcol].to_numpy(), method='fdr_bh')[1]
    elif scope == 'cluster':
        for _, idx in df.groupby('cluster', sort=False, observed=True).groups.items():
            p = df.loc[idx, pcol].to_numpy(dtype=float)
            ok = np.isfinite(p)
            if ok.any():
                df.loc[np.asarray(idx)[ok], name] = multipletests(p[ok], method='fdr_bh')[1]
    else:
        raise ValueError(scope)
    return df


def test_groups(summary, by):
    rows = []
    grouped = summary.groupby(by, observed=True, sort=True)
    for key, part in grouped:
        if not isinstance(key, tuple): key=(key,)
        row = dict(zip(by, key))
        x = part.loc[part['fmt_donor_group'].eq('ASD'), 'mean_proportion'].to_numpy(dtype=float)
        y = part.loc[part['fmt_donor_group'].eq('Neurotypical'), 'mean_proportion'].to_numpy(dtype=float)
        x = x[np.isfinite(x)]; y=y[np.isfinite(y)]
        row.update(n_ASD=len(x),n_NT=len(y),
                   mean_ASD=np.mean(x) if len(x) else np.nan,
                   mean_NT=np.mean(y) if len(y) else np.nan,
                   median_ASD=np.median(x) if len(x) else np.nan,
                   median_NT=np.median(y) if len(y) else np.nan)
        if len(x) and len(y):
            row['delta_percentage_points'] = 100 * (row['mean_ASD']-row['mean_NT'])
            row['relative_change_percent'] = 100*(row['mean_ASD']/row['mean_NT']-1) if row['mean_NT']>0 else np.nan
        else:
            row['delta_percentage_points']=np.nan
            row['relative_change_percent']=np.nan
        if min(len(x),len(y))>=MIN_GROUP_N and (np.var(x)+np.var(y))>0:
            row['p_welch']=ttest_ind(x,y,equal_var=False).pvalue
            row['p_mann_whitney']=mannwhitneyu(x,y,alternative='two-sided',method='auto').pvalue
        else:
            row['p_welch']=np.nan
            row['p_mann_whitney']=np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def validate_and_read():
    if not META.exists(): raise FileNotFoundError(META)
    metadata = pd.read_csv(META, sep='\t')
    if metadata['sample_ID'].duplicated().any():
        metadata=metadata.drop_duplicates('sample_ID')
    assert {'sample_ID','fmt_donor_group'}.issubset(metadata.columns)
    if not metadata['fmt_donor_group'].isin(['ASD','Neurotypical']).all():
        raise ValueError('Unknown groups in metadata')
    if metadata.set_index('sample_ID').loc['2_1Md','fmt_donor_group']!='ASD':
        raise ValueError('2_1Md must be ASD according to supplied corrected metadata')
    if metadata.set_index('sample_ID').loc['2_1M','fmt_donor_group']!='Neurotypical':
        raise ValueError('2_1M must be Neurotypical according to supplied corrected metadata')
    whole = pd.read_csv(SECTION/'allSamples_n16_cellTypeProportions_sectionSummary.tsv', sep='\t')
    cluster = pd.read_csv(CLUSTER/'allSamples_n16_cellTypeProportions_perCluster.tsv', sep='\t')
    needed={'sample_ID','cell_type','mean_proportion'}
    for d,label in [(whole,'whole'),(cluster,'cluster')]:
        if not needed.issubset(d.columns): raise ValueError(f'Missing required columns in {label}: {needed-set(d)}')
        if not np.isfinite(d['mean_proportion']).all():raise ValueError('Nonfinite proportions')
        if not d['mean_proportion'].between(0,1).all():raise ValueError('Proportions outside [0,1]')
    if 'cluster' not in cluster: raise ValueError('Missing cluster column')
    samples=set(whole.sample_ID.unique())
    if len(samples)!=16:raise ValueError(f'Expected 16 samples, found {len(samples)}')
    if not samples.issubset(set(metadata.sample_ID)): raise ValueError('Missing metadata for some sample IDs')
    if whole.duplicated(['sample_ID','cell_type']).any():raise ValueError('Duplicate sample/type entries in whole section summary')
    if cluster.duplicated(['sample_ID','cluster','cell_type']).any():raise ValueError('Duplicate sample/cluster/type entries')
    types=sorted(whole.cell_type.unique())
    if set(cluster.cell_type)!=set(types): raise ValueError('Different cell types in whole and clusters')
    if len(types)!=59: warnings.warn(f'Expected 59 types; got {len(types)}. Continuing with all available types')
    whole=whole.merge(metadata[['sample_ID','fmt_donor_group']],on='sample_ID',validate='many_to_one')
    cluster=cluster.merge(metadata[['sample_ID','fmt_donor_group']],on='sample_ID',validate='many_to_one')
    counts=whole[['sample_ID','fmt_donor_group']].drop_duplicates()['fmt_donor_group'].value_counts()
    if counts.get('ASD',0)!=10 or counts.get('Neurotypical',0)!=6:
        raise ValueError('Expected ASD=10, Neurotypical=6 based on 16 provided sections')
    return whole, cluster, types


def create_figures(whole, cluster, top):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not installed: skipping plots'); return
    fig,axes=plt.subplots(5,2,figsize=(11,15),sharex=True)
    groups=['Neurotypical','ASD']
    rng=np.random.default_rng(7)
    for ax,ct in zip(axes.ravel(),top):
        d=whole[whole.cell_type==ct]
        for pos, grp in enumerate(groups):
            vals=100*d.loc[d.fmt_donor_group==grp,'mean_proportion'].values
            ax.scatter(pos+rng.uniform(-.12,.12,len(vals)),vals,alpha=.7,s=22)
            if len(vals):ax.hlines(vals.mean(),pos-.18,pos+.18,lw=2.5)
        ax.set_title(ct,fontsize=9)
        ax.set_xticks([0,1],groups)
        ax.set_ylabel('Estimated proportion (%)',fontsize=8)
    fig.suptitle('Top 10 cell types by overall abundance (group-blind selection)')
    fig.tight_layout(rect=(0,0,1,.98))
    fig.savefig(SECTION/'top10_proportions_by_sample.pdf')
    plt.close(fig)
    # A heatmap of group-mean differences, NOT a significance heatmap.
    stats = pd.read_csv(CLUSTER/'ASD_vs_Neurotypical_perCluster_allTypes_stats.tsv',sep='\t')
    sub=stats[stats.cell_type.isin(top)]
    mat=sub.pivot(index='cluster',columns='cell_type',values='delta_percentage_points')
    mat=mat.reindex(columns=top)
    fig,ax=plt.subplots(figsize=(12,max(5,len(mat)*.35)))
    vmax=np.nanmax(np.abs(mat.to_numpy(dtype=float)))
    im=ax.imshow(mat.to_numpy(dtype=float),aspect='auto',cmap='coolwarm',vmin=-vmax,vmax=vmax)
    ax.set_xticks(range(len(mat.columns)),mat.columns,rotation=60,ha='right',fontsize=8)
    ax.set_yticks(range(len(mat.index)),[str(v) for v in mat.index])
    ax.set_ylabel('Leiden cluster')
    ax.set_title('ASD − Neurotypical (percentage points): top 10 abundant types')
    fig.colorbar(im,ax=ax,label='ASD − NT (percentage points)')
    fig.tight_layout();fig.savefig(CLUSTER/'top10_perCluster_delta_percentage_points.pdf');plt.close(fig)


def main():
    for directory in [SECTION,CLUSTER,QC]:directory.mkdir(parents=True,exist_ok=True)
    whole, cluster, types=validate_and_read()
    # Choose top 10 only for presentation, pooling equally across sections and ignoring group labels.
    abundance=whole.groupby('cell_type',sort=True)['mean_proportion'].mean().sort_values(ascending=False)
    ranking=abundance.rename('mean_proportion_across_all_16_sections').reset_index()
    ranking.insert(0,'rank',range(1,len(ranking)+1))
    ranking['selected_top10_for_figures']=ranking['rank']<=TOP_N
    ranking.to_csv(QC/'cellType_abundanceRanking_groupBlind.tsv',sep='\t',index=False)
    top=ranking.head(TOP_N)['cell_type'].tolist()

    ws=test_groups(whole,['cell_type'])
    for p in ['p_welch','p_mann_whitney']:
        bh(ws,p,'FDR_BH_allTypes_'+p)
    ws['selected_top10_for_figures']=ws.cell_type.isin(top)
    ws=ws.sort_values(['p_welch','cell_type'],na_position='last')
    ws.to_csv(SECTION/'ASD_vs_Neurotypical_allTypes_stats.tsv',sep='\t',index=False)
    ws[ws.selected_top10_for_figures].to_csv(SECTION/'ASD_vs_Neurotypical_top10_FOR_FIGURES.tsv',sep='\t',index=False)

    cs=test_groups(cluster,['cluster','cell_type'])
    # Global BH correction over ALL testable cluster x type pairs: primary per-cluster family.
    for p in ['p_welch','p_mann_whitney']:
        bh(cs,p,'FDR_BH_global_'+p)
        bh(cs,p,'FDR_BH_withinCluster_'+p,scope='cluster')  # secondary, less stringent
    cs['selected_top10_for_figures']=cs.cell_type.isin(top)
    cs=cs.sort_values(['p_welch','cluster','cell_type'],na_position='last').reset_index(drop=True)
    cs.to_csv(CLUSTER/'ASD_vs_Neurotypical_perCluster_allTypes_stats.tsv',sep='\t',index=False)
    cs[cs.selected_top10_for_figures].to_csv(CLUSTER/'ASD_vs_Neurotypical_perCluster_top10_FOR_FIGURES.tsv',sep='\t',index=False)

    create_figures(whole,cluster,top)
    lines=[
        'Analysis: all available cell types included in testing; NO abundance/prevalence cutoff.',
        'Top 10 selected by mean whole-section proportion across 16 sections, ignoring disease labels; DISPLAY ONLY.',
        'Welch t-test and two-sided Mann-Whitney tests use each section as an observation.',
        f'Minimum samples per group for testing: {MIN_GROUP_N}.',
        f'Whole section: {len(ws)} types, {ws.p_welch.notna().sum()} tested, {(ws.p_welch<.05).sum()} nominal Welch p<0.05, {(ws.FDR_BH_allTypes_p_welch<.05).sum()} FDR<0.05.',
        f'Per cluster: {len(cs)} candidate combinations, {cs.p_welch.notna().sum()} tested, {(cs.p_welch<.05).sum()} nominal Welch p<0.05, {(cs.FDR_BH_global_p_welch<.05).sum()} globally FDR<0.05.',
        'Important: global per-cluster FDR is BH across ALL testable cluster x cell-type comparisons, not just top 10.',
        'IMPORTANT LIMITATION: sections from the same animal are NOT independent biological replicates; these section-level p-values are exploratory.',
        'To make confirmatory animal-level claims, aggregate repeats by animal ID or fit a hierarchical model.',
        'Proportions are normalized posterior mean cell abundances; they are model estimates, not histological cell counts.',
    ]
    (QC/'statistical_analysis_summary.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print('\n'.join(lines))

if __name__=='__main__': main()
