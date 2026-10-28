#!/usr/bin/env python3
"""Published 59-category mouse brain cell2location reference, 16 slides, group comparisons.
Run from 00_run_all.sh. Requires raw Seurat count exports produced by 01_export_visium.R.
"""
import os, sys, logging, urllib.request, inspect
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.io import mmread
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import anndata as ad
import torch
import cell2location

ROOT = Path(os.environ['PROJECT_ROOT'])
HERE = Path(__file__).resolve().parent
OUT = ROOT/'results'/'maternalFMT_n16samples'/'cellTypeDeconvolution_cell2location'/'mouseBrain59'
DATA = OUT/'inputs'; REF = ROOT/'data'/'reference'/'cell2location_mouseBrain59'
META = HERE/'metadata_autismFMT_for_RCTD.tsv'
SOURCE = 'https://cell2location.cog.sanger.ac.uk/tutorial/mouse_brain_snrna/regression_model/RegressionGeneBackgroundCoverageTorch_65covariates_40532cells_12819genes/sc.h5ad'
N_CELLS_PER_LOCATION = float(os.getenv('N_CELLS_PER_LOCATION','8')) # Prior; adjust from histology nuclei if known
MAX_EPOCHS = int(os.getenv('C2L_MAX_EPOCHS','12000'))
BATCH_SIZE = int(os.getenv('C2L_BATCH_SIZE','1024'))
DEVICE = 'gpu' if torch.cuda.is_available() else 'cpu'
SEED = 7
np.random.seed(SEED); torch.manual_seed(SEED)
for folder in (OUT,DATA,REF,OUT/'perSpot',OUT/'perCluster',OUT/'wholeSection',OUT/'models',OUT/'qc'):folder.mkdir(parents=True,exist_ok=True)
logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s',handlers=[logging.StreamHandler(sys.stdout),logging.FileHandler(OUT/'run.log')])
log=logging.info

def fetch(url,dest):
 if dest.exists() and dest.stat().st_size>1000000:log('Using downloaded reference %s',dest);return
 tmp=dest.with_suffix('.download'); log('Downloading official mouse-brain reference: %s',url)
 try:
  urllib.request.urlretrieve(url,tmp)
  if tmp.stat().st_size<1000000:raise RuntimeError('Reference file suspiciously small')
  tmp.replace(dest)
 except Exception:
  tmp.unlink(missing_ok=True);raise

def get_signatures(ref_path):
 # From original cell2location tutorial, *already fitted* reference model.
 ref=ad.read_h5ad(ref_path,backed='r')
 try:
  if ref.raw is None:raise ValueError('Reference lacks .raw with fitted signature fields')
  var=ref.raw.var
  prefix='means_cov_effect_annotation_1_'
  cols=[c for c in var.columns if c.startswith(prefix)]
  if len(cols)!=59:raise ValueError(f'Expected 59 published signature columns; found {len(cols)}. Available: {list(var.columns[:15])}')
  s=var[cols].copy();s.columns=[c[len(prefix):] for c in cols]
  # Official mouse-brain tutorial: rescale signatures using reference regression model.
  scaling=np.asarray(ref.uns['regression_mod']['post_sample_means']['sample_scaling'],dtype=float).mean()
  if not np.isfinite(scaling) or scaling<=0:raise ValueError('Invalid reference scaling')
  s=s.apply(pd.to_numeric,errors='raise').astype('float32')*scaling
  s.index=s.index.astype(str)
  # Seurat commonly uses gene SYMBOLS; the published H5AD often indexes ENSEMBL IDs.
  # Choose gene namespace using only input names (not disease labels/statistical results).
  example=DATA/'1_1F'/'genes.tsv'
  if not example.exists():raise FileNotFoundError(f'Missing export for identifier check: {example}')
  input_genes=set(pd.read_csv(example,sep='\t')['gene'].astype(str))
  n_id=len(input_genes.intersection(s.index))
  if 'SYMBOL' in var.columns:
   candidate=var.loc[s.index,'SYMBOL'].astype(str)
   n_symbol=len(input_genes.intersection(candidate))
   if n_symbol>n_id:
    s.index=candidate.values
    s=s.loc[s.index.notna() & (s.index!='') & (s.index!='nan')]
    if s.index.duplicated().any():
     # Ambiguous symbols: exclude duplicate mappings rather than silently pick one gene.
     s=s.loc[~s.index.duplicated(keep=False)]
    log('Gene namespace: symbols; overlap before filtering=%d vs original IDs=%d',n_symbol,n_id)
   else:log('Gene namespace: reference IDs; overlap=%d vs symbols=%d',n_id,n_symbol)
  if s.index.duplicated().any():raise ValueError('Duplicated reference genes')
  if (s.to_numpy()<0).any() or not np.isfinite(s.to_numpy()).all():raise ValueError('Invalid reference signatures')
  log('Official regression signature scale factor: %.5g',scaling)
  s.to_csv(OUT/'qc'/'reference_signatures.tsv',sep='\t')
  return s
 finally:ref.file.close()

def read_sample(sample):
 d=DATA/sample
 if not (d/'matrix.mtx.gz').exists():raise FileNotFoundError('Missing R exported sample '+sample)
 x=mmread(d/'matrix.mtx.gz').tocsr().T
 genes=pd.read_csv(d/'genes.tsv',sep='\t')['gene'].astype(str).tolist()
 spots=pd.read_csv(d/'spots.tsv',sep='\t')
 if x.shape != (len(spots),len(genes)):raise ValueError('Matrix / gene / spot size mismatch for '+sample)
 if len(genes)!=len(set(genes)):raise ValueError('Duplicated genes for '+sample+'; fix gene identifiers first')
 a=ad.AnnData(X=x,obs=spots.set_index('spot'),var=pd.DataFrame(index=pd.Index(genes,name='gene')))
 a.obs_names=[f'{sample}__{v}' for v in a.obs_names]
 a.obs['sample_ID']=sample
 return a

def fit_one(sample,sig):
 dest=OUT/'perSpot'/f'{sample}_cell2location_proportions.tsv.gz'
 if dest.exists():log('SKIP COMPLETED %s',sample);return pd.read_csv(dest,sep='\t')
 a=read_sample(sample)
 common=a.var_names.intersection(sig.index)
 log('Sample %s: %d spots, %d overlapping genes, %s',sample,a.n_obs,len(common),DEVICE)
 if len(common)<3000:raise ValueError(f'Only {len(common)} overlapping genes for {sample}; check gene IDs (Ensembl vs symbols).')
 # Official tutorial recommends excluding mitochondrially encoded genes before spatial mapping.
 common=pd.Index([g for g in common if not str(g).lower().startswith('mt-')])
 a=a[:,common].copy()
 s=sig.loc[common].copy()
 cell2location.models.Cell2location.setup_anndata(a,batch_key='sample_ID')
 model=cell2location.models.Cell2location(a,cell_state_df=s,N_cells_per_location=N_CELLS_PER_LOCATION,detection_alpha=20)
 train_kwargs=dict(max_epochs=MAX_EPOCHS,batch_size=min(BATCH_SIZE,a.n_obs),train_size=1)
 if 'use_gpu' in inspect.signature(model.train).parameters:train_kwargs['use_gpu']=(DEVICE=='gpu')
 else:train_kwargs['accelerator']=('gpu' if DEVICE=='gpu' else 'cpu')
 model.train(**train_kwargs)
 posterior_kwargs={'num_samples':1000,'batch_size':min(BATCH_SIZE,a.n_obs)}
 # Use old/new hardware configuration without assuming a specific installed release.
 if 'use_gpu' in inspect.signature(model.export_posterior).parameters:posterior_kwargs['use_gpu']=(DEVICE=='gpu')
 a=model.export_posterior(a,sample_kwargs=posterior_kwargs)
 names=list(a.uns['mod']['factor_names'])
 if 'means_cell_abundance_w_sf' not in a.obsm:raise ValueError('Expected posterior mean cell abundance missing')
 abund=np.asarray(a.obsm['means_cell_abundance_w_sf'],dtype=float)
 if abund.shape!=(a.n_obs,len(names)):raise ValueError('Incorrect abundance matrix dimensions')
 totals=abund.sum(axis=1)
 if (totals<=0).any():raise ValueError('Zero total abundance for some spots')
 probs=abund/totals[:,None]
 df=pd.DataFrame(probs,columns=names);df.insert(0,'spot',a.obs_names);df.insert(1,'sample_ID',sample);df.insert(2,'cluster',a.obs['cluster'].values)
 df.to_csv(dest,sep='\t',index=False,compression='gzip')
 pd.DataFrame(abund,columns=names).assign(spot=a.obs_names).to_csv(OUT/'perSpot'/f'{sample}_estimated_cell_abundance.tsv.gz',sep='\t',index=False,compression='gzip')
 # save model for diagnosis and reproducibility
 model.save(str(OUT/'models'/sample),overwrite=True)
 log('COMPLETED %s',sample)
 return df

def summarise(df,group_cols,types):
 return df.melt(id_vars=['spot','sample_ID','cluster'],value_vars=types,var_name='cell_type',value_name='proportion').groupby(group_cols+['cell_type'],observed=True).agg(mean_proportion=('proportion','mean'),n_spots=('spot','nunique')).reset_index()

def main():
 if not META.exists():raise FileNotFoundError(f'Place provided metadata file next to script: {META}')
 metadata=pd.read_csv(META,sep='\t').drop_duplicates('sample_ID')
 if not {'sample_ID','fmt_donor_group'}.issubset(metadata):raise ValueError('Metadata requires sample_ID, fmt_donor_group')
 samples=['1_1F','1_1Fd','1_1M','2_1M','2_1Md','3_1F','3_1M','5_1M','5_3F','12_1M','13_1F','13_1M','15_1F','18_1F','18_1M','23_1F']
 if set(samples)-set(metadata.sample_ID):raise ValueError('Metadata missing samples')
 assert metadata.set_index('sample_ID').loc['2_1Md','fmt_donor_group']=='ASD'
 assert metadata.set_index('sample_ID').loc['2_1M','fmt_donor_group']=='Neurotypical'
 ref_path=REF/'sc.h5ad';fetch(SOURCE,ref_path);sig=get_signatures(ref_path)
 all_frames=[fit_one(sample,sig) for sample in samples]
 types=list(sig.columns)
 combined=pd.concat(all_frames,ignore_index=True)
 whole=summarise(combined,['sample_ID'],types)
 clusters=summarise(combined,['sample_ID','cluster'],types)
 whole.to_csv(OUT/'wholeSection'/'allSamples_n16_cellTypeProportions_sectionSummary.tsv',sep='\t',index=False)
 clusters.to_csv(OUT/'perCluster'/'allSamples_n16_cellTypeProportions_perCluster.tsv',sep='\t',index=False)
 w=whole.merge(metadata[['sample_ID','fmt_donor_group']],on='sample_ID',validate='many_to_one')
 results=[]
 for ct,g in w.groupby('cell_type'):
  x=g.loc[g.fmt_donor_group=='ASD','mean_proportion'].to_numpy();y=g.loc[g.fmt_donor_group=='Neurotypical','mean_proportion'].to_numpy()
  if len(x)!=10 or len(y)!=6:raise ValueError(f'Group assignment inconsistent for {ct}: ASD={len(x)} NT={len(y)}')
  p=float(ttest_ind(x,y,equal_var=False).pvalue)
  try:pm=float(mannwhitneyu(x,y,alternative='two-sided',method='auto').pvalue)
  except ValueError:pm=np.nan
  results.append(dict(cell_type=ct,n_asd=len(x),n_neurotypical=len(y),mean_asd=x.mean(),mean_neurotypical=y.mean(),delta_percentage_points=100*(x.mean()-y.mean()),relative_change_pct=100*(x.mean()/y.mean()-1) if y.mean()>0 else np.nan,p_welch=p,p_mann_whitney=pm))
 stats=pd.DataFrame(results)
 for pcol in ['p_welch','p_mann_whitney']:
  stats['FDR_'+pcol]=np.nan
  valid=stats[pcol].notna()
  if valid.any():stats.loc[valid,'FDR_'+pcol]=multipletests(stats.loc[valid,pcol],method='fdr_bh')[1]
 stats=stats.sort_values('p_welch');stats.to_csv(OUT/'wholeSection'/'ASD_vs_Neurotypical_stats.tsv',sep='\t',index=False)
 pd.Series(dict(reference_url=SOURCE,signature_count=len(types),N_cells_per_location=N_CELLS_PER_LOCATION,detection_alpha=20,max_epochs=MAX_EPOCHS,device=DEVICE,n_samples=16,n_asd=10,n_nt=6,warning='Exploratory statistical tests: repeat sections from the same animal treated as independent; potential pseudoreplication')).to_csv(OUT/'qc'/'run_parameters.tsv',sep='\t',header=False)
 log('DONE. Results: %s',OUT)
 log('Nominal Welch p<0.05: %d; FDR<0.05: %d',int((stats.p_welch<.05).sum()),int((stats.FDR_p_welch<.05).sum()))
if __name__=='__main__':main()
