#!/usr/bin/env python3
"""Score observed chemical responses against frozen external signed references."""
import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ.setdefault(key,'2')
from pathlib import Path
from datetime import datetime
import argparse,json,hashlib,traceback,zipfile
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
RESPONSE_SHA='24e24107a57a60d487a8960c7e85f5b70c0a04aa8c77bb522736847342420b03'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def cosine_reversal(x,d):
    den=np.linalg.norm(x,axis=1)*np.linalg.norm(d)
    return np.divide(-(x@d),den,out=np.full(len(x),np.nan),where=den>0)

def rank_reversal(x,d):
    # Average ties; negative Spearman over the fixed mapped gene universe.
    xr=pd.DataFrame(x).rank(axis=1,method='average').to_numpy()
    dr=pd.Series(d).rank(method='average').to_numpy()
    xr-=xr.mean(axis=1,keepdims=True);dr-=dr.mean()
    return cosine_reversal(xr,dr)

def aggregate_responses(r,identities):
    names,inv,count=np.unique(identities,return_inverse=True,return_counts=True)
    ans=np.zeros((2,2,len(names),r.shape[-1]))
    for c in range(2):
        for half in range(2):np.add.at(ans[c,half],inv,r[c,half])
    return names,ans/count[None,None,:,None]

def finite_or_none(value):
    return float(value) if np.isfinite(value) else None

def evaluate(r,groups,genes,refs,manifest):
    names,identity_r=aggregate_responses(r,groups)
    allrows=[];summary=[]
    for ref,q in refs.groupby('reference',sort=True):
        q=q.set_index('gene_id');mask=genes.pr_gene_id.isin(q.index).to_numpy()
        d=q.loc[genes.pr_gene_id[mask],'disease_effect'].to_numpy(float)
        cell=q.cell.iloc[0];ci=['MCF7','PC3'].index(cell)
        x=identity_r[ci][:,:,mask]
        scores=[rank_reversal(x[h],d) for h in range(2)]
        cosine=[cosine_reversal(x[h],d) for h in range(2)]
        meanrank=rank_reversal(x.mean(axis=0),d)
        # Fixed top 50 up and down genes by external reference effect, exact ID tie-break.
        ids=genes.pr_gene_id[mask].to_numpy()
        ordered=np.lexsort((ids,d));down=ordered[d[ordered]<0][:50];up=ordered[d[ordered]>0][-50:]
        if len(up)<25 or len(down)<25:raise ValueError('Insufficient signed-tail coverage')
        tail=[(x[h][:,down].mean(axis=1)-x[h][:,up].mean(axis=1)) for h in range(2)]
        for i,name in enumerate(names):
            m=manifest[manifest.identity_group.eq(name)]
            label=' | '.join(sorted(set(m.pert_inames.astype(str))))
            allrows.append(dict(reference=ref,cell=cell,identity_group=name,compound_names=label,
                member_ids='|'.join(m.pert_id),n_member_ids=len(m),n_genes=len(d),reference_role=q.role.iloc[0],
                reference_replication=q.replication.iloc[0],reversal_spearman_half0=scores[0][i],
                reversal_spearman_half1=scores[1][i],reversal_spearman_halfmean=meanrank[i],
                reversal_cosine_half0=cosine[0][i],reversal_cosine_half1=cosine[1][i],
                signed_tail_score_half0=tail[0][i],signed_tail_score_half1=tail[1][i],
                opposes_both_technical_halves=bool(scores[0][i]>0 and scores[1][i]>0),
                reinforces_both_technical_halves=bool(scores[0][i]<0 and scores[1][i]<0),
                response_rms_half0=float(np.sqrt(np.mean(x[0,i]**2))),response_rms_half1=float(np.sqrt(np.mean(x[1,i]**2)))))
        summary.append(dict(reference=ref,cell=cell,n_genes=len(d),n_identities=len(names),
            positive_both_halves=int(np.sum((scores[0]>0)&(scores[1]>0))),
            negative_both_halves=int(np.sum((scores[0]<0)&(scores[1]<0))),
            nonfinite_scores=int(np.sum(~np.isfinite(np.array(scores)))),
            top_tail_up_genes=len(up),top_tail_down_genes=len(down),
            half_score_agreement_descriptive=finite_or_none(pd.Series(scores[0]).corr(pd.Series(scores[1]),method='pearson'))))
    # Stable identity ordering, not efficacy/candidate ranking.
    return pd.DataFrame(allrows).sort_values(['reference','identity_group']),summary

def plot(scores,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    refs=list(scores.reference.unique());fig,axes=plt.subplots(1,len(refs),figsize=(14,4.6))
    for ax,ref in zip(np.atleast_1d(axes),refs):
        q=scores[scores.reference.eq(ref)]
        ax.scatter(q.reversal_spearman_half0,q.reversal_spearman_half1,s=5,alpha=.22,color='#1b6fa8',rasterized=True)
        ax.axhline(0,color='.5',lw=.7);ax.axvline(0,color='.5',lw=.7)
        ax.plot([-1,1],[-1,1],ls=':',color='.6',lw=.7)
        ax.set_xlim(-1,1);ax.set_ylim(-1,1)
        ax.set_xlabel('Reversal alignment: technical half 0')
        ax.set_ylabel('Reversal alignment: technical half 1')
        ax.set_title(ref.replace('_',' '),fontsize=9)
    fig.suptitle('Exploratory transcriptional opposition; phenotype reversal is not established',fontsize=11)
    fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(out/('reversal_score_halves.'+ext),dpi=180)
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--responses',type=Path,default=Path('/mnt/c/Wei/ICE-PRN/ICE_PRN_LINCS_DRUG_EXTRACTION_v1/drug_response_export_v1/responses.npz'))
    ap.add_argument('--out',type=Path,default=ROOT/('results_'+datetime.now().strftime('%Y%m%d_%H%M%S')))
    a=ap.parse_args()
    if a.out.exists():ap.error('Use a fresh output folder')
    a.out.mkdir(parents=True)
    status={'status':'RUNNING','models_fitted':False,'treatment_ranking_established':False,'phenotype_reversal_established':False}
    try:
        for name,h in json.loads((ROOT/'PACKAGE_HASHES.json').read_text()).items():
            if sha(ROOT/name)!=h:raise ValueError('Frozen package changed: '+name)
        if sha(a.responses)!=RESPONSE_SHA:raise ValueError('Response export hash mismatch')
        genes=pd.read_csv(ROOT/'measured_landmarks.tsv',sep='\t',dtype=str)
        manifest=pd.read_csv(ROOT/'compound_manifest.tsv',sep='\t',keep_default_na=False)
        refs=pd.read_csv(ROOT/'references/signed_references.tsv',sep='\t',dtype={'gene_id':str})
        with np.load(a.responses,allow_pickle=False) as z:
            for actual,expected in [(z['gene_ids'].astype(str),genes.pr_gene_id),(z['pert_ids'],manifest.pert_id),(z['identity_groups'],manifest.identity_group),(z['folds'],manifest.fold),(z['cells'],['MCF7','PC3'])]:
                if not np.array_equal(actual,np.asarray(expected)):raise ValueError('Frozen array/metadata axis mismatch')
            r=z['responses'].astype(float);groups=z['identity_groups']
            if z['profile_counts'].shape!=(2,2,4902) or z['profile_counts'].sum()!=53635 or np.any(z['profile_counts']<1):raise ValueError('Profile count mismatch')
        if r.shape!=(2,2,4902,978) or not np.isfinite(r).all():raise ValueError('Invalid response matrix')
        scores,summary=evaluate(r,groups,genes,refs,manifest)
        scores.to_csv(a.out/'identity_reversal_scores.tsv',sep='\t',index=False)
        pd.DataFrame(summary).to_csv(a.out/'reference_summary.tsv',sep='\t',index=False)
        refqc=json.loads((ROOT/'references/REFERENCE_QC.json').read_text())
        status.update(status='COMPLETE_EXPLORATORY_OPPOSITION_ONLY',chemical_identities=4185,source_response_sha256=RESPONSE_SHA,
            references=summary,reference_qc=refqc,
            limitations=['MCF7 versus MCF10A is a non-isogenic cell-line contrast, not a pure cancer effect',
                'PC3 versus RWPE1 has one deposited sample per line and no biological-replication uncertainty',
                'The prostate cell-line and patient-tissue signatures have essentially no agreement; scores remain separate',
                'L1000 z-scores cannot be added to external log fold changes to estimate distance to healthy',
                'Technical half agreement does not establish biological replication or clinical benefit',
                'Generic stress, toxicity and proliferation effects may generate reversal alignment; no independent viability endpoint tested'])
        plot(scores,a.out)
        from compare_contexts import compare
        paired=compare(scores,a.out)
        (a.out/'CONTEXT_COMPARISON.json').write_text(json.dumps(paired,indent=2,allow_nan=False)+'\n')
        status['context_comparison']=paired
        for n in ['ANALYSIS_PLAN.md','references/REFERENCE_QC.json']:(a.out/Path(n).name).write_bytes((ROOT/n).read_bytes())
    except Exception:status.update(status='FAILED_RETURN_FOR_REVIEW',error=traceback.format_exc())
    (a.out/'SUMMARY.json').write_text(json.dumps(status,indent=2,allow_nan=False)+'\n')
    with zipfile.ZipFile(a.out/'RETURN_TO_CHAT.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in a.out.iterdir():
            if p.suffix in ['.json','.tsv','.md','.png','.pdf','.svg']:z.write(p,p.name)
    print(json.dumps(status,indent=2));print('Upload:',a.out/'RETURN_TO_CHAT.zip',flush=True)
    if status['status'].startswith('FAILED'):raise SystemExit(1)
if __name__=='__main__':main()
