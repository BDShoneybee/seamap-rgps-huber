"""Reproducible audit, Monte Carlo and chronological archive replay. Python + NumPy + Matplotlib."""
from pathlib import Path
import csv, json, re, hashlib, math, time, platform
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'paper'
SRC=next((p for p in ROOT.iterdir() if p.name.startswith('0460')),None)
plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def savefig(name):
    plt.gcf().set_size_inches(6.69,4.6); plt.tight_layout()
    for ext in ['pdf','png']: plt.savefig(OUT/'figures'/f'{name}.{ext}',dpi=300,bbox_inches='tight')
    plt.close()
def dumpcsv(name,rows):
    if not rows:return
    with (OUT/'results'/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def mad(a):return float(np.median(np.abs(a-np.median(a))))
def wrap(a):return (a+np.pi)%(2*np.pi)-np.pi

def trend(t,y,alt=True,gate=True,huber=True,tau=10.,angle=False):
    """DLPPreProc Huber core, with explicit six-observation causal warmup.
    No polynomial wrapper; missing epochs use non-mutating deployed query rule.
    Angular ring residuals are an explicit extension, not claimed original behavior.
    """
    n=len(t);out=np.full(n,np.nan);reject=np.zeros(n,bool);prior=out.copy();scales=out.copy()
    buf=[];x=v=s=s0=a=0.;last=prev=0.;pdy=None;ready=False
    for i,(ti,z) in enumerate(zip(t,y)):
        if not ready:
            if np.isfinite(z):buf.append((ti,z))
            if len(buf)==6:
                bt,by=np.array(buf).T; dy=np.diff(by);dy=wrap(dy) if angle else dy
                s0=max(.35*1.4826*mad(dy),1e-4);s=s0
                v=float(np.median(dy/np.diff(bt)));x=z;last=ti;prev=z;ready=True
            out[i]=z;continue
        dt=max(ti-last,1e-3)
        vp=v*math.exp(-dt/12);xp=x+vp*dt
        prior[i]=xp
        if not np.isfinite(z):
            out[i]=x if dt>8 else xp
            if dt>100:out[i]=np.nan
            continue
        nu=wrap(z-xp) if angle else z-xp;dy=wrap(z-prev) if angle else z-prev
        if pdy is not None:
            a=min(1.,.88*a+.12) if dy*pdy<0 else max(0.,.90*a-.06)
        pdy=dy;prev=z
        sn=max(.97*s+.03*abs(nu),.35*s0)
        reject[i]=gate and abs(nu)>4*s
        w=min(1.,2.5*s/max(abs(nu),1e-12)) if huber else 1.
        aa=a if alt else 0.
        kx=(1-math.exp(-dt/tau))*(1-.45*aa)
        kv=(1-math.exp(-dt/24))/dt*(1-.60*aa)
        x=xp if reject[i] else xp+kx*w*nu
        v=vp if reject[i] else vp+kv*w*nu
        if angle:x=wrap(x)
        s=sn;last=ti;out[i]=x;scales[i]=s
    return out,reject,prior,scales

def kfilter(t,y,robust=False):
    out=np.full(len(y),np.nan);x=None;P=np.eye(2)*10;I=np.eye(2)
    for i,z in enumerate(y):
        if x is None:
            if np.isfinite(z):x=np.array([z,0.]);out[i]=z
            continue
        dt=t[i]-t[i-1];F=np.array([[1.,dt],[0.,1.]])
        Q=.0025*np.array([[dt**3/3,dt**2/2],[dt**2/2,dt]])
        x=F@x;P=F@P@F.T+Q
        if np.isfinite(z):
            nu=z-x[0];w=min(1.,2.5*math.sqrt(P[0,0]+1)/max(abs(nu),1e-12)) if robust else 1.
            R=1/max(w,1e-9);K=P[:,0]/(P[0,0]+R);x+=K*nu
            A=I-np.outer(K,[1.,0.]);P=A@P@A.T+np.outer(K,K)*R
        out[i]=x[0]
    return out

METHODS=['Raw/hold','Median-7','CV-KF','Huber-KF','Trend','No-alternation','No-gate']
def run(t,y):
    hold=[];med=[];hist=[];last=np.nan
    for z in y:
        if np.isfinite(z):last=z;hist.append(z)
        hold.append(last);med.append(np.median(hist[-7:]) if hist else np.nan)
    return {'Raw/hold':np.array(hold),'Median-7':np.array(med),'CV-KF':kfilter(t,y),
            'Huber-KF':kfilter(t,y,True),'Trend':trend(t,y)[0],
            'No-alternation':trend(t,y,alt=False)[0],'No-gate':trend(t,y,gate=False)[0]}

def scenario(seed,case,p=.05):
    r=np.random.default_rng(seed);t=np.arange(600,dtype=float)
    truth=100+2*np.sin(2*np.pi*t/300)+.002*t
    if case=='Maneuver':truth+=.10*np.maximum(t-260,0)-.10*np.maximum(t-340,0)
    y=truth+r.normal(size=len(t));bad=np.zeros(len(t),bool)
    if case in ['Alternating','Combined']:y+=1.5*(-1.)**np.arange(len(t))
    if case in ['Outliers','Combined']:
        bad=(r.random(len(t))<p)&(t>=30);y[bad]+=r.choice([-1,1],bad.sum())*12
    if case in ['Gaps','Combined']:y[180:190]=np.nan;y[390:420]=np.nan
    return t,truth,y,bad

def simulation():
    rows=[];cases=['Gaussian','Alternating','Outliers','Gaps','Combined','Maneuver']
    for case in cases:
        for seed in range(60):
            t,z,y,bad=scenario(81000+seed,case);res=run(t,y)
            for method,x in res.items():
                valid=(t>=30)&np.isfinite(x);e=x[valid]-z[valid]
                rows.append(dict(case=case,seed=seed,method=method,rmse=float(np.sqrt(np.mean(e**2))),mae=float(np.mean(abs(e))),p95=float(np.percentile(abs(e),95))))
        print('Finished Monte Carlo:',case,flush=True)
    dumpcsv('monte_carlo_trials.csv',rows)
    summary=[]
    for case in cases:
        for m in METHODS:
            a=np.array([r['rmse'] for r in rows if r['case']==case and r['method']==m]);se=a.std(ddof=1)/np.sqrt(len(a))
            summary.append(dict(case=case,method=m,mean=float(a.mean()),ci95=float(2.001*se)))
    dumpcsv('simulation_summary.csv',summary)
    fig,axs=plt.subplots(2,3)
    for ax,case in zip(axs.flat,cases):
        a=[s for s in summary if s['case']==case];ax.bar(np.arange(7),[x['mean'] for x in a],yerr=[x['ci95'] for x in a],color=['.7','.6','#73a9c2','#52a88d','#c45b47','#d5a25b','#9d89b5'])
        ax.set_xticks(np.arange(7));ax.set_xticklabels(['Raw','Med','KF','HKF','Tr','-Alt','-Gate'],rotation=45);ax.set_title(case);ax.set_ylabel('RMSE (normalized units)')
    savefig('simulation_comparison')
    t,z,y,bad=scenario(81000,'Combined');res=run(t,y)
    fig,axs=plt.subplots(2,1,sharex=True)
    axs[0].plot(t,y,'.',ms=2,color='.65',label='Observation');axs[0].plot(t,z,'k',lw=1.5,label='Truth')
    for m in ['Trend','Huber-KF']:axs[0].plot(t,res[m],lw=1,label=m);axs[1].plot(t,res[m]-z,label=m)
    for ax in axs:
        ax.axvspan(180,190,color='.8',alpha=.5);ax.axvspan(390,420,color='.8',alpha=.5);ax.legend(ncol=2,fontsize=8)
    axs[0].set_ylabel('Observation');axs[1].set_ylabel('Estimation error');axs[1].set_xlabel('Time (s)');savefig('simulation_timeseries')
    # Paired, run-level bootstrap comparison. No epoch-level pseudo-replication.
    paired=[];rng=np.random.default_rng(81244)
    for case in cases:
        tr=np.array([r['rmse'] for r in rows if r['case']==case and r['method']=='Trend'])
        for method in ['Raw/hold','Huber-KF','No-alternation','No-gate']:
            base=np.array([r['rmse'] for r in rows if r['case']==case and r['method']==method]);d=tr-base
            bs=np.mean(rng.choice(d,(3000,len(d))),axis=1)
            paired.append(dict(case=case,comparator=method,delta=float(d.mean()),lower=float(np.quantile(bs,.025)),upper=float(np.quantile(bs,.975))))
    dumpcsv('paired_bootstrap.csv',paired)
    sensitivity=[]
    for tau in [2,5,10,20]:
        for p in [0,.02,.05,.10,.20]:
            scores=[]
            for seed in range(20):
                t,z,y,_=scenario(92000+seed,'Combined',p);x=trend(t,y,tau=tau)[0];scores.append(np.sqrt(np.mean((x[30:]-z[30:])**2)))
            sensitivity.append(dict(tau=tau,p=p,rmse=float(np.mean(scores))))
    dumpcsv('sensitivity.csv',sensitivity)
    fig,ax=plt.subplots()
    for tau in [2,5,10,20]:
        a=[s for s in sensitivity if s['tau']==tau];ax.plot([s['p'] for s in a],[s['rmse'] for s in a],'o-',label=f'tau = {tau} s')
    ax.set_xlabel('Injected outlier probability');ax.set_ylabel('RMSE (normalized units)');ax.legend();savefig('sensitivity')
    return summary

def read_aligned(kind):
    if SRC is None:
        p=OUT/'data'/f'{kind}_archive.npz'
        with np.load(p) as a:t,y,old=a['t'],a['y'],a['archived_prediction']
        return p,t,y,old
    p=next(SRC.rglob(f'{kind}_aligned.csv'));rows=list(csv.DictReader(p.open(encoding='utf-8-sig')))
    t=[];ys=[];pred=[];offset=0;last=None
    for row in rows:
        a=list(map(float,row['Time'].split(':')));v=sum(x*60**j for j,x in enumerate(reversed(a)))
        # Bearing archive stores mm:ss; unwrap hour transitions before differencing.
        period=3600 if len(a)==2 else 86400
        if last is not None and v+offset<last-period/2:offset+=period
        v+=offset;t.append(v);last=v
        ys.append(float(row['Original']) if row['Original'] else np.nan);pred.append(float(row['Predicted']) if row['Predicted'] else np.nan)
    t=np.array(t);return p,t-t[0],np.array(ys),np.array(pred)

def replay():
    stats=[];inj=[];fig,axs=plt.subplots(2,2)
    for j,kind in enumerate(['range','bearing']):
        p,t,y,old=read_aligned(kind);finite=np.isfinite(y);warm=np.flatnonzero(finite & (np.arange(len(y))<300))
        assert warm[-1]<300, 'Calibration must precede evaluation'
        d=np.diff(y[warm]);scale=max(1.4826*mad(d)/np.sqrt(2),1e-4)
        center=float(np.median(y[warm]));yn=(y-center)/scale
        res=run(t,yn);mask=(np.arange(len(t))>=300)&finite
        for m,x in res.items():
            # Next recorded observation diagnostic; not independent truth accuracy.
            ix=np.flatnonzero(mask[:-1]&finite[1:]&(np.diff(t)>0)&(np.diff(t)<2))
            er=(x[ix]-yn[ix+1])*scale
            stats.append(dict(channel=kind,method=m,n=len(y),observed=int(finite.sum()),duration=float(t[-1]),dt_median=float(np.median(np.diff(t))),noise_scale=scale,next_rmse=float(np.sqrt(np.mean(er**2))),roughness=float(np.nanmedian(abs(np.diff(x[300:]))))*scale))
        for seed in range(30):
            rng=np.random.default_rng(130000+seed);yc=yn.copy();bad=(rng.random(len(y))<.02)&mask;yc[bad]+=rng.choice([-1,1],bad.sum())*12
            xc=run(t,yc)
            for m in METHODS:
                delta=(xc[m]-res[m])[300:];delta=delta[np.isfinite(delta)]
                inj.append(dict(channel=kind,seed=seed,method=m,divergence_rmse=float(np.sqrt(np.mean(delta**2)))))
        for m in ['Trend','Huber-KF']:
            axs[j,0].plot(t,(res[m]*scale+center),lw=.8,label=m)
        axs[j,0].plot(t,y,'.',color='.65',ms=1,label='Archived input',zorder=0)
        axs[j,0].set_ylabel(kind+(' (m)' if kind=='range' else ' (rad)'));axs[j,0].set_xlabel('Elapsed archive time (s)');axs[j,0].legend(fontsize=7)
        # Zoom chosen by fixed index, not by algorithm quality.
        sl=slice(1000,1250);axs[j,1].plot(t[sl],y[sl],'.',color='.6',ms=2,label='Input')
        for m in ['Trend','Huber-KF']:axs[j,1].plot(t[sl],(res[m]*scale+center)[sl],label=m,lw=1)
        axs[j,1].set_xlabel('Elapsed archive time (s)');axs[j,1].set_ylabel('Fixed-index detail');axs[j,1].legend(fontsize=7)
        outrows=[dict(time=float(t[i]),input=float(y[i]),trend=float(res['Trend'][i]*scale+center),huber_kf=float(res['Huber-KF'][i]*scale+center)) for i in range(len(t))]
        dumpcsv(f'{kind}_replay.csv',outrows)
        np.savez_compressed(OUT/'data'/f'{kind}_archive.npz',t=t,y=y,archived_prediction=old)
    savefig('archive_replay');dumpcsv('archive_metrics.csv',stats);dumpcsv('injection_trials.csv',inj)
    fig,axs=plt.subplots(1,2)
    for ax,kind in zip(axs,['range','bearing']):
        vals=[np.mean([r['divergence_rmse'] for r in inj if r['channel']==kind and r['method']==m]) for m in METHODS]
        ax.bar(range(7),vals);ax.set_xticks(range(7));ax.set_xticklabels(['Raw','Med','KF','HKF','Tr','-Alt','-Gate'],rotation=45);ax.set_title(kind);ax.set_ylabel('Paired output divergence / calibration scale')
    savefig('archive_injection');return stats

def log_audit():
    if SRC is None:
        events=list(csv.DictReader((OUT/'results/logged_event_counts.csv').open()))
        checks=list(csv.DictReader((OUT/'results/logged_equation_checks.csv').open()))
        for r in checks:
            for field,state,gain in [('x_error','x','kx'),('v_error','v','kv')]:
                error=abs(float(r[state+'_new'])-(float(r[state+'_pred'])+float(r[gain])*float(r['weight'])*float(r['innovation'])))
                assert abs(error-float(r[field]))<1e-12
        return events,checks
    # Full logs only: filtered extracts are duplicates and excluded.
    ps=[SRC/'huber算法/log-20260326-测试RGPS预处理器算法log.txt']+sorted((SRC/'huber算法').glob('20260309*/*log.txt'))
    num=r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)';field=re.compile(r'(\w+)=\s*'+num)
    rows=[];checks=[]
    for p in ps:
        if not p.exists():continue
        counter=Counter();pending={};times={}
        for line in p.read_text(encoding='utf-8',errors='replace').splitlines():
            match=re.search(r'\[(Huber|PreProc)\]\[([^]]+)\]',line)
            if not match:continue
            f={k:float(v) for k,v in field.findall(line)}
            if 'node' not in f or 'obsType' not in f:continue
            key=(int(f['node']),int(f['obsType']));event=match[2];counter[key+(event,)]+=1
            if event=='Obs':pending[key]=f;times.setdefault(key,[]).append(f['t'])
            if event=='Accept-Inlier' and key in pending:
                q=pending.pop(key)
                # Direct verification of state update from logged prior and logged gains.
                if all(k in f for k in ['Kx','Kv','w','xNew','vNew']):
                    checks.append(dict(session=p.parent.name,node=key[0],obs_type=key[1],x_pred=q['xPred'],v_pred=q['vPred'],innovation=q['nu'],kx=f['Kx'],kv=f['Kv'],weight=f['w'],x_new=f['xNew'],v_new=f['vNew'],x_error=abs(f['xNew']-(q['xPred']+f['Kx']*f['w']*q['nu'])),v_error=abs(f['vNew']-(q['vPred']+f['Kv']*f['w']*q['nu']))))
        for (node,obs,event),n in counter.items():rows.append(dict(session=p.parent.name if p.name=='log.txt' else p.stem,node=node,obs_type=obs,event=event,count=n))
    dumpcsv('logged_event_counts.csv',rows);dumpcsv('logged_equation_checks.csv',checks)
    return rows,checks

def checks():
    t,z,y,_=scenario(18,'Combined');a=trend(t,y)[0];yp=y.copy();yp[400:]+=10000;b=trend(t,yp)[0]
    assert np.allclose(a[:400],b[:400],equal_nan=True)
    th=3+.002*t;ya=wrap(th);aa=trend(t,ya,angle=True)[0];bb=trend(t,ya+2*np.pi,angle=True)[0]
    assert np.nanmax(abs(wrap(aa[30:]-bb[30:])))<1e-8
    assert np.isfinite(a).all()
    report={'prefix_causality':'pass','angle_2pi_invariance':'pass','gap_finite_outputs_under_100s':'pass','python':platform.python_version(),'numpy':np.__version__}
    (OUT/'results/verification.json').write_text(json.dumps(report,indent=2));return report

def inventory():
    if SRC is None:return
    rows=[]
    for p in sorted(SRC.rglob('*')):
        if p.is_file():rows.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest()))
    dumpcsv('source_manifest.csv',rows)

if __name__=='__main__':
    start=time.perf_counter();inventory();checks();logs,eq=log_audit();simulation();replay()
    (OUT/'results/run_info.json').write_text(json.dumps({'elapsed_s':time.perf_counter()-start,'log_events':len(logs),'equation_checks':len(eq),'base_mc_seed':81000,'mc_runs_per_case':60,'calibration_samples':300},indent=2))
    print('Completed all experiments',flush=True)
