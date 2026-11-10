from experiments import *

def extended():
    rows=[];stress=[]
    t=np.arange(600,dtype=float)
    for seed in range(60):
        rng=np.random.default_rng(150000+seed)
        rho=100+3*np.sin(2*np.pi*t/300);theta=3+.0006*t+.01*np.sin(2*np.pi*t/200)
        r=rho+.2*rng.normal(size=600);b=theta+.003*rng.normal(size=600)
        idx=(rng.random(600)<.05)&(t>30);r[idx]+=rng.choice([-1,1],idx.sum())*2.4
        idx=(rng.random(600)<.05)&(t>30);b[idx]+=rng.choice([-1,1],idx.sum())*.036
        rr=run(t,(r-100)/.2);bb=run(t,(b-3)/.003)
        for method in METHODS:
            rh=rr[method]*.2+100;th=bb[method]*.003+3
            # theta is clockwise from north, so east = range*sin(theta).
            er=np.sqrt((rh*np.sin(th)-rho*np.sin(theta))**2+(rh*np.cos(th)-rho*np.cos(theta))**2)
            rows.append(dict(seed=seed,method=method,position_rmse=float(np.sqrt(np.mean(er[30:]**2)))))
        if seed==0:
            fig,axs=plt.subplots(1,2)
            axs[0].plot(rho*np.sin(theta),rho*np.cos(theta),'k',label='Truth')
            for method in ['Raw/hold','Trend','Huber-KF']:
                rh=rr[method]*.2+100;th=bb[method]*.003+3
                axs[0].plot(rh*np.sin(th),rh*np.cos(th),lw=.7,alpha=.7,label=method)
            axs[0].set_xlabel('East (m)');axs[0].set_ylabel('North (m)');axs[0].axis('equal');axs[0].legend(fontsize=7)
            w=wrap(b);safe=trend(t,w,angle=True)[0];unsafe=trend(t,w)[0]
            axs[1].plot(t,wrap(safe-theta),label='Ring residual extension');axs[1].plot(t,wrap(unsafe-theta),label='Scalar residual core')
            axs[1].set_xlabel('Time (s)');axs[1].set_ylabel('Circular bearing error (rad)');axs[1].legend(fontsize=7);savefig('geometry_and_wrap')
    dumpcsv('synthetic_position_trials.csv',rows)
    for seed in range(30):
        rng=np.random.default_rng(160000+seed);truth=100+8*(t>=300);y=truth+rng.normal(size=600);res=run(t,y)
        for method,x in res.items():
            hit=np.flatnonzero(abs(x[300:]-108)<.8)
            stress.append(dict(seed=seed,method=method,post_step_rmse=float(np.sqrt(np.mean((x[300:360]-truth[300:360])**2))),first_entry_s=int(hit[0]) if len(hit) else 300))
        if seed==0:
            fig,ax=plt.subplots();ax.plot(t,y,'.',color='.75',ms=2,label='Observation');ax.plot(t,truth,'k',label='Truth')
            for method in ['Trend','No-gate','Huber-KF']:ax.plot(t,res[method],label=method)
            ax.set_xlim(280,380);ax.set_xlabel('Time (s)');ax.set_ylabel('Step response (normalized units)');ax.legend();savefig('step_response')
    dumpcsv('step_response_trials.csv',stress)
    # Detection is evaluated only where injected labels are defined, after warmup.
    detection=[]
    for seed in range(60):
        t,z,y,bad=scenario(81000+seed,'Outliers');_,rej,_,_=trend(t,y);use=t>=30
        tp=int(np.sum(rej&bad&use));fp=int(np.sum(rej&~bad&use));fn=int(np.sum(~rej&bad&use));tn=int(np.sum(~rej&~bad&use))
        detection.append(dict(seed=seed,tp=tp,fp=fp,fn=fn,tn=tn))
    dumpcsv('detection_trials.csv',detection)
    # Native equation validation accounts for six-digit debug output rounding.
    er=list(csv.DictReader((OUT/'results/logged_equation_checks.csv').open()))
    assert max(float(x['x_error']) for x in er)<2e-4
    assert max(float(x['v_error']) for x in er)<2e-6
    # Per-case and per-gap scores, missing epochs distinguished from observed epochs.
    gaps=[]
    for seed in range(60):
        t,z,y,_=scenario(81000+seed,'Combined');res=run(t,y)
        for method,x in res.items():
            for label,mask in [('observed',(t>=30)&np.isfinite(y)),('short_gap',(t>=180)&(t<190)),('long_gap',(t>=390)&(t<420))]:
                gaps.append(dict(seed=seed,method=method,subset=label,rmse=float(np.sqrt(np.mean((x[mask]-z[mask])**2)))))
    dumpcsv('gap_subset_trials.csv',gaps)
    a=np.ones(1200);tt=np.arange(1200,dtype=float);timing=[]
    for method in METHODS:
        elapsed=[]
        for _ in range(5):
            st=time.perf_counter()
            if method=='Trend':trend(tt,a)
            elif method=='No-alternation':trend(tt,a,alt=False)
            elif method=='No-gate':trend(tt,a,gate=False)
            elif method=='CV-KF':kfilter(tt,a)
            elif method=='Huber-KF':kfilter(tt,a,True)
            else:
                if method=='Median-7': [np.median(a[max(0,i-6):i+1]) for i in range(len(a))]
                else: a.copy()
            elapsed.append((time.perf_counter()-st)*1e6/len(a))
        timing.append(dict(method=method,median_us_per_epoch=float(np.median(elapsed))))
    dumpcsv('runtime.csv',timing)
    # Flow chart in vector and bitmap, all labels editable in this source.
    from matplotlib.patches import FancyBboxPatch
    fig,ax=plt.subplots(figsize=(6.69,5.8));ax.set_xlim(0,10);ax.set_ylim(0,10);ax.axis('off')
    boxes=[(5,9.3,'Range / bearing and timestamp'),(5,7.9,'Validate chronology; collect six past observations'),(5,6.5,'Initialize at latest epoch: median slope + MAD scale'),(5,5.1,'Damped trend prediction and innovation'),(2.5,3.6,'Observation missing\nQuery: extrapolate / freeze / unavailable'),(7.5,3.6,'Observation available\nAlternation score and pre-update scale gate'),(6,2,'Reject: prediction only\nAccept: Huber-weighted trend update'),(5,.5,'Update scale and timestamp; output trend and status')]
    for x,y,label in boxes:
        width=8.4 if x==5 else 4.5
        ax.add_patch(FancyBboxPatch((x-width/2,y-.42),width,.84,boxstyle='round,pad=.08',facecolor='#edf3f7',edgecolor='#44677c'))
        ax.text(x,y,label,ha='center',va='center',fontsize=8)
    def arrow(a,b):ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':'#44677c','lw':1.2})
    for y in [8.85,7.45,6.05]:arrow((5,y),(5,y-.5))
    arrow((4.5,4.65),(2.5,4.08));arrow((5.5,4.65),(7.5,4.08));arrow((7.5,3.15),(6,2.48));arrow((2.5,3.15),(2.5,1));arrow((6,1.55),(5,.98))
    plt.tight_layout()
    for ext in ['pdf','png']:fig.savefig(OUT/'figures'/f'flowchart.{ext}',dpi=300,bbox_inches='tight')
    plt.close()
    print('Extended experiments and equation checks passed')
if __name__=='__main__':extended()
