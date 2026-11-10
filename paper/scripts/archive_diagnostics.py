from experiments import *
fig,axs=plt.subplots(2,2)
rows=[]
for i,kind in enumerate(['range','bearing']):
    _,t,y,_=read_aligned(kind);available=np.isfinite(y);starts=[];lengths=[];start=None
    for k,v in enumerate(np.r_[available,True]):
        if not v and start is None:start=k
        if v and start is not None:starts.append(start);lengths.append(k-start);start=None
    axs[i,0].plot(t,available.astype(float),lw=.5);axs[i,0].axvspan(0,t[299],alpha=.3,color='#c45b47',label='Calibration epochs 0-299')
    axs[i,0].set_ylim(-.1,1.25);axs[i,0].set_ylabel(kind+' available');axs[i,0].set_xlabel('Elapsed time (s)');axs[i,0].legend(fontsize=7)
    axs[i,1].hist(lengths,bins=20,color='#73a9c2');axs[i,1].set_xlabel('Consecutive missing epochs');axs[i,1].set_ylabel('Gap count')
    rows.append(dict(channel=kind,calibration_valid=int(available[:300].sum()),calibration_last_index=299,evaluation_first_index=300,gap_count=len(lengths),longest_missing_run=max(lengths),observed=int(available.sum()),n=len(y)))
savefig('archive_diagnostics');dumpcsv('archive_protocol.csv',rows)
print(rows)
