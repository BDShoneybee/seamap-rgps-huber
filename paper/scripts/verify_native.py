from experiments import *
import subprocess,io,shutil
compiler=shutil.which('g++')
if compiler is None:raise SystemExit('g++ is required for native parity verification')
exe=OUT/'results'/('trend_core.exe' if platform.system()=='Windows' else 'trend_core')
subprocess.run([compiler,'-std=c++14','-O2',str(OUT/'scripts/trend_core.cpp'),'-o',str(exe)],check=True)
tests=[]
for seed in [81000,81009,81059]:
    t,z,y,_=scenario(seed,'Combined');ref,rej,_,_=trend(t,y)
    inp='time,value\n'+'\n'.join(f'{a:.17g},{b:.17g}' for a,b in zip(t,y))
    p=subprocess.run([str(exe)],input=inp,text=True,capture_output=True,check=True)
    rows=list(csv.DictReader(io.StringIO(p.stdout)));x=np.array([float(r['estimate']) for r in rows]);r=np.array([int(q['rejected']) for q in rows],bool)
    err=float(np.nanmax(abs(x-ref)));assert err<1e-10;assert np.array_equal(r,rej)
    tests.append(dict(seed=seed,max_absolute_error=err,rejection_match=True))
(OUT/'results/native_parity.json').write_text(json.dumps(tests,indent=2))
print('Native/Python parity passed:',len(tests),'sequences')
