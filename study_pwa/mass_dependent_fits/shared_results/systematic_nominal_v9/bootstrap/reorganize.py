import os

for t in ['010020','0200325','0325050','050075','075100']:
    base=f"bootstrap_t{t}"
    ofolder=f"results_t{t}"
    
    os.system(f'mkdir -p {ofolder}')
    #os.system(f'ln -snfr {base}/bootstrap_0/fit.log {ofolder}')
    j=0
    for i in range(100):
        if os.path.exists(f'{base}/bootstrap_{i}/etapi_result.fit'):  
            cmd=f'ln -snfr {base}/bootstrap_{i}/etapi_result.fit {ofolder}/etapi_result_{j}.fit'
            j+=1
            print(cmd)
            os.system(cmd)
