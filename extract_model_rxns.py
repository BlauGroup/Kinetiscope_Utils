from monty.serialization import loadfn

model_AutoTSs_path = "G:/My Drive/CRNs/AutoTS_output/model_TSs.json"
model_AutoTSs = loadfn(model_AutoTSs_path)

for i, d in enumerate(model_AutoTSs):
    ts = d['output']['molecule']
    name = d['name']
    ts.to(name + '.xyz', 'xyz')
    
