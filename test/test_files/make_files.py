
from lgdo import lh5, Table, Struct, Array
import numpy as np
import awkward as ak
import os
# make evt tier data
energy = ak.Array([[100],[200,300],[1000],[2000,3000]])
rawid = ak.Array([[1],[1,2],[3],[1,3]])
is_good = ak.Array([[True],[True,True],[True],[True,True]])

data = ak.Array({"geds":{"energy":energy,"rawid":rawid,"quality":{"is_good_channel":is_good}},"trigger":{"is_forced":[False,False,True,False]},
                "coincident":{"puls":[False,False,False,False]}})

data_table = Table(data)
lh5.write(data_table,"evt","evt_data.lh5",wo_mode="of")

det1 = Table(ak.Array({"energy":[1000,1020,900]}))
det2 = Table(ak.Array({"energy":[20]}))
det3 = Table(ak.Array({"energy":[1]}))

lh5.write(Struct({"det1":det1,"det2":det2,"det3":det3}),"hit",f"mc.lh5",wo_mode="of")

# make mc file
for phi in [0,1]:
    for z in [-1,1]:
        lh5.write(Struct({"det1":det1,"det2":det2,"det3":det3}),"hit",f"z_{z}_phi_{phi}/mc.lh5",wo_mode="of")
