
from lgdo import lh5, Table, Struct, Array
import numpy as np
import awkward as ak

# make evt tier data
energy = ak.Array([[100],[200,300],[1000],[2000,3000]])
rawid = ak.Array([[1],[1,2],[3],[1,3]])
data = ak.Array({"geds":{"energy":energy,"rawid":rawid},"trigger":{"is_forced":[False,False,True,False]},
                "coincident":{"pulser":[False,False,False,False]}})

data_table = Table(data)
lh5.write(data_table,"evt","evt_data.lh5",wo_mode="of")
