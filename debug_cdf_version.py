from pyspedas import mms_load_fpi


trange = ['2023-01-18','2023-01-19']
data = mms_load_fpi(trange=trange,probe=['1'], data_rate='fast', datatype='dis-moms')