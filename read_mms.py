from pyspedas import mms_load_fgm, mms_load_fpi, \
    mms_load_hpca, mms_load_feeps, mms_load_edp, \
    mms_load_mec,mms_part_getspec
from pyspedas import tplot_names,get_data,store_data
from pyspedas.projects.mms import mms_qcotrans 
from pytplot import tplot,options
import numpy as np
import matplotlib.pyplot as plt
from functions import *




t0 = '2019-3-31/18:50:00'
t1 = '2019-3-31/21:50:00'


mms_load_mec(trange=[t0,t1], probe='1',time_clip=True)

# hpca data
data = mms_load_hpca(trange=[t0,t1], probe=['1'],data_rate='srvy',time_clip=True)
print(data)
mms1_hpca_heplusplus_ion_bulk_velocity_meta = get_data('mms1_hpca_heplusplus_ion_bulk_velocity', metadata=True) # DBCS coordinates

data_cotrans = mms_qcotrans('mms1_hpca_heplusplus_ion_bulk_velocity','mms1_hpca_heplusplus_ion_bulk_velocity_gse',in_coord='dbcs',out_coord='gse',probe='1')



fpi_data = mms_load_fpi(trange=[t0,t1], probe=['1'],data_rate='fast',time_clip=True)
# print(fpi_data)

# # the partial data looks strange......
# mms1_dis_bulkv_part_dbcs_fast_meta = get_data('mms1_dis_bulkv_part_dbcs_fast', metadata=True)
# print(mms1_dis_bulkv_part_dbcs_fast_meta)

# mms1_dis_bulkv_part_gse_fast = get_data('mms1_dis_bulkv_part_gse_fast')
# mms1_dis_bulkv_part_gse_fast.y[mms1_dis_bulkv_part_gse_fast.y < 1] = np.nan
# store_data('mms1_dis_bulkv_part_gse_fast', data={'x':mms1_dis_bulkv_part_gse_fast.times, 'y':mms1_dis_bulkv_part_gse_fast.y})
# tplot(['mms1_dis_bulkv_part_dbcs_fast','mms1_dis_bulkv_part_gse_fast','mms1_dis_bulkv_gse_fast'])



options('mms1_hpca_heplusplus_ion_bulk_velocity_gse','legend_names',['$V_{x,GSE}$','$V_{y,GSE}$','$V_{z,GSE}$'])
options('mms1_hpca_heplusplus_ion_bulk_velocity_gse','ytitle','He++ V \ [km/s]')
tplot(['mms1_dis_bulkv_gse_fast','mms1_hpca_heplusplus_ion_bulk_velocity_gse','mms1_hpca_heplusplus_ion_bulk_velocity'])


# # test coordinate transform
# mms_load_mec(trange=[t0,t1], probe='1',time_clip=True) # necessary to transform coordinates
# data_cotrans = mms_qcotrans('mms1_dis_bulkv_dbcs_fast','mms1_dis_bulkv_dbcs_to_gse_fast',in_coord='dbcs',out_coord='gse',probe='1')
# print(data_cotrans)
# tplot(['mms1_dis_bulkv_dbcs_fast','mms1_dis_bulkv_dbcs_to_gse_fast','mms1_dis_bulkv_gse_fast'])

exit()



# read pitch angle data
# this command will generate a tplot variable 'mms1_des_dist_fast_pa'
data = mms_part_getspec(trange=[t0,t1], probe=['1'], instrument='fpi', species='e', 
                 data_rate='fast', output=['pa'], energy=[350,450]) 

# print(data)
# print(tplot_names())

mms1_des_dist_fast_pa_meta = get_data('mms1_des_dist_fast_pa', metadata=True)
print(mms1_des_dist_fast_pa_meta)
tplot(['mms1_des_dist_fast_pa'])

exit()

# read the location of MMS
mec_data = mms_load_mec(trange=[t0,t1], probe=['1','2','3','4'],time_clip=True)
# tplot(['mms1_mec_r_gsm', 'mms1_mec_v_gsm'])
# print(tplot_names())

mms1_mec_r_gse = get_data('mms1_mec_r_gse',units=True) # dt=True
# mms1_mec_r_gse_meta = get_data('mms1_mec_r_gse', metadata=True)
mms2_mec_r_gse = get_data('mms2_mec_r_gse',units=True) # dt=True
mms3_mec_r_gse = get_data('mms3_mec_r_gse',units=True) # dt=True
mms4_mec_r_gse = get_data('mms4_mec_r_gse',units=True) # dt=True

t_mec_mms1 = mms1_mec_r_gse.times
coord_mms1 = mms1_mec_r_gse.y
coord_mms1_RE = coord_mms1/Re

t_mec_mms2 = mms2_mec_r_gse.times
coord_mms2 = mms2_mec_r_gse.y
coord_mms2_RE = coord_mms2/Re

t_mec_mms3 = mms3_mec_r_gse.times
coord_mms3 = mms3_mec_r_gse.y
coord_mms3_RE = coord_mms3/Re

t_mec_mms4 = mms4_mec_r_gse.times
coord_mms4 = mms4_mec_r_gse.y
coord_mms4_RE = coord_mms4/Re


# tplot(['mms1_mec_r_gse','mms2_mec_r_gse','mms3_mec_r_gse','mms4_mec_r_gse'])


# fig = plt.figure()
# sub = fig.add_subplot(111)
# sub.plot(coord_mms1_RE[:,0],coord_mms1_RE[:,1],label='MMS1')
# sub.plot(coord_mms2_RE[:,0],coord_mms2_RE[:,1],label='MMS2')
# sub.plot(coord_mms3_RE[:,0],coord_mms3_RE[:,1],label='MMS3')
# sub.plot(coord_mms4_RE[:,0],coord_mms4_RE[:,1],label='MMS4')
# sub.legend()
# sub.set_xlabel('X (Re) / GSE')
# sub.set_ylabel('Y (Re) / GSE')
# fig.suptitle('MMS1-4 Trajectory: '+ t0 + ' -- ' + t1)
# plt.show()

# # calculate radius in Re and store in tplot
# radius = np.sqrt(coord_mms1_RE[:,0]**2 + coord_mms1_RE[:,1]**2 + coord_mms1_RE[:,2]**2)
# store_data("radius",data={'x':t_coord,'y':radius})


# read fgm data
fgm_data = mms_load_fgm(trange=[t0,t1], probe=['1','2','3','4'],data_rate='srvy',time_clip=True)

#mms1_fgm_b_gse = get_data('mms1_fgm_b_gse_srvy_l2',units=True)
# mms1_fgm_b_gse_meta = get_data('mms1_fgm_b_gse_srvy_l2', metadata=True)
# print(mms1_fgm_b_gse_meta)


# read fpi data
fpi_data = mms_load_fpi(trange=[t0,t1], probe=['1'],data_rate='fast',time_clip=True)


# print(fpi_data)
for i in range(len(fpi_data)):
    if 'des' in fpi_data[i]:
        print(fpi_data[i])


# keyword = '_dis_'
# fpi_data_dis = [s for s in fpi_data if keyword in s]
# print(fpi_data_dis)



mms1_des_pitchangdist_miden_fast_meta = get_data('mms1_des_pitchangdist_miden_fast', metadata=True)
print(mms1_des_pitchangdist_miden_fast_meta)


mms1_des_energy_fast = get_data('mms1_des_energy_fast',units=True)
mms1_des_energy_fast_meta = get_data('mms1_des_energy_fast',metadata=True)
print(mms1_des_energy_fast_meta)

print(mms1_des_energy_fast.y.shape)
print(mms1_des_energy_fast.y[0])

exit()



tplot(['mms1_fgm_b_gse_srvy_l2','mms1_dis_energyspectr_omni_fast',
    'mms1_dis_numberdensity_fast','mms1_des_numberdensity_fast',
    'mms1_des_pitchangdist_lowen_fast','mms1_des_pitchangdist_miden_fast',
    'mms1_des_pitchangdist_highen_fast','mms1_des_pitchangdist_avg'])
