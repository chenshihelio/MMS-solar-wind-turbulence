from functions import *
import matplotlib.pyplot as plt
import numpy as np
from pyspedas import mms_load_fgm, mms_load_fpi, mms_load_mec,mms_part_getspec,\
    mms_load_hpca
from pyspedas import tplot_names,get_data,store_data
from pyspedas.projects.mms import mms_qcotrans 
from pytplot import tplot,options
from pyspedas import minvar_matrix_make,tvector_rotate
from pyspedas import tlimit,options
from sys import exit



# pristine solar wind period for turbulence campaign
# trange = ['2019-2-24/14:40:00','2019-2-24/22:10:00']


trange_list = [['2019-2-24/14:40:00','2019-2-24/22:10:00'],
               ['2019-3-31/18:50:00','2019-3-31/21:50:00'],
               ['2019-4-7/19:00:00','2019-4-7/21:30:00'],
               ['2025-1-10/18:40:00','2025-1-11/02:00:00'],
               ['2025-1-28/8:30:00','2025-1-28/16:30:00']]


# for trange in trange_list:

# trange = trange_list[1]
trange = ['2019-3-31/20:25:00','2019-3-31/20:50:00']
trange_plot = ['2019-3-31/20:25:00','2019-3-31/20:40:00']

# read orbit data for MMS1
mec_data = mms_load_mec(trange=trange, probe=['1'],time_clip=True)
mms1_mec_r_gse = get_data('mms1_mec_r_gse',units=True) # dt=True

t_mec_mms1 = mms1_mec_r_gse.times
coord_mms1 = mms1_mec_r_gse.y
coord_mms1_RE = coord_mms1/Re




# read FGM data for MMS1
fgm_data = mms_load_fgm(trange=trange, probe=['1'],data_rate='srvy',time_clip=True)


# calculate MVA matrix 
mms1_fgm_b_gse_srvy_l2 = get_data('mms1_fgm_b_gse_srvy_l2')
store_data('mms1_fgm_bvec_gse_srvy_l2',data={'x':mms1_fgm_b_gse_srvy_l2.times,'y':mms1_fgm_b_gse_srvy_l2.y[:,0:3]})


minvar_matrix_make('mms1_fgm_bvec_gse_srvy_l2', tstart='2019-3-31/20:30:00', tstop='2019-3-31/20:36:00', 
                  newname= 'MVA_matrix')



# transform the B field
tvector_rotate('MVA_matrix','mms1_fgm_bvec_gse_srvy_l2','mms1_fgm_bvec_mva_srvy_l2')
tlimit(trange_plot)
options('mms1_fgm_bvec_mva_srvy_l2','legend_names',['$B_{l}$','$B_{n}$','$B_{m}$'])
# tplot(['mms1_fgm_bvec_gse_srvy_l2','mms1_fgm_bvec_mva_srvy_l2'])


# read fpi data
fpi_data = mms_load_fpi(trange=trange, probe=['1'],data_rate='fast',time_clip=True)
#print(fpi_data)
# tplot(['mms1_des_bulkv_gse_fast','mms1_des_bulkv_spintone_gse_fast'])

# read hpca data
hpca_data = mms_load_hpca(trange=trange, probe=['1'],data_rate='srvy',time_clip=True)
data_cotrans = mms_qcotrans('mms1_hpca_heplusplus_ion_bulk_velocity','mms1_hpca_heplusplus_ion_bulk_velocity_gse',
                            in_coord='dbcs',out_coord='gse',probe='1')



# get magnetic field data
mms1_fgm_b_gse = get_data('mms1_fgm_b_gse_srvy_l2',units=False)
mms1_fgm_b_gse_meta = get_data('mms1_fgm_b_gse_srvy_l2', metadata=True)

t_mms1_fgm_b_gse = mms1_fgm_b_gse.times
b_mms1_fgm_b_gse = mms1_fgm_b_gse.y[:,0:3]




# this command will generate a tplot variable 'mms1_des_dist_fast_pa'
PA_energy_range = [280,340]
data_pa = mms_part_getspec(trange=trange, probe=['1'], instrument='fpi', species='e', 
                data_rate='fast', output=['pa'], energy=PA_energy_range) 



# time step: very stable
dt_mms1_fgm_b_gse = np.diff(t_mms1_fgm_b_gse)
dt_ave = np.nanmean(dt_mms1_fgm_b_gse)
print('Average time step for MMS1_FGM_srvy = {:.3e} sec'.format(dt_ave))

# plt.plot(t_mms1_fgm_b_gse[0:-1],dt_mms1_fgm_b_gse)
# plt.xlabel('Time (s)')
# plt.ylabel('Time step (s)')
# plt.title('MMS1 FGM B GSE Time step')
# plt.show()


# select XX second for calculating PVI
tau = 60 # seconds
ndt = int(tau/dt_ave) # number of time steps to calculate the PVI

PVI_reult = calc_PVI(t_mms1_fgm_b_gse,b_mms1_fgm_b_gse,ndt)

t_PVI = PVI_reult['t']
PVI = PVI_reult['PVI']
delta_B_mean = PVI_reult['delta_arr_mean']


# store the PVI
store_data('mms1_fgm_b_gse_PVI',data={'x':t_PVI,'y':PVI})




# get date 
date0_str = trange[0].split('/')[0]


tlimit(trange_plot)
plot_objects = tplot(['mms1_fgm_b_gse_srvy_l2','mms1_fgm_bvec_mva_srvy_l2',
    'mms1_fgm_b_gse_PVI','mms1_des_dist_fast_pa',
    'mms1_dis_energyspectr_omni_fast',
    'mms1_dis_numberdensity_fast','mms1_des_numberdensity_fast',
    'mms1_dis_bulkv_gse_fast','mms1_hpca_heplusplus_ion_bulk_velocity_gse'],
    display=False,return_plot_objects=True) # ,save_png=date0_str + '.png'

fig = plot_objects[0]
subs = plot_objects[1]

subs[0].set_ylabel('$B_{GSE}$ [nT]')
subs[1].set_ylabel('$B_{MVA}$ [nT]')
subs[2].set_ylabel('PVI [60s]')
subs[3].set_ylabel('ePAD\n'+'{:d}-{:d}\n'.format(PA_energy_range[0],PA_energy_range[1]) + '[eV]')
subs[4].set_ylabel('DIS Energyspectrum OMNI')
subs[5].set_ylabel('DIS Density\n'+r'[cm$^{-3}$]')
subs[6].set_ylabel('DES Density\n'+r'[cm$^{-3}$]')
subs[7].set_ylabel('DIS Bulk Velocity GSE\n'+r'[km/s]')
subs[8].set_ylabel('He++ Bulk Velocity\n'+r'[km/s]')
# fig.savefig('./event_' + date0_str + '.png', dpi=300, bbox_inches='tight')
plt.show()


    # PAD
    # 'mms1_des_pitchangdist_lowen_fast','mms1_des_pitchangdist_miden_fast',
    #     'mms1_des_pitchangdist_highen_fast','mms1_des_pitchangdist_avg'