from functions import *
import matplotlib.pyplot as plt
import numpy as np
from pyspedas import mms_load_fgm, mms_load_fpi, mms_load_mec,mms_part_getspec,\
    mms_load_hpca
from pyspedas import tplot_names,get_data,store_data
from pyspedas.projects.mms import mms_qcotrans 
from pytplot import tplot,options
from pyspedas import minvar_matrix_make,tvector_rotate
from pyspedas import tlimit,options,time_string
from sys import exit
from pandas import read_csv
import datetime 

tau = 60 # seconds


# read the pristine solar wind periods provided by Yi Qi
data_pristine = read_csv('./input/MMS SW intervals - Sheet1.csv')

start_time_list = data_pristine['start'].to_list()
end_time_list = data_pristine['end'].to_list()


ninterval = len(start_time_list)

for iinter in range(ninterval):

# iinter = 2 # as an example

    trange = [start_time_list[iinter], end_time_list[iinter]]


    # Read data
    # read FGM data for MMS1
    fgm_data = mms_load_fgm(trange=trange, probe=['1'],data_rate='srvy',time_clip=True)


    # get magnetic field data
    mms1_fgm_b_gse = get_data('mms1_fgm_b_gse_srvy_l2',units=False)

    t_mms1_fgm_b_gse = mms1_fgm_b_gse.times
    b_mms1_fgm_b_gse = mms1_fgm_b_gse.y[:,0:3]

    # time step: very stable
    dt_mms1_fgm_b_gse = np.diff(t_mms1_fgm_b_gse)
    dt_ave = np.nanmean(dt_mms1_fgm_b_gse)


    # plot magnetic field only
    # fig, sub = plt.subplots(1, 1, figsize=(10, 3))
    # sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,0], label=r'$B_x$')
    # sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,1], label=r'$B_y$')
    # sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,2], label=r'$B_z$')
    # sub.plot(t_mms1_fgm_b_gse, mms1_fgm_b_gse.y[:,3], label=r'$|B|$',color='k')
    # sub.set_ylabel(r'$B_{GSE}$ [nT]')
    # sub.legend(ncols=4,loc='lower right')
    # sub.set_xlabel('Time', fontsize=13)
    # plt.show()

    tplot(['mms1_fgm_b_gse_srvy_l2'])


    # # select XX second for calculating PVI
    # ndt = int(tau/dt_ave) # number of time steps to calculate the PVI


    # PVI_reult = calc_PVI(t_mms1_fgm_b_gse,b_mms1_fgm_b_gse,ndt)

    # t_PVI = PVI_reult['t']
    # PVI = PVI_reult['PVI']
    # delta_B_mean = PVI_reult['delta_arr_mean']


    # # store the PVI
    # store_data('mms1_fgm_b_gse_PVI',data={'x':t_PVI,'y':PVI})

    # date0_str = trange[0].split('/')[0]





    # # identify discontinuities
    # PVI_threshold = 3.0
    # PVI_threshold_exit = 1.0

    # nrec = len(PVI)
    # n_discon = 0
    # stat = False 
    # ind_arr = np.zeros([nrec,2]) + np.nan

    # print('Begin identify discontinuities......')
    # for i in range(nrec):
    #     progress_bar(i+1, nrec, bar_length=100)
    #     if PVI[i] >= PVI_threshold_exit and not stat:
    #         stat = not stat
    #         ind_arr[n_discon, 0] = i
    #     if PVI[i] < PVI_threshold_exit and stat:
    #         stat = not stat
    #         ind_arr[n_discon, 1] = i

    #         # find the peak
    #         peak_PVI = np.nanmax(PVI[int(ind_arr[n_discon, 0]): (int(ind_arr[n_discon, 1]) + 1)])
    #         # if the peak is low, discard
    #         if peak_PVI < PVI_threshold:
    #             ind_arr[n_discon, 0] = np.nan
    #             ind_arr[n_discon, 1] = np.nan
    #             continue

    #         n_discon += 1

    # print('\n')

    # ind_arr_discon = np.copy(ind_arr[0:n_discon,:]).astype(int)

    # # discontinuity start times [n_discon]
    # t_discon_start = t_PVI[ind_arr_discon[:,0]]
    # # discontinuity end times [n_discon]
    # t_discon_end = t_PVI[ind_arr_discon[:,1]]


    # fig, subs = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

    # # calculate the xticks and xticklabels
    # # t_xticks = sub.get_xticks()
    # dt_tick_0 = datetime.datetime(2023,1,6,18,0,0)
    # dt_tick_1 = datetime.datetime(2023,1,9,6,0,0)

    # dt_ticks = np.arange(dt_tick_0, dt_tick_1 + datetime.timedelta(seconds=10), 
    #                     datetime.timedelta(hours=6)).astype(datetime.datetime)


    # dt_ticks_double = np.zeros(len(dt_ticks))
    # for i in range(len(dt_ticks)):
    #     dt_ticks_double[i] = time_double(dt_ticks[i].strftime('%Y-%m-%d/%H:%M:00'))

    # xticklabels = []
    # for tick in dt_ticks_double:
    #     xticklabels.append(time_string(tick, fmt='%H:%M\n%m/%d'))


    # sub = subs[0]
    # sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,0], label=r'$B_x$')
    # sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,1], label=r'$B_y$')
    # sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,2], label=r'$B_z$')
    # sub.plot(t_mms1_fgm_b_gse, mms1_fgm_b_gse.y[:,3], label=r'$|B|$',color='k')
    # sub.set_ylabel(r'$B_{GSE}$ [nT]')
    # sub.legend(ncols=4,loc='lower right')
    # sub.set_xticks(dt_ticks_double)
    # sub.set_xticklabels(xticklabels)

    # sub = subs[1]
    # sub.plot(t_PVI, PVI, label='PVI (tau={:d}s)'.format(tau))
    # sub.set_ylabel(r'PVI($\tau=$60s)')

    # for i in range(n_discon):
    #     t_PVI_discon = t_PVI[ind_arr_discon[i,0]: ind_arr_discon[i,1]+1]
    #     PVI_discon = PVI[ind_arr_discon[i,0]: ind_arr_discon[i,1]+1]
    #     ind_max_PVI_discon = np.nanargmax(PVI_discon)

    #     sub.scatter(t_PVI_discon[ind_max_PVI_discon], PVI_discon[ind_max_PVI_discon], color='r', s=10)


    # sub.set_xticks(dt_ticks_double)
    # sub.set_xticklabels(xticklabels)

    # sub.set_xlabel('Time', fontsize=13)

    # fig.suptitle(f'{n_discon} discontinuities identified with PVI>{PVI_threshold:.1f}', fontsize=13)

    # fig.tight_layout(rect=[0,0,1,0.99])


    # fig.savefig('./discontinuities_search_' + date0_str + '.png',dpi=200)
    # plt.show()