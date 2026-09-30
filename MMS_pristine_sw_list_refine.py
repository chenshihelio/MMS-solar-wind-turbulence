import os
from functions import *
import matplotlib.pyplot as plt
import numpy as np
from pyspedas import mms_load_fgm, mms_load_fpi, mms_load_mec,mms_part_getspec,\
    mms_load_hpca
from pyspedas import tplot_names,get_data,store_data, avg_data
from pyspedas.projects.mms import mms_qcotrans 
from pytplot import tplot,options
from pyspedas import minvar_matrix_make,tvector_rotate
from pyspedas import tlimit,options,time_string
from sys import exit
from pandas import read_csv
import pandas as pd
import datetime 


'''
# read the pristine solar wind periods provided by Yi Qi
data_pristine = read_csv('./input/MMS SW intervals - Sheet1.csv')

start_time_list = data_pristine['start'].to_list()
end_time_list = data_pristine['end'].to_list()


ninterval = len(start_time_list)

for iinter in range(ninterval):
    trange = [start_time_list[iinter], end_time_list[iinter]]


    # Read data
    # read FGM data for MMS1
    fgm_data = mms_load_fgm(trange=trange, probe=['1'],data_rate='srvy',time_clip=True)


    # get magnetic field data
    mms1_fgm_b_gse = get_data('mms1_fgm_b_gse_srvy_l2',units=False)

    t_mms1_fgm_b_gse = mms1_fgm_b_gse.times
    b_mms1_fgm_b_gse = mms1_fgm_b_gse.y[:,0:3]


    babs_mms1_fgm_b_gse = mms1_fgm_b_gse.y[:,3] # absolute value of the magnetic field






    plot_objects = tplot(['mms1_fgm_b_gse_srvy_l2'],display=False,return_plot_objects=True)

    fig, subs = plot_objects[0], plot_objects[1]

    fig.suptitle('Interval {}: {} - {}'.format(iinter+1, time_string(t_mms1_fgm_b_gse[0],fmt='%m-%d %H:%M'), 
                            time_string(t_mms1_fgm_b_gse[-1],fmt='%m-%d %H:%M')))

    plt.show()



    continue


    # time step: very stable
    dt_mms1_fgm_b_gse = np.diff(t_mms1_fgm_b_gse)
    dt_ave = np.nanmean(dt_mms1_fgm_b_gse)


    
    # total time duration of the interval
    T_interval = t_mms1_fgm_b_gse[-1] - t_mms1_fgm_b_gse[0]

    
    # find median of all babs values
    babs_median = np.nanmedian(babs_mms1_fgm_b_gse)

    print('Median(|B|) = ', babs_median, 'nT')


    # average the magnetic field data with a window of 10 minutes (600 seconds)
    avg_data('mms1_fgm_b_gse_srvy_l2', res=600)


    # get the averaged |B|
    t_avg = get_data('mms1_fgm_b_gse_srvy_l2-avg').times
    babs_avg = get_data('mms1_fgm_b_gse_srvy_l2-avg').y[:,3]


    # in case that the begining is still inside the magnetosheath,
    # we scan from the beginning, if |B| is very large,
    # we determine the boundary at which |B| drops to a value close to the median value of |B| in the pristine solar wind 
    # (e.g., 2 times of the median value).
    
    babs_sheath_ratio_threshold = 1.5
    it_boundary_begin = 0
    for it in range(len(t_avg)//4):
        if babs_avg[it] < babs_sheath_ratio_threshold*babs_median:
            # t_boundary = t_avg[it]
            # print('Boundary time = ', time_string(t_boundary))
            it_boundary_begin = it
            break

    t_boundary_begin = t_avg[it_boundary_begin]

    # scan the end of the interval
    it_boundary_begin = len(t_avg)-1
    for it in range(len(t_avg)-1, len(t_avg)*3//4, -1):
        if babs_avg[it] < babs_sheath_ratio_threshold*babs_median:
            # t_boundary = t_avg[it]
            # print('Boundary time = ', time_string(t_boundary))
            it_boundary_end = it
            break

    t_boundary_end = t_avg[it_boundary_end]



    # plot_objects = tplot(['mms1_fgm_b_gse_srvy_l2','mms1_fgm_b_gse_srvy_l2-avg'],display=False,return_plot_objects=True)
    # fig = plot_objects[0]
    # subs = plot_objects[1]


    xticks = np.linspace(t_mms1_fgm_b_gse[0], t_mms1_fgm_b_gse[-1], 5)
    xticklabels = [time_string(t,fmt='%H:%M\n%m-%d') for t in xticks]

    xticks_minor = np.linspace(t_mms1_fgm_b_gse[0], t_mms1_fgm_b_gse[-1], 21)

    fig, subs = plt.subplots(2,1, figsize=[10,5], sharex=True)

    sub = subs[0]
    sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,0], label=r'$B_x$')
    sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,1], label=r'$B_y$')
    sub.plot(t_mms1_fgm_b_gse, b_mms1_fgm_b_gse[:,2], label=r'$B_z$')
    sub.plot(t_mms1_fgm_b_gse, babs_mms1_fgm_b_gse, label=r'$|B|$',color='k')
    sub.set_ylabel('B [nT]')
    sub.legend(loc='upper right')
    sub.axvline(t_boundary_begin, color='r', linestyle='--', label='Boundary (begin)', lw=2)
    sub.axvline(t_boundary_end, color='g', linestyle='--', label='Boundary (end)', lw=2)
    sub.set_xticks(xticks)
    sub.set_xticklabels(xticklabels)
    sub.set_xticks(xticks_minor, minor=True)
    sub.set_xlim(t_mms1_fgm_b_gse[0], t_mms1_fgm_b_gse[-1])


    sub = subs[1]
    sub.plot(t_avg, babs_avg, label=r'|B| (averaged)', color='k')
    sub.set_ylabel('Averaged |B| [nT]')
    sub.legend(loc='upper right')
    sub.axvline(t_boundary_begin, color='r', linestyle='--', label='Boundary (begin)', lw=2)
    sub.axvline(t_boundary_end, color='g', linestyle='--', label='Boundary (end)', lw=2)
    sub.set_xticks(xticks)
    sub.set_xticklabels(xticklabels)
    sub.set_xticks(xticks_minor, minor=True)
    sub.set_xlim(t_mms1_fgm_b_gse[0], t_mms1_fgm_b_gse[-1])

    fig.suptitle('Interval {}: {} - {}'.format(iinter+1, time_string(t_mms1_fgm_b_gse[0],fmt='%m-%d %H:%M'), 
                            time_string(t_mms1_fgm_b_gse[-1],fmt='%m-%d %H:%M')))

    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig('./pristine_sw_interval_{}.png'.format(iinter+1))


    data_pristine['start'][iinter] = time_string(t_boundary_begin,fmt='%Y-%m-%d/%H:%M:%S') 
    data_pristine['end'][iinter] = time_string(t_boundary_end,fmt='%Y-%m-%d/%H:%M:%S')

    

# data_pristine.to_csv('./MMS_pristine_sw_intervals_refined.csv', index=False)
'''



# read the eye-selected intervals 
data = read_csv('./input/MMS SW intervals - Sheet1_eye_selected.csv')

print(data.columns)


print('Total number of intervals before selection = ', len(data))

# select rows with "adjusted" == 1
data_adjusted = data[data['adjusted'] == 1]

print('Total number of intervals after selection = ', len(data_adjusted))


# some rows have second intervals, we split them into separate rows
for i in range(len(data_adjusted)):
    if pd.isna(data_adjusted['interval 2 start'].iloc[i]):
        continue
    else:
        new_row = data_adjusted.iloc[i].copy()
        new_row['start'] = data_adjusted['interval 2 start'].iloc[i]
        new_row['end'] = data_adjusted['interval 2 end'].iloc[i]


        data_adjusted = pd.concat([data_adjusted, pd.DataFrame([new_row])], ignore_index=True)



print('Total number of intervals after splitting = ', len(data_adjusted))



if not os.path.exists('./figure/selected_intervals'):
    os.makedirs('./figure/selected_intervals')

# sort according to the start time
data_adjusted['start_datetime'] = pd.to_datetime(data_adjusted['start'])
data_adjusted = data_adjusted.sort_values(by='start_datetime').reset_index(drop=True)


# save the adjusted intervals to a new csv file, 
# delete the "interval 2 start" and "interval 2 end" columns
data_adjusted = data_adjusted.drop(columns=['interval 2 start', 'interval 2 end', 'start_datetime',
                                            'adjusted'])
data_adjusted.to_csv('./input/MMS_pristine_sw_intervals_selected.csv', index=False)

exit()

start_time_list = data_adjusted['start'].to_list()
end_time_list = data_adjusted['end'].to_list()

# plot the intervals
ninterval = len(start_time_list)

for iinter in range(ninterval):
    print('Processing interval {}/{}: {} - {}'.format(iinter+1, ninterval, 
                            start_time_list[iinter], end_time_list[iinter]))
    trange = [start_time_list[iinter], end_time_list[iinter]]


    # Read data
    # read FGM data for MMS1
    fgm_data = mms_load_fgm(trange=trange, probe=['1'],data_rate='srvy',time_clip=True)


    # get magnetic field data
    mms1_fgm_b_gse = get_data('mms1_fgm_b_gse_srvy_l2',units=False)

    t_mms1_fgm_b_gse = mms1_fgm_b_gse.times
    b_mms1_fgm_b_gse = mms1_fgm_b_gse.y[:,0:3]


    babs_mms1_fgm_b_gse = mms1_fgm_b_gse.y[:,3] # absolute value of the magnetic field






    plot_objects = tplot(['mms1_fgm_b_gse_srvy_l2'],display=False,return_plot_objects=True)

    fig, subs = plot_objects[0], plot_objects[1]

    fig.suptitle('Interval {}: {} - {}'.format(iinter+1, time_string(t_mms1_fgm_b_gse[0],fmt='%m-%d %H:%M'), 
                            time_string(t_mms1_fgm_b_gse[-1],fmt='%m-%d %H:%M')))

    fig.savefig('./figure/selected_intervals/interval_{}.png'.format(iinter+1),dpi=200)
    plt.close(fig)