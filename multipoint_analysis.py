from pyspedas import mms_load_fgm, mms_load_fpi, \
    mms_load_hpca, mms_load_feeps, mms_load_edp, \
    mms_load_mec,mms_part_getspec
from pyspedas import tplot_names,get_data,store_data, \
    tsmooth,avg_data,time_datetime,time_double,subtract,\
    tinterpol
from pyspedas.projects.mms import mms_qcotrans 
from pytplot import tplot,options
from pyspedas import subtract_average
import numpy as np
import matplotlib.pyplot as plt
from functions import *
from scipy.constants import mu_0, m_p
import os
from itertools import combinations


mms_id_list = ['1','2','3','4']  # list of MMS probes

Nsat = 4  # number of satellites


# trange = ['2019-2-24/12:10:00','2019-2-24/22:10:00'] # string-configuration
trange = ['2019-3-31/18:50:00','2019-3-31/21:50:00']  # tetrahedron configuration
# trange = ['2019-3-31/18:50:00','2019-3-31/20:50:00']  # tetrahedron configuration -- shorter for comparison
# trange = ['2017-1-18/00:45:53','2017-1-18/00:49:43'] # (Bandyopadhyay2018 ApJ 866 106) # sheath

# trange = ['2017-11-24/01:10:03','2017-11-24/02:10:03'] # (Bandyopadhyay2018 ApJ 866 81) # sw
# trange = ['2017-11-23/22:40:00','2017-11-24/02:30:00'] # extended (Bandyopadhyay2018 ApJ 866 81) # sw

trange0_str = trange[0].replace(':', '-').replace('/', '-')
trange1_str = trange[1].replace(':', '-').replace('/', '-')



# # the analyzed data--------------------
# file_output = './output/Yaglom_law_' + trange0_str + '.npz' # '_shorter_interval.npz'

# if os.path.exists(file_output):
#     data = np.load(file_output)

#     trange = data['trange']
#     tau_arr = data['tau_arr']
#     ndens_avg = data['ndens_avg']
#     bulkv_avg = data['bulkv_avg']
#     B_avg = data['B_avg']
#     coord_diff_mms_avg = data['coord_diff_mms_avg']
#     dX_center = data['dX_center']
#     Yp_arr = data['Yp_arr']
#     Ym_arr = data['Ym_arr']


#     dX_abs = np.sqrt(np.sum(dX_center**2, axis=1))  # absolute distance-increment of the center of the tetrahedron


#     ntau = len(tau_arr)  # number of time increments

#     dX_direction = np.zeros((ntau, 3))  # direction of the center of the tetrahedron
#     for itau in range(ntau):
#         dX_direction[itau, :] = dX_center[itau, :] / np.linalg.norm(dX_center[itau, :])  # direction of the center of the tetrahedron

#     dX_direction = dX_direction[0,:]



#     tau_start_analysis = 60 #s
#     tau_end_analysis = 3600 #s

#     ind_tau_start = np.where(tau_arr >= tau_start_analysis)[0][0]  # index of the start time increment
#     ind_tau_end = np.where(tau_arr <= tau_end_analysis)[0][-1]  # index of the end time increment



#     # For each dX, we have 4*3+1=13 legs
#     # from the 13 legs, we can form C_{13}^4 = 715 tetrahedrons
#     nleg = 13
#     ncomb = 715  # number of combinations of tetrahedrons
    
#     # the leg corresponding to 11 or 22 or 33 or 44 is marked as 00. But we will average (11-44).

#     # fisrt, list all possible 13 legs
#     leg_list = np.zeros((nleg,2),dtype=int)
#     leg_list[0,0] = 0
#     leg_list[0,1] = 0
#     ind = 1
#     for i in range(1,5):
#         for j in range(1,5):
#             if j == i:
#                 continue

#             leg_list[ind,0] = i 
#             leg_list[ind,1] = j
#             ind += 1

#     # print(leg_list)

#     # second, find all combinations of 4 legs from the 13 legs
#     tetra_info = np.zeros((ncomb,4),dtype=int) # (tetrahedron index, leg index)
#     comb = combinations(range(nleg), 4)  # combinations of 4 legs from the 13 legs
#     ind_tetra = 0
#     for c in comb:
#         for ic in range(4):
#             tetra_info[ind_tetra, ic] = c[ic]

#         ind_tetra += 1


#     # print(tetra_info)

#     # exit()

#     # third, calculate the relative spatial difference, w.r.t. dX_center, of these legs
#     dx_list = np.zeros((nleg,3))
#     for ileg in range(1,nleg):
#         leg = leg_list[ileg,:]
#         if leg[0]==leg[1]:
#             print('error: leg[0] == leg[1], should not happen')

#         dx_list[ileg,:] = (coord_diff_mms_avg[leg[1]-1,:] - coord_diff_mms_avg[leg[0]-1,:]) 

#     # for i in range(nleg):
#     #     print(leg_list[i,:], dx_list[i,:])

#     # exit()


#     dX_subarr = dX_center[ind_tau_start:ind_tau_end+1]
#     tau_subarr = tau_arr[ind_tau_start:ind_tau_end+1]
#     dX_abs_subarr = dX_abs[ind_tau_start:ind_tau_end+1] 


#     divYp_subarr = np.zeros(len(tau_subarr)) # mean of all tetrahedrons
#     divYm_subarr = np.zeros(len(tau_subarr)) # mean of all tetrahedrons
#     divYp_std_subarr = np.zeros(len(tau_subarr)) # std of all tetrahedrons
#     divYm_std_subarr = np.zeros(len(tau_subarr)) # std of the mean of all tetrahedrons

#     print('Calculating the divergence of Y+ and Y- for each tau...')
#     for itau in range(ind_tau_start, ind_tau_end+1):
#         progress_bar(itau-ind_tau_start, ind_tau_end-ind_tau_start+1, 100)
#         tau = tau_arr[itau]
#         dX0 = dX_center[itau]

#         # print('Calculating for tau = {:.1f} s, dX0 = {:.1f} km'.format(tau, np.linalg.norm(dX0)))

#         # first, calculate Y+ and Y- for the 13 legs
#         Yp = np.zeros((nleg,3))
#         Ym = np.zeros((nleg,3))
#         for ileg in range(nleg):
#             leg = leg_list[ileg,:]

#             if leg[0]==0 and leg[1]==0: # leg 00, average 
#                 for isat in range(Nsat):
#                     Yp[ileg,:] = Yp[ileg,:] + Yp_arr[itau,isat,isat,:] 
#                     Ym[ileg,:] = Ym[ileg,:] + Ym_arr[itau,isat,isat,:]
#                 Yp[ileg,:] = Yp[ileg,:] / Nsat
#                 Ym[ileg,:] = Ym[ileg,:] / Nsat
#             else:
#                 if leg[0] == leg[1]:
#                     print('error: leg[0] == leg[1], should not happen')

#                 Yp[ileg,:] = Yp_arr[itau,leg[0]-1,leg[1]-1,:]
#                 Ym[ileg,:] = Ym_arr[itau,leg[0]-1,leg[1]-1,:]

#         # print(Yp,Ym)
#         # exit()


#         # second, get the corresponding satellite separation -- not needed
#         dX = np.zeros((nleg,3))
#         for ileg in range(nleg):
#             # dX = dX0 + (r_isat1 - r_isat0)
#             dX[ileg,:] = dX0[:] + dx_list[ileg,:]



#         # functions used to study why a tetrahedron is invalid
#         def if_repeat_pairs(legs):
#             n = len(legs)

#             for i in range(n-1):
#                 leg0 = legs[i]

#                 for j in range(i+1,n):
#                     leg1 = legs[j]
#                     if leg1[0] == leg0[1] and leg1[1] == leg0[0]:
#                         return True 
#             return False
        
#         def if_two_repeat_pairs(legs):
#             n = len(legs)

#             n_pair = 0
#             for i in range(n-1):
#                 leg0 = legs[i]

#                 for j in range(i+1,n):
#                     leg1 = legs[j]
#                     if leg1[0] == leg0[1] and leg1[1] == leg0[0]:
#                         n_pair += 1
#             if n_pair == 2:
#                 return True
#             else:
#                 return False

#         def if_has_00(legs):
#             for leg in legs:
#                 if leg[0] == 0 and leg[1] == 0:
#                     return True
#             return False

#         def if_same_tetra(tetra0_legs, tetra1_legs, if_exact=False):
#             num_same_legs = 0
#             for i in range(4):
#                 leg1 = tetra1_legs[i]

#                 for j in range(4):
#                     leg0 = tetra0_legs[j]

#                     if if_exact:
#                         if leg1[0] == leg0[0] and leg1[1] == leg0[1]:
#                             num_same_legs += 1
#                     else:
#                         if (leg1[0] == leg0[1] and leg1[1] == leg0[0]) \
#                             or (leg1[0] == leg0[0] and leg1[1] == leg0[1]):
#                             num_same_legs += 1

#             if num_same_legs == 4:
#                 return True
#             else:
#                 return False



#         # third, take each tetrahedron
#         divYp_tetra = []
#         divYm_tetra = []

#         n_valid_tetra = 0
#         n_invalid_tetra = 0
#         n_invalid_with_00 = 0
#         n_invalid_with_two_pairs = 0
#         n_invalid_with_00_and_one_pair = 0
#         n_invalid_with_00_and_no_pair = 0 
#         n_invalid_tetra_without_repeats = 0
#         n_invalid_tetra_without_repeats_and_00 = 0
#         n_invalid_tetra_one_pair_no_00 = 0
#         for itetra in range(ncomb):
#             legs_tetra = tetra_info[itetra,:]  # indices of the 4 legs
#             legs_this_tetra = leg_list[legs_tetra,:] # get the legs of this tetrahedron

#             dx_vec = dx_list[legs_tetra,:]  # get the dx vector of the legs of this tetrahedron, 4 legs * 3 components


#             # determine whether the tetrahedron is valid
#             v1 = dx_vec[1,:] - dx_vec[0,:]
#             v2 = dx_vec[2,:] - dx_vec[0,:]
#             v3 = dx_vec[3,:] - dx_vec[0,:]

#             vol_norm = np.linalg.norm(v1) * np.linalg.norm(v2) * np.linalg.norm(v3)
#             vol_tetra = np.abs(np.dot(v1, np.cross(v2, v3))) # volume of the tetrahedron

            

#             if vol_tetra/vol_norm < 1e-8:
#                 n_invalid_tetra += 1
#                 # print('Tetrahedron #{:03d} is invalid. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

#                 # see play_with_tetrahedron.py for the details
#                 if if_has_00(legs_this_tetra):
#                     n_invalid_with_00 += 1
#                     # print('Tetrahedron #{:03d} is invalid with 00. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

#                 if if_two_repeat_pairs(legs_this_tetra):
#                     n_invalid_with_two_pairs += 1
#                     # print('Tetrahedron #{:03d} is invalid with two repeat pairs. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

#                 if if_has_00(legs_this_tetra) and if_repeat_pairs(legs_this_tetra):
#                     n_invalid_with_00_and_one_pair += 1
#                     # print('Tetrahedron #{:03d} is invalid with 00 and one repeat pair. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

#                 if if_has_00(legs_this_tetra) and not if_repeat_pairs(legs_this_tetra):
#                     n_invalid_with_00_and_no_pair += 1
#                     # print('Tetrahedron #{:03d} is invalid with 00 and no repeat pairs. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

#                 if not if_repeat_pairs(legs_this_tetra):
#                     n_invalid_tetra_without_repeats += 1
#                     # print('Tetrahedron #{:03d} is invalid without repeat pairs. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

#                 if if_repeat_pairs(legs_this_tetra) and not if_has_00(legs_this_tetra) and not if_two_repeat_pairs(legs_this_tetra) :
#                     n_invalid_tetra_one_pair_no_00 += 1
#                     # print('Tetrahedron #{:03d} is invalid with one pair and without 00. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

#                 if not if_repeat_pairs(legs_this_tetra) and not if_has_00(legs_this_tetra):
#                     n_invalid_tetra_without_repeats_and_00 += 1
#                     # print('Tetrahedron #{:03d} is invalid without repeat pairs and without 00. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))
#             else:
#                 # calculate divergence div_l(Y+(l)) and div_l(Y-(l))
#                 n_valid_tetra += 1

#                 Yp_vec = Yp[legs_tetra,:]  # get the Y+ vector of the tetrahedron
#                 divYp = calculate_divF_tetrahedron(dx_vec, Yp_vec)

#                 Ym_vec = Ym[legs_tetra,:]  # get the Y- vector of the tetrahedron
#                 divYm = calculate_divF_tetrahedron(dx_vec, Ym_vec)

#                 divYp_tetra.append(divYp)
#                 divYm_tetra.append(divYm)


#         divYp_tetra = np.array(divYp_tetra)
#         divYm_tetra = np.array(divYm_tetra)

#         divYp_subarr[itau-ind_tau_start] = np.mean(divYp_tetra)  # mean of all tetrahedrons
#         divYm_subarr[itau-ind_tau_start] = np.mean(divYm_tetra)  # mean of all tetrahedrons

#         divYp_std_subarr[itau-ind_tau_start] = np.std(divYp_tetra)  # std of all tetrahedrons
#         divYm_std_subarr[itau-ind_tau_start] = np.std(divYm_tetra)  # std of all tetrahedrons


#         # # plot a histogram of the divergence
#         # fig, subs = plt.subplots(2, 1, figsize=(8, 6))
#         # sub_p = subs[0]
#         # sub_m = subs[1]

#         # sub_p.hist(divYp_tetra, bins=20, alpha=0.5, label='divY+')
#         # sub_m.hist(divYm_tetra, bins=20, alpha=0.5, label='divY-')

#         # sub_p.set_xlabel('divY+')
#         # sub_m.set_xlabel('divY-')

#         # sub_p.set_ylabel('Count')
#         # sub_m.set_ylabel('Count')

#         # sub_p.legend()
#         # sub_m.legend()
#         # plt.tight_layout()
#         # plt.show()


#         # print('Number of valid tetrahedrons: {}'.format(n_valid_tetra))
#         # print('Number of invalid tetrahedrons: {}'.format(n_invalid_tetra))
#         # print('Number of invalid tetrahedrons with 00: {}'.format(n_invalid_with_00))
#         # print('Number of invalid tetrahedrons with 00 and one repeat pair: {}'.format(n_invalid_with_00_and_one_pair))
#         # print('Number of invalid tetrahedrons with 00 and no repeat pairs: {}'.format(n_invalid_with_00_and_no_pair))
#         # print('Number of invalid tetrahedrons without repeat pairs: {}'.format(n_invalid_tetra_without_repeats))
#         # print('Number of invalid tetrahedrons without repeat pairs and without 00: {}'.format(n_invalid_tetra_without_repeats_and_00))
#         # print('Number of invalid tetrahedrons with two repeat pairs: {}'.format(n_invalid_with_two_pairs))
#         # print('Number of invalid tetrahedrons with one pair and without 00: {}'.format(n_invalid_tetra_one_pair_no_00))

#         # exit()

        
#     print('\n')


#     # unit of div Y is (km/s)^3 / km = km^2 / s^3 = 10^6 m^2 / s^3 = 10^6 J / (kg s)

#     # according to div(Y) = -4 * epsilon, we can calculate epsilon:
#     # epsilon_p = -0.25 * divYp
#     # epsilon_m = -0.25 * divYm

#     epsilon_p_div = -0.25 * divYp_subarr
#     epsilon_m_div = -0.25 * divYm_subarr

#     epsilon_p_div_std = -0.25 * divYp_std_subarr
#     epsilon_m_div_std = -0.25 * divYm_std_subarr



#     '''
#     # plot in two panels
#     fig,subs = plt.subplots(2,1,figsize=(8, 8))

#     sub_p = subs[0]
#     sub_m = subs[1]

#     sub_p.plot(dX_abs_subarr, divYp_subarr, label='divY+')
#     # plot a shade 
#     low_bound = divYp_subarr - divYp_std_subarr
#     high_bound = divYp_subarr + divYp_std_subarr
#     sub_p.fill_between(dX_abs_subarr, low_bound, high_bound, alpha=0.2, color='C1')

#     sub_m.plot(dX_abs_subarr, divYm_subarr, label='divY-')
#     low_bound = divYm_subarr - divYm_std_subarr
#     high_bound = divYm_subarr + divYm_std_subarr
#     sub_m.fill_between(dX_abs_subarr, low_bound, high_bound, alpha=0.2, color='C2')

#     sub_p.set_xlabel(r'$\Delta X$ km')
#     sub_m.set_xlabel(r'$\Delta X$ km')
#     sub_p.set_ylabel(r'$\nabla_{l} \cdot Y^+$')
#     sub_m.set_ylabel(r'$\nabla_{l} \cdot Y^-$')

#     plt.tight_layout()

#     # fig.savefig('./figure/divYpYm_' + trange0_str + '_shorter_interval.png', dpi=300)

#     plt.show()
#     plt.close(fig)
#     '''



#     # # plot Y_{1,1}, Y_{1,2}, Y_{1,3}, Y_{1,4} ------------------------------
#     # fig, subs = plt.subplots(2, 1, figsize=(8,6))

#     # sub_p = subs[0]
#     # sub_m = subs[1]

#     # isat0 = 0
#     # for isat in range(Nsat):
#     #     Yp_vec = Yp_arr[:,isat0, isat, :]
#     #     Yp_para = np.dot(Yp_vec, dX_direction)

#     #     Ym_vec = Ym_arr[:,isat0, isat, :]
#     #     Ym_para = np.dot(Ym_vec, dX_direction)


#     #     sub_p.plot(dX_abs, np.abs(Yp_para), label=r'$Y^+_{%1d,%1d}$' % (isat0+1, isat+1)) 
#     #     sub_m.plot(dX_abs, np.abs(Ym_para), label=r'$Y^-_{%1d,%1d}$' % (isat0+1, isat+1))

#     # sub_p.set_ylabel(r'$Y^+_{\parallel}$')
#     # sub_m.set_ylabel(r'$Y^-_{\parallel}$')

#     # sub_p.set_yscale('log',base=10)
#     # sub_p.set_xscale('log',base=10)
#     # sub_m.set_yscale('log',base=10)
#     # sub_m.set_xscale('log',base=10)


#     # sub_m.set_xlabel(r'$\Delta X$ km')
#     # # sub_m.set_xlabel(r'$\tau$ (s)')

#     # sub_p.legend()
#     # sub_m.legend()
#     # plt.tight_layout()
#     # plt.show()



    
#     # plot the self-increment Yp_parallel 
#     epsilon_p_para_list = []
#     epsilon_m_para_list = []

#     fig, subs = plt.subplots(2, 1, figsize=(8,6))

#     sub_p = subs[0]
#     sub_m = subs[1]


#     for isat in range(Nsat):
#         # print('MMS', mms_id_list[isat])
#         Yp_vec = Yp_arr[:,isat, isat, :]
#         Yp_para = np.dot(Yp_vec, dX_direction)

#         # print(Yp_para)

#         Ym_vec = Ym_arr[:,isat, isat, :]
#         Ym_para = np.dot(Ym_vec, dX_direction)

#         # according to Y_parallel = -4/3 * l * epsilon, we can calculate epsilon:
#         # epsilon = -3/4 * d Y_parallel / d l
#         epsilon_p_para = np.zeros(Yp_para.shape)
#         epsilon_m_para = np.zeros(Ym_para.shape)

#         for itau in range(len(Yp_para)-1):
#             epsilon_p_para[itau] = -3/4 * (Yp_para[itau+1] - Yp_para[itau]) / (dX_abs[itau+1] - dX_abs[itau])
#         epsilon_p_para[-1] = -3/4 * (Yp_para[-1] - Yp_para[-2]) / (dX_abs[-1] - dX_abs[-2])  # last point

#         for itau in range(len(Ym_para)-1):
#             epsilon_m_para[itau] = -3/4 * (Ym_para[itau+1] - Ym_para[itau]) / (dX_abs[itau+1] - dX_abs[itau])
#         epsilon_m_para[-1] = -3/4 * (Ym_para[-1] - Ym_para[-2]) / (dX_abs[-1] - dX_abs[-2])  # last point

#         epsilon_p_para_list.append(epsilon_p_para)
#         epsilon_m_para_list.append(epsilon_m_para)


#         # plot 
#         sub_p.plot(dX_abs, Yp_para, label=f'MMS{mms_id_list[isat]}') # dX_abs
#         sub_m.plot(dX_abs, Ym_para, label=f'MMS{mms_id_list[isat]}')

#         # sub_p.plot(tau_arr, Yp_para, label=f'MMS{mms_id_list[isat]}') # dX_abs
#         # sub_m.plot(tau_arr, Ym_para, label=f'MMS{mms_id_list[isat]}')


#         # x_ref = np.linspace(dX_abs[0], dX_abs[-1], 100)
#         # coeff = np.polyfit(dX_abs, np.abs(Yp_para), 1)
#         # slope = coeff[0]
#         # print('Yp_parallel slope:',slope)

#         # slope_arr = np.linspace(slope/5,slope*5,5)
#         # for islp in range(len(slope_arr)):
#         #     y_ref = x_ref * slope_arr[islp] 
#         #     sub_p.plot(x_ref, y_ref, linestyle='--', color='k', alpha=0.5)

#         # coeff = np.polyfit(dX_abs, np.abs(Ym_para), 1)
#         # slope = coeff[0]
#         # print('Ym_parallel slope:',slope)

#         # slope_arr = np.linspace(slope/5,slope*5,5)
#         # for islp in range(len(slope_arr)):
#         #     y_ref = x_ref * slope_arr[islp] 
#         #     sub_m.plot(x_ref, y_ref, linestyle='--', color='k', alpha=0.5)

#     sub_p.set_ylabel(r'$Y^+_{\parallel}$')
#     sub_m.set_ylabel(r'$Y^-_{\parallel}$')

#     # sub_p.set_yscale('log',base=10)
#     # sub_p.set_xscale('log',base=10)
#     # sub_m.set_yscale('log',base=10)
#     # sub_m.set_xscale('log',base=10)


#     sub_m.set_xlabel(r'$\Delta X$ km')
#     # sub_m.set_xlabel(r'$\tau$ (s)')

#     sub_p.legend()
#     sub_m.legend()
#     plt.tight_layout()

#     #fig.savefig('./figure/YpYm_parallel_' + trange0_str + '_shorter_interval.png', dpi=300)
#     plt.show()
#     plt.close(fig)
    




#     # plot the div Y and Y_parallel into one figure
#     fig, subs = plt.subplots(2,1, figsize=(8,6))

#     sub = subs[0]
#     sub.plot(dX_abs_subarr, divYp_subarr, label=r'$ \nabla_{l} \cdot Y^+$')
#     # plot a shade 
#     low_bound = divYp_subarr - divYp_std_subarr
#     high_bound = divYp_subarr + divYp_std_subarr
#     sub.fill_between(dX_abs_subarr, low_bound, high_bound, alpha=0.2, color='gray')

#     sub.plot(dX_abs_subarr, divYm_subarr, label=r'$ \nabla_{l} \cdot Y^-$')
#     low_bound = divYm_subarr - divYm_std_subarr
#     high_bound = divYm_subarr + divYm_std_subarr
#     sub.fill_between(dX_abs_subarr, low_bound, high_bound, alpha=0.2, color='gray')

#     sub.legend()

#     sub.set_xlabel(r'$\Delta X$ km')
#     sub.set_ylabel(r'$\nabla_{l} \cdot Y^\pm$ km$^2$/s$^3$')



#     sub = subs[1]
#     epsilon_p_para_list = []
#     epsilon_m_para_list = []

#     for isat in range(Nsat):
#         # print('MMS', mms_id_list[isat])
#         Yp_vec = Yp_arr[:,isat, isat, :]
#         Yp_para = np.dot(Yp_vec, dX_direction)

#         # print(Yp_para)

#         Ym_vec = Ym_arr[:,isat, isat, :]
#         Ym_para = np.dot(Ym_vec, dX_direction)

#         # according to Y_parallel = -4/3 * l * epsilon, we can calculate epsilon:
#         # epsilon = -3/4 * d Y_parallel / d l
#         epsilon_p_para = np.zeros(Yp_para.shape)
#         epsilon_m_para = np.zeros(Ym_para.shape)

#         for itau in range(len(Yp_para)-1):
#             epsilon_p_para[itau] = -3/4 * (Yp_para[itau+1] - Yp_para[itau]) / (dX_abs[itau+1] - dX_abs[itau])
#         epsilon_p_para[-1] = -3/4 * (Yp_para[-1] - Yp_para[-2]) / (dX_abs[-1] - dX_abs[-2])  # last point

#         for itau in range(len(Ym_para)-1):
#             epsilon_m_para[itau] = -3/4 * (Ym_para[itau+1] - Ym_para[itau]) / (dX_abs[itau+1] - dX_abs[itau])
#         epsilon_m_para[-1] = -3/4 * (Ym_para[-1] - Ym_para[-2]) / (dX_abs[-1] - dX_abs[-2])  # last point

#         epsilon_p_para_list.append(epsilon_p_para)
#         epsilon_m_para_list.append(epsilon_m_para)


#         # plot 
#         sub.plot(dX_abs, Yp_para, label=f'MMS{mms_id_list[isat]}',color=f'C{isat}') # dX_abs
#         sub.plot(dX_abs, Ym_para ,color=f'C{isat}',linestyle='--')


#     sub.set_ylabel(r'$Y^\pm_{\parallel}$ km$^3$/s$^3$')
#     # sub.set_yscale('log',base=10)
#     # sub.set_xscale('log',base=10)


#     sub.set_xlabel(r'$\Delta X$ km')
#     leg = sub.legend(loc='upper left')

#     sub.add_artist(leg)

#     line1, = sub.plot(0,0,color='k',label=r'$Y^+_{\parallel}$')
#     line2, = sub.plot(0,0,color='k',linestyle='--',label=r'$Y^-_{\parallel}$')
#     sub.legend(handles=[line1, line2], loc='lower right')


#     plt.tight_layout()

#     fig.savefig('./figure/divY_and_Y_parallel_' + trange0_str + '.png', dpi=200)
#     plt.show()







#     # plot epsilon calculated using different methods
#     epsilon_p_para_avg = np.mean(np.array(epsilon_p_para_list), axis=0)
#     epsilon_m_para_avg = np.mean(np.array(epsilon_m_para_list), axis=0)

#     fig, subs = plt.subplots(2,1, figsize=(8,6))

#     sub = subs[0]

#     sub.plot(dX_abs_subarr, epsilon_p_div, label=r'-$\nabla \cdot Y^+$ / 4')
#     sub.plot(dX_abs_subarr, epsilon_m_div, label=r'-$\nabla \cdot Y^-$ / 4')

#     sub.fill_between(dX_abs_subarr, epsilon_p_div - epsilon_p_div_std, epsilon_p_div + epsilon_p_div_std, alpha=0.2)
#     sub.fill_between(dX_abs_subarr, epsilon_m_div - epsilon_m_div_std, epsilon_m_div + epsilon_m_div_std, alpha=0.2)

#     sub.set_xlim([dX_abs[0], dX_abs[-1]])
#     sub.legend()


#     sub = subs[1]
#     sub.plot(dX_abs, epsilon_p_para_avg, label=r'-$\frac{3}{4} \frac{d Y^+_{\parallel}}{d \Delta X}$')
#     sub.plot(dX_abs, epsilon_m_para_avg, label=r'-$\frac{3}{4} \frac{d Y^-_{\parallel}}{d \Delta X}$')

#     sub.set_xlim([dX_abs[0], dX_abs[-1]])

#     sub.legend()

#     plt.tight_layout()
#     # fig.savefig('./figure/epsilon_comparison_' + trange0_str + '_shorter_interval.png', dpi=300)
#     plt.show()
#     plt.close(fig)


#     # print(Yp_arr.shape)

# exit()





# read orbit data for all MMS probes
mms_load_mec(trange=trange, probe=mms_id_list, time_clip=True)

mms_mec_r_gse = []
for mms_id in mms_id_list:
    mms_mec_r_gse.append(get_data(f'mms{mms_id}_mec_r_gse', units=True))  # dt=True


# get the coordinates of the four MMS probes in GSE
t_mec_mms = []
coord_mms = []
coord_mms_RE = []
for i, mms_id in enumerate(mms_id_list):
    t_mec_mms.append(mms_mec_r_gse[i].times)
    coord_mms.append(mms_mec_r_gse[i].y)
    coord_mms_RE.append(coord_mms[i] / Re)  # convert to RE




# get velocity of the four MMS probes
mms_mec_v_gse = []
for mms_id in mms_id_list:
    mms_mec_v_gse.append(get_data(f'mms{mms_id}_mec_v_gse', units=True))  # dt=True


v_mms = []
for i, mms_id in enumerate(mms_id_list):
    v_mms.append(mms_mec_v_gse[i].y)  # velocity in GSE coordinates


# # plot velocity of the four MMS probes
# fig, subs = plt.subplots(4,1, figsize=(10, 10))
# for i, mms_id in enumerate(mms_id_list):
#     sub = subs[i]
#     sub.plot(t_mec_mms[i], v_mms[i][:,0], label=r'$V_x$')
#     sub.plot(t_mec_mms[i], v_mms[i][:,1], label=r'$V_y$')
#     sub.plot(t_mec_mms[i], v_mms[i][:,2], label=r'$V_z$')
#     sub.set_xlabel('Time (s)')
#     sub.set_ylabel('Velocity (km/s)')
#     sub.set_title(f'MMS{mms_id} Velocity in GSE')
#     sub.legend()

# fig.suptitle('MMS1-4 Velocity in GSE: ' + trange[0] + ' -- ' + trange[1])
# plt.show()
# plt.close(fig)




# judge whether the four time series are identical
allclose_mec_12 = np.allclose(t_mec_mms[0], t_mec_mms[1])
allclose_mec_13 = np.allclose(t_mec_mms[0], t_mec_mms[2])
allclose_mec_14 = np.allclose(t_mec_mms[0], t_mec_mms[3])
if allclose_mec_12 and allclose_mec_13 and allclose_mec_14:
    print("All four MMS probes have identical time series for MEC data.")
else:
    print("The time series for MEC data are not identical across all MMS probes.")
    if not allclose_mec_12:
        print("MMS1 and MMS2 have different time series.")
    if not allclose_mec_13:
        print("MMS1 and MMS3 have different time series.")
    if not allclose_mec_14:
        print("MMS1 and MMS4 have different time series.")

    exit()



t_mec_center = np.copy(t_mec_mms[0])  # use the time of MMS1 as the reference time for all satellites

coord_center_mms = (coord_mms[0] + coord_mms[1] + coord_mms[2] + coord_mms[3]) / Nsat


coord_idff_mms = []
for i in range(Nsat):
    coord_idff_mms.append(coord_mms[i] - coord_center_mms)








# read the DIS data for all MMS probes
fpi_data = mms_load_fpi(trange=trange, probe=mms_id_list, data_rate='fast', time_clip=True)
fgm_data = mms_load_fgm(trange=trange, probe=mms_id_list, data_rate='srvy', time_clip=True)
mec_data = mms_load_mec(trange=trange, probe=mms_id_list, time_clip=True)


# despun ion velocity
for mms_id in mms_id_list:
    subtract(f'mms{mms_id}_dis_bulkv_gse_fast', f'mms{mms_id}_dis_bulkv_spintone_gse_fast')



'''
# check time steps----------------------------------------------------------------------
# get time steps for the four MMS probes - bulkv
t_mms_dis_bulkv_gse_fast = []
for mms_id in mms_id_list:
    t_mms_dis_bulkv_gse_fast.append(get_data(f'mms{mms_id}_dis_bulkv_gse_fast', units=True).times)


dt_mms_dis_bulkv_gse_fast = []
for i in range(Nsat):
    dt_mms_dis_bulkv_gse_fast.append(np.diff(t_mms_dis_bulkv_gse_fast[i]))


dt_ave_mms_dis_bulkv_gse_fast = []
dt_std_mms_dis_bulkv_gse_fast = []
dt_std_norm_mms_dis_bulkv_gse_fast = []
for i in range(Nsat):
    dt_ave_mms_dis_bulkv_gse_fast.append(np.nanmean(dt_mms_dis_bulkv_gse_fast[i]))
    dt_std_mms_dis_bulkv_gse_fast.append(np.nanstd(dt_mms_dis_bulkv_gse_fast[i]))
    dt_std_norm_mms_dis_bulkv_gse_fast.append(dt_std_mms_dis_bulkv_gse_fast[i] / dt_ave_mms_dis_bulkv_gse_fast[i])

for i, mms_id in enumerate(mms_id_list):
    print(f'dt_ave_mms_dis_bulkv_gse_fast MMS{mms_id} = {dt_ave_mms_dis_bulkv_gse_fast[i]:.4f} s')


for i, mms_id in enumerate(mms_id_list):
    print(f"variation(dt)/average(dt) MMS{mms_id}_dis_bulkv_gse_fast = {dt_std_norm_mms_dis_bulkv_gse_fast[i]:.4f}")



# get time steps for the four MMS probes - number density
t_mms_dis_ndens_fast = []
for mms_id in mms_id_list:
    t_mms_dis_ndens_fast.append(get_data(f'mms{mms_id}_dis_numberdensity_fast', units=True).times)

dt_mms_dis_ndens_fast = []
for i in range(Nsat):
    dt_mms_dis_ndens_fast.append(np.diff(t_mms_dis_ndens_fast[i]))

dt_ave_mms_dis_ndens_fast = []
dt_std_mms_dis_ndens_fast = []
dt_std_norm_mms_dis_ndens_fast = []
for i in range(Nsat):
    dt_ave_mms_dis_ndens_fast.append(np.nanmean(dt_mms_dis_ndens_fast[i]))
    dt_std_mms_dis_ndens_fast.append(np.nanstd(dt_mms_dis_ndens_fast[i]))
    dt_std_norm_mms_dis_ndens_fast.append(dt_std_mms_dis_ndens_fast[i] / dt_ave_mms_dis_ndens_fast[i])

for i, mms_id in enumerate(mms_id_list):
    print(f'dt_ave_mms_dis_numberdensity_fast MMS{mms_id} = {dt_ave_mms_dis_ndens_fast[i]:.4f} s')

for i, mms_id in enumerate(mms_id_list):
    print(f"variation(dt)/average(dt) MMS{mms_id}_dis_numberdensity_fast = {dt_std_norm_mms_dis_ndens_fast[i]:.4f}")



# check if the time series of bulkv and number density are identical
allcose_bulkv_ndens_1 = np.allclose(t_mms_dis_bulkv_gse_fast[0], t_mms_dis_ndens_fast[0])
allcose_bulkv_ndens_2 = np.allclose(t_mms_dis_bulkv_gse_fast[1], t_mms_dis_ndens_fast[1])
allcose_bulkv_ndens_3 = np.allclose(t_mms_dis_bulkv_gse_fast[2], t_mms_dis_ndens_fast[2])
allcose_bulkv_ndens_4 = np.allclose(t_mms_dis_bulkv_gse_fast[3], t_mms_dis_ndens_fast[3])

if allcose_bulkv_ndens_1 and allcose_bulkv_ndens_2 and allcose_bulkv_ndens_3 and allcose_bulkv_ndens_4:
    print("All four MMS probes have identical time series for dis_bulkv_gse_fast and dis_numberdensity_fast data.")
else:
    print("The time series for dis_bulkv_gse_fast and dis_numberdensity_fast data are not identical across all MMS probes.")
    if not allcose_bulkv_ndens_1:
        print("MMS1 has different time series.")
    if not allcose_bulkv_ndens_2:
        print("MMS2 has different time series.")
    if not allcose_bulkv_ndens_3:
        print("MMS3 has different time series.")
    if not allcose_bulkv_ndens_4:
        print("MMS4 has different time series.")
    print("Exiting the script due to non-identical time series.")
    exit()



# check time steps for magnetic field
t_mms_fgm_b_gse_srvy_l2 = []
for mms_id in mms_id_list:
    t_mms_fgm_b_gse_srvy_l2.append(get_data(f'mms{mms_id}_fgm_b_gse_srvy_l2', units=True).times)

dt_mms_fgm_b_gse_srvy_l2 = []
for i in range(Nsat):
    dt_mms_fgm_b_gse_srvy_l2.append(np.diff(t_mms_fgm_b_gse_srvy_l2[i]))

dt_ave_mms_fgm_b_gse_srvy_l2 = []
dt_std_mms_fgm_b_gse_srvy_l2 = []
dt_std_norm_mms_fgm_b_gse_srvy_l2 = []
for i in range(Nsat):
    dt_ave_mms_fgm_b_gse_srvy_l2.append(np.nanmean(dt_mms_fgm_b_gse_srvy_l2[i]))
    dt_std_mms_fgm_b_gse_srvy_l2.append(np.nanstd(dt_mms_fgm_b_gse_srvy_l2[i]))
    dt_std_norm_mms_fgm_b_gse_srvy_l2.append(dt_std_mms_fgm_b_gse_srvy_l2[i] / dt_ave_mms_fgm_b_gse_srvy_l2[i])



# the magnetic field data have similar average time steps, but the variations can be large.
for i, mms_id in enumerate(mms_id_list):
    print(f'dt_ave_mms_fgm_b_gse_srvy_l2 MMS{mms_id} = {dt_ave_mms_fgm_b_gse_srvy_l2[i]:.4f} s')
for i, mms_id in enumerate(mms_id_list):
    print(f"variation(dt)/average(dt) MMS{mms_id}_fgm_b_gse_srvy_l2 = {dt_std_norm_mms_fgm_b_gse_srvy_l2[i]:.4f}")

exit()
#---------------------------------------------------------------------
'''



# interplolate all fields to the MMS1 DIS
tinterpol(['mms?_dis_numberdensity_fast','mms?_dis_bulkv_gse_fast','mms?_fgm_b_gse_srvy_l2',
           'mms?_mec_r_gse'], 
          'mms1_dis_numberdensity_fast', suffix='_interp')


# tplot(['mms1_dis_numberdensity_fast', 'mms1_dis_numberdensity_fast_interp',
#        'mms4_fgm_b_gse_srvy_l2', 'mms4_fgm_b_gse_srvy_l2_interp'])  # plot the interpolated data


# subtract average and get v1
for mms_id in mms_id_list:
    # subtract average from velocity
    subtract_average(f'mms{mms_id}_dis_bulkv_gse_fast_interp', f'mms{mms_id}_dis_bulkv1_gse_fast_interp')

    # print average velocity
    v_bulkv = get_data(f'mms{mms_id}_dis_bulkv_gse_fast_interp', units=False).y
    v_bulkv_ave = np.nanmean(v_bulkv, axis=0)
    print(f'Average velocity for MMS{mms_id}: {v_bulkv_ave} km/s')


# exit()


# tplot(['mms1_dis_numberdensity_fast_interp','mms1_dis_bulkv1_gse_fast_interp',
#        'mms1_fgm_b_gse_srvy_l2_interp'], save_png= './figure/' +  trange0_str + '.png')


# exit()





mms_mec_r_gse = []
for mms_id in mms_id_list:
    mms_mec_r_gse.append(get_data(f'mms{mms_id}_mec_r_gse_interp', units=False))  # dt=True

# get the coordinates of the four MMS probes in GSE
t_mec_mms = []
coord_mms = []
for i, mms_id in enumerate(mms_id_list):
    t_mec_mms.append(mms_mec_r_gse[i].times)
    coord_mms.append(mms_mec_r_gse[i].y)

coord_center_mms = (coord_mms[0] + coord_mms[1] + coord_mms[2] + coord_mms[3]) / Nsat


coord_diff_mms = []
for i in range(Nsat):
    coord_diff_mms.append(coord_mms[i] - coord_center_mms)

coord_diff_mms_abs = []
for i in range(Nsat):
    coord_diff_mms_abs.append(np.sqrt(np.sum(coord_diff_mms[i]**2, axis=1)))

norm_variation_coord_diff_mms = []
for i in range(Nsat):
    norm_variation_coord_diff_mms.append((np.nanmax(coord_diff_mms_abs[i]) - np.nanmin(coord_diff_mms_abs[i])) / np.nanmean(coord_diff_mms_abs[i]))

for i, mms_id in enumerate(mms_id_list):
    print(f"Normalized variation of distance from center for MMS{mms_id}: {norm_variation_coord_diff_mms[i]:.4f}")

t_mec_center = np.copy(t_mec_mms[0])



# plot the absolute distances of the four MMS satellites from the center of the tetrahedron
fig, sub = plt.subplots(1,1)
for i, mms_id in enumerate(mms_id_list):
    sub.plot(t_mec_center, coord_diff_mms_abs[i], label=f'MMS{mms_id}')
sub.set_xlabel('Time (s)')
sub.set_ylabel('Distance from Center (km)')
sub.set_title('Distance of MMS Satellites from Center of Tetrahedron')
sub.legend()
fig.suptitle('MMS1-4 Relative Distances: ' + trange[0] + ' -- ' + trange[1])
plt.tight_layout()
plt.show()


# plot the relative positions of the four MMS satellites
fig, subs = plt.subplots(1,3, figsize=(15, 5))

sub = subs[0]
# X-Y plane
sub.plot(coord_diff_mms[0][:,0], coord_diff_mms[0][:,1], label='MMS1')
sub.plot(coord_diff_mms[1][:,0], coord_diff_mms[1][:,1], label='MMS2')
sub.plot(coord_diff_mms[2][:,0], coord_diff_mms[2][:,1], label='MMS3')
sub.plot(coord_diff_mms[3][:,0], coord_diff_mms[3][:,1], label='MMS4')
sub.set_xlabel('X (km) / GSE')
sub.set_ylabel('Y (km) / GSE')
sub.set_title('Relative Positions in X-Y Plane')
sub.legend()

sub = subs[1]
# X-Z plane
sub.plot(coord_diff_mms[0][:,0], coord_diff_mms[0][:,2], label='MMS1')
sub.plot(coord_diff_mms[1][:,0], coord_diff_mms[1][:,2], label='MMS2')
sub.plot(coord_diff_mms[2][:,0], coord_diff_mms[2][:,2], label='MMS3')
sub.plot(coord_diff_mms[3][:,0], coord_diff_mms[3][:,2], label='MMS4')
sub.set_xlabel('X (km) / GSE')
sub.set_ylabel('Z (km) / GSE')
sub.set_title('Relative Positions in X-Z Plane')

sub = subs[2]
# Y-Z plane
sub.plot(coord_diff_mms[0][:,1], coord_diff_mms[0][:,2], label='MMS1')
sub.plot(coord_diff_mms[1][:,1], coord_diff_mms[1][:,2], label='MMS2')
sub.plot(coord_diff_mms[2][:,1], coord_diff_mms[2][:,2], label='MMS3')
sub.plot(coord_diff_mms[3][:,1], coord_diff_mms[3][:,2], label='MMS4')
sub.set_xlabel('Y (km) / GSE')
sub.set_ylabel('Z (km) / GSE')
sub.set_title('Relative Positions in Y-Z Plane')
fig.suptitle('Relative Positions of MMS1-4: ' + trange[0] + ' -- ' + trange[1])
plt.tight_layout()
plt.show()


# plot the trajectory of the four MMS satellites in all three planes
fig,subs = plt.subplots(1,3, figsize=(15, 5))

sub = subs[0]
sub.plot(coord_mms[0][:,0],coord_mms[0][:,1],label='MMS1')
sub.plot(coord_mms[1][:,0],coord_mms[1][:,1],label='MMS2')
sub.plot(coord_mms[2][:,0],coord_mms[2][:,1],label='MMS3')
sub.plot(coord_mms[3][:,0],coord_mms[3][:,1],label='MMS4')
sub.legend()
sub.set_xlabel('X (km) / GSE')
sub.set_ylabel('Y (km) / GSE')

sub = subs[1]
sub.plot(coord_mms[0][:,0],coord_mms[0][:,2],label='MMS1')
sub.plot(coord_mms[1][:,0],coord_mms[1][:,2],label='MMS2')
sub.plot(coord_mms[2][:,0],coord_mms[2][:,2],label='MMS3')
sub.plot(coord_mms[3][:,0],coord_mms[3][:,2],label='MMS4')
sub.set_xlabel('X (km) / GSE')
sub.set_ylabel('Z (km) / GSE')

sub = subs[2]
sub.plot(coord_mms[0][:,1],coord_mms[0][:,2],label='MMS1')
sub.plot(coord_mms[1][:,1],coord_mms[1][:,2],label='MMS2')
sub.plot(coord_mms[2][:,1],coord_mms[2][:,2],label='MMS3')
sub.plot(coord_mms[3][:,1],coord_mms[3][:,2],label='MMS4')
sub.set_xlabel('Y (km) / GSE')
sub.set_ylabel('Z (km) / GSE')

fig.suptitle('MMS1-4 Trajectory: '+ trange[0] + ' -- ' + trange[1])
plt.show()

exit()






# calculate average solar wind speed and density
ndens_avg = []
bulkv_avg = []
B_avg = []

for mms_id in mms_id_list:
    ndens = get_data(f'mms{mms_id}_dis_numberdensity_fast_interp')
    bulkv = get_data(f'mms{mms_id}_dis_bulkv_gse_fast_interp')
    B = get_data(f'mms{mms_id}_fgm_b_gse_srvy_l2_interp')

    t = ndens.times  # use the time of the number density data


    # take the time range 
    t_range = np.where((t >= time_double(trange[0])) & (t <= time_double(trange[1])))[0]

    ndens_avg.append(np.nanmean(ndens.y[t_range]))
    bulkv_avg.append(np.nanmean(bulkv.y[t_range], axis=0))  # average over the time range
    B_avg.append(np.nanmean(B.y[t_range], axis=0))  # average over the time range


    # get average coord_diff_mms in time



# for i, mms_id in enumerate(mms_id_list):
#     print(f'MMS{mms_id} average number density: {ndens_avg[i]:.4f} cm^-3')
#     print(f'MMS{mms_id} average bulk velocity: {bulkv_avg[i][0]:.4f} km/s, {bulkv_avg[i][1]:.4f} km/s, {bulkv_avg[i][2]:.4f} km/s')
#     print(f'MMS{mms_id} average magnetic field: {B_avg[i][0]:.4f} nT, {B_avg[i][1]:.4f} nT, {B_avg[i][2]:.4f} nT')



ndens_avg = np.mean(np.array(ndens_avg))
bulkv_avg = np.mean(np.array(bulkv_avg), axis=0)  # average over the four MMS probes
B_avg = np.mean(np.array(B_avg), axis=0)  # average over the four MMS

print(f'Average number density: {ndens_avg:.4f} cm^-3')
print(f'Average bulk velocity: {bulkv_avg[0]:.4f} km/s, {bulkv_avg[1]:.4f} km/s, {bulkv_avg[2]:.4f} km/s')
print(f'Average magnetic field: {B_avg[0]:.4f} nT, {B_avg[1]:.4f} nT, {B_avg[2]:.4f} nT')


coord_diff_mms_avg = []
for i in range(Nsat):
    # print(coord_diff_mms[i].shape)
    coord_diff_mms_avg.append(np.nanmean(coord_diff_mms[i], axis=0))  # average over the time range





polarity_B0 = np.sign(B_avg[0])



# calculate z+ and z- for each MMS probe 
zp_list = []
zm_list = []
for mms_id in mms_id_list:
    V = get_data(f'mms{mms_id}_dis_bulkv_gse_fast_interp').y
    B = get_data(f'mms{mms_id}_fgm_b_gse_srvy_l2_interp').y[:,0:3]

    t = get_data(f'mms{mms_id}_dis_bulkv_gse_fast_interp').times  # use the time of the bulk velocity data

    t_range = np.where((t >= time_double(trange[0])) & (t <= time_double(trange[1])))[0]

    V = V[t_range]
    B = B[t_range]

    v1 = np.zeros(V.shape)
    b1 = np.zeros(B.shape)
    v1[:,0] = V[:,0] - bulkv_avg[0]  # subtract the average bulk velocity
    v1[:,1] = V[:,1] - bulkv_avg[1]  # subtract the average bulk velocity
    v1[:,2] = V[:,2] - bulkv_avg[2]  # subtract the average bulk velocity

    b1[:,0] = B[:,0] - B_avg[0]  # subtract the average magnetic field
    b1[:,1] = B[:,1] - B_avg[1]  # subtract the average magnetic field
    b1[:,2] = B[:,2] - B_avg[2]  # subtract the average magnetic field

    va1 = b1 * 1e-9 /np.sqrt(ndens_avg * 1e6 * m_p * mu_0) * 1e-3 # km/s

    if polarity_B0 > 0: # sunward magnetic field: outward wave is z+ = u + va
        zp = v1 + va1
        zm = v1 - va1
    else: # anti-sunward magnetic field: outward wave is z+ = u - va
        zp = v1 - va1
        zm = v1 + va1

    
    zp_list.append(zp)
    zm_list.append(zm)


    store_data(f'mms{mms_id}_zp_gse', data={'x': t, 'y': zp})
    store_data(f'mms{mms_id}_zm_gse', data={'x': t, 'y': zm})

    # print(np.nanstd(zp),np.nanstd(zm))  # print the standard deviation of z+ and z- for each MMS probe


# exit()

t_arr = np.copy(t)


tplot(['mms1_zp_gse', 'mms2_zp_gse', 'mms3_zp_gse', 'mms4_zp_gse'])
tplot(['mms1_zm_gse', 'mms2_zm_gse', 'mms3_zm_gse', 'mms4_zm_gse'])


# Begin evaluate the structure function 


# Step 1: give series of tau: time increment 
T_tot = t_arr[-1] - t_arr[0]  # total time duration
dt_arr = t_arr - t_arr[0]

tau_min = 60 # minimum time increment, e.g. 60s
tau_max = 3600 # maximum time increment, e.g. 3600s

ind_tau_min = np.where(dt_arr >= tau_min)[0][0]  # index of the minimum time increment
ind_tau_max = np.where(dt_arr <= tau_max)[0][-1]  # index of the maximum time increment

tau_arr = dt_arr[ind_tau_min:ind_tau_max+1]  # time increments from tau_min to tau_max


ntau = len(tau_arr)  # number of time increments

dX_center = np.zeros([ntau,3])  # distance-increment of the center of the tetrahedron in GSE coordinates



Yp_arr = np.zeros((ntau, Nsat, Nsat, 3))  # Yp_arr[tau, mms_id1, mms_id2, 3] for (dz+^2) * dz- structure function
Ym_arr = np.zeros((ntau, Nsat, Nsat, 3))  # Ym_arr[tau, mms_id1, mms_id2, 3] for (dz-^2) * dz+ structure function




# Step 2: for each tau, we will scan the whole time series and calculate the structure function
for itau in range(len(tau_arr)):
    progress_bar(itau+1,len(tau_arr), bar_length=50)  # print progress bar
    tau = tau_arr[itau]
    dindex = int(tau/dt_arr[1])
    

    dX_center[itau,:] = - bulkv_avg * tau

    # determine the start and end indices for the array
    start_index = 0
    end_index = len(t_arr) - dindex - 1

    if end_index <= start_index:
        print(f"Error: end_index ({end_index}) is not greater than start_index ({start_index}) for tau = {tau}.")
        continue


    if end_index + dindex >= len(t_arr):
        print(f"Error: end_index + dindex ({end_index + dindex}) exceeds the length of t_arr ({len(t_arr)}) for tau = {tau}.")
        continue


    for isat in range(Nsat):
        for jsat in range(Nsat):
            # between MMS(isat+1) and MMS(jsat+1)
            Yp_tmp = np.zeros(3)
            Ym_tmp = np.zeros(3)

            n_int = 0
            for i in range(start_index, end_index + 1):
                index0 = i 
                index1 = i + dindex 

                dzp = zp_list[jsat][index1,:] - zp_list[isat][index0,:]
                dzm = zm_list[jsat][index1,:] - zm_list[isat][index0,:]

                if np.isnan(dzp).any() or np.isnan(dzm).any():
                    continue


                Yp_tmp += (dzp[0]*dzp[0] + dzp[1]*dzp[1] + dzp[2]*dzp[2]) * dzm
                Ym_tmp += (dzm[0]*dzm[0] + dzm[1]*dzm[1] + dzm[2]*dzm[2]) * dzp

                n_int += 1  # count the number of intervals

            Yp_arr[itau, isat, jsat, :] = Yp_tmp / n_int
            Ym_arr[itau, isat, jsat, :] = Ym_tmp / n_int


np.savez('./output/Yaglom_law_' + trange0_str + '_shorter_interval.npz', trange = trange, 
         tau_arr = tau_arr, ndens_avg = ndens_avg, bulkv_avg = bulkv_avg, 
         B_avg = B_avg, coord_diff_mms_avg = coord_diff_mms_avg, 
         dX_center = dX_center, Yp_arr = Yp_arr, Ym_arr = Ym_arr)




exit()
















# t_window to remove spin tones
t_window_de_spin = 20


# dt_in_double = time_double('2018-1-1 00:00:20.000000') - time_double('2018-1-1 00:00:00')  # 20 seconds in double format

# # seems to be in seconds, i.e. equal to t_window_de_spin
# print(dt_in_double)


# name for dis_density : 'mms1_dis_numberdensity_fast'
mms_dis_ndens_fast_names = []
for mms_id in mms_id_list:
    mms_dis_ndens_fast_names.append(f'mms{mms_id}_dis_numberdensity_fast')


# name for velocity : 'mms1_dis_bulkv_gse_fast'

# # this avg_data function cannot give identical time arrays because the time stamps of the four MMS probes are not exactly identical
# mms_dis_ndens_fast_avg_names = avg_data(mms_dis_ndens_fast_names, trange=trange,res=t_window_de_spin,
#                                   suffix='_avg')


# average density data for all MMS probes
mms_dis_ndens_fast_avg = []
mms_dis_ndens_fast_avg_names = []
for i, mms_id in enumerate(mms_id_list):
    avg_data_name = tvariable_smooth(mms_dis_ndens_fast_names[i],trange,t_window_de_spin,suffix='_avg')

    mms_dis_ndens_fast_avg.append(get_data(avg_data_name, units=True))  # dt=True
    mms_dis_ndens_fast_avg_names.append(avg_data_name)  # store the name of the averaged data

# average velocity data for all MMS probes
mms_dis_bulkv_gse_fast_avg = []
mms_dis_bulkv_gse_fast_avg_names = []
for i, mms_id in enumerate(mms_id_list):
    avg_data_name = tvariable_smooth(f'mms{mms_id}_dis_bulkv_gse_fast',trange,t_window_de_spin,suffix='_avg')

    mms_dis_bulkv_gse_fast_avg.append(get_data(avg_data_name, units=True))  # dt=True
    mms_dis_bulkv_gse_fast_avg_names.append(avg_data_name)  # store the name of the averaged data


# average magnetic field data for all MMS probes
mms_fgm_b_gse_srvy_l2_avg = []
mms_fgm_b_gse_srvy_l2_avg_names = []
for i, mms_id in enumerate(mms_id_list):
    avg_data_name = tvariable_smooth(f'mms{mms_id}_fgm_b_gse_srvy_l2',trange,t_window_de_spin,suffix='_avg')

    mms_fgm_b_gse_srvy_l2_avg.append(get_data(avg_data_name, units=True))  # dt=True
    mms_fgm_b_gse_srvy_l2_avg_names.append(avg_data_name)  # store the name of the averaged data




tplot([mms_dis_ndens_fast_avg_names[3],mms_dis_bulkv_gse_fast_avg_names[3],
       mms_fgm_b_gse_srvy_l2_avg_names[3]])  # plot the averaged data



exit()



mms_dis_numberdensity_fast = []
for mms_id in mms_id_list:
    mms_dis_numberdensity_fast.append(get_data(f'mms{mms_id}_dis_numberdensity_fast', units=True))


t_mms_dis_ndens_fast = []
for i in range(Nsat):
    t_mms_dis_ndens_fast.append(mms_dis_numberdensity_fast[i].times)


# print(t_mms1_dis_numberdensity_fast[2348:])
# print(t_mms2_dis_numberdensity_fast[2348:])
# print(t_mms3_dis_numberdensity_fast[2348:])
# print(t_mms4_dis_numberdensity_fast[2348:])

# calculate dt 
dt_mms_dis_ndens_fast = []
for i in range(Nsat):
    dt_mms_dis_ndens_fast.append(np.diff(t_mms_dis_ndens_fast[i]))

dt_ave_mms_dis_ndens_fast = []
dt_std_mms_dis_ndens_fast = []
dt_std_norm_mms_dis_ndens_fast = []
for i in range(Nsat):
    dt_ave_mms_dis_ndens_fast.append(np.nanmean(dt_mms_dis_ndens_fast[i]))
    dt_std_mms_dis_ndens_fast.append(np.nanstd(dt_mms_dis_ndens_fast[i]))
    dt_std_norm_mms_dis_ndens_fast.append(dt_std_mms_dis_ndens_fast[i] / dt_ave_mms_dis_ndens_fast[i])

# # whether the dt_ave are identical for the four MMS probes  -- not exactly identical, but very close
# if dt_ave_mms1_dis_numberdensity_fast == dt_ave_mms2_dis_numberdensity_fast and \
#     dt_ave_mms1_dis_numberdensity_fast == dt_ave_mms3_dis_numberdensity_fast and \
#     dt_ave_mms1_dis_numberdensity_fast == dt_ave_mms4_dis_numberdensity_fast:
#     print("The average time steps for dis_numberdensity_fast data are identical across all MMS probes.")
# else:
#     print("The average time steps for dis_numberdensity_fast data are not identical across all MMS probes.")
#     print(f"MMS1: {dt_ave_mms1_dis_numberdensity_fast:.4f} s")
#     print(f"MMS2: {dt_ave_mms2_dis_numberdensity_fast:.4f} s")
#     print(f"MMS3: {dt_ave_mms3_dis_numberdensity_fast:.4f} s")
#     print(f"MMS4: {dt_ave_mms4_dis_numberdensity_fast:.4f} s")
#     exit()

for i, mms_id in enumerate(mms_id_list):
    print(f"variation(dt)/average(dt) MMS{mms_id}_dis_numberdensity_fast = {dt_std_norm_mms_dis_ndens_fast[i]:.4f}")


if dt_std_norm_mms_dis_ndens_fast[0] > 0.001 or \
   dt_std_norm_mms_dis_ndens_fast[1] > 0.001 or \
   dt_std_norm_mms_dis_ndens_fast[2] > 0.001 or \
   dt_std_norm_mms_dis_ndens_fast[3] > 0.001:
    print("The time series of dis_numberdensity_fast data have large variations in time steps.")
    exit()
else:
    print("The time series of dis_numberdensity_fast data have stable time steps.")

exit()














npoint_width1 = int(t_window_de_spin / dt_ave_mms1_dis_numberdensity_fast)  # number of points in the window
tsmooth('mms1_dis_numberdensity_fast', width=npoint_width1, suffix='_smooth')

npoint_width2 = int(t_window_de_spin / dt_ave_mms2_dis_numberdensity_fast)  # number of points in the window
tsmooth('mms2_dis_numberdensity_fast', width=npoint_width2, suffix='_smooth')

npoint_width3 = int(t_window_de_spin / dt_ave_mms3_dis_numberdensity_fast)  # number of points in the window
tsmooth('mms3_dis_numberdensity_fast', width=npoint_width3, suffix='_smooth')

npoint_width4 = int(t_window_de_spin / dt_ave_mms4_dis_numberdensity_fast)  # number of points in the window
tsmooth('mms4_dis_numberdensity_fast', width=npoint_width4, suffix='_smooth')


# get the smoothed data
mms1_dis_numberdensity_fast_despin = get_data('mms1_dis_numberdensity_fast_smooth', units=True)
mms2_dis_numberdensity_fast_despin = get_data('mms2_dis_numberdensity_fast_smooth', units=True)
mms3_dis_numberdensity_fast_despin = get_data('mms3_dis_numberdensity_fast_smooth', units=True)
mms4_dis_numberdensity_fast_despin = get_data('mms4_dis_numberdensity_fast_smooth', units=True)





# plot the dis_numberdensity_fast data for all four MMS probes

fig, subs = plt.subplots(4, 1, figsize=(10, 8), sharex=True)

sub = subs[0]
sub.plot(t_mms1_dis_numberdensity_fast, mms1_dis_numberdensity_fast.y, label='MMS1')
sub.plot(mms1_dis_numberdensity_fast_despin.times,
         mms1_dis_numberdensity_fast_despin.y, label='MMS1 (smoothed)', linestyle='--')
sub.set_xlabel('Time (s)')
sub.set_ylabel('Density (cm$^{-3}$)')
sub.set_title('MMS1 DIS Number Density Fast')
sub.legend()


sub = subs[1]
sub.plot(t_mms2_dis_numberdensity_fast, mms2_dis_numberdensity_fast.y, label='MMS2')
sub.plot(mms2_dis_numberdensity_fast_despin.times,
            mms2_dis_numberdensity_fast_despin.y, label='MMS2 (smoothed)', linestyle='--')
sub.set_xlabel('Time (s)')
sub.set_ylabel('Density (cm$^{-3}$)')
sub.set_title('MMS2 DIS Number Density Fast')
sub.legend()

sub = subs[2]
sub.plot(t_mms3_dis_numberdensity_fast, mms3_dis_numberdensity_fast.y, label='MMS3')
sub.plot(mms3_dis_numberdensity_fast_despin.times,
            mms3_dis_numberdensity_fast_despin.y, label='MMS3 (smoothed)', linestyle='--')
sub.set_xlabel('Time (s)')
sub.set_ylabel('Density (cm$^{-3}$)')
sub.set_title('MMS3 DIS Number Density Fast')
sub.legend()

sub = subs[3]
sub.plot(t_mms4_dis_numberdensity_fast, mms4_dis_numberdensity_fast.y, label='MMS4')
sub.plot(mms4_dis_numberdensity_fast_despin.times,
            mms4_dis_numberdensity_fast_despin.y, label='MMS4 (smoothed)', linestyle='--')
sub.set_xlabel('Time (s)')
sub.set_ylabel('Density (cm$^{-3}$)')  
sub.set_title('MMS4 DIS Number Density Fast')
sub.legend()

fig.suptitle('MMS1-4 DIS Number Density Fast: ' + trange[0] + ' -- ' + trange[1])
plt.tight_layout()
plt.show()




y = mms1_dis_numberdensity_fast.y 
y_fft = np.fft.rfft(y)
freqs = np.fft.rfftfreq(len(y), d=dt_ave_mms1_dis_numberdensity_fast)

y_smooth = mms1_dis_numberdensity_fast_despin.y
y_smooth_fft = np.fft.rfft(y_smooth)
freqs_smooth = np.fft.rfftfreq(len(y_smooth), d=dt_ave_mms1_dis_numberdensity_fast)


# plot the FFT of the dis_numberdensity_fast data for MMS1
fig, subs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
sub = subs[0]
sub.plot(freqs, np.abs(y_fft)**2, label='MMS1 DIS Number Density Fast FFT')
sub.set_xlabel('Frequency (Hz)')
sub.set_ylabel('Power')
sub.set_title('FFT of MMS1 DIS Number Density Fast')
sub.legend()
sub.set_xscale('log',base=10)
sub.set_yscale('log',base=10)


sub = subs[1]
sub.plot(freqs_smooth, np.abs(y_smooth_fft)**2, label='MMS1 DIS Number Density Fast Smoothed FFT', color='orange')
sub.set_xlabel('Frequency (Hz)')
sub.set_ylabel('Power')
sub.set_title('FFT of MMS1 DIS Number Density Fast Smoothed')
sub.legend()
sub.set_xscale('log',base=10)
sub.set_yscale('log',base=10)

fig.suptitle('MMS1 DIS Number Density Fast FFT: ' + trange[0] + ' -- ' + trange[1])

plt.tight_layout()
plt.show()




# # check if the time series of dis_numberdensity_fast are identical
# # somehow the four times series have differnet lengths, so we cannot use np.allclose directly
# # but they all start at the exact same time, so we do not need to worry too much
# allclose_dis_12 = np.allclose(t_mms1_dis_numberdensity_fast, t_mms2_dis_numberdensity_fast)
# allclose_dis_13 = np.allclose(t_mms1_dis_numberdensity_fast, t_mms3_dis_numberdensity_fast)
# allclose_dis_14 = np.allclose(t_mms1_dis_numberdensity_fast, t_mms4_dis_numberdensity_fast)

# if allclose_dis_12 and allclose_dis_13 and allclose_dis_14:
#     print("All four MMS probes have identical time series for dis_numberdensity_fast data.")
# else:
#     print("The time series for dis_numberdensity_fast data are not identical across all MMS probes.")
#     if not allclose_dis_12:
#         print("MMS1 and MMS2 have different time series.")
#     if not allclose_dis_13:
#         print("MMS1 and MMS3 have different time series.")
#     if not allclose_dis_14:
#         print("MMS1 and MMS4 have different time series.")

#     exit()