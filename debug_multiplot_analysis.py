import numpy as np 
from itertools import combinations
from functions import *
import matplotlib.pyplot as plt

mms_id_list = ['1','2','3','4']  # list of MMS probes

Nsat = 4  # number of satellites
trange = ['2019-3-31/18:50:00','2019-3-31/21:50:00']  # tetrahedron configuration



trange0_str = trange[0].replace(':', '-').replace('/', '-')


# the analyzed data--------------------
file_output = './output/Yaglom_law_' + trange0_str + '.npz' # '_shorter_interval.npz'


data = np.load(file_output)

trange = data['trange']
tau_arr = data['tau_arr']
ndens_avg = data['ndens_avg']
bulkv_avg = data['bulkv_avg']
B_avg = data['B_avg']
coord_diff_mms_avg = data['coord_diff_mms_avg']
dX_center = data['dX_center']
Yp_arr = data['Yp_arr']
Ym_arr = data['Ym_arr']




dX_abs = np.sqrt(np.sum(dX_center**2, axis=1))  # absolute distance-increment of the center of the tetrahedron


ntau = len(tau_arr)  # number of time increments

dX_direction = np.zeros((ntau, 3))  # direction of the center of the tetrahedron
for itau in range(ntau):
    dX_direction[itau, :] = dX_center[itau, :] / np.linalg.norm(dX_center[itau, :])  # direction of the center of the tetrahedron

dX_direction = dX_direction[0,:]




# print(dX_direction)






tau_start_analysis = 60 #s
tau_end_analysis = 3600 #s

ind_tau_start = np.where(tau_arr >= tau_start_analysis)[0][0]  # index of the start time increment
ind_tau_end = np.where(tau_arr <= tau_end_analysis)[0][-1]  # index of the end time increment



# For each dX, we have 4*3+1=13 legs
# from the 13 legs, we can form C_{13}^4 = 715 tetrahedrons
nleg = 13
ncomb = 715  # number of combinations of tetrahedrons

# the leg corresponding to 11 or 22 or 33 or 44 is marked as 00. But we will average (11-44).

# fisrt, list all possible 13 legs
leg_list = np.zeros((nleg,2),dtype=int)
leg_list[0,0] = 0
leg_list[0,1] = 0
ind = 1
for i in range(1,5):
    for j in range(1,5):
        if j == i:
            continue

        leg_list[ind,0] = i 
        leg_list[ind,1] = j
        ind += 1

# print(leg_list)
# exit()

# second, find all combinations of 4 legs from the 13 legs
tetra_info = np.zeros((ncomb,4),dtype=int) # (tetrahedron index, leg index)
comb = combinations(range(nleg), 4)  # combinations of 4 legs from the 13 legs
ind_tetra = 0
for c in comb:
    for ic in range(4):
        tetra_info[ind_tetra, ic] = c[ic]

    ind_tetra += 1

# print(tetra_info)


# third, calculate the relative spatial difference, w.r.t. dX_center, of these legs
dx_list = np.zeros((nleg,3))
for ileg in range(1,nleg):
    leg = leg_list[ileg,:]
    if leg[0]==leg[1]:
        print('error: leg[0] == leg[1], should not happen')

    dx_list[ileg,:] = (coord_diff_mms_avg[leg[1]-1,:] - coord_diff_mms_avg[leg[0]-1,:]) 



dx_abs_list = np.sqrt(np.sum(dx_list**2, axis=1))  # absolute distance-increment of the legs

# print(dx_abs_list)



# begin debugging------------------------------

# Method 1: calculate dY/dl_parallel, roughly 
tau_for_estimate_dY_dl_para = np.arange(600,1800+1,100)
dX_abs_for_estimate_dY_dl_para = np.zeros(len(tau_for_estimate_dY_dl_para))  
Yp_para_arr = np.zeros((nleg,len(tau_for_estimate_dY_dl_para)))
Ym_para_arr = np.zeros((nleg,len(tau_for_estimate_dY_dl_para)))

for i in range(len(tau_for_estimate_dY_dl_para)):
    tau_target = tau_for_estimate_dY_dl_para[i]  # target time increment
    itau = np.abs(tau_arr - tau_target).argmin()  # index of the target time increment

    dX_abs_for_estimate_dY_dl_para[i] = dX_abs[itau]  # absolute distance-increment of the center of the tetrahedron
    print('|dX| = {:.3f} km'.format(dX_abs_for_estimate_dY_dl_para[i])) 

    # first, calculate Y+ and Y- for the 13 legs
    Yp = np.zeros((nleg,3))
    Ym = np.zeros((nleg,3))
    for ileg in range(nleg):
        leg = leg_list[ileg,:]

        if leg[0]==0 and leg[1]==0: # leg 00, average 
            for isat in range(Nsat):
                Yp[ileg,:] = Yp[ileg,:] + Yp_arr[itau,isat,isat,:] 
                Ym[ileg,:] = Ym[ileg,:] + Ym_arr[itau,isat,isat,:]
            Yp[ileg,:] = Yp[ileg,:] / Nsat
            Ym[ileg,:] = Ym[ileg,:] / Nsat
        else:
            if leg[0] == leg[1]:
                print('error: leg[0] == leg[1], should not happen')

            Yp[ileg,:] = Yp_arr[itau,leg[0]-1,leg[1]-1,:]
            Ym[ileg,:] = Ym_arr[itau,leg[0]-1,leg[1]-1,:]


        Yp_para_arr[ileg, i] = np.dot(Yp[ileg, :], dX_direction)  # parallel component of Yp
        Ym_para_arr[ileg, i] = np.dot(Ym[ileg, :], dX_direction)  # parallel component of Ym


Yp_para_1234_arr = np.zeros((4, len(tau_for_estimate_dY_dl_para)))
Ym_para_1234_arr = np.zeros((4, len(tau_for_estimate_dY_dl_para)))
for i in range(len(tau_for_estimate_dY_dl_para)):
    tau_target = tau_for_estimate_dY_dl_para[i]  # target time increment
    itau = np.abs(tau_arr - tau_target).argmin()  # index of the target time increment

    for isat in range(Nsat):
        Yp = np.zeros(3)
        Ym = np.zeros(3)
    


        Yp = Yp_arr[itau,isat,isat,:] 
        Ym = Ym_arr[itau,isat,isat,:]

        Yp_para_1234_arr[isat, i] = np.dot(Yp, dX_direction)  # parallel component of Yp
        Ym_para_1234_arr[isat, i] = np.dot(Ym, dX_direction)  # parallel component of Ym


    



fig, subs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

sub = subs[0]

legs_plot = [5,8,6,11,9,12] # [1,4,2,7,3,10]

for ileg in legs_plot:
    sub.plot(dX_abs_for_estimate_dY_dl_para, Yp_para_arr[ileg, :], label='Leg {}{}'.format(leg_list[ileg, 0], leg_list[ileg, 1]))

# auto-calculation is similar for the four satellites
for isat in range(Nsat):
    sub.plot(dX_abs_for_estimate_dY_dl_para, Yp_para_1234_arr[isat, :], linestyle='--', lw=2, label='MMS{}'.format(isat+1))

sub.set_ylabel('Y+ (parallel) [nT]')
sub.legend()

sub = subs[1]
for ileg in legs_plot:
    sub.plot(dX_abs_for_estimate_dY_dl_para, Ym_para_arr[ileg, :], label='Leg {}{}'.format(leg_list[ileg, 0], leg_list[ileg, 1]))

for isat in range(Nsat):
    sub.plot(dX_abs_for_estimate_dY_dl_para, Ym_para_1234_arr[isat, :], linestyle='--', lw=2, label='MMS{}'.format(isat+1))

sub.set_ylabel('Y- (parallel) [nT]')
sub.set_xlabel('|dX| [km]')

sub.legend()

plt.tight_layout()
plt.show()

exit()








# Method 2: use the tetrahedron to calculate div Y+ and div Y-
# tau_target
tau_target = 1800 # seconds
itau = np.abs(tau_arr - tau_target).argmin()  # index of the target time increment


dX_abs = dX_abs[itau]  # absolute distance-increment of the center of the tetrahedron
print('|dX| = {:.3f} km'.format(dX_abs)) 

# first, calculate Y+ and Y- for the 13 legs
Yp = np.zeros((nleg,3))
Ym = np.zeros((nleg,3))
for ileg in range(nleg):
    leg = leg_list[ileg,:]

    if leg[0]==0 and leg[1]==0: # leg 00, average 
        for isat in range(Nsat):
            Yp[ileg,:] = Yp[ileg,:] + Yp_arr[itau,isat,isat,:] 
            Ym[ileg,:] = Ym[ileg,:] + Ym_arr[itau,isat,isat,:]
        Yp[ileg,:] = Yp[ileg,:] / Nsat
        Ym[ileg,:] = Ym[ileg,:] / Nsat
    else:
        if leg[0] == leg[1]:
            print('error: leg[0] == leg[1], should not happen')

        Yp[ileg,:] = Yp_arr[itau,leg[0]-1,leg[1]-1,:]
        Ym[ileg,:] = Ym_arr[itau,leg[0]-1,leg[1]-1,:]


Yp_para = np.dot(Yp, dX_direction)  # parallel component of Yp
Ym_para = np.dot(Ym, dX_direction)  # parallel component of Ym

for ileg in range(nleg):
    print(leg_list[ileg,:], Yp_para[ileg], Ym_para[ileg])


# third, take each tetrahedron
divYp_tetra = []
divYm_tetra = []


for itetra in range(ncomb):
    legs_tetra = tetra_info[itetra,:]  # indices of the 4 legs
    legs_this_tetra = leg_list[legs_tetra,:] # get the legs of this tetrahedron

    dx_vec = dx_list[legs_tetra,:]  # get the dx vector of the legs of this tetrahedron, 4 legs * 3 components


    # determine whether the tetrahedron is valid
    v1 = dx_vec[1,:] - dx_vec[0,:]
    v2 = dx_vec[2,:] - dx_vec[0,:]
    v3 = dx_vec[3,:] - dx_vec[0,:]

    vol_norm = np.linalg.norm(v1) * np.linalg.norm(v2) * np.linalg.norm(v3)
    vol_tetra = np.abs(np.dot(v1, np.cross(v2, v3))) # volume of the tetrahedron

    

    if vol_tetra/vol_norm < 1e-8:
        n_invalid_tetra += 1
        # print('Tetrahedron #{:03d} is invalid. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))
    else:
        # calculate divergence div_l(Y+(l)) and div_l(Y-(l))
        n_valid_tetra += 1

        Yp_vec = Yp[legs_tetra,:]  # get the Y+ vector of the tetrahedron
        divYp = calculate_divF_tetrahedron(dx_vec, Yp_vec)

        Ym_vec = Ym[legs_tetra,:]  # get the Y- vector of the tetrahedron
        divYm = calculate_divF_tetrahedron(dx_vec, Ym_vec)

        divYp_tetra.append(divYp)
        divYm_tetra.append(divYm)


divYp_tetra = np.array(divYp_tetra)
divYm_tetra = np.array(divYm_tetra)


print()