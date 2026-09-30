# after a random-based test,
# we find that the number of invalid tetrahedrons is stable
# (1) (ij,ji,kl,lk) -> 15
# (2) (00,ij,ji,any) -> 60
# (3) (ij,ji,ik,jk) -> 48
# (4) (00,ij,ik,jk) -> 32
# (5) (1i,1j,ik,jk)  -> 3*6 = 18
# total: 173

import numpy as np 
from itertools import combinations
from functions import calculate_divF_tetrahedron


def if_repeat_pairs(legs):
    n = len(legs)

    for i in range(n-1):
        leg0 = legs[i]

        for j in range(i+1,n):
            leg1 = legs[j]
            if leg1[0] == leg0[1] and leg1[1] == leg0[0]:
                return True 
    return False

def if_two_repeat_pairs(legs):
    n = len(legs)

    n_pair = 0
    for i in range(n-1):
        leg0 = legs[i]

        for j in range(i+1,n):
            leg1 = legs[j]
            if leg1[0] == leg0[1] and leg1[1] == leg0[0]:
                n_pair += 1
    if n_pair == 2:
        return True
    else:
        return False

def if_has_00(legs):
    for leg in legs:
        if leg[0] == 0 and leg[1] == 0:
            return True
    return False

def if_same_tetra(tetra0_legs, tetra1_legs, if_exact=False):
    num_same_legs = 0
    for i in range(4):
        leg1 = tetra1_legs[i]

        for j in range(4):
            leg0 = tetra0_legs[j]

            if if_exact:
                if leg1[0] == leg0[0] and leg1[1] == leg0[1]:
                    num_same_legs += 1
            else:
                if (leg1[0] == leg0[1] and leg1[1] == leg0[0]) \
                    or (leg1[0] == leg0[0] and leg1[1] == leg0[1]):
                    num_same_legs += 1

    if num_same_legs == 4:
        return True
    else:
        return False



# Known linear vector field: F(x) = A x + b
A = np.array([[1.0, 2.0, 0.0],
              [0.0, -2.0, 1.0],
              [0.0, 0.0, 2.0]])
b = np.array([0.5, -0.5, 1.0])

# Analytical divergence is the trace of A
div_exact = np.trace(A)

# # Define tetrahedron vertices
# x0 = np.array([0.0, 0.0, 0.0])
# x1 = np.array([1.0, 0.0, 0.0])
# x2 = np.array([0.0, 1.0, 0.0])
# x3 = np.array([0.0, 0.0, 1.0])
# x = np.array([x0, x1, x2, x3])  # Shape (4, 3)

x_center = np.random.rand(3) * 2 - 1 # np.array([0.13, 0.5, -0.1])

x = np.zeros((4, 3))  # Shape (4, 3)
for i in range(4):
    rand_deviate= np.random.rand(3) - 0.5
    x[i,:] = x_center + rand_deviate

print("Tetrahedron vertices (x):")
print(x)

# Compute F at each vertex
f = np.array([A @ xi + b for xi in x])  # Shape (4, 3)

# Evaluate your divergence function
div_computed = calculate_divF_tetrahedron(x, f)

# Compare with exact divergence
print(f"Computed divergence: {div_computed}")
print(f"Exact divergence:    {div_exact}")
print(f"Error:               {abs(div_computed - div_exact)}")


exit()









# test the invalid tetrahedron------------------------------------------
# four points : corrrespoinding to four satellites 
# generate random array of shape (4,3)
x_sat = (np.random.rand(4, 3) - 0.5)*20  # random positions in 3D space

print(x_sat)

x_cent = np.mean(x_sat, axis=0)  # center of the four points
coord_diff = x_sat - x_cent  # relative position of each point w.r.t



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

# second, find all combinations of 4 legs from the 13 legs
tetra_info = np.zeros((ncomb,4),dtype=int) # (tetrahedron index, leg index)
comb = combinations(range(nleg), 4)  # combinations of 4 legs from the 13 legs
ind_tetra = 0
for c in comb:
    for ic in range(4):
        tetra_info[ind_tetra, ic] = c[ic]

    ind_tetra += 1


# print(tetra_info)
# exit()

# third, calculate the relative spatial difference, w.r.t. x_center, of these legs
dx_list = np.zeros((nleg,3))
for ileg in range(1,nleg):
    leg = leg_list[ileg,:]
    if leg[0]==leg[1]:
        print('error: leg[0] == leg[1], should not happen')

    dx_list[ileg,:] = (coord_diff[leg[1]-1,:] - coord_diff[leg[0]-1,:]) 



# print(dx_list)





# third, take each tetrahedron
n_invalid_tetra = 0
n_invalid_with_00 = 0
n_invalid_with_two_pairs = 0
n_invalid_with_00_and_one_pair = 0
n_invalid_with_00_and_no_pair = 0 
n_invalid_tetra_without_repeats = 0
n_invalid_tetra_without_repeats_and_00 = 0
n_invalid_tetra_one_pair_no_00 = 0

tetra_no_repeats_and_00 = []
count_tetra_no_repeats_and_00 = []


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


    # if if_same_tetra(tetra_tmp, legs_this_tetra, if_exact=True):
    #     print('Find tetrahedron #{:03d}, legs: {}'.format(itetra, legs_this_tetra))
    #     print(dx_vec)
    #     print(vol_tetra, vol_norm, vol_tetra/vol_norm)
    #     break 




    if vol_tetra/vol_norm < 1e-8:
        n_invalid_tetra += 1
        # print('Tetrahedron #{:03d} is invalid. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

        if if_has_00(legs_this_tetra):
            n_invalid_with_00 += 1
            # print('Tetrahedron #{:03d} is invalid with 00. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

        if if_two_repeat_pairs(legs_this_tetra):
            n_invalid_with_two_pairs += 1
            # print('Tetrahedron #{:03d} is invalid with two repeat pairs. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

        if if_has_00(legs_this_tetra) and if_repeat_pairs(legs_this_tetra):
            n_invalid_with_00_and_one_pair += 1
            # print('Tetrahedron #{:03d} is invalid with 00 and one repeat pair. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

        if if_has_00(legs_this_tetra) and not if_repeat_pairs(legs_this_tetra):
            n_invalid_with_00_and_no_pair += 1
            # print('Tetrahedron #{:03d} is invalid with 00 and no repeat pairs. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

        if not if_repeat_pairs(legs_this_tetra):
            n_invalid_tetra_without_repeats += 1
            # print('Tetrahedron #{:03d} is invalid without repeat pairs. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

        if if_repeat_pairs(legs_this_tetra) and not if_has_00(legs_this_tetra) and not if_two_repeat_pairs(legs_this_tetra) :
            n_invalid_tetra_one_pair_no_00 += 1
            # print('Tetrahedron #{:03d} is invalid with one pair and without 00. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))

        if not if_repeat_pairs(legs_this_tetra) and not if_has_00(legs_this_tetra):
            n_invalid_tetra_without_repeats_and_00 += 1
            # print('Tetrahedron #{:03d} is invalid without repeat pairs and without 00. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))
            
            

            # if if_same_tetra(tetra_tmp, legs_this_tetra):
            #     print(legs_this_tetra)
            #     continue


            # find whether there are already same tetrahedron in tetra_no_repeats_and_00
            if_same = False
            for itetra, tetra in enumerate(tetra_no_repeats_and_00):
                if if_same_tetra(tetra, legs_this_tetra):
                    if_same = True
                    count_tetra_no_repeats_and_00[itetra] += 1
                    break
            if not if_same:
                tetra_no_repeats_and_00.append(legs_this_tetra)
                count_tetra_no_repeats_and_00.append(1)

            # # draw the tetrahedron
            # fig = plt.figure(figsize=(8,6))
            # sub = fig.add_subplot(111,projection='3d')
            # sub.set_title('Tetrahedron #{:03d} is invalid without repeat pairs and without 00. Legs: {}'.format(itetra, leg_list[legs_tetra,:]))
            
            # sub.scatter(dx_vec[:,0], dx_vec[:,1], dx_vec[:,2], c='r', marker='o')
            # sub.legend()
            # plt.show()

    # print(legs_tetra, dx_vec)

    # exit()

print('Number of invalid tetrahedrons: {}'.format(n_invalid_tetra))
print('Number of invalid tetrahedrons with 00: {}'.format(n_invalid_with_00))
print('Number of invalid tetrahedrons with 00 and one repeat pair: {}'.format(n_invalid_with_00_and_one_pair))
print('Number of invalid tetrahedrons with 00 and no repeat pairs: {}'.format(n_invalid_with_00_and_no_pair))
print('Number of invalid tetrahedrons without repeat pairs: {}'.format(n_invalid_tetra_without_repeats))
print('Number of invalid tetrahedrons without repeat pairs and without 00: {}'.format(n_invalid_tetra_without_repeats_and_00))
print('Number of invalid tetrahedrons with two repeat pairs: {}'.format(n_invalid_with_two_pairs))
print('Number of invalid tetrahedrons with one pair and without 00: {}'.format(n_invalid_tetra_one_pair_no_00))

for i in range(len(tetra_no_repeats_and_00)):
    print(tetra_no_repeats_and_00[i], count_tetra_no_repeats_and_00[i])

exit()