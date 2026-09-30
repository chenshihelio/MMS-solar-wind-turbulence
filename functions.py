import numpy as np
from pyspedas import get_data, store_data,time_datetime,time_double
import sys

Re = 6371.2  #km


def calc_PVI(t,arr,ndt):
    # t: 1D array of time
    # arr: array of data
    # ndt: number of time steps to calculate the PVI
    # Note: we have assumed that t is uniformly sampled
    # and that the time step is constant


    if len(arr.shape) == 1: # arr is a 1D array
        # make it a 2D array for compatibility
        arr = np.expand_dims(arr, axis=1)

    nt = arr.shape[0]
    nvar = arr.shape[1]

    # some checks
    if nt != len(t):
        print("Error: t and arr must have the same length")
        return None
    
    if ndt > nt:
        print("Error: ndt must be less than nt")
        return None
    
    if ndt < 1:
        print("Error: ndt must be greater than 0")
        return None
    

    # calculate the PVI - the last ndt points will be lost 
    t_pvi = t[:nt-ndt]
    delta_arr = np.zeros(nt - ndt)

    # first step: calculate |delta arr| throughout the time series
    for i in range(nt-ndt):
        delta_arr[i] = np.linalg.norm(arr[i+ndt,:] - arr[i,:])

    delta_arr_mean = np.nanmean(delta_arr)

    # second step: calculate the PVI
    pvi = delta_arr/delta_arr_mean
    
    return {'t': t_pvi, 'PVI': pvi, 'delta_arr_mean': delta_arr_mean}





def tvariable_smooth(original_tplot_name,trange_4_smooth,delta_t_seconds,suffix='_avg',newname=None):
    # Chen Shi, 2025-06-16
    # this function averages the data exactly based on the given trange_4_smooth and delta_t_seconds
    # useful for averaging different data products with different time resolutions or timestamps

    if len(trange_4_smooth) != 2:
        raise ValueError("trange_4_smooth must be a list or tuple with two elements: [start_time, end_time]")

    # time arrays for averaging data
    t0_double = time_double(trange_4_smooth[0])  # start time in double format
    t1_double = time_double(trange_4_smooth[1])  # end time in double format
    t_bound = np.arange(t0_double, t1_double, delta_t_seconds)  # time bounds in double format
    t_center = t_bound[:-1] + delta_t_seconds / 2  # time centers in double format


    orig_data = get_data(original_tplot_name, units=True)  # dt=True

    # unit = orig_data.unit

    orig_time = orig_data.times
    orig_y = orig_data.y


    if orig_y.ndim == 1:  # if the data is 1D, reshape it to 2D
        orig_y = np.expand_dims(orig_y, axis=1)
    elif orig_y.ndim > 2:  # if the data is more than 2D, raise an error
        raise ValueError(f"Data for {original_tplot_name} has more than 2 dimensions, which is not supported.")
    

    avg_y = np.zeros([len(t_center),orig_y.shape[1]])  # initialize the averaged data array


    for j in range(len(t_center)):
        # find the indices of the original data that fall within the current time window
        mask = (orig_time >= t_bound[j]) & (orig_time < t_bound[j+1])
        
        if np.any(mask):
            avg_y[j,:] = np.nanmean(orig_y[mask,:], axis=0)  # average the data in the current time window
        else:
            avg_y[j,:] = np.nan

        
    # create a new tplot variable with the averaged data
    if newname is not None:
        avg_data_name = newname
    else:
        avg_data_name = f'{original_tplot_name}{suffix}'

    success = store_data(avg_data_name, data={'x': t_center, 'y': avg_y})

    if not success:
        print(f"Error: Failed to store data for {avg_data_name}.")
        return None

    return avg_data_name



def calculate_divF_tetrahedron(X_vertices, F_vertices):
    # X_vertices: 4x3 array of the vertices of the tetrahedron
    # F_vertices: 4x3 array of the face centers of the tetrahedron

    if X_vertices.shape != (4, 3) or F_vertices.shape != (4, 3):
        raise ValueError("X_vertices and F_vertices must be 4x3 arrays.")

    X_center = np.average(X_vertices, axis=0)  # center of the tetrahedron

    ind_list = [0,1,2,3]

    v1 = X_vertices[1,:] - X_vertices[0,:]
    v2 = X_vertices[2,:] - X_vertices[0,:]
    v3 = X_vertices[3,:] - X_vertices[0,:]

    volume = np.abs(np.dot(v1, np.cross(v2, v3))) / 6  # volume of the tetrahedron
    FdotA = 0
    for i in range(4):
        ind_surf = ind_list.copy()
        ind_surf.remove(i)  # remove the index of the vertex that is not on

        F_average = (F_vertices[ind_surf[0],:] + F_vertices[ind_surf[1],:] + F_vertices[ind_surf[2],:]) / 3

        v1 = X_vertices[ind_surf[1],:] - X_vertices[ind_surf[0],:]
        v2 = X_vertices[ind_surf[2],:] - X_vertices[ind_surf[0],:]
        area_vector = np.cross(v1, v2) / 2

        v3 = X_vertices[i,:] - X_center # vector from center of the tetrahedron to the vertex i
        # judge whether the area vector points outward or inward
        if np.dot(area_vector, v3) > 0:
            area_vector = -area_vector

        FdotA = FdotA + np.dot(F_average, area_vector)  # dot product of the face center and the area vector

    return FdotA/volume  # return the divergence of the tetrahedron

def progress_bar(current, total, bar_length=50):
    fraction = current / total
    arrow = '=' * int(fraction * bar_length)
    padding = ' ' * (bar_length - len(arrow))
    percent = int(fraction * 100)
    sys.stdout.write(f'\r[{arrow}{padding}] {percent}%')
    sys.stdout.flush()


if __name__ == "__main__":
    print("You are running the functions.py file directly.")
    print("This is not intended to be run directly.")
    print("Please run the main script instead.")