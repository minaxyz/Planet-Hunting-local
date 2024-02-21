from math import floor, ceil

def binarySearch(array, element):
    lb = 0
    ub = len(array) - 1
    while ub > lb + 1:
        mid = (lb + ub)>>1
        if array[mid] > element:
            ub = mid - 1
        else:
            lb = mid
    return ub if array[ub] <= element else lb

def interpolationSearch(array, element):
    if element < array[0]:
        return 0
    elif element > array[-1]:
        return len(array) - 1
    
    lb = 0
    ub = len(array) - 1
    f = 0
    while ub > lb + 1:
        mid = (ub - lb)*(element - array[f := floor(lb)])/(array[ceil(ub)] - array[f]) + lb
        if array[floor(mid)] > element:
            ub = mid - 1
        else:
            lb = mid
    return ub if array[ub := floor(ub)] <= element else floor(lb)

def estimatePeriodicSignal(xSum, xiSum, N, s0=0, s1=0, s2=0):
    """A least squares algorithm for the period and phase of a signal. Allows for data to have missing points. 
    
    The points of the periodic signal must conform to the model:
    `ith point = ip + c`,
    where `i` is the index of a point, `p` is the period of the signal, `c` the phase of the signal.

    ----------
    Arguments:
        xSum (float) - The sum of periodic point in the signal (`sum(ith point for each point)`).
        
        xiSum (float) - The sum of the product of each point of the signal and its corresponding index index (`sum(i*ith point for each point)`).
        
        N (int) - The number of points in the periodic signal.
    
    Optional Arguments:
        s0 (int) - The number of missing points.
        
        s1 (int) - The sum of the indexes of the missing points (`sum(i for each missing point)`).
        
        s2 (int) - The sum of the square of indexes of the missing points (`sum(i^2 for each missing point)`).

    ----------
    Returns:
        tuple - (period, phase):
            period (float) - The period of the signal.
            
            phase (float) - The phase of the signal.
    """
    midN = (N-1)/2
    i, j = N - s0, s1 - N*midN
    d1, d2 = midN*s0 - s1, midN*s1 - s2
    period = (12*(i*xiSum +  j*xSum))/(i*(N*(N*N - 1) + 12*d2) + 12*j*d1)
    phase = ((xSum - period*d1)/i) - midN*period
    return period, phase
