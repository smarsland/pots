
# Resample
def resample(points,npoints):
    npots = np.shape(points)[0]
    newpots = np.zeros((npots,npoints,2))
    newpots_closed = np.zeros((npots,2*npoints-2,2))
    for i in range(npots):
        ind = np.max(np.where(points[i,:,0] != 0)[0])
        p = np.squeeze(points[i,:ind+1,:])
        dp = np.diff(p,axis=0)
        pts = np.zeros(len(dp)+1)
        pts[1:] = np.cumsum(np.sqrt(np.sum(dp*dp,axis=1)))
        newpts = np.linspace(0,pts[-1],npoints)
        newpots[i,:,0] = np.interp(newpts,pts,p[:,0])
        newpots[i,:,1] = np.interp(newpts,pts,p[:,1])
        newpots_closed[i,:npoints,:] = np.copy(newpots[i,:,:]) - newpots[i,0,:]
        extra = np.squeeze(np.copy(newpots[i,-2:0:-1,:])- newpots[i,0,:])
        extra[:,0] = -extra[:,0]
        newpots_closed[i,npoints:,:] = np.copy(extra)
    return newpots, newpots_closed

# ++++ Plotting things +++

# Read the metadata, get some colours sorted

import ast
import matplotlib.colors as mcolors

def metadata():
    f = open('GoodPots/PotIndex.txt','r')
    index = f.readlines()
    f.close()
    a = []
    names = np.zeros((325),dtype=int)
    count = 0
    for x in range(len(index)):
        a.append(ast.literal_eval(index[x]))
        #print(ast.literal_eval(index[x]))
        names[count] = int(ast.literal_eval(index[x])[2])
        count+=1
    c = list(mcolors.CSS4_COLORS)
    return names, c
