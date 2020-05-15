
# Procrustes alignment, basic ASM/Eigenshape

# Procrustes alignment 
from scipy.spatial import procrustes

def procrustes_align(pots):
    newpots = np.zeros(np.shape(pots))
    #disparity = np.zeros((npots,npots))
    for j in range(npots):
        #mtx1, newpots[i+1,:,:], disparity[i,j] = procrustes(pots[i,:,:],pots[j,:,:])
        mtx1, newpots[j,:,:], _ = procrustes(pots[0,:,:],pots[j,:,:])
    newpots[0,:,:] = mtx1
    # Scale to -1:1
    newpots[:,:,0] /= np.max(newpots[:,:,0])
    newpots[:,:,1] /= np.max(newpots[:,:,1])
    return newpots

# Simple Active Shape Model (Eigenfaces)

def asm(file):
    #pots = np.loadtxt('pots_closed.csv',delimiter=',')
    pots = np.loadtxt(file,delimiter=',')
    npoints = np.shape(pots)[1]//2
    # Compute the covariance matrix
    C = np.cov(pots.T)

    # Get the eigenvalues and eigenvectors
    evals,evecs = np.linalg.eig(C)

    # Now need to sort them into descending order
    indices = np.argsort(evals)
    indices = indices[::-1]
    evecs = evecs[:,indices]
    evals = evals[indices]
    evecs = np.real(evecs)  
    evals = np.real(evals)  
    m = np.mean(pots,axis=0).reshape((npoints,2))
    #evecs = evecs.reshape((npoints,2,npoints*2))
    return m, evals, evecs
