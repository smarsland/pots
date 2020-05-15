
# +++ VTK Format ++++
import vtk
from vtk.util import numpy_support

def writeVTK(newpots,filenames,ext):
    npots,npoints,_ = np.shape(newpots)
    for i in range(npots):
        a = np.zeros((npoints,3))
        a[:,:2] = np.squeeze(newpots[i,:,:])
        VTK_data = numpy_support.numpy_to_vtk(num_array=a, deep=True, array_type=vtk.VTK_FLOAT)

        vtk_points = vtk.vtkPoints()
        vtk_points.SetData(VTK_data)

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(vtk_points)

        writer = vtk.vtkPolyDataWriter()
        writer.SetInputData(polydata)
        writer.SetFileName(filenames[i]+ext+'.vtk')
        writer.Write()

# This reads VTK
def readVTK(filename):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName('filename')
    #reader.SetFileName('output/BayesianAtlas__Reconstruction__skull__subject_australopithecus.vtk')
    reader.Update()
    data = reader.GetOutput()

    return numpy_support.vtk_to_numpy(data.GetPoints().GetData()).astype('float64')
