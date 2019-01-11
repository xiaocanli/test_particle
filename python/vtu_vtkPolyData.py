pdi = self.GetPolyDataInput()
pdo =  self.GetPolyDataOutput()
newPoints = vtk.vtkPoints()
numPoints = pdi.GetNumberOfPoints()

for i in range(0, numPoints):
    coord = pdi.GetPoint(i)
    x, y, z = coord[:3]
    newPoints.InsertPoint(i, x, y, z)

pdo.SetPoints(newPoints)
aPolyLine = vtk.vtkPolyLine()
aPolyLine.GetPointIds().SetNumberOfIds(numPoints)
for i in range(0,numPoints):
    aPolyLine.GetPointIds().SetId(i, i)

pdo.Allocate(1, 1)
pdo.InsertNextCell(aPolyLine.GetCellType(), aPolyLine.GetPointIds())
