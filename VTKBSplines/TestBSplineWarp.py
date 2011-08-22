#!/usr/bin/env python
import vtk
import bsplines

# first, create an image that looks like
# graph paper by combining two image grid
# sources via vtkImageBlend
imageGrid1 = vtk.vtkImageGridSource()
imageGrid1.SetGridSpacing(4,4,0)
imageGrid1.SetGridOrigin(0,0,0)
imageGrid1.SetDataExtent(0,255,0,255,0,0)
imageGrid1.SetDataScalarTypeToUnsignedChar()

imageGrid2 = vtk.vtkImageGridSource()
imageGrid2.SetGridSpacing(16,16,0)
imageGrid2.SetGridOrigin(0,0,0)
imageGrid2.SetDataExtent(0,255,0,255,0,0)
imageGrid2.SetDataScalarTypeToUnsignedChar()

table1 = vtk.vtkLookupTable()
table1.SetTableRange(0,1)
table1.SetValueRange(1.0,0.7)
table1.SetSaturationRange(0.0,1.0)
table1.SetHueRange(0.12,0.12)
table1.SetAlphaRange(1.0,1.0)
table1.Build()

table2 = vtk.vtkLookupTable()
table2.SetTableRange(0,1)
table2.SetValueRange(1.0,0.0)
table2.SetSaturationRange(0.0,0.0)
table2.SetHueRange(0.0,0.0)
table2.SetAlphaRange(0.0,1.0)
table2.Build()

map1 = vtk.vtkImageMapToColors()
map1.SetInputConnection(imageGrid1.GetOutputPort())
map1.SetLookupTable(table1)

map2 = vtk.vtkImageMapToColors()
map2.SetInputConnection(imageGrid2.GetOutputPort())
map2.SetLookupTable(table2)

blend = vtk.vtkImageBlend()
blend.AddInput(map1.GetOutput())
blend.AddInput(map2.GetOutput())

# next, create a ThinPlateSpline transform, which
# will then be used to create the B-spline transform
p1 = vtk.vtkPoints()
p1.SetNumberOfPoints(8)
p1.SetPoint(0,0,0,0)
p1.SetPoint(1,0,255,0)
p1.SetPoint(2,255,0,0)
p1.SetPoint(3,255,255,0)
p1.SetPoint(4,96,96,0)
p1.SetPoint(5,96,159,0)
p1.SetPoint(6,159,159,0)
p1.SetPoint(7,159,96,0)

p2 = vtk.vtkPoints()
p2.SetNumberOfPoints(8)
p2.SetPoint(0,0,0,0)
p2.SetPoint(1,0,255,0)
p2.SetPoint(2,255,0,0)
p2.SetPoint(3,255,255,0)
p2.SetPoint(4,96,159,0)
p2.SetPoint(5,159,159,0)
p2.SetPoint(6,159,96,0)
p2.SetPoint(7,96,96,0)

thinPlate = vtk.vtkThinPlateSplineTransform()
thinPlate.SetSourceLandmarks(p2)
thinPlate.SetTargetLandmarks(p1)
thinPlate.SetBasisToR2LogR()

# convert the thin plate spline into a B-spline, by
# sampling it onto a grid and then computing the
# B-spline coefficients
transformToGrid = vtk.vtkTransformToGrid()
transformToGrid.SetInput(thinPlate)
transformToGrid.SetGridSpacing(16.0, 16.0, 1.0)
transformToGrid.SetGridOrigin(0.0, 0.0, 0.0)
transformToGrid.SetGridExtent(0,16, 0,16, 0,0)

grid = vtk.vtkImageBSplineCoefficients()
grid.SetInputConnection(transformToGrid.GetOutputPort())
grid.Update()

# create the B-spline transform, scale the deformation by
# half to demonstrate how deformation scaling works
transform = vtk.vtkBSplineTransform()
transform.SetCoefficients(grid.GetOutput())
transform.SetDisplacementScale(0.5)
transform.SetBorderModeToZero()

# invert the transform before passing it to vtkImageReslice
transform.Inverse()

# reslice the image through the B-spline transform,
# using B-spline interpolation and the "Repeat"
# boundary condition
prefilter = vtk.vtkImageBSplineCoefficients()
prefilter.SetInputConnection(blend.GetOutputPort())
prefilter.SetBorderModeToRepeat()

reslice = vtk.vtkImageBSplineReslice()
reslice.SetInputConnection(prefilter.GetOutputPort())
reslice.SetResliceTransform(transform)
reslice.WrapOn()
reslice.SetInterpolationModeToBSpline()
reslice.SetOutputSpacing(1.0, 1.0, 1.0)
reslice.SetOutputOrigin(-32.0, -32.0, 0.0)
reslice.SetOutputExtent(0, 319, 0, 319, 0, 0)

# set the window/level to 255.0/127.5 to view full range
iren = vtk.vtkRenderWindowInteractor()
viewer = vtk.vtkImageViewer()
viewer.SetupInteractor(iren)
viewer.SetInputConnection(reslice.GetOutputPort())
viewer.SetColorWindow(255.0)
viewer.SetColorLevel(127.5)
viewer.SetZSlice(0)
viewer.Render()
iren.Start()
# --- end of script --
