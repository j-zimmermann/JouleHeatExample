# state file generated using paraview version 5.5.2

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [829, 630]
renderView1.InteractionMode = '2D'
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 0
renderView1.CameraPosition = [0.004489566864676492, 0.0, 2.7320508075688776]
renderView1.CameraFocalPoint = [0.004489566864676492, 0.0, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476
renderView1.Background = [0.32, 0.34, 0.43]
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.XTitleFontFile = ''
renderView1.AxesGrid.YTitleFontFile = ''
renderView1.AxesGrid.ZTitleFontFile = ''
renderView1.AxesGrid.XLabelFontFile = ''
renderView1.AxesGrid.YLabelFontFile = ''
renderView1.AxesGrid.ZLabelFontFile = ''

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'PVD Reader'
jouleHeatingTpvd = PVDReader(FileName='/home/julius/Documents/gitlab/Bio-EM/fenics-examples/JouleHeatingFenics/FenicsStudy/JouleHeatingT.pvd')
jouleHeatingTpvd.PointArrays = ['f_29']

# create a new 'Annotate Time Filter'
annotateTimeFilter1 = AnnotateTimeFilter(Input=jouleHeatingTpvd)
annotateTimeFilter1.Format = 'Time: %.0f'
annotateTimeFilter1.Scale = 50.0

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from jouleHeatingTpvd
jouleHeatingTpvdDisplay = Show(jouleHeatingTpvd, renderView1)

# get color transfer function/color map for 'f_29'
f_29LUT = GetColorTransferFunction('f_29')
f_29LUT.AutomaticRescaleRangeMode = 'Never'
f_29LUT.RGBPoints = [293.0, 0.231373, 0.298039, 0.752941, 342.87910473450575, 0.865003, 0.865003, 0.865003, 392.75820946901155, 0.705882, 0.0156863, 0.14902]
f_29LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'f_29'
f_29PWF = GetOpacityTransferFunction('f_29')
f_29PWF.Points = [293.0, 0.0, 0.5, 0.0, 392.75820946901155, 1.0, 0.5, 0.0]
f_29PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
jouleHeatingTpvdDisplay.Representation = 'Surface'
jouleHeatingTpvdDisplay.ColorArrayName = ['POINTS', 'f_29']
jouleHeatingTpvdDisplay.LookupTable = f_29LUT
jouleHeatingTpvdDisplay.OSPRayScaleArray = 'f_29'
jouleHeatingTpvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
jouleHeatingTpvdDisplay.SelectOrientationVectors = 'None'
jouleHeatingTpvdDisplay.ScaleFactor = 0.1
jouleHeatingTpvdDisplay.SelectScaleArray = 'f_29'
jouleHeatingTpvdDisplay.GlyphType = 'Arrow'
jouleHeatingTpvdDisplay.GlyphTableIndexArray = 'f_29'
jouleHeatingTpvdDisplay.GaussianRadius = 0.005
jouleHeatingTpvdDisplay.SetScaleArray = ['POINTS', 'f_29']
jouleHeatingTpvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
jouleHeatingTpvdDisplay.OpacityArray = ['POINTS', 'f_29']
jouleHeatingTpvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
jouleHeatingTpvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
jouleHeatingTpvdDisplay.SelectionCellLabelFontFile = ''
jouleHeatingTpvdDisplay.SelectionPointLabelFontFile = ''
jouleHeatingTpvdDisplay.PolarAxes = 'PolarAxesRepresentation'
jouleHeatingTpvdDisplay.ScalarOpacityFunction = f_29PWF
jouleHeatingTpvdDisplay.ScalarOpacityUnitDistance = 0.1068434939318827

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
jouleHeatingTpvdDisplay.ScaleTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, 392.75820946901155, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
jouleHeatingTpvdDisplay.OpacityTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, 392.75820946901155, 1.0, 0.5, 0.0]

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
jouleHeatingTpvdDisplay.DataAxesGrid.XTitleFontFile = ''
jouleHeatingTpvdDisplay.DataAxesGrid.YTitleFontFile = ''
jouleHeatingTpvdDisplay.DataAxesGrid.ZTitleFontFile = ''
jouleHeatingTpvdDisplay.DataAxesGrid.XLabelFontFile = ''
jouleHeatingTpvdDisplay.DataAxesGrid.YLabelFontFile = ''
jouleHeatingTpvdDisplay.DataAxesGrid.ZLabelFontFile = ''

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
jouleHeatingTpvdDisplay.PolarAxes.PolarAxisTitleFontFile = ''
jouleHeatingTpvdDisplay.PolarAxes.PolarAxisLabelFontFile = ''
jouleHeatingTpvdDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
jouleHeatingTpvdDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

# show data from annotateTimeFilter1
annotateTimeFilter1Display = Show(annotateTimeFilter1, renderView1)

# trace defaults for the display properties.
annotateTimeFilter1Display.FontFile = ''
annotateTimeFilter1Display.Position = [0.012412545235223163, 0.9344444444444444]

# setup the color legend parameters for each legend in this view

# get color legend/bar for f_29LUT in view renderView1
f_29LUTColorBar = GetScalarBar(f_29LUT, renderView1)
f_29LUTColorBar.Orientation = 'Horizontal'
f_29LUTColorBar.WindowLocation = 'AnyLocation'
f_29LUTColorBar.Position = [0.3349577804583837, 0.8720634920634921]
f_29LUTColorBar.Title = 'T [K]'
f_29LUTColorBar.ComponentTitle = ''
f_29LUTColorBar.TitleFontFile = ''
f_29LUTColorBar.LabelFontFile = ''
f_29LUTColorBar.AutomaticLabelFormat = 0
f_29LUTColorBar.LabelFormat = '%-#6.1f'
f_29LUTColorBar.RangeLabelFormat = '%-#6.1f'
f_29LUTColorBar.ScalarBarLength = 0.3300000000000009

# set color bar visibility
f_29LUTColorBar.Visibility = 1

# show color legend
jouleHeatingTpvdDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(annotateTimeFilter1)
# ----------------------------------------------------------------