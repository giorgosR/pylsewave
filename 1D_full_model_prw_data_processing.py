#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML MultiBlock Data Reader'
multi_ = XMLMultiBlockDataReader(FileName=['D:\\gitRepositories\\pulseWavePy\\multi_00.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_01.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_02.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_03.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_04.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_05.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_06.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_07.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_08.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_09.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_10.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_11.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_12.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_13.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_14.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_15.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_16.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_17.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_18.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_19.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_20.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_21.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_22.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_23.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_24.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_25.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_26.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_27.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_28.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_29.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_30.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_31.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_32.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_33.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_34.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_35.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_36.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_37.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_38.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_39.vtm'])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1565, 859]

# show data in view
multi_Display = Show(multi_, renderView1)

# get color transfer function/color map for 'radius'
radiusLUT = GetColorTransferFunction('radius')

# trace defaults for the display properties.
multi_Display.Representation = 'Surface'
multi_Display.ColorArrayName = ['POINTS', 'radius']
multi_Display.LookupTable = radiusLUT
multi_Display.OSPRayScaleArray = 'radius'
multi_Display.OSPRayScaleFunction = 'PiecewiseFunction'
multi_Display.SelectOrientationVectors = 'None'
multi_Display.ScaleFactor = 69.74477589074523
multi_Display.SelectScaleArray = 'radius'
multi_Display.GlyphType = 'Arrow'
multi_Display.GlyphTableIndexArray = 'radius'
multi_Display.DataAxesGrid = 'GridAxesRepresentation'
multi_Display.PolarAxes = 'PolarAxesRepresentation'
multi_Display.GaussianRadius = 34.872387945372616
multi_Display.SetScaleArray = ['POINTS', 'radius']
multi_Display.ScaleTransferFunction = 'PiecewiseFunction'
multi_Display.OpacityArray = ['POINTS', 'radius']
multi_Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
multi_Display.ScaleTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 19.75566907693728, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
multi_Display.OpacityTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 19.75566907693728, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [240.21042251586914, 348.69646966736764, 10000.0]

# show color bar/color legend
multi_Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Tube'
tube1 = Tube(Input=multi_)
tube1.Scalars = ['POINTS', 'radius']
tube1.Vectors = [None, '1']
tube1.Radius = 6.974477589074523

# Properties modified on tube1
tube1.Vectors = [None, '']
tube1.NumberofSides = 20
tube1.Radius = 0.0
tube1.VaryRadius = 'By Absolute Scalar'
tube1.RadiusFactor = 0.0

# show data in view
tube1Display = Show(tube1, renderView1)

# trace defaults for the display properties.
tube1Display.Representation = 'Surface'
tube1Display.ColorArrayName = ['POINTS', 'radius']
tube1Display.LookupTable = radiusLUT
tube1Display.OSPRayScaleArray = 'radius'
tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tube1Display.SelectOrientationVectors = 'None'
tube1Display.ScaleFactor = 69.78099973350764
tube1Display.SelectScaleArray = 'radius'
tube1Display.GlyphType = 'Arrow'
tube1Display.GlyphTableIndexArray = 'radius'
tube1Display.DataAxesGrid = 'GridAxesRepresentation'
tube1Display.PolarAxes = 'PolarAxesRepresentation'
tube1Display.GaussianRadius = 34.89049986675382
tube1Display.SetScaleArray = ['POINTS', 'radius']
tube1Display.ScaleTransferFunction = 'PiecewiseFunction'
tube1Display.OpacityArray = ['POINTS', 'radius']
tube1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube1Display.ScaleTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 19.75566907693728, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube1Display.OpacityTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 19.75566907693728, 1.0, 0.5, 0.0]

# hide data in view
Hide(multi_, renderView1)

# show color bar/color legend
tube1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Extract Block'
extractBlock1 = ExtractBlock(Input=tube1)

# Properties modified on extractBlock1
extractBlock1.BlockIndices = [8, 9, 14, 12, 13, 10, 11]

# show data in view
extractBlock1Display = Show(extractBlock1, renderView1)

# trace defaults for the display properties.
extractBlock1Display.Representation = 'Surface'
extractBlock1Display.ColorArrayName = ['POINTS', 'radius']
extractBlock1Display.LookupTable = radiusLUT
extractBlock1Display.OSPRayScaleArray = 'radius'
extractBlock1Display.OSPRayScaleFunction = 'PiecewiseFunction'
extractBlock1Display.SelectOrientationVectors = 'None'
extractBlock1Display.ScaleFactor = 20.494680786132815
extractBlock1Display.SelectScaleArray = 'radius'
extractBlock1Display.GlyphType = 'Arrow'
extractBlock1Display.GlyphTableIndexArray = 'radius'
extractBlock1Display.DataAxesGrid = 'GridAxesRepresentation'
extractBlock1Display.PolarAxes = 'PolarAxesRepresentation'
extractBlock1Display.GaussianRadius = 10.247340393066407
extractBlock1Display.SetScaleArray = ['POINTS', 'radius']
extractBlock1Display.ScaleTransferFunction = 'PiecewiseFunction'
extractBlock1Display.OpacityArray = ['POINTS', 'radius']
extractBlock1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
extractBlock1Display.ScaleTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 2.541784675028596, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
extractBlock1Display.OpacityTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 2.541784675028596, 1.0, 0.5, 0.0]

# hide data in view
Hide(tube1, renderView1)

# show color bar/color legend
extractBlock1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(extractBlock1Display, ('POINTS', 'flow'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(radiusLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
extractBlock1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
extractBlock1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'flow'
flowLUT = GetColorTransferFunction('flow')

# rescale color and/or opacity maps used to exactly fit the current data range
extractBlock1Display.RescaleTransferFunctionToDataRange(False, True)

# Rescale transfer function
flowLUT.RescaleTransferFunction(320.206660306, 3362.31913572)

# get opacity transfer function/opacity map for 'flow'
flowPWF = GetOpacityTransferFunction('flow')

# Rescale transfer function
flowPWF.RescaleTransferFunction(320.206660306, 3362.31913572)

# Rescale transfer function
flowLUT.RescaleTransferFunction(242.259727269, 10304.7038051)

# Rescale transfer function
flowPWF.RescaleTransferFunction(242.259727269, 10304.7038051)

# rescale color and/or opacity maps used to exactly fit the current data range
extractBlock1Display.RescaleTransferFunctionToDataRange(False, True)

# hide color bar/color legend
extractBlock1Display.SetScalarBarVisibility(renderView1, False)

# Rescale transfer function
flowLUT.RescaleTransferFunction(242.259727269, 10304.7038051)

# Rescale transfer function
flowPWF.RescaleTransferFunction(242.259727269, 10304.7038051)

# show color bar/color legend
extractBlock1Display.SetScalarBarVisibility(renderView1, True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [240.21042251586914, 348.69646966736764, 10000.0]
renderView1.CameraFocalPoint = [240.21042251586914, 348.69646966736764, 0.0]
renderView1.CameraParallelScale = 382.82537413540354

Show()
Render()
Interact()

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).