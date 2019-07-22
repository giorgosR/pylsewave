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
tube1.NumberofSides = 50
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

# # Hide the scalar bar for this color map if no visible data is colored by it.
# HideScalarBarIfNotNeeded(radiusLUT, renderView1)

# # Properties modified on tube1Display
# tube1Display.Opacity = 0.11

# # Properties modified on tube1Display
# tube1Display.Opacity = 0.36

# # Properties modified on tube1Display
# tube1Display.Opacity = 0.27

# # set active source
# SetActiveSource(multi_)

# # show data in view
# multi_Display = Show(multi_, renderView1)

# # show color bar/color legend
# multi_Display.SetScalarBarVisibility(renderView1, True)

# # get layout
# layout1 = GetLayout()

# # split cell
# layout1.SplitHorizontal(0, 0.5)

# # set active view
# SetActiveView(None)

# # place view in the layout
# layout1.AssignView(2, renderView2)

# # show data in view
# multi_Display_1 = Show(multi_, renderView2)

# # trace defaults for the display properties.
# multi_Display_1.Representation = 'Surface'
# multi_Display_1.ColorArrayName = ['POINTS', 'radius']
# multi_Display_1.LookupTable = radiusLUT
# multi_Display_1.OSPRayScaleArray = 'radius'
# multi_Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
# multi_Display_1.SelectOrientationVectors = 'None'
# multi_Display_1.ScaleFactor = 69.74477589074523
# multi_Display_1.SelectScaleArray = 'radius'
# multi_Display_1.GlyphType = 'Arrow'
# multi_Display_1.GlyphTableIndexArray = 'radius'
# multi_Display_1.DataAxesGrid = 'GridAxesRepresentation'
# multi_Display_1.PolarAxes = 'PolarAxesRepresentation'
# multi_Display_1.GaussianRadius = 34.872387945372616
# multi_Display_1.SetScaleArray = ['POINTS', 'radius']
# multi_Display_1.ScaleTransferFunction = 'PiecewiseFunction'
# multi_Display_1.OpacityArray = ['POINTS', 'radius']
# multi_Display_1.OpacityTransferFunction = 'PiecewiseFunction'

# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# multi_Display_1.ScaleTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 19.75566907693728, 1.0, 0.5, 0.0]

# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# multi_Display_1.OpacityTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 19.75566907693728, 1.0, 0.5, 0.0]

# # show color bar/color legend
# multi_Display_1.SetScalarBarVisibility(renderView2, True)

# # reset view to fit data
# renderView2.ResetCamera()

# # hide data in view
# Hide(multi_, renderView2)

# # set active source
# SetActiveSource(tube1)

# # show data in view
# tube1Display_1 = Show(tube1, renderView2)

# # trace defaults for the display properties.
# tube1Display_1.Representation = 'Surface'
# tube1Display_1.ColorArrayName = ['POINTS', 'radius']
# tube1Display_1.LookupTable = radiusLUT
# tube1Display_1.OSPRayScaleArray = 'radius'
# tube1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
# tube1Display_1.SelectOrientationVectors = 'None'
# tube1Display_1.ScaleFactor = 69.78099973350764
# tube1Display_1.SelectScaleArray = 'radius'
# tube1Display_1.GlyphType = 'Arrow'
# tube1Display_1.GlyphTableIndexArray = 'radius'
# tube1Display_1.DataAxesGrid = 'GridAxesRepresentation'
# tube1Display_1.PolarAxes = 'PolarAxesRepresentation'
# tube1Display_1.GaussianRadius = 34.89049986675382
# tube1Display_1.SetScaleArray = ['POINTS', 'radius']
# tube1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
# tube1Display_1.OpacityArray = ['POINTS', 'radius']
# tube1Display_1.OpacityTransferFunction = 'PiecewiseFunction'

# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# tube1Display_1.ScaleTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 19.75566907693728, 1.0, 0.5, 0.0]

# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# tube1Display_1.OpacityTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 19.75566907693728, 1.0, 0.5, 0.0]

# # show color bar/color legend
# tube1Display_1.SetScalarBarVisibility(renderView2, True)

# # reset view to fit data
# renderView2.ResetCamera()

# # Hide the scalar bar for this color map if no visible data is colored by it.
# HideScalarBarIfNotNeeded(radiusLUT, renderView2)

# # Properties modified on tube1Display_1
# tube1Display_1.Opacity = 0.42

# # Properties modified on tube1Display_1
# tube1Display_1.Opacity = 0.8

# # reset view to fit data
# renderView2.ResetCamera()

# # set active view
# SetActiveView(renderView1)

# # reset view to fit data
# renderView1.ResetCamera()

# # set active view
# SetActiveView(renderView2)

# # reset view to fit data
# renderView2.ResetCamera()

# # create a new 'Extract Block'
# extractBlock1 = ExtractBlock(Input=tube1)

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [14, 13, 12, 11, 10, 9, 8]

# # show data in view
# extractBlock1Display = Show(extractBlock1, renderView2)

# # trace defaults for the display properties.
# extractBlock1Display.Representation = 'Surface'
# extractBlock1Display.ColorArrayName = ['POINTS', 'radius']
# extractBlock1Display.LookupTable = radiusLUT
# extractBlock1Display.OSPRayScaleArray = 'radius'
# extractBlock1Display.OSPRayScaleFunction = 'PiecewiseFunction'
# extractBlock1Display.SelectOrientationVectors = 'None'
# extractBlock1Display.ScaleFactor = 20.494680786132815
# extractBlock1Display.SelectScaleArray = 'radius'
# extractBlock1Display.GlyphType = 'Arrow'
# extractBlock1Display.GlyphTableIndexArray = 'radius'
# extractBlock1Display.DataAxesGrid = 'GridAxesRepresentation'
# extractBlock1Display.PolarAxes = 'PolarAxesRepresentation'
# extractBlock1Display.GaussianRadius = 10.247340393066407
# extractBlock1Display.SetScaleArray = ['POINTS', 'radius']
# extractBlock1Display.ScaleTransferFunction = 'PiecewiseFunction'
# extractBlock1Display.OpacityArray = ['POINTS', 'radius']
# extractBlock1Display.OpacityTransferFunction = 'PiecewiseFunction'

# # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
# extractBlock1Display.ScaleTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 2.541784675028596, 1.0, 0.5, 0.0]

# # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
# extractBlock1Display.OpacityTransferFunction.Points = [0.7034934881058851, 0.0, 0.5, 0.0, 2.541784675028596, 1.0, 0.5, 0.0]

# # hide data in view
# Hide(tube1, renderView2)

# # show color bar/color legend
# extractBlock1Display.SetScalarBarVisibility(renderView2, True)

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [14, 13, 12, 11, 10, 9, 8, 7]

# # update the view to ensure updated data information
# renderView2.Update()

# # set active source
# SetActiveSource(multi_)

# # show data in view
# multi_Display_1 = Show(multi_, renderView2)

# # show color bar/color legend
# multi_Display_1.SetScalarBarVisibility(renderView2, True)

# # set active source
# SetActiveSource(tube1)

# # show data in view
# tube1Display_1 = Show(tube1, renderView2)

# # set active source
# SetActiveSource(extractBlock1)

# # rescale color and/or opacity maps used to exactly fit the current data range
# extractBlock1Display.RescaleTransferFunctionToDataRange(False, True)

# # set scalar coloring
# ColorBy(extractBlock1Display, ('POINTS', 'velocity'))

# # Hide the scalar bar for this color map if no visible data is colored by it.
# HideScalarBarIfNotNeeded(radiusLUT, renderView2)

# # rescale color and/or opacity maps used to include current data range
# extractBlock1Display.RescaleTransferFunctionToDataRange(True, False)

# # show color bar/color legend
# extractBlock1Display.SetScalarBarVisibility(renderView2, True)

# # get color transfer function/color map for 'velocity'
# velocityLUT = GetColorTransferFunction('velocity')

# # get color legend/bar for radiusLUT in view renderView2
# radiusLUTColorBar = GetScalarBar(radiusLUT, renderView2)

# # hide data in view
# Hide(multi_, renderView2)

# # rescale color and/or opacity maps used to exactly fit the current data range
# extractBlock1Display.RescaleTransferFunctionToDataRange(False, True)

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [14, 13, 12, 11, 10, 9, 8, 7, 3, 2, 1]

# # update the view to ensure updated data information
# renderView1.Update()

# # update the view to ensure updated data information
# renderView2.Update()

# # Rescale transfer function
# radiusLUT.RescaleTransferFunction(0.703493488106, 19.7556690769)

# # get opacity transfer function/opacity map for 'radius'
# radiusPWF = GetOpacityTransferFunction('radius')

# # Rescale transfer function
# radiusPWF.RescaleTransferFunction(0.703493488106, 19.7556690769)

# # Rescale transfer function
# velocityLUT.RescaleTransferFunction(-4.91880015688, 273.51826859)

# # get opacity transfer function/opacity map for 'velocity'
# velocityPWF = GetOpacityTransferFunction('velocity')

# # Rescale transfer function
# velocityPWF.RescaleTransferFunction(-4.91880015688, 273.51826859)

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [14, 13, 12, 11, 10, 9, 8, 7, 4, 3, 2, 1]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [14, 13, 12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1]

# # update the view to ensure updated data information
# renderView2.Update()

# # hide data in view
# Hide(tube1, renderView2)

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [14, 13, 12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 16, 15]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [14, 13, 12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [14, 13, 12, 11, 10, 9, 8, 7, 5, 4, 19, 3, 18, 2, 1]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 22, 21, 20, 19, 18, 14, 13]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 26, 25, 24, 23, 22, 21, 20, 19, 18, 14, 13]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 14, 13]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 11, 10, 9, 8, 7, 5, 32, 4, 3, 2, 1, 31, 24, 23, 22, 21, 20, 19, 18, 14, 13]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 11, 10, 9, 8, 35, 7, 34, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 14, 13]

# # update the view to ensure updated data information
# renderView2.Update()

# # Rescale transfer function
# velocityLUT.RescaleTransferFunction(-6.402216607, 273.51826859)

# # Rescale transfer function
# velocityPWF.RescaleTransferFunction(-6.402216607, 273.51826859)

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 39, 11, 38, 10, 37, 9, 36, 8, 35, 7, 34, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 14, 13, 40]

# # update the view to ensure updated data information
# renderView2.Update()

# # Rescale transfer function
# velocityLUT.RescaleTransferFunction(-23.8918270107, 273.51826859)

# # Rescale transfer function
# velocityPWF.RescaleTransferFunction(-23.8918270107, 273.51826859)

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 39, 11, 38, 10, 37, 9, 36, 8, 35, 7, 34, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 43, 42, 14, 41, 13, 40]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 39, 11, 38, 10, 37, 9, 36, 8, 35, 7, 34, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 46, 45, 44, 43, 42, 14, 41, 13, 40]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 39, 11, 38, 10, 37, 9, 36, 8, 35, 7, 34, 5, 4, 3, 2, 1, 24, 23, 22, 21, 48, 20, 19, 18, 47, 46, 45, 44, 43, 42, 14, 41, 13, 40]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [12, 39, 11, 38, 10, 37, 9, 36, 8, 35, 7, 34, 5, 4, 3, 2, 1, 24, 23, 50, 22, 49, 21, 48, 20, 19, 18, 47, 46, 45, 44, 43, 42, 14, 41, 13, 40]

# # update the view to ensure updated data information
# renderView2.Update()

# # Rescale transfer function
# velocityLUT.RescaleTransferFunction(-23.8918270107, 291.291716547)

# # Rescale transfer function
# velocityPWF.RescaleTransferFunction(-23.8918270107, 291.291716547)

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [42, 41, 40, 39, 38, 37, 36, 35, 34, 52, 51, 50, 49, 48, 14, 13, 12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 47, 46, 45, 44, 43]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [42, 41, 40, 39, 38, 37, 36, 35, 34, 54, 53, 52, 51, 50, 49, 48, 14, 13, 12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 47, 46, 45, 44, 43]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [42, 41, 40, 39, 38, 37, 36, 35, 34, 56, 55, 54, 53, 52, 51, 50, 49, 48, 14, 13, 12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 47, 46, 45, 44, 43]

# # update the view to ensure updated data information
# renderView2.Update()

# # Properties modified on extractBlock1
# extractBlock1.BlockIndices = [42, 41, 40, 39, 38, 37, 36, 35, 34, 55, 54, 53, 52, 51, 50, 49, 48, 14, 13, 12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 24, 23, 22, 21, 20, 19, 18, 47, 46, 45, 44, 43]

# # update the view to ensure updated data information
# renderView2.Update()

# # set active source
# SetActiveSource(tube1)

# # show data in view
# tube1Display_1 = Show(tube1, renderView2)

# # set active source
# SetActiveSource(tube1)

# # Properties modified on tube1Display_1
# tube1Display_1.Opacity = 0.49

# # Properties modified on tube1Display_1
# tube1Display_1.Opacity = 0.2

# # set active view
# SetActiveView(renderView1)

# # set active view
# SetActiveView(renderView2)

# #### saving camera placements for all active views

# # current camera placement for renderView1
# renderView1.InteractionMode = '2D'
# renderView1.CameraPosition = [240.21420669555664, 348.7760926410556, 1628.6484414398767]
# renderView1.CameraFocalPoint = [240.21420669555664, 348.7760926410556, 0.0]
# renderView1.CameraParallelScale = 424.08544749575185

# # current camera placement for renderView2
# renderView2.CameraPosition = [253.36714849417365, 348.7760926410556, 1628.5953290877949]
# renderView2.CameraFocalPoint = [240.21420669555664, 348.7760926410556, 0.0]
# renderView2.CameraParallelScale = 424.08544749575185

Show()
Render()
Interact()

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).