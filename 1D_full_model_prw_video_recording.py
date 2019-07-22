#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1565, 859]

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitHorizontal(0, 0.5)

# set active view
SetActiveView(None)

# place view in the layout
layout1.AssignView(2, renderView2)

# set active view
SetActiveView(renderView1)

# split cell
layout1.SplitHorizontal(1, 0.5)

# set active view
SetActiveView(None)

# place view in the layout
layout1.AssignView(4, renderView3)

# set active view
SetActiveView(renderView2)

# split cell
layout1.SplitHorizontal(2, 0.5)

# set active view
SetActiveView(None)

# place view in the layout
layout1.AssignView(6, renderView4)

# set active view
SetActiveView(renderView1)

# create a new 'XML MultiBlock Data Reader'
multi_ = XMLMultiBlockDataReader(FileName=['D:\\gitRepositories\\pulseWavePy\\multi_00.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_01.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_02.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_03.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_04.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_05.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_06.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_07.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_08.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_09.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_10.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_11.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_12.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_13.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_14.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_15.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_16.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_17.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_18.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_19.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_20.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_21.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_22.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_23.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_24.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_25.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_26.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_27.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_28.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_29.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_30.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_31.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_32.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_33.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_34.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_35.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_36.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_37.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_38.vtm', 'D:\\gitRepositories\\pulseWavePy\\multi_39.vtm'])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# show data in view
multi_Display = Show(multi_, renderView1)

# get color transfer function/color map for 'radius'
radiusLUT = GetColorTransferFunction('radius')
radiusLUT.RGBPoints = [0.702254450766134, 0.231373, 0.298039, 0.752941, 10.2289617638517, 0.865003, 0.865003, 0.865003, 19.7556690769373, 0.705882, 0.0156863, 0.14902]
radiusLUT.ScalarRangeInitialized = 1.0

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

# show color bar/color legend
multi_Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
renderView4.Update()

# create a new 'Tube'
tube1 = Tube(Input=multi_)
tube1.Scalars = ['POINTS', 'radius']
tube1.Vectors = [None, '1']
tube1.Radius = 6.974477589074523

# Properties modified on tube1
tube1.Vectors = [None, '']
tube1.NumberofSides = 40
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

# Rescale transfer function
radiusLUT.RescaleTransferFunction(0.701850270547, 21.5403208384)

# get opacity transfer function/opacity map for 'radius'
radiusPWF = GetOpacityTransferFunction('radius')
radiusPWF.Points = [0.702254450766134, 0.0, 0.5, 0.0, 19.7556690769373, 1.0, 0.5, 0.0]
radiusPWF.ScalarRangeInitialized = 1

# Rescale transfer function
radiusPWF.RescaleTransferFunction(0.701850270547, 21.5403208384)

# set active view
SetActiveView(renderView3)

# set active source
SetActiveSource(tube1)

# show data in view
tube1Display_1 = Show(tube1, renderView3)

# trace defaults for the display properties.
tube1Display_1.Representation = 'Surface'
tube1Display_1.ColorArrayName = ['POINTS', 'radius']
tube1Display_1.LookupTable = radiusLUT
tube1Display_1.OSPRayScaleArray = 'radius'
tube1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
tube1Display_1.SelectOrientationVectors = 'None'
tube1Display_1.ScaleFactor = 69.78091971874237
tube1Display_1.SelectScaleArray = 'radius'
tube1Display_1.GlyphType = 'Arrow'
tube1Display_1.GlyphTableIndexArray = 'radius'
tube1Display_1.DataAxesGrid = 'GridAxesRepresentation'
tube1Display_1.PolarAxes = 'PolarAxesRepresentation'
tube1Display_1.GaussianRadius = 34.890459859371184
tube1Display_1.SetScaleArray = ['POINTS', 'radius']
tube1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
tube1Display_1.OpacityArray = ['POINTS', 'radius']
tube1Display_1.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube1Display_1.ScaleTransferFunction.Points = [0.7022544507661337, 0.0, 0.5, 0.0, 19.68237082607396, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube1Display_1.OpacityTransferFunction.Points = [0.7022544507661337, 0.0, 0.5, 0.0, 19.68237082607396, 1.0, 0.5, 0.0]

# show color bar/color legend
tube1Display_1.SetScalarBarVisibility(renderView3, True)

# reset view to fit data
renderView3.ResetCamera()

# set scalar coloring
ColorBy(tube1Display_1, ('POINTS', 'pressure'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(radiusLUT, renderView3)

# rescale color and/or opacity maps used to include current data range
tube1Display_1.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
tube1Display_1.SetScalarBarVisibility(renderView3, True)

# get color transfer function/color map for 'pressure'
pressureLUT = GetColorTransferFunction('pressure')
pressureLUT.RGBPoints = [0.00856925408043476, 0.231373, 0.298039, 0.752941, 0.009542100934403846, 0.865003, 0.865003, 0.865003, 0.010514947788372907, 0.705882, 0.0156863, 0.14902]
pressureLUT.ScalarRangeInitialized = 1.0

# Rescale transfer function
pressureLUT.RescaleTransferFunction(0.00786017101309, 0.0214987041236)

# get opacity transfer function/opacity map for 'pressure'
pressurePWF = GetOpacityTransferFunction('pressure')
pressurePWF.Points = [0.00856925408043476, 0.0, 0.5, 0.0, 0.010514947788372907, 1.0, 0.5, 0.0]
pressurePWF.ScalarRangeInitialized = 1

# Rescale transfer function
pressurePWF.RescaleTransferFunction(0.00786017101309, 0.0214987041236)

# set active view
SetActiveView(renderView2)

# show data in view
tube1Display_2 = Show(tube1, renderView2)

# trace defaults for the display properties.
tube1Display_2.Representation = 'Surface'
tube1Display_2.ColorArrayName = ['POINTS', 'radius']
tube1Display_2.LookupTable = radiusLUT
tube1Display_2.OSPRayScaleArray = 'radius'
tube1Display_2.OSPRayScaleFunction = 'PiecewiseFunction'
tube1Display_2.SelectOrientationVectors = 'None'
tube1Display_2.ScaleFactor = 69.78096587955952
tube1Display_2.SelectScaleArray = 'radius'
tube1Display_2.GlyphType = 'Arrow'
tube1Display_2.GlyphTableIndexArray = 'radius'
tube1Display_2.DataAxesGrid = 'GridAxesRepresentation'
tube1Display_2.PolarAxes = 'PolarAxesRepresentation'
tube1Display_2.GaussianRadius = 34.89048293977976
tube1Display_2.SetScaleArray = ['POINTS', 'radius']
tube1Display_2.ScaleTransferFunction = 'PiecewiseFunction'
tube1Display_2.OpacityArray = ['POINTS', 'radius']
tube1Display_2.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube1Display_2.ScaleTransferFunction.Points = [0.7026503860269266, 0.0, 0.5, 0.0, 19.729271585784527, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube1Display_2.OpacityTransferFunction.Points = [0.7026503860269266, 0.0, 0.5, 0.0, 19.729271585784527, 1.0, 0.5, 0.0]

# show color bar/color legend
tube1Display_2.SetScalarBarVisibility(renderView2, True)

# reset view to fit data
renderView2.ResetCamera()

# set scalar coloring
ColorBy(tube1Display_2, ('POINTS', 'flow'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(radiusLUT, renderView2)

# rescale color and/or opacity maps used to include current data range
tube1Display_2.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
tube1Display_2.SetScalarBarVisibility(renderView2, True)

# get color transfer function/color map for 'flow'
flowLUT = GetColorTransferFunction('flow')
flowLUT.AutomaticRescaleRangeMode = 'Never'
flowLUT.RGBPoints = [-153930.62598506358, 0.231373, 0.298039, 0.752941, 208717.69433148086, 0.865003, 0.865003, 0.865003, 571366.0146480246, 0.705882, 0.0156863, 0.14902]
flowLUT.ScalarRangeInitialized = 1.0

# set active view
SetActiveView(renderView4)

# show data in view
tube1Display_3 = Show(tube1, renderView4)

# trace defaults for the display properties.
tube1Display_3.Representation = 'Surface'
tube1Display_3.ColorArrayName = ['POINTS', 'radius']
tube1Display_3.LookupTable = radiusLUT
tube1Display_3.OSPRayScaleArray = 'radius'
tube1Display_3.OSPRayScaleFunction = 'PiecewiseFunction'
tube1Display_3.SelectOrientationVectors = 'None'
tube1Display_3.ScaleFactor = 69.78091971874237
tube1Display_3.SelectScaleArray = 'radius'
tube1Display_3.GlyphType = 'Arrow'
tube1Display_3.GlyphTableIndexArray = 'radius'
tube1Display_3.DataAxesGrid = 'GridAxesRepresentation'
tube1Display_3.PolarAxes = 'PolarAxesRepresentation'
tube1Display_3.GaussianRadius = 34.890459859371184
tube1Display_3.SetScaleArray = ['POINTS', 'radius']
tube1Display_3.ScaleTransferFunction = 'PiecewiseFunction'
tube1Display_3.OpacityArray = ['POINTS', 'radius']
tube1Display_3.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tube1Display_3.ScaleTransferFunction.Points = [0.7022544507661337, 0.0, 0.5, 0.0, 19.68237082607396, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tube1Display_3.OpacityTransferFunction.Points = [0.7022544507661337, 0.0, 0.5, 0.0, 19.68237082607396, 1.0, 0.5, 0.0]

# show color bar/color legend
tube1Display_3.SetScalarBarVisibility(renderView4, True)

# reset view to fit data
renderView4.ResetCamera()

# set scalar coloring
ColorBy(tube1Display_3, ('POINTS', 'velocity'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(radiusLUT, renderView4)

# rescale color and/or opacity maps used to include current data range
tube1Display_3.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
tube1Display_3.SetScalarBarVisibility(renderView4, True)

# get color transfer function/color map for 'velocity'
velocityLUT = GetColorTransferFunction('velocity')
velocityLUT.RGBPoints = [-114.03544090308917, 0.231373, 0.298039, 0.752941, 70.25035429821698, 0.865003, 0.865003, 0.865003, 254.5361494995231, 0.705882, 0.0156863, 0.14902]
velocityLUT.ScalarRangeInitialized = 1.0

# Rescale transfer function
velocityLUT.RescaleTransferFunction(-452.885300208, 1085.27785439)

# get opacity transfer function/opacity map for 'velocity'
velocityPWF = GetOpacityTransferFunction('velocity')
velocityPWF.Points = [-114.03544090308917, 0.0, 0.5, 0.0, 254.5361494995231, 1.0, 0.5, 0.0]
velocityPWF.ScalarRangeInitialized = 1

# Rescale transfer function
velocityPWF.RescaleTransferFunction(-452.885300208, 1085.27785439)

# set active view
SetActiveView(renderView1)

# Rescale transfer function
tube1Display.ScaleTransferFunction.RescaleTransferFunction(0.702254450766, 19.6823708261)

# Rescale transfer function
tube1Display.OpacityTransferFunction.RescaleTransferFunction(0.702254450766, 19.6823708261)

# set active view
SetActiveView(renderView3)

# set active view
SetActiveView(renderView1)

# set active view
SetActiveView(renderView2)

# Rescale transfer function
tube1Display_2.ScaleTransferFunction.RescaleTransferFunction(0.702254450766, 19.6823708261)

# Rescale transfer function
tube1Display_2.OpacityTransferFunction.RescaleTransferFunction(0.702254450766, 19.6823708261)

# set active view
SetActiveView(renderView1)

# set active view
SetActiveView(renderView4)

# set active view
SetActiveView(renderView1)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# set active view
SetActiveView(renderView3)

# set active view
SetActiveView(renderView1)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# set active view
SetActiveView(renderView3)

# set active view
SetActiveView(renderView2)

# set active view
SetActiveView(renderView3)

# set active view
SetActiveView(renderView1)

# set active view
SetActiveView(renderView3)

# reset view to fit data
renderView3.ResetCamera()

# set active view
SetActiveView(renderView1)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# current camera placement for renderView4
renderView4.CameraPosition = [255.59792065249823, 315.87271830090157, 1502.7612971210963]
renderView4.CameraFocalPoint = [255.59792065249823, 315.87271830090157, 0.0]
renderView4.CameraParallelScale = 399.79100714962914

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [255.59792065249823, 315.87271830090157, 1502.7612971210963]
renderView1.CameraFocalPoint = [255.59792065249823, 315.87271830090157, 0.0]
renderView1.CameraParallelScale = 399.79100714962914

# current camera placement for renderView3
renderView3.CameraPosition = [255.59792065249823, 315.87271830090157, 1502.7612971210963]
renderView3.CameraFocalPoint = [255.59792065249823, 315.87271830090157, 0.0]
renderView3.CameraParallelScale = 399.79100714962914

# current camera placement for renderView2
renderView2.CameraPosition = [255.59792065249823, 315.87271830090157, 1502.7612971210963]
renderView2.CameraFocalPoint = [255.59792065249823, 315.87271830090157, 0.0]
renderView2.CameraParallelScale = 399.79100714962914

# save animation
SaveAnimation('D:\gitRepositories\pulseWavePy\Adan_56.avi', layout1, SaveAllViews=1,
    ImageResolution=[1536, 856],
    SeparatorWidth=0,
    OverrideColorPalette='WhiteBackground',
    FrameRate=10,
    FrameWindow=[0, 39])

#### saving camera placements for all active views

# current camera placement for renderView4
renderView4.CameraPosition = [255.59792065249823, 315.87271830090157, 1502.7612971210963]
renderView4.CameraFocalPoint = [255.59792065249823, 315.87271830090157, 0.0]
renderView4.CameraParallelScale = 399.79100714962914

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [255.59792065249823, 315.87271830090157, 1502.7612971210963]
renderView1.CameraFocalPoint = [255.59792065249823, 315.87271830090157, 0.0]
renderView1.CameraParallelScale = 399.79100714962914

# current camera placement for renderView3
renderView3.CameraPosition = [255.59792065249823, 315.87271830090157, 1502.7612971210963]
renderView3.CameraFocalPoint = [255.59792065249823, 315.87271830090157, 0.0]
renderView3.CameraParallelScale = 399.79100714962914

# current camera placement for renderView2
renderView2.CameraPosition = [255.59792065249823, 315.87271830090157, 1502.7612971210963]
renderView2.CameraFocalPoint = [255.59792065249823, 315.87271830090157, 0.0]
renderView2.CameraParallelScale = 399.79100714962914
#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).