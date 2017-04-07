#Definition of inputs and outputs
#==================================
##BC=group
##PG04_WaterQualityParameters_CoastColour_L1P=name
##ParameterSelection|Icol|Icol|false; true
##ParameterSelection|Calibration|Calibration|false; true
##ParameterSelection|Smile|Smile|true; false
##ParameterSelection|Equalization|Equalization|true; false
##ParameterSelection|IgnoreSeaIceClim|Ignore sea ice climatology|false; true
##ParameterNumber|CloudBufferWidth|lnsert size of cloud buffer in pixel|0|None|2
##ParameterNumber|CloudScreeningAmbiguous|Cloud screening ambiguous threshold|0|None|1.4
##ParameterNumber|CloudScreeningSure|Cloud screening sure threshold|0|None|1.8
##ParameterNumber|GlintCloudThresholdAddition|Glint cloud screening addition|0|None|0.1
##ParameterSelection|OutputCloudProbabilityFeatureValue|Write cloud probability feature value|false; true
