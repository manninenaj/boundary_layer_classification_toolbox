# boundary_layer_classification_toolbox
Boundary layer classification toolbox

This toolbox includes all the functions that are required for boundary layer classfication.


Script_BL_classification
- preprocessHalo
	- loadHaloVert
	- correctHaloRipples
		- calculateBKG_v2
	- correctBackground
		- my_robustfit
	- correct_focus
	- calculate_dl_SNR
	- write_nc_silent
- combineHaloDBSnVAD
	- loadHaloWinds
	- medianfilter
- loadModel
- LLJ_detection
- calcWindQuantities
	- calcTKE
		- medianfilter
	- windowSlider
		- my_skewness
- decimal2daten
- suncycle % OR pvl_ephemeris
- createBitfield
- create_BL_categorization
- write_nc_silent
