# boundary_layer_classification_toolbox
Boundary layer classification toolbox

This toolbox includes all the functions that are required for boundary layer classfication.


Script_boundary_layer_classification.m
- preprocessHalo.m
	- loadHaloVert.m
	- correctHaloRipples.m
		- calculateBKG_v2.m
	- correctBackground.m
		- my_robustfit.m
	- correct_focus.m
	- calculate_dl_SNR.m
	- write_nc_silent.m
- combineHaloDBSnVAD.m
	- loadHaloWinds.m
	- medianfilter.m
- loadModel.m
- LLJ_detection.m
- calcWindQuantities.m
	- calcTKE.m
		- medianfilter.m
	- windowSlider.m
		- my_skewness.m
- decimal2daten.m
- suncycle.m % OR pvl_ephemeris.m
- createBitfield.m
- create_BL_categorization.m
- write_nc_silent.m
