# boundary_layer_classification_toolbox
Boundary layer classification toolbox

This toolbox includes all the functions that are required for boundary layer classfication.


Script_boundary_layer_classification.m
- preprocessHalo.m
	- loadHaloVert.m
		- jd2date.m
	- correctHaloRipples.m
		- calculateBKG.m
	- correctBackground.m
		- my_robustfit.m
			- statrobustwfun.m
			- statremovenan.m
			- statrobustfit.m
	- correct_focus.m
	- calculate_dl_SNR.m
	- write_nc_silent.m
- combineHaloDBSnVAD.m
	- loadHaloWinds.m
	- medianfilter.m
- loadModel.m
- detectLLJ.m
- calcWindQuantities.m
	- calcTKE.m
		- medianfilter.m
	- windowSlider.m
		- my_skewness.m
- decimal2daten.m
- TKE_connected_with.m
- suncycle.m
- createBitfield.m
- create_BL_categorization.m
- write_nc_silent.m
