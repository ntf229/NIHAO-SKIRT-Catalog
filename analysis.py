# Make plots from spatially integrated fits file

import numpy as np
import matplotlib.pyplot as plt
import fitsio
import matplotlib.patches as mpatches
import os

def dustPedia():
	# Load data
	dp_file = np.loadtxt('DustPediaData/DustPedia_Aperture_Photometry_2.2.csv',dtype=str)
	dp_file2 = np.loadtxt('DustPediaData/DustPedia_HyperLEDA_Herschel.csv',dtype=str) 
	dp_file3 = np.loadtxt('DustPediaData/dustpedia_cigale_results_final_version.csv',dtype=str)
	dp_num = len(dp_file) - 1 # first line is header
	dp_flux_index = np.asarray([7,10,13,16,19,22,25,28,31,34,37,40,43,46,70,73,76,79,82,85]) # DP indicies corresponding to band_names
	dp_err_index = dp_flux_index + 1 # index of DustPedia flux errors corresponding to band_names
	# Initialize DustPedia arrays
	dp_flux = np.zeros((dp_num, len(bands))) # first index specifies galaxy, second index specifies band 
	dp_err = np.zeros((dp_num, len(bands)))
	dp_bool = np.zeros((dp_num, len(bands)), dtype=bool) # True means flux consistent with 0 (to be colored red) 
	dp_axisRatio = np.zeros(dp_num)
	dp_disk = np.zeros(dp_num) # 0 means not disk, 1 means disk (t >= 0)
	dp_stellarMass = np.zeros(dp_num)
	dp_dustMass = np.zeros(dp_num)
	dp_SFR = np.zeros(dp_num)
	dp3_names = [] 
	for i in range(dp_num):
		g3 = i+1 # skip first line of file (headers)
		dp3_names.append(dp_file3[g3].split(',')[0])
	dp3_names = np.asarray(dp3_names) # doesn't include first line of dp3_file
	# Fill DustPedia arrays
	for i in range(dp_num):
		g = i+1 # index of current galaxy
		params = dp_file[g].split(',')
		params2 = dp_file2[g].split(',')
		g3 = np.where(dp3_names == params[0])[0][0] + 1 
		params3 = dp_file3[g3].split(',')
		t = float(params2[3])
		if t >= 0:
			dp_disk[i] = 1
		dp_axisRatio[i] = 1./float(params[4]) # change from a/b to b/a 
		dp_stellarMass[i] = params3[3] # in solar masses (CIGALE)
		dp_dustMass[i] = params3[17] # in solar masses (CIGALE)
		dp_SFR[i] = params3[1] # in solar masses per year (CIGALE)
		for j in range(len(bands)):
			if params[int(dp_flux_index[j])]:
				dp_flux[i,j] = float(params[int(dp_flux_index[j])])
				dp_err[i,j] = float(params[int(dp_err_index[j])])
				if dp_flux[i,j] - 2*dp_err[i,j] <= 0:
					dp_flux[i,j] = 2 * dp_err[i,j] # if flux is consistent with 0, set to 2*sigma 
					dp_bool[i,j] = True # signals flux is consistent with 0 (upper limit)
				else:
					dp_bool[i,j] = False
			else:
				dp_flux[i,j] = float("NaN") # set flux and error to NaN if photometry not available 
				dp_err[i,j] = float("NaN")
				dp_bool[i,j] = True 
	diskMask = dp_disk == 1
	faceOnMask = dp_axisRatio > 0.85 
	return dp_flux, dp_err, dp_bool, dp_axisRatio, dp_stellarMass, dp_dustMass, dp_SFR, diskMask, faceOnMask

def dustPediaPlots():
	plt.figure(figsize=(10,8))
	massMask = np.log10(dp_stellarMass_og) > 7.
	plt.scatter(np.log10(dp_stellarMass_og[massMask]), np.log10(dp_SFR_og[massMask]/dp_stellarMass_og[massMask]), 
				color='k', alpha=0.5, s=15, linewidth=0)
	# before mass cut
	left = np.amin(np.log10(stellarMass))
	bottom =  np.amin(np.log10(sSFR))
	width = np.amax(np.log10(stellarMass)) - left
	height =  np.amax(np.log10(sSFR)) - bottom
	rect = mpatches.Rectangle((left, bottom), width, height, fill=False, color="red", linewidth=2)
	plt.gca().add_patch(rect)
	# after mass cut
	left = np.amin(np.log10(stellarMass[sphMassMask]))
	bottom =  np.amin(np.log10(sSFR[sphMassMask]))
	width = np.amax(np.log10(stellarMass[sphMassMask])) - left
	height =  np.amax(np.log10(sSFR[sphMassMask])) - bottom
	rect = mpatches.Rectangle((left, bottom), width, height, fill=False, color="k", linewidth=2)
	plt.gca().add_patch(rect)
	plt.xlabel(r'$\log_{10}(Stellar \; Mass \, / \, M_{\odot})$', fontsize=28)
	plt.ylabel(r'$\log_{10}(sSFR \, / \, yr^{-1})$',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/dp_sSFR_stellarMass.png', dpi=300, bbox_inches='tight', pad_inches=0.25)
	plt.close()

def plotAvAxisRatio():
	plt.figure(figsize=(10,8))
	for i in range(len(axisRatios)):
		if ~sphMassMask[i]:
			plt.scatter(axisRatios[i], Av[i], color='red', linewidth=0)
	for i in range(len(axisRatios)): # separate loop so black points always on top
		if sphMassMask[i]:
			plt.scatter(axisRatios[i], Av[i], color='k', linewidth=0)
	plt.xlabel('Axis Ratio', fontsize=16)
	plt.ylabel(r'$A_{V}$',fontsize=16)
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig('Plots/AvAxisRatio.png', dpi=300)
	plt.close()

def plotAvStellarMass():
	plt.figure(figsize=(10,8))
	for i in range(len(axisRatios)):
		if ~sphMassMask[i]:
			plt.scatter(np.log10(stellarMass[i]), Av[i], color='red', linewidth=0)
		else:
			plt.scatter(np.log10(stellarMass[i]), Av[i], color='k', linewidth=0)
	plt.xlabel(r'$\log_{10}(Stellar \; Mass \, / \, M_{\odot})$', fontsize=24)
	plt.ylabel(r'$A_{V}$',fontsize=24)
	plt.xticks(fontsize=24)
	plt.yticks(fontsize=24)
	plt.savefig('Plots/AvStellarMass.png', dpi=300, bbox_inches='tight', pad_inches=0.25)
	plt.close()

def plotAttenuation():
	plt.figure(figsize=(10,8))
	for i in range(len(axisRatios)):
		if sphMassMask[i]:
			if AvMask[i]:
				plt.plot(np.log10(attenuation_wave), attenuation_mags[i], 
						 color='k', alpha=0.5, linewidth=0.1)
			else:
				plt.plot(np.log10(attenuation_wave), attenuation_mags[i], 
						 color='blue', alpha=0.5, linewidth=0.1)
		else:
			plt.plot(np.log10(attenuation_wave), attenuation_mags[i], 
					 color='red', alpha=0.5, linewidth=0.1)
	plt.xlabel(r'$\log_{10}(\lambda \, / \, \AA)$', fontsize=28)
	plt.ylabel(r'$A_{\lambda}$',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/attenuation.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def plotAttenuationNoColors():
	plt.figure(figsize=(10,8))
	for i in range(len(axisRatios)):
		plt.plot(np.log10(attenuation_wave), attenuation_mags[i], color='k', alpha=0.5, linewidth=0.1)
	plt.xlabel('log('+r'$\lambda$'+' / '+r'$\AA$'+')', fontsize=24)
	plt.ylabel(r'$A_{\lambda}$',fontsize=24)
	plt.xticks(fontsize=24)
	plt.yticks(fontsize=24)
	plt.savefig('Plots/attenuationNoColors.png', dpi=600)
	plt.close()

def plotAttenuationNorm():
	plt.figure(figsize=(10,8))
	for i in range(len(axisRatios)):
		if sphMassMask[i] & AvMask[i]:
			plt.plot(np.log10(attenuation_wave), attenuation_mags[i]/Av[i], 
					 color='k', alpha=0.5, linewidth=0.1)
	plt.xlabel(r'$\log_{10}(\lambda \, / \, \AA)$', fontsize=28)
	plt.ylabel(r'$A_{\lambda} \, / \, A_{V}$',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/attenuation_norm.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def plotAttenuationNormSeparate(galaxy):
	os.system('mkdir -p Plots/SeparateAttenuationNorm/')
	nameMask = names == galaxy
	plt.figure(figsize=(10,8))
	for i in range(10):
		if axisRatios[nameMask][i] == np.amin(axisRatios[nameMask]):
			color = 'red'
			linewidth = 1.5
		elif axisRatios[nameMask][i] == np.amax(axisRatios[nameMask]):
			color = 'blue'
			linewidth = 1.5
		else:
			color = 'k'
			linewidth = 0.3
		plt.plot(np.log10(attenuation_wave), attenuation_mags[nameMask][i]/Av[nameMask][i], 
				 color=color, alpha=1, linewidth=linewidth)
	plt.xlabel(r'$\log_{10}(\lambda \, / \, \AA)$', fontsize=28)
	plt.ylabel(r'$A_{\lambda} \, / \, A_{V}$',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/SeparateAttenuationNorm/'+galaxy+'.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def plotEnergyBalance():
	plt.figure(figsize=(10,8))
	for i in range(len(axisRatios)):
		if sphMassMask[i]:
			plt.scatter(axisRatios[i], np.log10(attEnergy[i]/emitEnergy[i]), color='k', linewidth=0)
	plt.axhline(y=0, color='k')
	plt.xlim((0,1))
	plt.xlabel(r'$Axis \; Ratio$', fontsize=28)
	plt.ylabel(r'$\log_{10}(Attenuated \, / \, Emitted \; Energy)$',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/energy_balance.png', dpi=300, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def plotSpectra(galaxy):
	# Plot edge-on, face-on, and no-dust spectra
	os.system('mkdir -p Plots/Spectra/')
	plt.figure(figsize=(10,8))
	waveMask = (wave >= 1e3) & (wave <= 1e7)
	nameMask = names == galaxy
	edgeIndex = np.argmin(axisRatios[nameMask])
	faceIndex = np.argmax(axisRatios[nameMask])
	plt.plot(np.log10(wave[waveMask]), np.log10(spectrum[nameMask][edgeIndex][waveMask]), 
			 color='red', label='edge-on', linewidth=1.5)
	plt.plot(np.log10(wave[waveMask]), np.log10(spectrum[nameMask][faceIndex][waveMask]), 
			 color='blue', label='face-on', linewidth=1.5)
	plt.plot(np.log10(wave[waveMask]), np.log10(spectrum_nodust[nameMask][faceIndex][waveMask]), 
			 color='k', label='no-dust', linewidth=1.5)
	plt.legend(fontsize=28)
	plt.xlabel(r'$\log_{10}(\lambda \, / \, \AA)$', fontsize=28)
	plt.ylabel(r'$\log_{10}(f_{\nu} \, / \, Jy)$',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/Spectra/'+galaxy+'_SEDs.png', dpi=300, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def plotSpectraNoIR(galaxy):
	# Make plots of spectra zoomed in to attenuation region with all axis ratios
	os.system('mkdir -p Plots/SpectraNoIR/')
	plt.figure(figsize=(10,8))
	waveMask = (wave >= 1e3) & (wave <= 2e4)
	nameMask = names == galaxy
	plt.plot(np.log10(wave[waveMask]), np.log10(spectrum_nodust[nameMask][0][waveMask]), 
			 color='k', label='no-dust', linewidth=1.5)
	for i in range(10):
		if axisRatios[nameMask][i] == np.amin(axisRatios[nameMask]):
			color = 'red'
			linewidth = 1.5
			label = 'edge-on'
		elif axisRatios[nameMask][i] == np.amax(axisRatios[nameMask]):
			color = 'blue'
			linewidth = 1.5
			label = 'face-on'
		else:
			color = 'k'
			linewidth = 0.3
			label = None
		plt.plot(np.log10(wave[waveMask]), np.log10(spectrum[nameMask][i][waveMask]), 
				 color=color, label=label, linewidth=linewidth)
	plt.legend(fontsize=28)
	plt.xlabel(r'$\log_{10}(\lambda \, / \, \AA)$', fontsize=28)
	plt.ylabel(r'$\log_{10}(f_{\nu} \, / \, Jy)$',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/SpectraNoIR/'+galaxy+'_SEDsNoIR.png', dpi=600, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def plotSpectraSeparate(galaxy):
	# Make separate plots of edge-on and face-on spectra
	os.system('mkdir -p Plots/SpectraSeparate/')
	# edge-on
	plt.figure(figsize=(10,8))
	waveMask = (wave >= 1e3) & (wave <= 1e7)
	nameMask = names == galaxy
	edgeIndex = np.argmin(axisRatios[nameMask])
	faceIndex = np.argmax(axisRatios[nameMask])
	plt.plot(np.log10(wave[waveMask]), np.log10(spectrum[nameMask][edgeIndex][waveMask]), color='k')
	plt.xlabel('log('+r'$\lambda$'+' / '+r'$\AA$'+')', fontsize=28)
	plt.ylabel('log('+r'$f_{\nu}$'+' / Janskys)',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	ylim = plt.gca().get_ylim()
	plt.savefig('Plots/SpectraSeparate/'+galaxy+'_edge_SEDs.png', dpi=300, bbox_inches='tight', pad_inches=0.5)
	plt.close()
	# face-on
	plt.figure(figsize=(10,8))
	plt.plot(np.log10(wave[waveMask]), np.log10(spectrum[nameMask][faceIndex][waveMask]), color='k')	
	plt.xlabel('log('+r'$\lambda$'+' / '+r'$\AA$'+')', fontsize=28)
	plt.ylabel('log('+r'$f_{\nu}$'+' / Janskys)',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.gca().set_ylim(ylim)
	plt.savefig('Plots/SpectraSeparate/'+galaxy+'_face_SEDs.png', dpi=300, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def plotSpectraSeparateAllARs(galaxy):
	# Make plots of spectra with all axis ratios
	nameMask = names == galaxy
	waveMask = (wave >= 1e3) & (wave <= 1e7)
	yMin = []
	yMax = []
	for i in range(10):
		yMin.append(np.amin(np.log10(spectrum[nameMask][i][waveMask])))
		yMax.append(np.amax(np.log10(spectrum[nameMask][i][waveMask])))
	padding = (np.amax(yMax) - np.amin(yMin))*0.05
	yRange = [np.amin(yMin)-padding,np.amax(yMax)+padding]
	for i in range(10):
		os.system('mkdir -p Plots/SpectraAllAxisRatios/'+galaxy+'_SEDs/')
		plt.figure(figsize=(10,8))
		plt.plot(np.log10(wave[waveMask]), np.log10(spectrum[nameMask][i][waveMask]), color='k')
		plt.xlabel('log('+r'$\lambda$'+' / '+r'$\AA$'+')', fontsize=28)
		plt.ylabel('log('+r'$f_{\nu}$'+' / Janskys)',fontsize=28)
		plt.xticks(fontsize=28)
		plt.yticks(fontsize=28)
		plt.gca().set_ylim(yRange)
		plt.savefig('Plots/SpectraAllAxisRatios/'+galaxy+'_SEDs/axisRatio'+
					str(axisRatios[nameMask][i])+'_SED.png', dpi=300, bbox_inches='tight', pad_inches=0.5)
		plt.close()

def colorColorPlots():
	os.system('mkdir -p Plots/colorPlots/')
	plt.figure(figsize=(10,8))
	dpLabelBool = True # DustPedia label (triggers only once)
	# DustPedia propogation of errors (log scale)
	xflux1 = dp_flux[:,ratio_indicies[i][0]]
	xflux2 = dp_flux[:,ratio_indicies[i][1]]
	yflux1 = dp_flux[:,ratio_indicies[i][2]]
	yflux2 = dp_flux[:,ratio_indicies[i][3]]
	xerr1 = dp_err[:,ratio_indicies[i][0]]
	xerr2 = dp_err[:,ratio_indicies[i][1]]
	yerr1 = dp_err[:,ratio_indicies[i][2]]
	yerr2 = dp_err[:,ratio_indicies[i][3]]
	x = xflux1 / xflux2
	y = yflux1 / yflux2
	xerr = x * np.sqrt((xerr1/xflux1)**2 + (xerr2/xflux2)**2)
	yerr = y * np.sqrt((yerr1/yflux1)**2 + (yerr2/yflux2)**2)
	# put on log scale
	xerrs = np.asarray([np.log10(x) - np.log10(x - xerr), np.log10(x + xerr) - np.log10(x)])
	yerrs = np.asarray([np.log10(y) - np.log10(y - yerr), np.log10(y + yerr) - np.log10(y)])
	colors = np.empty(len(dp_flux[:,0]),dtype=str)
	for d in range(len(dp_flux[:,0])):
		# plot arguments
		bad = False
		xuplims = False
		xlolims = False
		uplims = False
		lolims = False
		colors[d] = 'k'
		if dp_bool[d,ratio_indicies[i][0]] and dp_bool[d,ratio_indicies[i][1]]: # bad data point (0/0)
			bad = True # don't plot this point
		elif dp_bool[d,ratio_indicies[i][0]]: # x is an upper limit (numerator is upper limit)
			xuplims = True
			colors[d] = 'red'
		elif dp_bool[d,ratio_indicies[i][1]]: # x is a lower limit (denominator is upper limit)
			xlolims = True
			colors[d] = 'red'
		if dp_bool[d,ratio_indicies[i][2]] and dp_bool[d,ratio_indicies[i][3]]: # bad data point (0/0)
			bad = True # don't plot this point
		elif dp_bool[d,ratio_indicies[i][2]]: # y is an upper limit (numerator is upper limit)
			uplims = True
			colors[d] = 'red'
		elif dp_bool[d,ratio_indicies[i][3]]: # y is a lower limit (denominator is upper limit)
			lolims = True
			colors[d] = 'red'
		if ~bad:
			if (colors[d] == 'k') & dpLabelBool:
				dpLabel = 'DustPedia'
				dpLabelBool = False
			else:
				dpLabel = None
			plt.errorbar(np.log10(x[d]), np.log10(y[d]), 
						 xerr=xerrs[:,[d]], yerr=yerrs[:,[d]], elinewidth=1, marker='o',
						 markersize=12, linewidth=0, color=colors[d], zorder=0, alpha=0.2,
						 xuplims=xuplims, xlolims=xlolims, uplims=uplims, lolims=lolims, 
						 label=dpLabel, markeredgewidth=0.0)
	plt.scatter(np.log10(flux[sphMassMask,ratio_indicies[i][0]] / flux[sphMassMask,ratio_indicies[i][1]]), 
				np.log10(flux[sphMassMask,ratio_indicies[i][2]] / flux[sphMassMask,ratio_indicies[i][3]]), 
				marker='D', s=120, zorder=10, alpha=0.2, c='green', label='NIHAO', linewidth=0)
	plt.xlabel(r'$\log_{10}[f_{\nu}$('+band_names[ratio_indicies[i][0]]+') / '+
			   r'$f_{\nu}($'+band_names[ratio_indicies[i][1]]+r'$)]$', fontsize=28)
	plt.ylabel(r'$\log_{10}[f_{\nu}$('+band_names[ratio_indicies[i][2]]+') / '+
			   r'$f_{\nu}($'+band_names[ratio_indicies[i][3]]+r'$)]$', fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.legend(fontsize=28)
	plt.savefig('Plots/colorPlots/'+plot_names[i]+'.png',dpi=300, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def colorColorPlotsMassCut():
	# compare color-color plots before and after mass cut
	os.system('mkdir -p Plots/colorPlotsMassCut/')
	fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)
	fig.subplots_adjust(wspace=0)
	dpLabelBool = True # DustPedia label (triggers only once)
	# DustPedia propogation of errors (log scale)
	xflux1 = dp_flux_og[:,ratio_indicies[i][0]]
	xflux2 = dp_flux_og[:,ratio_indicies[i][1]]
	yflux1 = dp_flux_og[:,ratio_indicies[i][2]]
	yflux2 = dp_flux_og[:,ratio_indicies[i][3]]
	xerr1 = dp_err_og[:,ratio_indicies[i][0]]
	xerr2 = dp_err_og[:,ratio_indicies[i][1]]
	yerr1 = dp_err_og[:,ratio_indicies[i][2]]
	yerr2 = dp_err_og[:,ratio_indicies[i][3]]
	x = xflux1 / xflux2
	y = yflux1 / yflux2
	xerr = x * np.sqrt((xerr1/xflux1)**2 + (xerr2/xflux2)**2)
	yerr = y * np.sqrt((yerr1/yflux1)**2 + (yerr2/yflux2)**2)
	# put on log scale
	xerrs = np.asarray([np.log10(x) - np.log10(x - xerr), np.log10(x + xerr) - np.log10(x)])
	yerrs = np.asarray([np.log10(y) - np.log10(y - yerr), np.log10(y + yerr) - np.log10(y)])
	colors = np.empty(len(dp_flux_og[:,0]),dtype=str)
	for d in range(len(dp_flux_og[:,0])):
		bad = False
		xuplims = False
		xlolims = False
		uplims = False
		lolims = False
		colors[d] = 'k'
		if dp_bool_og[d,ratio_indicies[i][0]] and dp_bool_og[d,ratio_indicies[i][1]]: # bad data point (0/0)
			bad = True # don't plot this point
		elif dp_bool_og[d,ratio_indicies[i][0]]: # x is an upper limit (numerator is upper limit)
			xuplims = True
			colors[d] = 'red'
		elif dp_bool_og[d,ratio_indicies[i][1]]: # x is a lower limit (denominator is upper limit)
			xlolims = True
			colors[d] = 'red'
		if dp_bool_og[d,ratio_indicies[i][2]] and dp_bool_og[d,ratio_indicies[i][3]]: # bad data point (0/0)
			bad = True # don't plot this point
		elif dp_bool_og[d,ratio_indicies[i][2]]: # y is an upper limit (numerator is upper limit)
			uplims = True
			colors[d] = 'red'
		elif dp_bool_og[d,ratio_indicies[i][3]]: # y is a lower limit (denominator is upper limit)
			lolims = True
			colors[d] = 'red'
		if ~bad:
			# before sph mass cut
			if dpMassMask_og[d] & dpSFRMask_separate_og[d]:
				ax1.errorbar(np.log10(x[d]), np.log10(y[d]), 
							 xerr=xerrs[:,[d]], yerr=yerrs[:,[d]], elinewidth=0.5, marker='o',
							 markersize=6, linewidth=0, color=colors[d], zorder=0, alpha=0.1,
							 xuplims=xuplims, xlolims=xlolims, uplims=uplims, lolims=lolims, markeredgewidth=0.0)
			# after sph mass cut
			if dpMassMask[d] & dpSFRMask_separate[d]:
				if (colors[d] == 'k') & dpLabelBool:
					dpLabel = 'DustPedia'
					dpLabelBool = False
				else:
					dpLabel = None
				ax2.errorbar(np.log10(x[d]), np.log10(y[d]), 
							 xerr=xerrs[:,[d]], yerr=yerrs[:,[d]], elinewidth=0.5, marker='o',
							 markersize=6, linewidth=0, color=colors[d], zorder=0, alpha=0.1,
							 xuplims=xuplims, xlolims=xlolims, uplims=uplims, lolims=lolims, 
							 label=dpLabel, markeredgewidth=0.0)
	ax1.scatter(np.log10(flux[:,ratio_indicies[i][0]] / flux[:,ratio_indicies[i][1]]), 
				np.log10(flux[:,ratio_indicies[i][2]] / flux[:,ratio_indicies[i][3]]), 
				marker='D', s=40, zorder=10, alpha=0.1, c='green', linewidth=0)
	ax2.scatter(np.log10(flux[sphMassMask,ratio_indicies[i][0]] / flux[sphMassMask,ratio_indicies[i][1]]), 
				np.log10(flux[sphMassMask,ratio_indicies[i][2]] / flux[sphMassMask,ratio_indicies[i][3]]), 
				marker='D', s=40, zorder=10, alpha=0.1, c='green', label='NIHAO', linewidth=0)
	ax1.tick_params(axis='both', labelsize=14, left=True, bottom=True)
	ax2.tick_params(axis='both', labelsize=14, left=True, bottom=True)
	ax1.set_xticks([-3, -2, -1])
	ax2.set_xticks([-3, -2, -1])
	ax1.set_title('Before Mass Cut', fontsize=14)
	ax2.set_title('After Mass Cut', fontsize=14)
	ax2.legend(fontsize=14)
	ax = fig.add_subplot(111, frameon=False)
	ax.tick_params('both', labelbottom=False, labelleft=False, left=False, bottom=False)
	ax.set_xlabel(r'$\log_{10}[f_{\nu}$('+band_names[ratio_indicies[i][0]]+') / '+
				  r'$f_{\nu}($'+band_names[ratio_indicies[i][1]]+r'$)]$', fontsize=14)
	ax.set_ylabel(r'$\log_{10}[f_{\nu}$('+band_names[ratio_indicies[i][2]]+') / '+
				  r'$f_{\nu}($'+band_names[ratio_indicies[i][3]]+r'$)]$', fontsize=14)
	ax.xaxis.set_label_coords(0.5, -0.09)
	ax.yaxis.set_label_coords(-0.12, 0.5)
	plt.savefig('Plots/colorPlotsMassCut/'+plot_names[i]+'.png',dpi=300, bbox_inches='tight', pad_inches=0.25)
	plt.close()

def plotSFH(galaxy):
	os.system('mkdir -p Plots/SFHs/')
	plt.figure(figsize=(10,8))
	nameMask = singleNames == galaxy
	bins = ages[nameMask, :][0]
	weights = SFH[nameMask, :][0]
	bins /= 1e9 # convert to Gyrs
	plt.hist(bins, bins=bins, weights=weights)
	plt.xlabel('Age [Gyrs]', fontsize=28)
	plt.ylabel(r'$Mass \, [M_{\odot}]$',fontsize=28)
	plt.yscale('log')
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/SFHs/'+galaxy+'_SFH.png', dpi=300, bbox_inches='tight', pad_inches=0.5)
	plt.close()

def plotCEH(galaxy):
	os.system('mkdir -p Plots/CEHs/')
	plt.figure(figsize=(10,8))
	nameMask = singleNames == galaxy
	xbins = ages[nameMask, :][0]
	ybins = metals[nameMask, :][0]
	ageBins, metalBins = np.meshgrid(xbins, ybins)
	weights = CEH[nameMask, :, :][0]
	ageBins /= 1e9 # convert to Gyrs
	plt.hist2d(ageBins.flatten(), metalBins.flatten(), bins=[ageBins[0,:],metalBins[:,0]], 
			   weights=weights.flatten(), cmap='Greys')
	plt.xlabel(r'$Age \, [Gyrs]$', fontsize=28)
	plt.ylabel(r'$Metallicity$',fontsize=28)
	plt.xticks(fontsize=28)
	plt.yticks(fontsize=28)
	plt.savefig('Plots/CEHs/'+galaxy+'_CEH.png', dpi=300, bbox_inches='tight', pad_inches=0.5)
	plt.close()

# Load fits file
galaxies = fitsio.read('nihao-integrated-seds.fits', ext='GALAXIES') # one per galaxy
summary = fitsio.read('nihao-integrated-seds.fits', ext='SUMMARY') # 10 per galaxy (one per viewing orientation)
wave = fitsio.read('nihao-integrated-seds.fits', ext='WAVE')
spectrum = fitsio.read('nihao-integrated-seds.fits', ext='SPEC')
spectrum_nodust = fitsio.read('nihao-integrated-seds.fits', ext='SPECNODUST')
attenuation_wave = fitsio.read('nihao-integrated-seds.fits', ext='ATTENUATION_WAVE')
attenuation_mags = fitsio.read('nihao-integrated-seds.fits', ext='ATTENUATION_MAGS')

singleNames = galaxies['name'] # one per galaxy
singleStellarMass = galaxies['stellar_mass'] # one per galaxy, in M_sun
names = summary['name']
stellarMass = summary['stellar_mass'] # in solar mass
SFR = summary['sfr'] # in M_sun per year
dustMass = summary['dust_mass'] # in M_sun
axisRatios = summary['axis_ratio']
Av = summary['Av'] # in magnitudes
bands = summary['bands'][0]
flux = summary['flux'] # in Jansky
flux_noDust = summary['flux_nodust'] # in Jansky
attEnergy = summary['attenuated_energy'] # in 10^-23 erg s^-1 cm^-2
emitEnergy = summary['emitted_energy'] # in 10^-23 erg s^-1 cm^-2
SFH = galaxies['sfh'] # in M_sun (as a function of ages)
CEH = galaxies['ceh'] # in M_sun  (as a function of ages and metals)
ages = galaxies['ages'] # age bins in yeras for SFH and CEH
metals = galaxies['metals'] # metallicity bins for CEH

# Broadband photometric bandpass names
band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', 'J', 'H', 'K', 'W1', 
              'W2', 'W3', 'W4', r'$70 \mu m$', r'$100 \mu m$', r'$160 \mu m$', 
              r'$250 \mu m$', r'$350 \mu m$', r'$500 \mu m$']

sSFR = SFR / stellarMass

# create masks
sphMassMask = np.log10(stellarMass) > 9.5 
AvMask = Av > 0.2
singleSphMassMask = np.log10(singleStellarMass) > 9.5 
minAvMask = [] # only include galaxies with minimum Av > 0.2 
for i in range(len(singleNames)):
	nameMask = names == singleNames[i]
	if (np.amin(Av[nameMask]) > 0.2):
		minAvMask.append(True)
	else:
		minAvMask.append(False)

# Original DustPedia data, not matched to NIHAO mass and sSFR ranges
dp_flux_og, dp_err_og, dp_bool_og, dp_axisRatio_og, dp_stellarMass_og, dp_dustMass_og, dp_SFR_og, dp_diskMask_og, dp_faceOnMask_og = dustPedia()
dp_sSFR_og = dp_SFR_og / dp_stellarMass_og

# Match DustPedia stellar mass and sSFR to simulations (before mass cut)
dpMassMask_og = (dp_stellarMass_og >= np.amin(stellarMass)) & (dp_stellarMass_og <= np.amax(stellarMass))
dpSFRMask_og = (dp_SFR_og[dpMassMask_og]/dp_stellarMass_og[dpMassMask_og] >= np.amin(sSFR)) & (dp_SFR_og[dpMassMask_og]/dp_stellarMass_og[dpMassMask_og] <= np.amax(sSFR))
dpSFRMask_separate_og = (dp_SFR_og/dp_stellarMass_og >= np.amin(sSFR)) & (dp_SFR_og/dp_stellarMass_og <= np.amax(sSFR))

# Match DustPedia stellar mass and sSFR to simulations (after mass cut)
dpMassMask = (dp_stellarMass_og >= np.amin(stellarMass[sphMassMask])) & (dp_stellarMass_og <= np.amax(stellarMass[sphMassMask]))
dpSFRMask = (dp_SFR_og[dpMassMask]/dp_stellarMass_og[dpMassMask] >= np.amin(sSFR[sphMassMask])) & (dp_SFR_og[dpMassMask]/dp_stellarMass_og[dpMassMask] <= np.amax(sSFR[sphMassMask]))
dpSFRMask_separate = (dp_SFR_og/dp_stellarMass_og >= np.amin(sSFR[sphMassMask])) & (dp_SFR_og/dp_stellarMass_og <= np.amax(sSFR[sphMassMask]))
dp_diskMask = dp_diskMask_og[dpMassMask][dpSFRMask]
dp_faceOnMask = dp_faceOnMask_og[dpMassMask][dpSFRMask][dp_diskMask]
dp_flux = dp_flux_og[dpMassMask][dpSFRMask]
dp_err = dp_err_og[dpMassMask][dpSFRMask]
dp_bool = dp_bool_og[dpMassMask][dpSFRMask]
dp_axisRatio = dp_axisRatio_og[dpMassMask][dpSFRMask]
dp_stellarMass = dp_stellarMass_og[dpMassMask][dpSFRMask]
dp_dustMass = dp_dustMass_og[dpMassMask][dpSFRMask]
dp_SFR = dp_SFR_og[dpMassMask][dpSFRMask]
dp_sSFR = dp_SFR / dp_stellarMass

os.system('mkdir -p Plots/')

# Make plots 
dustPediaPlots()
plotAvAxisRatio()
plotAvStellarMass()
plotAttenuation()
plotAttenuationNoColors()
plotAttenuationNorm()
plotEnergyBalance()
for i in range(len(singleNames)): # loop over all galaxies
	plotSpectra(singleNames[i])
	plotSpectraNoIR(singleNames[i])
	plotSpectraSeparate(singleNames[i])
	plotSpectraSeparateAllARs(singleNames[i])
	plotSFH(singleNames[i])
	plotCEH(singleNames[i])
for i in range(len(singleNames[minAvMask])): # only loop over galaxies with minimum Av > 0.2
	plotAttenuationNormSeparate(singleNames[minAvMask][i])

# color-color plots (0 = FUV, 6 = z, 7-19 given by plot_names)
ratio_indicies = [[0,6,7,6], [0,6,8,6], [0,6,9,6], [0,6,10,6], [0,6,11,6], [0,6,12,6], [0,6,13,6], 
				  [0,6,14,6], [0,6,15,6], [0,6,16,6], [0,6,17,6], [0,6,18,6], [0,6,19,6]]
plot_names = ["J","H","K","W1","W2","W3","W4","pacs70","pacs100","pacs160","spire250","spire350","spire500"]
for i in range(len(ratio_indicies)): 
	colorColorPlots()
	colorColorPlotsMassCut()




