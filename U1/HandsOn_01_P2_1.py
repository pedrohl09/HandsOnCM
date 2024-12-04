# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 19:54:21 2024

@author: pedro
"""

import numpy as np
import matplotlib.pyplot as plt

# Input parameters
dR = 150                          # Radius of the hexagon
dShad = 50                        # Shadowing decorrelation distance
dPasso = 7                        # Distance between measurement points
dSigmaShad = 8                    # Standard deviation of log-normal shadowing

# Calculation of other variables depending on input parameters
dDimXOri = 5 * dR                 # X dimension of the grid
dDimYOri = 6 * np.sqrt(3 / 4) * dR  # Y dimension of the grid

# Reference matrix with the position of each grid point (relative to bottom-left corner)
dDimY = int(np.ceil(dDimYOri + (dDimYOri % dPasso)))   # Adjust Y dimension to cover the entire grid
dDimX = int(np.ceil(dDimXOri + (dDimXOri % dPasso)))   # Adjust X dimension to cover the entire grid
mtPosx, mtPosy = np.meshgrid(np.arange(0, dDimX + 1, dPasso), np.arange(0, dDimY + 1, dPasso))
mtPontosMedicao = mtPosx + 1j * mtPosy

# Equidistant points matrix with dShad spacing
dDimYS = int(np.ceil(dDimYOri + (dDimYOri % dShad)))
dDimXS = int(np.ceil(dDimXOri + (dDimXOri % dShad)))
mtPosxShad, mtPosyShad = np.meshgrid(np.arange(0, dDimXS + 1, dShad), np.arange(0, dDimYS + 1, dShad))
mtPosShad = mtPosxShad + 1j * mtPosyShad

# Shadowing samples for grid points
mtShadowingSamples = dSigmaShad * np.random.randn(*mtPosyShad.shape)

# Initialization of the shadowing correction matrix
mtShadowingCorr = np.zeros(mtPosx.shape)

# Main calculation loop
for il in range(mtPosx.shape[0]):
    for ic in range(mtPosx.shape[1]):
        dshadPoint = mtPontosMedicao[il, ic]

        dXIndexP1 = np.real(dshadPoint) / dShad
        dYIndexP1 = np.imag(dshadPoint) / dShad

        if dXIndexP1.is_integer() and dYIndexP1.is_integer():
            dXIndexP1 = int(np.floor(dXIndexP1))
            dYIndexP1 = int(np.floor(dYIndexP1))
            plt.plot(np.real(mtPosShad[dYIndexP1, dXIndexP1]), np.imag(mtPosShad[dYIndexP1, dXIndexP1]), 'g*')
            mtShadowingCorr[il, ic] = mtShadowingSamples[dYIndexP1, dXIndexP1]
        else:
            dXIndexP1 = int(np.floor(dXIndexP1))
            dYIndexP1 = int(np.floor(dYIndexP1))

            if dXIndexP1 == mtPosyShad.shape[1] - 1 and dYIndexP1 == mtPosyShad.shape[0] - 1:
                dXIndexP2, dYIndexP2 = dXIndexP1 - 1, dYIndexP1
                dXIndexP4, dYIndexP4 = dXIndexP1 - 1, dYIndexP1 - 1
                dXIndexP3, dYIndexP3 = dXIndexP1, dYIndexP1 - 1
            elif dXIndexP1 == mtPosyShad.shape[1] - 1:
                dXIndexP2, dYIndexP2 = dXIndexP1 - 1, dYIndexP1
                dXIndexP4, dYIndexP4 = dXIndexP1 - 1, dYIndexP1 + 1
                dXIndexP3, dYIndexP3 = dXIndexP1, dYIndexP1 + 1
            elif dYIndexP1 == mtPosyShad.shape[0] - 1:
                dXIndexP2, dYIndexP2 = dXIndexP1 + 1, dYIndexP1
                dXIndexP4, dYIndexP4 = dXIndexP1 + 1, dYIndexP1 - 1
                dXIndexP3, dYIndexP3 = dXIndexP1, dYIndexP1 - 1
            else:
                dXIndexP2, dYIndexP2 = dXIndexP1 + 1, dYIndexP1
                dXIndexP4, dYIndexP4 = dXIndexP1 + 1, dYIndexP1 + 1
                dXIndexP3, dYIndexP3 = dXIndexP1, dYIndexP1 + 1

            dDistX = (np.mod(np.real(dshadPoint), dShad)) / dShad
            dDistY = (np.mod(np.imag(dshadPoint), dShad)) / dShad
            dStdNormFactor = np.sqrt((1 - 2 * dDistY + 2 * dDistY ** 2) * (1 - 2 * dDistX + 2 * dDistX ** 2))

            dSample1 = mtShadowingSamples[dYIndexP1, dXIndexP1]
            dSample2 = mtShadowingSamples[dYIndexP2, dXIndexP2]
            dSample3 = mtShadowingSamples[dYIndexP3, dXIndexP3]
            dSample4 = mtShadowingSamples[dYIndexP4, dXIndexP4]

            mtShadowingCorr[il, ic] = (
                ((1 - dDistY) * (dSample1 * (1 - dDistX) + dSample2 * dDistX) +
                dDistY * (dSample3 * (1 - dDistX) + dSample4 * dDistX))
                / dStdNormFactor
            )

# Plotting the shadowing attenuation surface
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(mtPosx, mtPosy, mtShadowingCorr, cmap='viridis')
plt.show()