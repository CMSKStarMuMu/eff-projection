# Instructions

This repository contains all the macros used to parametrise the **signal efficiency** for the Bd->KstMuMu analysis in Run2, using the **projection method**.
The workflow is structured as a chain of macros, each one taking as input a root file produced in the previous step.
The initial inputs are two ntuples, produced with the code contained in the
[B0KstMuMuNtuple repo](https://github.com/CMSKStarMuMu/B0KstMuMuNtuple)
( [snapshot](https://github.com/CMSKStarMuMu/B0KstMuMuNtuple/commit/d898840ee78df072b0d3862e3c141df79f0aeb5b) ),
 while the final output is a RooAbsReal object containing the 3D function that describes the efficiency.

## Macro list

#### Macros to produce the official efficiency description

```c++
void createEffHist(int q2Bin, int tagFlag=1, int xbins=25, int ybins = 0, int zbins = 0)
```
* `createEffHist.cc` to create binned efficiency from the ntuples
  * **q2Bin** is the bin number in the [0;8] range. Use `-1` to run recursively on all the bins
  * **tagFlag** is the flavor-tag status. Use `1` for correctly tagged, use `0` for wrongly tagged, use `2` to run recursively on both
  * **xbins** is the number of bins to use in the cos(theta_K) axis
  * **ybins** is the number of bins to use in the cos(theta_L) axis. Use `0` to copy the `xbins` value
  * **zbins** is the number of bins to use in the phi axis. Use `0` to copy the `xbins` value
  * Input: ntuple paths, to be specified in the code
  * Output: `./effHist_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins].root` containing two TH3D objects: `denHistb[q2Bin][ct|wt]__ctK_ctL_phi` and `numHistb[q2Bin][ct|wt]__ctK_ctL_phi`
  * Plots (to be activated/dis-activated with a hard-coded boolean `plot` ): three canvases with the efficiency histogram sliced for each of the three variables ( `./effHist_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins]_[CTK|CTL|PHI]slices_dp[depth].pdf`, where `depth` is one hundred times the dimension of the slices in the hidden variables (hard-coded variable)

```c++
void projEff_spHarm_fromHist(int q2Bin, int tagFlag, int maxOrder = 5, int xbins=25, int ybins = 0, int zbins = 0)
```
* `projEff_spHarm_fromHist.cc` to project the binned efficiency on spherical harmonics and save the RooAbsReal function
  * **maxOrder** is the maximum order of spherical harmonic functions to use
  * Input: `./effHist_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins].root` containing two TH3D objects: `denHistb[q2Bin][ct|wt]__ctK_ctL_phi` and `numHistb[q2Bin][ct|wt]__ctK_ctL_phi`
  * Output: `./effProjection_sh[maxOrder]o_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins].root`, containing a RooWorkspace `ws`, containing a RooAbsReal object `projectedFunc` and all its dependencies, included the three RooRealVar objects `ctK`, `ctL`, `phi`
  * Plots (to be activated/dis-activated with a hard-coded boolean `plot` ): three canvases with both efficiency function and histogram sliced for each of the three variables ( `./effProj_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins]_[CTK|CTL|PHI]slices_sh[maxOrder]o_dp[depth].pdf`, where `depth` is one hundred times the dimension of the slices in the hidden variables (hard-coded variable)

```c++
void plotEff(int q2Bin, int tagFlag, int maxOrder, int xbins=25, int ybins = 0, int zbins = 0)
```
* `plotEff.cc` to plot the efficiency and run the closure test
  * Input:`./effProjection_sh[maxOrder]o_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins].root`, containing a RooWorkspace `ws`, containing a RooAbsReal object `projectedFunc\
` and the three RooRealVar objects `ctK`, `ctL`, `phi`
  * Plots: `./effProj_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins]_1DProj_sh[maxOrder]o.pdf` and `./effProj_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins]_2DProj_sh[maxOrder]o.pdf` containing 1D and 2D projections of the efficiency function, respectively. `./closure_b[q2Bin][ct|wt]_[xbins]_[ybins]_[zbins]_sh[maxOrder]o.pdf` containing the closure test with comparison of RECO and efficiency-corrected GEN distributions

#### Macros under development, outdated, or to test alternative methods

* `test*.cc` set of macros using toy pseudo-experiments originated with analytical function, to test efficiency parameterisation
* `projEff_spHarm_kernel.cc` projection method applied to the output of a simple adaptive KDE sampling
