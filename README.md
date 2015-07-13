# JetTrackCorrelations

Implementation of this code for jet-track analysis:

(1) HT_make_ntuples contains code used to skim the HiForest to produce mini-ntuples with only the quantities relevant for this analysis.  JFF-dependent JEC are also applied to jtpt at this analysis step, and written into a branch called "corrpt"

(2) HT_Analyzer_All/HT_Analyzer_All_JFFCorr2.C is a unified analyzer that produces correlations from PbPb and pp data as well as pythia and pythia+hydjet simulation.

The output from this analyzer is then (locally) taken through the analysis procedure described in AN-14-082 and AN-15-159:

(3) Mixed event correction (me_correct/me_correct3.cxx runs on data, me_correct_mc/me_correct_mc3.cxx runs on *merged* monte carlo simulation)

(4) Background subtraction (bg_fit/bg_fit3.cxx)

(5) Determination of the background fluctuation bias correction (spill_over/spill_over_hydjet_only.cxx)

(6) Determination of the residual JFF-JEC/swapping correction (jff_residual/jff_residual.cxx)

(7) Projection, correction, and application systematic uncertainties (study_yield/study_yield.cxx)

(8) Final plotting of dEta and dPhi correlations and integrated excess yield (final_plots/PAS_Plots8.cxx)

(9) Width determination (width_determination/width_determination2.cxx)
