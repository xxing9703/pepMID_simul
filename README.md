# pepMID_simul
peptide MID simulation and turnover rate fittings in batch mode
This code is for quantification of peptide turnover rate (r) using isotope tracing. peptide_mid.m simulates the isotope distribution (I) of a peptide, for both unlabeled and labeled. The measured isotope distribution is fitted to the equation below by least square fitting.
I_measured=I_unlabeled * r +I_labeled * (1-r)
