# MR_Tools
A package for processing and analyzing of MRI data and MR spectrsocopy research data.

The package contains functionality for:

  1)FileIO -- reading propretary Siemiens MRI and Bruker NMR raw data file formats. Also contains functionality for processing the raw time domain data into spatial / freqencey domain images / spectra. Outputs data to hdf5 for storage.
  
  2)MRI   -- contains tools for manipulating MRI images, masking / seqmenting images and selecting regions on intrest. Supports MRI - spectroscopy research workflow.
  
  3)SIM  -- is a python package for preforming Bloch Mcconnell spin-physics based simulations of MRI pulse sequences. Contains functionality for chaining the physcs based simulations to constrained nonlinear least squares optimizers in scipy and fitting experimental MRI data. Heavy numerical calculations are done in Numpy and Cython for speed.  Unit tests exists for this package. 
