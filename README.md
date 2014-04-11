# PySurvey

## Overview

**PySurvey** is a python package suitable for analysing radio survey data released in FITS format. Features of the analysis include mainly calulation of column density distributions. The analysis can be performed for mosaics of different surveys (i.e., CGPS, SGPS, VGPS, LAB, Dame) and for two different species, namely atomic hydrogen (HI) and carbon monoxide molecular (CO), a tracer of molecular hydrogen.

To speed up the analysis process, each mosaic is divided into smaller sub-mosaics. These sub-mosaics can then be processed in parallel if more than one cpu is available. The global variable *glob_ncpu* in **SurveyUtils.py** allows to set the number of cpus that determines the number of sub-mosaics being processed in parrallel. At the end of the analysis all sub-mosaics are brought back together.
