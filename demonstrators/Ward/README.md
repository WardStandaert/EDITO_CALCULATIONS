Scripts are named according to the numbering of the DATA folder:
Script 1-2 resamples the downloaded environmental layers (1.DOWNLOAD) to the preprocessed files (2.PREPROCESSED)
Scripts 2-3 start from the preprocessed files (2.PREPROCESSED) and create model outcomes (3.MODEL_OUTPUT) in different steps:
* 2-3_1 checks assumptions, removes outliers, filters occurrences and converts the preprocessed files (2.PREPROCESSED) to the desired formats for modelling
* 2-3_2 tests different sampling and model settings
* 2-3_3 creates models using best sampling and model settings and generates output (3.MODEL_OUTPUT)