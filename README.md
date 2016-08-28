# Frame-Analysis

This is C code for analyzing things such as G87 format trajectory files.


  The program takes one command line argument, say, "xxx.g87", and outputs to "xxx.g87.yyyy.avg",
    where yyyy is the type of analysis run.


  If CHOP is nonzero, it will analyze NUMBER initial frames, beginning with
    the first, utilizing every initial frame until NUMBER have been analyzed,
    never extending the analysis beyond NUMBER frames from the beginning.

  If CHOP is zero, it will analyze NUMBER initial frames, beginning with 
    the first in the file, spreading the rest evenly throughout the file,
    unless the largest delta is greater than half the file, in which case
    the initial frames are spread evenly throughout the first half-ish of
    the file.


  At each initial frame, it will analyze final frames with

    delta = (final-initial) = 1, 2, 4, 8, 16, ..., (2^T), 2*(2^T), 3*(2^T), ..., N*(2^T)

    where T = TCYCLE and N = NCYCLE.
