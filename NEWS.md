# betanorm 0.1.1

* Hardened the internal C++ kernels against malformed arguments, none of
  which the exported functions can produce: the BMIQ block gather/scatter no
  longer overflow `int` on an `NA` or very large block start, the scatter
  refuses a destination it would otherwise silently coerce, and the quantile
  kernel sorts `NaN` with a defined order and clamps every target index.
  Results are bit-identical; the numeric snapshots are unchanged.

# betanorm 0.1.0

* Relicensed from MIT to GPL (>= 3). The quantile normalizer is derived from
  preprocessCore and the beta-mixture EM from RPMM, the latter via Steve
  Horvath's BMIQcalibration, which also underlies the BMIQ calibration
  (itself adapted from Andrew Teschendorff's BMIQ). Their authors are now
  credited as copyright holders; see `inst/COPYRIGHTS`. No code change.
