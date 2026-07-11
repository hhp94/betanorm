For nL = 2, use truncated component quantile mapping joined at the gold threshold.

Let:

    (t_s): sample U/M threshold
    (t_g): gold-standard U/M threshold
    (F_{sU},F_{sM}): fitted sample component CDFs
    (F_{gU},F_{gM}): fitted gold component CDFs

U class: (x \le t_s)

[
u=\frac{F_{sU}(x)}{F_{sU}(t_s)}
]

[
g(x)=F_{gU}^{-1}\left(uF_{gU}(t_g)\right)
]
M class: (x > t_s)

Using upper tails for numerical stability:

[
r=\frac{1-F_{sM}(x)}{1-F_{sM}(t_s)}
]

[
g(x)=F_{gM}^{-1}\left(
1-r\left[1-F_{gM}(t_g)\right]
\right)
]

In R terms:

# U
u <- pbeta(x, sample.aU, sample.bU) /
     pbeta(sample.threshold, sample.aU, sample.bU)

out <- qbeta(
  u * pbeta(gold.threshold, gold.aU, gold.bU),
  gold.aU,
  gold.bU
)

# M, using upper tails
r <- pbeta(
  x, sample.aM, sample.bM,
  lower.tail = FALSE
) / pbeta(
  sample.threshold, sample.aM, sample.bM,
  lower.tail = FALSE
)

out <- qbeta(
  r * pbeta(
    gold.threshold, gold.aM, gold.bM,
    lower.tail = FALSE
  ),
  gold.aM,
  gold.bM,
  lower.tail = FALSE
)

This guarantees:

    both sides map the sample cut to the gold cut;
    continuity at the U/M boundary;
    monotonicity and global rank preservation;
    endpoints remain 0 and 1;
    conditional quantiles within each hard class are preserved.

It may retain a slope kink, but not a jump. Removing the slope kink would require additional smoothing, such as a monotone spline bridge, and would further depart from component quantile normalization.

A global mixture-CDF map would also be continuous, but it can normalize away differences in U/M mixture proportions. The truncated-component approach better preserves BMIQ’s class-wise semantics.