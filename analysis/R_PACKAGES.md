# R package dependencies

The analysis scripts use packages from the following groups. Pin exact versions
in the publication environment (for example with `renv`) before final release.

Core analysis: `tidyverse`, `readr`, `readxl`, `lubridate`, `stringr`, `purrr`,
`ggplot2`, `lme4`, `lmerTest`, `emmeans`, `visreg`, `ggeffects`, `ggpubr`,
`svglite`, `plotly`, `htmlwidgets`.

QC/plots additionally use: `patchwork`, `flexplot`, `ggseg`, `gghalves`.

Gradient workflow additionally uses: `broom`, `SCORPIUS`, `conflicted`, `magick`,
`sf`, `ggpmisc`, `cowplot`, `pbapply`, `RcppCNPy`, `reticulate`, `ggsegSchaefer`,
`ggseg.formats`, plus packages referenced by the upstream-derived utility
functions when those functions are used.
