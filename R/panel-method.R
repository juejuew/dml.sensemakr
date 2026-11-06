# Normalize public method names and metadata from earlier adapter releases.
.panel_adapter_method <- function(method = "dml") {
  if (is.null(method)) return("dml")
  if (identical(method, c("dml", "replicate", "dml_refit")) ||
      identical(method, c("dml", "replicate")) ||
      identical(method, c("crossfit", "parametric"))) method <- method[1L]
  method <- match.arg(method, c("dml", "replicate", "dml_refit", "crossfit", "parametric"))
  switch(method, crossfit = "dml", parametric = "replicate", method)
}
