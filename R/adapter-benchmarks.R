# Benchmarking for the did::att_gt()/DRDID::drdid_panel() adapters
# (R/did-adapter.R, R/drdid-adapter.R), analogous to dml_benchmark()/
# bench_fun() (R/benchmarks.R) for ordinary plm/npm dml() fits -- see there
# for the general gain.Y/gain.D/rho/delta definitions and their
# interpretation, all reused UNCHANGED here (normalize_benchmark_groups(),
# benchmark_gain_y(), benchmark_gain_d(), benchmark_rho(),
# warn_if_negative_gain()).
#
# THE KEY DIFFERENCE FROM plm/npm: there, benchmarking refits dml() itself
# with a covariate dropped from x, and theta.s/sigma2.s/nu2.s ALL shift
# together as a consequence (they're all estimated by the one refit). Here,
# theta.s is EXTERNAL (read from fit$ATT/mp$att[k], never re-estimated by
# this package) -- dropping a covariate from OUR OWN cross-fitted sigma2.s/
# nu2.s nuisances alone would never move theta.s, forcing Bias = theta.sj -
# theta.s = 0 identically and making rho/delta meaningless. So a
# leave-one-covariate-out comparison here has to refit the EXTERNAL
# estimator too (DRDID::drdid_panel()/did::att_gt() with the covariate
# dropped), not just did_cell_scale_nuisances()/did_cell_sigma2_nu2().
#
# Output shape: unlike bench_fun(), which produces one row per cf.reps
# repetition (dml() cross-fits multiple times), there is exactly one
# "repetition" here (did_cell_scale_nuisances()'s own L-fold cross-fitting
# is internal to a single sigma2.s/nu2.s estimate, not repeated). Both
# functions below therefore return $benchmarks[[covar]] as a 1-row
# data.frame and $benchmarks_psis[[covar]]$psi.* as length-1 lists --
# deliberately matching the shape bench_fun() would produce for cf.reps = 1
# -- so the EXISTING summary.dml_benchmark()/print.dml_benchmark() methods
# (R/benchmarks.R) work on this output completely unchanged.
#
# psi.GY/psi.GD/psi.rho below are the same formulas as bench_fun()'s (not
# factored into a shared helper, to avoid touching the working plm/npm
# code path in R/benchmarks.R) -- if those formulas are ever revised there,
# revise them here too.
#
# VERIFIED, AND A KNOWN CAVEAT (see project notes -- both checked by
# Monte Carlo, and against the EXISTING plm dml_benchmark() pathway under a
# comparable DGP, not just this file's new code):
#  - theta.sj/Bias/delta (the "without" refit and its difference from the
#    "with" fit) reproduce an independently-built "drop this covariate"
#    fit exactly.
#  - gain.Y calibrates well (analytical SE / Monte Carlo SD close to 1).
#  - gain.D calibrates well AWAY from nu.sq.wo near 0, but -- like nu2.s's
#    own SE (R/did-adapter.R) -- is fragile with a small denominator: near
#    nu.sq.wo = 0, both the point estimate and its SE blow up together.
#    Not a formula bug; a finite-sample feature of a ratio estimator with a
#    near-zero denominator.
#  - rho frequently saturates at its +/-1 cap (benchmark_rho()'s
#    pmin(1, ...)), and its SE noticeably understates the true sampling SD
#    even excluding capped draws. This is NOT specific to this adapter:
#    the SAME formula, exercised against the EXISTING plm dml_benchmark()
#    pathway with linear (yreg = dreg = "lm") nuisances and a single omitted
#    regressor, hits the +/-1 cap on effectively every draw too -- a known
#    exact identity in that linear/single-omitted-regressor case (Bias^2 =
#    V.g*V.a exactly), not a numerical defect. rho's point estimate should
#    be trusted; its SE should be treated with caution, here and for
#    plm/npm alike -- this is a pre-existing property of the general
#    benchmarking framework, not something introduced by this file.
#
# THEORETICAL-CONSISTENCY AUDIT (six-component hybrid estimator: theta.s/
# theta.s.wo external, sigma2.s/sigma2.s.wo/nu2.s/nu2.s.wo paper-native
# DML, combined via the manuscript's joint delta-method formulas -- see
# test-15-adapter-benchmark-audit.R for the tests this produced):
#  - SAME ANALYSIS POPULATION: for drdid, guaranteed by construction (D/
#    deltaY/w are the SAME vectors for both fits; dropping X's columns
#    cannot change row count; DRDID::drdid_panel() has no NA-driven row-
#    dropping of its own -- confirmed against its source). For did, this
#    was a REAL, CONFIRMED BUG before this audit: did::att_gt() does
#    complete-case exclusion across ALL xformla covariates during
#    preprocessing, so a covariate with missing values can cause att_gt()
#    to exclude a unit from the FULL fit that then silently RE-ENTERS the
#    REDUCED fit once that covariate is dropped (empirically confirmed:
#    dropping a covariate with a single NA changed a synthetic panel's
#    fitted population from 229 to 230 units, and the affected unit's own
#    (g,t) cell membership along with it). Fixed by check_same_did_
#    population() (whole-panel AND per-cell), called before any Bias/
#    gain.Y/gain.D/rho is computed for a dropped term -- errors clearly
#    rather than silently mixing two populations.
#  - ACTIVE TRIMMING checked on BOTH fits: already held, for free -- both
#    did_cell_short_results() and drdid_cell_short_results() call
#    validate_no_active_trimming() internally on whatever sample/
#    covariates they're given, and both the "with" and "without" short.
#    results below are built by calling those same functions on the full
#    and reduced covariate sets respectively. Verified with a constructed
#    DGP where dropping a covariate genuinely activates trimming that
#    was NOT active with it present (not just a same-trim.level check).
#  - EXTERNAL ATT IFs: already correctly aligned for drdid (D/deltaY/w
#    identical between fits, so fit$att.inf.func/fit.wo$att.inf.func are
#    automatically observation-aligned). For did, alignment relies on
#    mp$DIDparams$time_invariant_data's row order being determined solely
#    by (tname, gname, idname) sorting -- NOT by xformla -- so it is
#    identical between mp and mp.wo whenever check_same_did_population()
#    passes; that check is therefore also what makes the positional
#    alignment valid, not a separate re-indexing step.
#  - SCALE IFs (sigma2.s/nu2.s and their "wo" counterparts): already
#    correct -- did_cell_scale_nuisances()/did_cell_sigma2_nu2() (R/did-
#    adapter.R) return genuine Neyman-orthogonal influence functions, not
#    raw scores, unchanged by this audit.
#  - COVARIANCE PRESERVED: already correct BY CONSTRUCTION -- psi.GY/
#    psi.GD/psi.rho are linear combinations of the six components' own
#    (correctly-aligned) psis, so psi.sd() on those combinations captures
#    every pairwise covariance automatically (Var(sum) = sum(Var) +
#    2*sum(Cov)); no independence is assumed anywhere. This is exactly
#    bench_fun()'s own construction (R/benchmarks.R), reused unchanged.
#  - COVARIATE-DROPPING SEMANTICS: did drops an exact xformla TERM (not a
#    raw column) -- this does NOT cascade to interactions/transformations
#    built from the same variable (confirmed empirically: dropping "x1"
#    from ~x1 + x1:x2 leaves x1:x2 in the reduced formula, understating
#    x1's true removal effect). Fixed by check_no_term_leakage(): for a
#    SIMPLE (single-variable) dropped term, errors if any remaining term
#    still involves that variable, naming the offending term(s) and
#    suggesting they be dropped together. Interaction/transformed terms
#    (e.g. "x1:x2", "I(x1^2)") CAN be benchmarked directly -- the whole,
#    exact term is unambiguously removed, nothing further to check. drdid
#    benchmarks raw model.matrix COLUMNS (dml_benchmark()'s own existing
#    interface), not terms -- a factor's dummy columns must be grouped
#    explicitly via a named list (list(region = c("regionA","regionB")))
#    to be dropped together; the default, single-column-per-name behavior
#    benchmarks one dummy at a time. This is a real interface difference
#    from did (documented, not changed -- unifying it would need drdid_
#    dml_benchmark() to accept a formula it currently has no use for).
#  - OUTPUT SHAPE: combine.median()/combine.mean() (R/dml.R) were checked
#    to reduce to an exact no-op for a single "repetition" (verified
#    algebraically: ss = 0 when there is only one value, so se stays
#    exactly the input se) -- no inappropriate variance inflation/
#    deflation from treating the adapter's single cross-fit as if it were
#    several independent ones.

# Benchmarks a single DRDID::drdid_panel() comparison. `fit`/`D`/`deltaY`/
# `X`/`w` are the same arguments drdid_cell_short_results() takes (the
# ORIGINAL, "with every covariate" comparison); `benchmark_covariates` is
# the same argument dml_benchmark() takes (a character vector of columns in
# `X`, each benchmarked on its own, or a named list of column groups to
# drop together).
#
# For each covariate group: refits DRDID::drdid_panel() on the SAME
# (D, deltaY, w) with that group's columns dropped from X (same trim.level
# as `fit`, read from fit$argu$trim.level), then drdid_cell_short_results()
# on that refit -- which reruns the active-trimming check (for BOTH the
# full and reduced fit -- see file header) and the cross-fitted sigma2.s/
# nu2.s estimator on the reduced covariate set.
#
# Same-population invariant (see did_dml_benchmark()'s equivalent checks):
# for drdid this is guaranteed BY CONSTRUCTION, not checked -- D/deltaY/w
# are literally the same vectors reused for both the full and reduced fit
# (only X's COLUMNS change via X[, -index.o, drop=FALSE], which cannot
# change the ROW count), and DRDID::drdid_panel() has no complete-case row-
# dropping of its own (confirmed against its source: no NA handling at
# all -- it errors immediately on NA covariates via a downstream fastglm
# failure, rather than silently excluding rows). The anyNA() check below
# turns that eventual low-level crash into a clear, immediate one instead
# -- it is defensive/UX, not closing a real population-mismatch risk the
# way did_dml_benchmark()'s checks do (did::att_gt() DOES silently drop
# incomplete rows during preprocessing; drdid_panel() does not).
drdid_dml_benchmark <- function(fit, D, deltaY, X, w = NULL, benchmark_covariates,
                                tol = 1e-6, cf.folds = 5, cf.seed = 1) {
  validate_drdid_fit(fit)

  D <- as.integer(D)
  deltaY <- as.numeric(deltaY)
  X <- as.matrix(X)
  n <- length(D)
  if (length(deltaY) != n || nrow(X) != n) {
    stop("D, deltaY, and X must all have the same number of observations.")
  }
  if (anyNA(X)) {
    stop("X contains missing values. DRDID::drdid_panel() does not handle missing ",
        "covariates (it errors, rather than excluding incomplete rows) -- remove or ",
        "impute missing values before benchmarking.")
  }
  if (is.null(w)) w <- rep(1, n)
  w <- as.numeric(w) / mean(as.numeric(w))

  trim.level <- if (!is.null(fit$argu$trim.level)) fit$argu$trim.level else 0.995

  covariate.groups <- normalize_benchmark_groups(benchmark_covariates, X)

  sr <- drdid_cell_short_results(fit, D = D, deltaY = deltaY, X = X, w = w,
                                 tol = tol, cf.folds = cf.folds, cf.seed = cf.seed)

  benchmarks <- list()
  benchmarks_psis <- list()

  for (covar in names(covariate.groups)) {
    cols <- covariate.groups[[covar]]
    cat("\n=== Computing benchmarks using covariate:", covar, " ===\n\n")
    index.o <- which(colnames(X) %in% cols)
    Xo <- X[, -index.o, drop = FALSE]

    fit.wo <- DRDID::drdid_panel(y1 = deltaY, y0 = rep(0, n), D = D,
                                 covariates = Xo, i.weights = w,
                                 inffunc = TRUE, trim.level = trim.level)
    sr.wo <- drdid_cell_short_results(fit.wo, D = D, deltaY = deltaY, X = Xo, w = w,
                                      tol = tol, cf.folds = cf.folds, cf.seed = cf.seed)

    theta.s    <- sr$estimates$theta.s;    theta.s.wo <- sr.wo$estimates$theta.s
    sigma.sq   <- sr$estimates$sigma2.s;   sigma.sq.wo <- sr.wo$estimates$sigma2.s
    nu.sq      <- sr$estimates$nu2.s;      nu.sq.wo   <- sr.wo$estimates$nu2.s

    # Bias = theta.sj - theta.s (without minus with) -- see bench_fun()'s
    # comment (R/benchmarks.R) for why this direction, not its negation, is
    # the one that makes rho's sign match the manuscript's convention.
    Bias  <- theta.s.wo - theta.s
    V.g   <- sigma.sq.wo - sigma.sq
    V.a   <- nu.sq - nu.sq.wo
    valid <- V.g > 0 && V.a > 0

    Cor    <- benchmark_rho(Bias, V.g, V.a, valid)
    Gain.Y <- benchmark_gain_y(sigma.sq, V.g)
    Gain.D <- benchmark_gain_d(nu.sq.wo, V.a)
    warn_if_negative_gain(covar, Gain.Y, Gain.D, suggest_solutions = FALSE)

    benchmarks[[covar]] <- data.frame(
      gain.Y = Gain.Y, gain.D = Gain.D, rho = Cor,
      theta.s = theta.s, theta.sj = theta.s.wo, delta = Bias
    )

    psi.theta.s <- sr$psis$psi.theta.s;    psi.theta.s.wo <- sr.wo$psis$psi.theta.s
    psi.sigma2.s <- sr$psis$psi.sigma2.s;  psi.sigma2.s.wo <- sr.wo$psis$psi.sigma2.s
    psi.nu2.s <- sr$psis$psi.nu2.s;        psi.nu2.s.wo <- sr.wo$psis$psi.nu2.s

    psi.GY <- (sigma.sq * psi.sigma2.s.wo - sigma.sq.wo * psi.sigma2.s) / (sigma.sq^2)
    psi.GD <- (nu.sq.wo * psi.nu2.s - nu.sq * psi.nu2.s.wo) / (nu.sq.wo^2)
    psi.rho <- if (!valid) {
      rep(0, n)
    } else {
      ((psi.theta.s.wo - psi.theta.s) / sqrt(V.g * V.a)) -
        ((Bias * (psi.sigma2.s.wo - psi.sigma2.s)) / (2 * V.g^(3 / 2) * sqrt(V.a))) -
        ((Bias * (psi.nu2.s - psi.nu2.s.wo)) / (2 * V.a^(3 / 2) * sqrt(V.g)))
    }

    benchmarks_psis[[covar]] <- list(
      psi.theta.s = list(psi.theta.s), psi.sigma2.s = list(psi.sigma2.s), psi.nu2.s = list(psi.nu2.s),
      psi.theta.s.wo = list(psi.theta.s.wo), psi.sigma2.s.wo = list(psi.sigma2.s.wo), psi.nu2.s.wo = list(psi.nu2.s.wo),
      psi.GY = list(psi.GY), psi.GD = list(psi.GD), psi.rho = list(psi.rho)
    )
  }

  out <- list(benchmarks = benchmarks, benchmarks_psis = benchmarks_psis)
  class(out) <- "dml_benchmark"
  out
}

# Errors if dropping `dropped_term` (a SIMPLE variable, e.g. "x1" -- not
# itself an interaction/transformation) leaves that variable still
# influencing the model through some OTHER remaining term (e.g. "x1:x2",
# "I(x1^2)"). update(xformla, ~ . - term) only removes the EXACT term
# named -- confirmed against R's own formula algebra: dropping "x1" from
# ~x1 + x1:x2 leaves x1:x2 in place (R does not cascade-remove derived
# terms), so a naive "benchmark x1" would silently leave x1's interactive
# effect on x2 in BOTH the "with" and "without" fits, understating its true
# removal effect. Only checked for simple (single-variable) dropped terms:
# if the dropped term is ITSELF an interaction/transformation (e.g. the
# user explicitly benchmarks "x1:x2" or "I(x1^2)"), that exact term is
# fully and unambiguously removed, so there is nothing further to check.
check_no_term_leakage <- function(dropped_term, remaining_term_labels) {
  dropped_vars <- all.vars(stats::reformulate(dropped_term))
  if (length(dropped_vars) != 1L || !identical(dropped_vars, dropped_term)) {
    return(invisible(NULL))
  }
  leaking <- remaining_term_labels[vapply(remaining_term_labels, function(rt) {
    dropped_term %in% all.vars(stats::reformulate(rt))
  }, logical(1))]
  if (length(leaking) > 0L) {
    stop("Benchmarking covariate '", dropped_term, "': dropping this term still leaves '",
        dropped_term, "' influencing the model through ",
        paste0("'", leaking, "'", collapse = ", "), " (formula algebra does not cascade-remove ",
        "interaction/transformed terms derived from a dropped main effect). Drop those term(s) ",
        "together with '", dropped_term, "', e.g. benchmark_covariates = c(\"", dropped_term,
        "\", ", paste0("\"", leaking, "\"", collapse = ", "), "), or benchmark a different, ",
        "self-contained term.")
  }
  invisible(NULL)
}

# Errors if dropping a covariate changed att_gt()'s analysis population --
# att_gt() does complete-case exclusion across ALL xformla covariates
# during preprocessing (confirmed empirically: a covariate with a single
# NA causes did_standardization() to drop that unit's rows entirely when
# the covariate is IN the model, but not when it is dropped -- the unit
# silently "re-enters" the reduced fit). Checked twice: once for the whole
# panel (cheap, catches the common NA-driven case in one shot) and once
# per (group, t) cell actually being benchmarked (the invariant that
# actually matters for THIS cell's Bias/V.g/V.a comparison, and would also
# catch any other hypothetical full/reduced asymmetry the panel-level
# check doesn't).
check_same_did_population <- function(mp, mp.wo, term) {
  idname <- mp$DIDparams$idname
  ids.full <- mp$DIDparams$time_invariant_data[[idname]]
  ids.wo   <- mp.wo$DIDparams$time_invariant_data[[mp.wo$DIDparams$idname]]
  if (!setequal(ids.full, ids.wo)) {
    stop("Benchmarking covariate '", term, "': dropping it changes att_gt()'s analysis ",
        "population (", length(ids.full), " units with '", term, "' present vs. ",
        length(ids.wo), " units without it). This is most likely because '", term,
        "' has missing values that att_gt() excludes when the covariate is present but not ",
        "when it is dropped -- benchmarking requires the SAME analysis population for both ",
        "fits. Remove rows with missing '", term, "' from the data before calling att_gt(), ",
        "or choose a different benchmark covariate.")
  }
  invisible(NULL)
}

# Benchmarks a did::att_gt() "MP" object (est_method = "dr") by dropping one
# variable at a time from its xformla and refitting att_gt() ENTIRELY (every
# (g,t) cell at once, via mp$call -- the same match.call()-and-eval() trick
# bench_fun() uses on model$call, so the same caveat applies: any object
# mp$call references -- e.g. the original `data` -- must be reachable from
# wherever did_dml_benchmark() is CALLED from, not just from inside this
# package). See file header for why refitting the EXTERNAL estimator, not
# just this package's own scale nuisances, is necessary here.
#
# `benchmark_covariates` is a character vector of variable names as they
# appear in mp$DIDparams$xformla (e.g. "lpop"), each dropped on its own by
# formula update (update(xformla, ~ . - term)) -- deliberately a formula
# term, not a raw model.matrix column name (unlike drdid_dml_benchmark()):
# dropping a whole term lets a factor's dummy columns all get dropped
# together automatically via formula algebra, which a column-name-based
# interface (normalize_benchmark_groups(), as drdid_dml_benchmark() uses)
# cannot do losslessly in general (contrasts/interactions). See
# check_no_term_leakage() for what IS and ISN'T handled for interactions/
# transformations built from the dropped variable.
#
# Two invariants are checked before any Bias/gain.Y/gain.D/rho is computed
# for a dropped term, and error clearly (abort that term's benchmarking
# entirely, not a per-cell skip) if violated -- see check_same_did_
# population()/check_no_term_leakage() above:
#  (1) the full and reduced fits must share the exact same analysis
#      population (whole panel AND the specific cell being benchmarked);
#  (2) the dropped variable must not still influence the model through
#      some other remaining term.
# Active propensity-score trimming is checked too, for BOTH the full and
# reduced fit, but that happens for free: did_cell_short_results() (called
# for every cell of both `model` and `model.wo` below) already calls
# validate_no_active_trimming() internally on whatever sample/covariates
# it's given -- so a cell whose REDUCED fit newly trims a control (even
# though the full fit didn't) fails inside did_to_dml(mp.wo, ...) itself,
# which skips that cell with a warning identifying why (the same skip-on-
# failure path used for e.g. overlap/rank-condition failures).
#
# Returns a NAMED LIST of dml_benchmark-classed objects, one per (g,t) cell
# that both the original and every "without" fit computed successfully --
# named the same way did_to_dml()'s results$groups is (did_slot_name(),
# e.g. "g2004_t2006"). Each element has the exact same shape
# drdid_dml_benchmark() returns for a single comparison, so
# summary()/print() (R/benchmarks.R) work on each cell's result unchanged:
# e.g. print(did_dml_benchmark(mp, "lpop")$g2004_t2006). Cells where the
# "without" refit itself failed for a reason OTHER than a population
# mismatch (e.g. a new overlap/rank-condition failure, or newly-active
# trimming, after dropping the covariate) are silently absent from a given
# term's contribution rather than aborting the whole benchmark -- mirroring
# did_to_dml()'s own per-cell skip-on-failure behavior.
did_dml_benchmark <- function(mp, benchmark_covariates, cf.folds = 5, cf.seed = 1) {
  validate_did_scope(mp)

  xformla <- mp$DIDparams$xformla
  term.labels <- attr(stats::terms(xformla), "term.labels")
  missing.vars <- setdiff(benchmark_covariates, term.labels)
  if (length(missing.vars) > 0) {
    stop("Variable(s) not found in this MP object's xformla (", deparse(xformla),
        "): ", paste(missing.vars, collapse = ", "), ".")
  }

  model <- did_to_dml(mp, cf.folds = cf.folds, cf.seed = cf.seed)

  out <- list()
  for (term in benchmark_covariates) {
    cat("\n=== Computing benchmarks using covariate:", term, " ===\n\n")

    reduced.xformla <- stats::update(xformla, stats::as.formula(paste("~ . -", term)))
    check_no_term_leakage(term, attr(stats::terms(reduced.xformla), "term.labels"))

    mp.call <- mp$call
    mp.call$xformla <- reduced.xformla
    mp.wo <- eval(mp.call)
    check_same_did_population(mp, mp.wo, term)

    model.wo <- did_to_dml(mp.wo, cf.folds = cf.folds, cf.seed = cf.seed)

    for (k in seq_along(mp$group)) {
      group <- mp$group[k]; t <- mp$t[k]
      slot <- did_slot_name(group, t)
      if (is.null(model$results$groups[[slot]]) || is.null(model.wo$results$groups[[slot]])) next

      # Defense in depth: the whole-panel check above already guarantees
      # this in the common (NA-driven) case, but verify THIS cell's own
      # reconstructed sample matches too -- see check_same_did_population().
      ids.cell.full <- did_cell_sample(mp, group, t)$ids
      ids.cell.wo   <- did_cell_sample(mp.wo, group, t)$ids
      if (!setequal(ids.cell.full, ids.cell.wo)) {
        stop("Benchmarking covariate '", term, "': (group, t) = (", group, ", ", t, ")'s own ",
            "reconstructed analysis sample differs between the full and reduced fits even ",
            "though the whole panel's unit set matched -- this should not happen for this ",
            "adapter's supported scope; please report this as a bug.")
      }

      sr    <- model$results$groups[[slot]][[1]]
      sr.wo <- model.wo$results$groups[[slot]][[1]]

      theta.s <- sr$estimates$theta.s;    theta.s.wo <- sr.wo$estimates$theta.s
      sigma.sq <- sr$estimates$sigma2.s;  sigma.sq.wo <- sr.wo$estimates$sigma2.s
      nu.sq <- sr$estimates$nu2.s;        nu.sq.wo <- sr.wo$estimates$nu2.s

      Bias  <- theta.s.wo - theta.s
      V.g   <- sigma.sq.wo - sigma.sq
      V.a   <- nu.sq - nu.sq.wo
      valid <- V.g > 0 && V.a > 0

      Cor    <- benchmark_rho(Bias, V.g, V.a, valid)
      Gain.Y <- benchmark_gain_y(sigma.sq, V.g)
      Gain.D <- benchmark_gain_d(nu.sq.wo, V.a)
      warn_if_negative_gain(paste0(term, " (", slot, ")"), Gain.Y, Gain.D, suggest_solutions = FALSE)

      bench.slot <- if (is.null(out[[slot]])) list(benchmarks = list(), benchmarks_psis = list()) else out[[slot]]

      bench.slot$benchmarks[[term]] <- data.frame(
        gain.Y = Gain.Y, gain.D = Gain.D, rho = Cor,
        theta.s = theta.s, theta.sj = theta.s.wo, delta = Bias
      )

      psi.theta.s <- sr$psis$psi.theta.s;    psi.theta.s.wo <- sr.wo$psis$psi.theta.s
      psi.sigma2.s <- sr$psis$psi.sigma2.s;  psi.sigma2.s.wo <- sr.wo$psis$psi.sigma2.s
      psi.nu2.s <- sr$psis$psi.nu2.s;        psi.nu2.s.wo <- sr.wo$psis$psi.nu2.s

      psi.GY <- (sigma.sq * psi.sigma2.s.wo - sigma.sq.wo * psi.sigma2.s) / (sigma.sq^2)
      psi.GD <- (nu.sq.wo * psi.nu2.s - nu.sq * psi.nu2.s.wo) / (nu.sq.wo^2)
      psi.rho <- if (!valid) {
        rep(0, length(psi.theta.s))
      } else {
        ((psi.theta.s.wo - psi.theta.s) / sqrt(V.g * V.a)) -
          ((Bias * (psi.sigma2.s.wo - psi.sigma2.s)) / (2 * V.g^(3 / 2) * sqrt(V.a))) -
          ((Bias * (psi.nu2.s - psi.nu2.s.wo)) / (2 * V.a^(3 / 2) * sqrt(V.g)))
      }

      bench.slot$benchmarks_psis[[term]] <- list(
        psi.theta.s = list(psi.theta.s), psi.sigma2.s = list(psi.sigma2.s), psi.nu2.s = list(psi.nu2.s),
        psi.theta.s.wo = list(psi.theta.s.wo), psi.sigma2.s.wo = list(psi.sigma2.s.wo), psi.nu2.s.wo = list(psi.nu2.s.wo),
        psi.GY = list(psi.GY), psi.GD = list(psi.GD), psi.rho = list(psi.rho)
      )

      out[[slot]] <- bench.slot
    }
  }

  for (slot in names(out)) class(out[[slot]]) <- "dml_benchmark"
  out
}

# print.dml_benchmark()/summary.dml_benchmark() (R/benchmarks.R) expect a
# single dml_benchmark object -- did_dml_benchmark()/drdid_dml_benchmark()
# above instead return a NAMED LIST of them, one per (g,t) cell / comparison
# group, since there's no single whole-sample theta.s the way plm/npm have.
# print.dml.sensemakr()/summary.dml.sensemakr() (R/sensemakr.R) call this
# helper instead of print.dml_benchmark() directly, so they handle both
# shapes without needing to know which adapter (if any) produced `bb`.
print_bench_bounds <- function(bb, digits) {
  if (inherits(bb, "dml_benchmark")) {
    print.dml_benchmark(bb, digits = digits)
  } else {
    for (nm in names(bb)) {
      cat("\n--", nm, "--\n")
      print.dml_benchmark(bb[[nm]], digits = digits)
    }
  }
  invisible(bb)
}
