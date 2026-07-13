# LEVAiM 0.1.0

* Added a frozen parse-specify-build-compile-fit API for spline and GP
  longitudinal trajectory models.
* Added exact and Hilbert-space approximate GP fitting with RBF, exponential,
  Matern-3/2, and Matern-5/2 covariance kernels.
* Added Gaussian, binomial, and beta response families across spline and GP
  engines with explicit response-boundary validation.
* Added subject-grouped validation, posterior GP derivatives, trajectory
  moments, plotting, assay-wide spline testing, trajectory comparison, and
  distance-matrix workflows.
* Added schema-versioned LEVAiM serialization with legacy raw-RDS compatibility.
