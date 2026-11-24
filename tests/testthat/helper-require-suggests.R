# With LAMB_REQUIRE_SUGGESTS=true (set in CI) a missing suggested package fails
# the test instead of skipping it, so reader tests cannot drop out unnoticed.
skip_if_not_installed <- function(pkg, minimum_version = NULL) {
    if (identical(Sys.getenv("LAMB_REQUIRE_SUGGESTS"), "true") && !requireNamespace(pkg, quietly = TRUE)) {
        testthat::fail(paste0("Suggested package '", pkg, "' is not installed (LAMB_REQUIRE_SUGGESTS=true)."))
    }
    testthat::skip_if_not_installed(pkg, minimum_version)
}
