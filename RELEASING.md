# Releasing to CRAN

Every push to `master` and every pull request runs `R CMD check --as-cran` on
Linux (R devel, release, oldrel), macOS and Windows, plus an ASAN/UBSAN build
(`.github/workflows/R-CMD-check.yaml`). The release workflow
(`.github/workflows/cran-release.yaml`) uses those checks as its gate.

## Steps

1. **Prepare the release on a branch and merge it to `master`:**
   - Bump `Version:` in `DESCRIPTION`. CRAN only accepts a version greater
     than the one it already has.
   - Add a `# AhoCorasickTrie <version>` section at the top of `NEWS.md`.
   - Update `cran-comments.md`: say what changed, and keep the check results
     and reverse dependency sections accurate. CRAN reviewers read this file.
   - If you changed roxygen comments or Rcpp exports, regenerate the docs and
     exports first: `Rcpp::compileAttributes(); roxygen2::roxygenise()`.

2. **Dry run.** In GitHub, open Actions, choose "CRAN release" and click
   "Run workflow", or run:

   ```sh
   gh workflow run cran-release.yaml
   ```

   This checks that the version is newer than CRAN's and that `NEWS.md` has
   the heading. It then runs the full check matrix, checks the reverse
   dependencies against the new version, and builds and checks the exact
   tarball that would be submitted. The tarball is uploaded as a run artifact.

3. **Submit.** Run the same workflow from `master` with `submit` checked:

   ```sh
   gh workflow run cran-release.yaml -f submit=true
   ```

   This repeats the dry run. It then uploads the tarball and `cran-comments.md`
   to CRAN and creates a **draft** GitHub release `v<version>`.

4. **Confirm the submission.** CRAN emails a confirmation link to the
   maintainer address in `DESCRIPTION`. The submission is not queued until
   that link is clicked.

5. **Wait for CRAN.** Usually you get an automated "pretest" email within a few
   hours, then a human review. That can take from a day to a week or more.
   - If CRAN accepts, publish the draft release. This creates the git tag:
     `gh release edit v<version> --draft=false`
   - If CRAN asks for changes, fix them on `master`, update `cran-comments.md`
     with a "Resubmission" section that answers each point, then rerun step 3.
     You can reuse the same version number if this version was never accepted.
     The draft release is updated in place.

## Notes

- The submission step calls `devtools:::upload_cran()`. That is the
  non-interactive core of `devtools::submit_cran()`, so it uses CRAN's real
  web form. If CRAN changes the form, updated devtools releases follow it.
- To submit by hand instead, run `devtools::submit_cran()` locally, or upload
  the tarball from the workflow artifact at
  <https://cran.r-project.org/submit.html>.
- For an extra pre-submission check on CRAN's own Windows R-devel machine, run
  `devtools::check_win_devel()`. Results are emailed to the maintainer.
- After acceptance, watch
  <https://cran.r-project.org/web/checks/check_results_AhoCorasickTrie.html>.
  CRAN sometimes reports sanitizer or new-compiler problems days later.
