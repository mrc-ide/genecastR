## Render the quickstart vignette to GitHub-flavoured Markdown, so that it is
## readable in the repository browser without installing the package.
##
## The Rmd is copied into docs/ before rendering rather than rendered with
## `output_dir = "docs"`.  With output_dir, knitr cannot express the figure
## directory relative to the Rmd and writes absolute paths into the Markdown,
## which break on GitHub and leak the local home directory.  Rendering from
## inside docs/ keeps the paths relative.
##
## Needs pandoc.  Outside RStudio, point RSTUDIO_PANDOC at the bundled copy:
##   export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64
##
## Run from the package root:  Rscript data-raw/render_quickstart.R

src <- "vignettes/GeneCastRQuickstart.Rmd"
stopifnot(file.exists(src))

unlink("docs", recursive = TRUE)
dir.create("docs", showWarnings = FALSE)

tmp <- file.path("docs", "quickstart.Rmd")
invisible(file.copy(src, tmp, overwrite = TRUE))

## `on.exit` fires on function exit, not at the top level of a script, so the
## copy is removed explicitly -- and before the check below, so a failure
## cannot leave a stray Rmd in docs/.
ok <- tryCatch({
  rmarkdown::render(tmp,
                    output_format = rmarkdown::github_document(
                      html_preview = FALSE),
                    quiet = TRUE)
  TRUE
}, error = function(e) {unlink(tmp); stop(e)})
unlink(tmp)

md <- readLines(file.path("docs", "quickstart.md"), warn = FALSE)
if (any(grepl(normalizePath(".", winslash = "/"), md, fixed = TRUE)))
  stop("absolute paths leaked into docs/quickstart.md", call. = FALSE)

message("wrote docs/quickstart.md (",
        length(list.files("docs", "[.]png$", recursive = TRUE)),
        " figures)")
