# Contributing to taxdiv

Thank you for your interest in contributing to taxdiv! This document provides
guidelines for contributing to this package.

## How to Contribute

### Reporting Bugs

If you find a bug, please open a
[GitHub issue](https://github.com/mgorgoz/taxonomic-diversity-r/issues/new?template=bug_report.md)
with:

- A minimal reproducible example
- The output of `sessionInfo()`
- Expected vs actual behavior

### Suggesting Features

Feature requests are welcome! Please open a
[GitHub issue](https://github.com/mgorgoz/taxonomic-diversity-r/issues/new?template=feature_request.md)
describing:

- What you would like to see
- Why it would be useful
- Example usage if possible

### Pull Requests

1. Fork the repository and create a new branch from `main`
2. Make your changes
3. Add or update tests in `tests/testthat/`
4. Run `devtools::check()` and ensure 0 errors, 0 warnings, 0 notes
5. Update documentation if needed (`devtools::document()`)
6. Submit a pull request with a clear description of the changes

## Development Setup

```r
# Install development dependencies
install.packages(c("devtools", "testthat", "roxygen2", "knitr", "rmarkdown"))

# Clone and install
devtools::install_dev_deps()
devtools::load_all()

# Run tests
devtools::test()

# Full check
devtools::check()
```

## Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use roxygen2 for documentation
- All exported functions must have examples and tests
- Write commit messages in English

## Code of Conduct

Please note that this project follows a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating, you agree
to abide by its terms.

## Questions?

For questions about usage or development, open an issue on GitHub or contact
the maintainer at muratgorgoz350@gmail.com.
