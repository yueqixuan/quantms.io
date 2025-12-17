# Contributing to QPX

Thank you for your interest in contributing to QPX! This document provides guidelines and instructions for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
  - [Fork and Clone the Repository](#fork-and-clone-the-repository)
  - [Setting Up the Development Environment](#setting-up-the-development-environment)
- [Development Workflow](#development-workflow)
  - [Creating Feature Branches](#creating-feature-branches)
  - [Code Style and Formatting](#code-style-and-formatting)
  - [Running Tests](#running-tests)
  - [Running Linting and Code Quality Checks](#running-linting-and-code-quality-checks)
- [Submitting Changes](#submitting-changes)
  - [Creating a Pull Request](#creating-a-pull-request)
  - [Pull Request Description Best Practices](#pull-request-description-best-practices)
  - [Code Review Process](#code-review-process)
- [Issue Guidelines](#issue-guidelines)
  - [Reporting Bugs](#reporting-bugs)
  - [Requesting Features](#requesting-features)
- [Additional Information](#additional-information)
  - [Building and Previewing Documentation](#building-and-previewing-documentation)
  - [Contact Information](#contact-information)

## Code of Conduct

As part of our efforts toward delivering open and inclusive science, we follow the [Contributor Covenant Code of Conduct for Open Source Projects](https://www.contributor-covenant.org/version/2/0/code_of_conduct/). Please read and adhere to this code of conduct in all your interactions with the project.

## Getting Started

### Fork and Clone the Repository

1. Fork the repository on GitHub by clicking the "Fork" button on the [QPX repository page](https://github.com/bigbio/qpx).

2. Clone your fork locally:
   ```bash
   git clone https://github.com/your-username/qpx.git
   cd qpx
   ```

3. Add the upstream repository as a remote:
   ```bash
   git remote add upstream https://github.com/bigbio/qpx.git
   ```

### Setting Up the Development Environment

QPX uses [Poetry](https://python-poetry.org/) for dependency management. Make sure you have Python 3.10 or higher installed.

1. Install Poetry if you haven't already:
   ```bash
   pip install poetry
   ```

2. Install all dependencies including development dependencies:
   ```bash
   poetry install
   ```

   This will install:
   - All package dependencies
   - Development dependencies (pytest, pytest-timeout, hypothesis)
   - Documentation dependencies (mkdocs, mkdocs-material)

3. Verify the installation:
   ```bash
   poetry run qpxc --version
   ```

## Development Workflow

### Creating Feature Branches

Always create a new branch for your work:

```bash
# Update your local main branch
git checkout main
git pull upstream main

# Create a new feature branch
git checkout -b feature/your-feature-name
```

Use descriptive branch names:
- `feature/add-new-converter` for new features
- `fix/bug-description` for bug fixes
- `docs/update-readme` for documentation changes

### Code Style and Formatting

QPX follows Python code style standards enforced by Black and flake8.

#### Black Formatting

All Python code must be formatted with [Black](https://github.com/psf/black):

```bash
# Format all code
poetry run black .

# Check formatting without making changes
poetry run black . --check
```

#### Flake8 Linting

Code should pass flake8 checks:

```bash
# Run flake8 with project settings
poetry run flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
poetry run flake8 . --count --exit-zero --max-complexity=30 --max-line-length=130 --statistics --ignore F401,W503,W504
```

### Running Tests

QPX uses [pytest](https://pytest.org/) for testing. Tests are organized with markers for different test categories:

#### Test Markers

- **`integration`**: Integration tests that may require external resources or take longer to run
- **`large_data`**: Tests that require large data files

#### Running Tests

Run all tests except integration tests (recommended for local development):
```bash
poetry run pytest -vv -m "not integration"
```

Run all tests including integration tests:
```bash
poetry run pytest -vv
```

Run only integration tests:
```bash
poetry run pytest -vv -m "integration"
```

Exclude both integration and large data tests:
```bash
poetry run pytest -vv -m "not integration and not large_data"
```

Run a specific test file:
```bash
poetry run pytest -vv tests/test_fragment_annotations.py
```

Run with verbose output and show test durations:
```bash
poetry run pytest -vv --durations=10
```

### Running Linting and Code Quality Checks

Before submitting your changes, ensure all quality checks pass:

```bash
# Format code with Black
poetry run black .

# Run flake8 linting
poetry run flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics
poetry run flake8 . --count --exit-zero --max-complexity=30 --max-line-length=130 --statistics --ignore F401,W503,W504

# Run tests
poetry run pytest -vv -m "not integration"
```

## Submitting Changes

### Creating a Pull Request

1. Commit your changes with clear, descriptive commit messages:
   ```bash
   git add .
   git commit -m "Add feature: brief description of changes"
   ```

2. Push your branch to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

3. Go to the [QPX repository](https://github.com/bigbio/qpx) and click "New Pull Request".

4. Select your fork and branch, then click "Create Pull Request".

### Pull Request Description Best Practices

Write clear and comprehensive PR descriptions:

- **Title**: Use a concise, descriptive title
- **Description**: Include:
  - What changes were made and why
  - Any related issue numbers (e.g., "Fixes #123")
  - Steps to test the changes
  - Screenshots or examples if applicable
  - Any breaking changes or migration notes

Example PR template:
```markdown
## Description
Brief description of what this PR does.

## Related Issue
Fixes #123

## Changes Made
- Added new converter for XYZ format
- Updated documentation
- Added tests

## Testing
- [ ] All tests pass
- [ ] Code is formatted with Black
- [ ] Flake8 checks pass
- [ ] Manually tested with sample data

## Screenshots (if applicable)
```

### Code Review Process

1. All PRs require review before merging
2. Address review comments by pushing new commits to your branch
3. Once approved, a maintainer will merge your PR
4. Your contribution will be acknowledged in the project

## Issue Guidelines

### Reporting Bugs

When reporting bugs, please include:

1. **Clear title**: Summarize the issue in one line
2. **Description**: Detailed description of the bug
3. **Steps to reproduce**: Step-by-step instructions to reproduce the issue
4. **Expected behavior**: What you expected to happen
5. **Actual behavior**: What actually happened
6. **Environment**:
   - QPX version (`qpxc --version`)
   - Python version
   - Operating system
7. **Error messages**: Include full error messages and stack traces
8. **Sample data** (if applicable): Provide sample files or data to reproduce the issue

### Requesting Features

When requesting new features:

1. **Clear title**: Summarize the feature request
2. **Use case**: Describe the problem or use case this feature would solve
3. **Proposed solution**: Suggest how the feature might work
4. **Alternatives**: Mention any alternative solutions you've considered
5. **Additional context**: Any other relevant information

## Additional Information

### Building and Previewing Documentation

QPX uses [MkDocs](https://www.mkdocs.org/) with the Material theme for documentation.

#### Install Documentation Dependencies

Documentation dependencies are included when you run `poetry install`.

#### Build and Serve Documentation Locally

```bash
# Serve documentation locally with live reload
poetry run mkdocs serve
```

The documentation will be available at `http://127.0.0.1:8000/`.

#### Build Documentation for Production

```bash
poetry run mkdocs build
```

This generates static HTML files in the `site/` directory.

### Contact Information

For questions or assistance, you can contact the maintainers:

**PRIDE Team at EMBL-EBI:**
- Yasset Perez-Riverol (European Bioinformatics Institute - EMBL-EBI, U.K.)

**Collaborators:**
- Ping Zheng (Chongqing Key Laboratory of Big Data for Bio Intelligence, Chongqing University of Posts and Telecommunications, Chongqing, China)

You can also:
- Open an issue on [GitHub Issues](https://github.com/bigbio/qpx/issues)
- Start a discussion on [GitHub Discussions](https://github.com/bigbio/qpx/discussions) (if enabled)

## License

By contributing to QPX, you agree that your contributions will be licensed under the [Apache License 2.0](LICENSE).

## Adding Yourself as a Contributor

If you contribute to this project, please add your name to the list of contributors in the README.md file.

---

Thank you for contributing to QPX! Your efforts help make mass spectrometry data analysis more accessible to the scientific community.
