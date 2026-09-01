# SounderPy pytest suite

Copy these files into the repository's `tests/` directory:

- `conftest.py`
- `test_retrieval_and_plotting.py`
- `test_file_io.py` (replace the stale existing file)

Merge the settings from `pyproject_pytest_snippet.toml` into the existing
`[tool.pytest.ini_options]` section of `pyproject.toml`.

## Install test tooling

```bash
python -m pip install -U pytest pytest-cov
```

## Normal test run

This runs offline API tests, file-I/O tests, and plot rendering tests. Live
retrieval tests are skipped unless explicitly enabled.

```bash
pytest -v
```

## Fast run without computationally heavier plots

```bash
pytest -m "not slow" -v
```

## Plotting tests only

```bash
pytest -m slow -v
```

## Live data retrieval tests

These contact remote NOAA/NCEI/Siphon/IEM data services and can fail when an
upstream service is unavailable even if SounderPy itself is correct.

```bash
SOUNDERPY_RUN_NETWORK_TESTS=1 pytest -m network -v
```

Run the RAP case alone:

```bash
SOUNDERPY_RUN_NETWORK_TESTS=1 pytest -m network -k rap -v
```

Run the RAOB case alone:

```bash
SOUNDERPY_RUN_NETWORK_TESTS=1 pytest -m network -k raob -v
```

Run the BUFKIT case alone:

```bash
SOUNDERPY_RUN_NETWORK_TESTS=1 pytest -m network -k bufkit -v
```

## Coverage

```bash
pytest --cov=sounderpy --cov-report=term-missing
```

For offline/fast coverage only:

```bash
pytest -m "not network and not slow" --cov=sounderpy --cov-report=term-missing
```
