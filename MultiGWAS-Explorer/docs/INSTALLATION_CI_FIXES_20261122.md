# Installation CI corrections

Two distinct installation failures were investigated.

1. The [earlier Windows run](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/runs/45033718485)
   rejected the recorded Cygwin Python in the isolated SASPy fixture and selected
   an unrelated `python.exe`. The shared discovery probe required a system
   `purelib` directory even when SASPy was importable through `PYTHONPATH`.
   Discovery now returns the actual import root from the SASPy specification.
   The regression explicitly makes the system `purelib` path absent and still
   verifies interpreter selection, CRLF handling, submitted HTML and remote
   directory-listing results without contacting SAS.
2. The [subsequent Docker run](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/runs/45100152221)
   omitted `tools/sync_macroparas_bundle.py` because of the existing
   `tools/*` Docker exclusion. The normal checkout tests passed, but the image
   could not import the bundle checker. The required tool is now explicitly
   included. A fast COPY-only Docker build verifies every file in the runtime
   manifest survives context filtering before the full dependency build starts.

The manifest hashes are refreshed for the changed discovery source and test.
The original port revision remains the source baseline; these installation
amendments are recorded separately. No research datasets, plotting macros,
SAS target-selection rules or private configurations are changed.

The shared discovery correction is also in MCP4SAS commit
`3098227cda343e4fec1e0342afb8bd88c80b1118`.
Focused bootstrap, discovery and portable-query regressions passed locally for
both repositories. After manifest refresh, the complete Explorer runtime suite
passed again (620 Perl assertions and 50 Python unit tests, plus source/manifest
checks).
The Docker daemon is unavailable on this local machine, so actual Docker
context/build verification runs in GitHub Actions rather than being claimed
from a mock build. Fresh installation, synthetic rendering and Apptainer
results must be read from the corresponding workflow jobs; offline checks do
not establish SAS authentication or scientific validity.
