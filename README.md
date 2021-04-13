# predict_loading_concentration

## Dependencies:
* `python3`
* `pipenv`

## Installation:
`./install.sh`
Note this has the following output as of right now:
```
pipenv.patched.notpip._internal.exceptions.InstallationError: Command errored out with exit status 1: /Users/jdmccauley/.local/share/virtualenvs/library_loading_concentration_prediction-8ryN5LMN/bin/python /usr/local/lib/python3.9/site-packages/pipenv/patched/notpip install --ignore-installed --no-user --prefix /private/var/folders/22/d1_56bp93jq_z9m0qy5ql0z80000gq/T/pip-build-env-tdi9o8gr/overlay --no-warn-script-location --no-binary :none: --only-binary :none: -i https://pypi.org/simple -- wheel setuptools Cython 'numpy==1.9.3; python_version=='"'"'3.5'"'"'' 'numpy==1.12.1; python_version=='"'"'3.6'"'"'' 'numpy==1.13.1; python_version>='"'"'3.7'"'"'' Check the logs for full command output.
```
This is fine, just move on.

## Running:
`./run.sh sample_file.csv ladder_file.csv output_file.pdf`
