#! /bin/sh

# As of August 19th, 2025
twine upload --repository testpypi dist/*

# Either way, test it
firefox https://test.pypi.org/project/corestone/
