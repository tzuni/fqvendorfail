#!/usr/bin/env python

import importlib
import os
import subprocess
from textwrap import dedent
from types import SimpleNamespace

from setuptools import Command, find_packages, setup

PACKAGE = "fqvendorfail"
REPO = "fqvendorfail"

GITHUB_REPO = "https://github.com/tzuni/{}".format(REPO)

INSTALL_REQUIRES = []

DEV_REQUIRES = [
    'detect-secrets==0.13.1',
    'flake8',
    'isort',
    'pre-commit',
]

TESTS_REQUIRE = [
    'mock',
    'pytest',
    'pytest-cov',
]

setup(
    name=PACKAGE,
    description="Filter to remove FASTQ reads flagged as 'fail' by vendor software.",
    url=GITHUB_REPO,
    version='1.0.0',
    python_requires=">=3.6",
    packages=['fqvendorfail'],
    install_requires=INSTALL_REQUIRES,
    package_data={"": ["*.gz"]},
    entry_points = {
        'console_scripts': ['fqvendorfail=fqvendorfail.__main__:main']
    }
)

# __END__
