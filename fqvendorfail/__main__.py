#!/usr/bin/env python
"""
Python Project Template
"""

import argparse
import collections
import datetime
import logging
import os
import sys
from typing import List, Optional

from fqvendorfail.vfail import vendorfail

try:
    from python_project import __version__
except Exception:
    __version__ = '0.0.0'

log = logging.getLogger(__name__)


def setup_logger(args):
    """Apply logging config from CLI args."""

    logging.basicConfig(
        level=args.log_level, format=args.log_format,
    )


def setup_parser():
    parser = argparse.ArgumentParser()

    logging_group = parser.add_argument_group("Logging")
    logging_group.add_argument(
        '--log-level',
        choices=[
            logging.INFO,
            logging.DEBUG,
            logging.CRITICAL,
            logging.WARNING,
            logging.ERROR,
        ],
        type=int,
        default=logging.INFO,
    )
    logging_group.add_argument(
        '--log-format',
        type=str,
        metavar='STR',
        default="%(asctime)s %(name)s:%(lineno)s %(levelname)s | %(message)s",
    )

    parser.add_argument('--data-dir', default=None)
    parser.add_argument('--version', action='version', version=__version__)

    parser.add_argument(
        '--output-prefix',
        '-o',
        dest='output_prefix',
        help='output prefix for filtered files',
        default='vf.',
        required=True,
    )
    parser.add_argument(
        'fastq_files',
        nargs='+',
        help='one or more correlated FASTQ files to be filtered together',
    )

    return parser


def process_args(argv: Optional[List] = None) -> collections.namedtuple:
    """Process args to NamedTuple.
    """

    parser = setup_parser()
    argv = argv or sys.argv[1:]

    if argv:
        args, unknown_args = parser.parse_known_args(argv)
    else:
        args, unknown_args = parser.parse_known_args()

    if args.data_dir:
        os.makedirs(args.data_dir)

    args_dict = vars(args)

    # Process extras list
    args_dict['extras'] = unknown_args

    # Recast to immutable namedtuple
    run_args = collections.namedtuple('RunArgs', list(args_dict.keys()))
    return run_args(**args_dict)


def run(run_args) -> int:
    """Method for running script logic."""

    ret_val = 0

    start_time = datetime.datetime.now()

    log.info("Running process...")
    vendorfail(run_args.fastq_files, run_args.output_prefix, log)
    # Log runtime info
    end_time = datetime.datetime.now()
    run_time = end_time - start_time
    log.info("Run time: %d seconds", run_time.seconds)
    return ret_val


def main(argv=None) -> int:
    """Main Entrypoint."""
    exit_code = 0

    args = process_args(argv)

    setup_logger(args)

    log.info("Process called with %s", args)
    log.info("ARGS {}".format(sys.argv))

    try:
        exit_code = run(args)
    except Exception as e:
        log.exception(e)
        exit_code = 1
    return exit_code


if __name__ == "__main__":
    """CLI Entrypoint"""

    status_code = 0
    try:
        status_code = main()
    except Exception as e:
        log.exception(e)
        sys.exit(1)
    sys.exit(status_code)


# __END__