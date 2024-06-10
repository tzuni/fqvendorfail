import gzip
import os
import re
from collections import namedtuple
from enum import Enum
from io import SEEK_CUR, BytesIO
from itertools import islice
from logging import Logger
from typing import Callable, Dict, Iterator, List, Optional, Text

Block = namedtuple('FastqBlock', ['position', 'line0', 'length'])
Range = namedtuple('Range', ['position', 'length'])


class Strategy(Enum):
    SKIP = 'Non-valid sequence identifier line'
    LINE_SCAN = 'Scanning Lines'
    SEEK_SCAN = 'Scanning ID Lines Only'


class Error(Exception):
    """Base class for vfail exceptions"""

    pass


class IrregularRecordLengthException(Error):
    """
    Exception Raised if irregularities are detected in the size of FASTQ
    records.
    """

    def __init__(self, message):
        self.message = message


def flexopen(filename: Text) -> Callable:
    """
    Flexible open function to deal with reading gzip
    """
    if filename.endswith(('.fq', '.fastq')):
        return open
    elif filename.endswith('.gz'):
        return gzip.open


def iter_linescan(filename: Text) -> Iterator[Block]:
    """
    FASTQ record iterator
    Uses readline() to read all 4 lines of a FASTQ record

    It's able to deal with variable length sequences but slower
    """
    openfun = flexopen(filename)

    with openfun(filename, 'rb') as fq:
        lines = []
        pos = fq.tell()
        while True:
            # read position at start of line
            line = fq.readline()
            if len(line) == 0:
                # maybe add check for set of lines < 4
                break
            lines.append(line)
            if len(lines) == 4:
                yield Block(pos, lines[0], sum(map(len, lines)))
                lines = []
                pos = fq.tell()


def iter_seekscan(filename: Text) -> Iterator[Range]:
    """
    FASTQ record iterator
    Uses readline() once then seek() to jump past 2nd, 3rd, and 4th lines

    Faster but depends on uniform sequence lengths and nothing on the 'plus'
    line.
    """
    openfun = flexopen(filename)
    with openfun(filename, 'rb') as fq:
        # determine length of seq, plus, and qual lines to use as offset
        line0 = fq.readline()
        lines = [fq.readline() for i in range(3)]
        offset = sum(map(len, lines))
        # reset to top of file
        fq.seek(0)
        # prepare to scan
        line = None
        pos = 0
        while True:
            line = fq.readline()
            # check for end of file
            if len(line) == 0:
                # maybe add check for set of lines < 4
                break
            # check for irregular record size
            # binary string slices as an int; int value of b'@' is 64
            if line[0] != 64:
                # somewhere the expected sequence length or plus line length changed
                # give warning and fail back to line_scan strategy
                raise IrregularRecordLengthException(
                    'FASTQ record sizes have become unpredictable'
                )
            fq.seek(offset, SEEK_CUR)
            # yield record
            yield Block(pos, line, len(line) + offset)
            pos = fq.tell()


def block_fail(block: Block) -> bool:
    """
    Checks whether the vendor fail flag has been set
    """
    FAIL = ':Y:'
    return FAIL in block.line0.decode("utf-8")


def bad_range(block: Block) -> Range:
    """
    Converts a Block to a Range
    """
    return Range(block.position, block.length)


def block_scan(
    fq_list: List, fq_iterfunc: Callable[[str], Iterator[Block]]
) -> Dict[str, List[Block]]:
    """
    Iterate over Blocks in FASTQ files saving any flagged blocks in bad_ranges
    """
    # save list of (position, length) tuples per fastq file
    bad_ranges = {fq: [] for fq in fq_list}
    # list of fastq generators (one per file)
    fq_gen_list = [fq_iterfunc(fq) for fq in fq_list]
    # print("reading from {}".format(fq_list))
    for blocks in zip(*fq_gen_list):
        # check if any blocks in this iteration were failed
        if any(map(block_fail, blocks)):
            # save all blocks in lists by file so that the corresponding read
            # is removed from all files
            for fq, block in zip(fq_list, blocks):
                bad_ranges[fq].append(bad_range(block))
    return bad_ranges


def validate_casava_header(line: str) -> bool:
    """
    Checks whether a FASTQ ID line matches the Casava 1.8+ format needed to
    look for vendor-flagged reads
    """
    # regex to match casava v1.8 fastq header
    # start = r'@'
    # instrument = r'[-_A-Z0-9]+'
    # run = r'\d+'
    # flowcell = r'[A-Z0-9]+'
    # lane = r'\d{1,2}'
    # tile = r'\d+'
    # x = r'\d+'
    # y = r'\d+'
    # member = r'[123]'
    # filtered = r'[YN]'
    # control_bits = r'\d+'
    # index = r'[ATGC0-9]+'

    casava_18_line_re = re.compile(
        r'@[-_A-Z0-9]+:'
        r'\d+:'
        r'[A-Z0-9]+:'
        r'\d{1,2}:'
        r'\d+:'
        r'\d+:'
        r'\d+'
        r' '
        r'[123]:'
        r'[YN]:'
        r'\d+:'
        r'[ATGC0-9]+'
    )
    match = casava_18_line_re.match(line)
    return match


def check_format(filename: Text, log: Logger) -> Strategy:
    """
    Quickly Analyze input files to determine best scanning strategy
    """
    openfun = flexopen(filename)
    with openfun(filename, 'rb') as fq:
        l1 = fq.readline()
        l2 = fq.readline()
        l3 = fq.readline()
        l4 = fq.readline()

    # -- Checks for conditions that would prevent any filtering

    # check to see that we are dealing with a casava v 1.8 file
    if not validate_casava_header(l1.decode('utf-8')):
        log.warn("Could not parse ID line: {}".format(l1.decode('utf-8')))
        log.warn("Unable to check this file.")
        return Strategy.SKIP

    # -- Checks for conditions that would mean a slower LINE_SCAN is required

    # check that '+' line is otherwise empty
    if l3.decode('utf-8') != '+\n':
        log.warn(
            "FASTQ record has non-empty 'plus' line: {}".format(l3.decode('utf-8'))
        )
        log.warn("Using Line-Scan method due to unpredictable record lengths.")
        return Strategy.LINE_SCAN

    # check that sequence lengths are uniform in first N records
    peek_records = 10000
    last_length = None
    count = 0
    for block in iter_linescan(filename):
        count += 1
        if count > peek_records:
            break
        if last_length is None:
            last_length = block.length - len(block.line0)
        elif last_length != block.length - len(block.line0):
            log.warn("Irregular FASTQ entry lengths detected")
            log.warn("Using Line-Scan method due to unpredictable record lengths.")
            return Strategy.LINE_SCAN

    # -- Passed all checks
    return Strategy.SEEK_SCAN


def merge_range_set(ranges: List[Range]) -> List[Range]:
    """
    Merge adjacent 'bad' ranges into continuous blocks that can be skipped
    when copying data
    """
    start = end = None
    merged_ranges = []
    for rng in ranges:
        # initialize with new range
        # print(rng)
        if start is None:
            start = rng.position
            end = rng.position + rng.length
        elif end == rng.position:
            # overlap
            end = rng.position + rng.length
        else:
            # complete range
            merged_ranges.append(Range(position=start, length=end - start))
            # start new range
            start = rng.position
            end = rng.position + rng.length
    # complete range
    merged_ranges.append(Range(position=start, length=end - start))
    return merged_ranges


def merge_ranges(
    bad_ranges: dict, fq_list: List, log: Logger
) -> Dict[str, List[Range]]:
    """
    Accept dictionary of bad ranges
    If the range lists are empty return them as-is
    Otherwise merge each individual list
    """
    if len(bad_ranges[fq_list[0]]) == 0:
        log.info("All reads passed. No merging required.")
        merged_ranges = bad_ranges
    else:
        merged_ranges = {fq: merge_range_set(bad_ranges[fq]) for fq in fq_list}
    return merged_ranges


def copy_chunks(
    src: BytesIO, dest: BytesIO, length: int, chunk_size: int = 64 * 2 ** 20
) -> None:
    """
    Copy data from source to destincation in large chunks.
    """
    nbytes = 0
    while length > 0:
        if chunk_size <= length:
            nbytes = chunk_size
        else:
            nbytes = length
        length -= nbytes
        # copy lines
        dest.write(src.read(nbytes))


def copy_last_chunk(
    src: BytesIO, dest: BytesIO, chunk_size: int = 64 * 2 ** 20
) -> None:
    """
    Copy the remaining good data in source file in a way that doesn't depend
    on knowing where the file ends since this isn't always possible in a
    gzipped file without decompressing it entirely.
    """
    while True:
        chunk = src.read(chunk_size)
        if len(chunk) == 0:
            break
        dest.write(chunk)


def symlink(filename: str, outfile: str) -> None:
    """
    Create symbolic link to source file.
    Useful in cases where FASTQ ID lines have no flagging information or where
    no bad reads were found.
    """
    os.symlink(os.path.basename(filename), outfile)


def get_outfile(out_prefix: str, ordinal: int) -> str:
    """
    Compose output filename from output prefix and pair ordinal.
    """
    out_filename = f'{out_prefix}_R{ordinal}.fq.gz'
    return out_filename


def get_outfile_dict(
    fq_list: List[str], out_prefix: str, pair_ordinals: List[str] = ['1', '2']
):
    """
    Compose output file names using a constant output prefix and a list of
    possible fastq ordinals.
    """
    outfile_dict = {
        fq: get_outfile(out_prefix, ordinal)
        for fq, ordinal in zip(fq_list, pair_ordinals)
    }
    return outfile_dict


def write_pass(fq: str, ranges: List[Range], out_filename: str) -> None:
    """
    Write good data from source `fq` file while skipping bad data in `ranges`.
    """
    # special case where bad_ranges is empty -> make symlink
    if len(ranges) == 0:
        symlink(fq, out_filename)
    else:
        openfun = flexopen(fq)
        # write good data in large chunks
        with openfun(fq, 'rb') as src, openfun(out_filename, 'wb') as dest:
            for range in ranges:
                # find distance from current position to beginning of bad range
                length = range.position - src.tell()
                # copy chunk of good data from src to dest
                copy_chunks(src, dest, length)
                # move to end of bad range
                src.seek(range.length, SEEK_CUR)
            # copy remaining good data
            copy_last_chunk(src, dest)


def scan(fq_list: List[str], strategy: Strategy, log: Logger) -> Dict[str, List[Block]]:
    """
    Executes the scanning strategy and returns a dictionary of bad ranges
    """

    if strategy == Strategy.SKIP:
        log.info("Unable to scan for fail-flagged reads: creating symlink")
        return {fq: [] for fq in fq_list}

    if strategy == Strategy.SEEK_SCAN:
        # try to execute SEEK_SCAN with fallback to LINE_SCAN
        try:
            log.info("Fast-scanning reads")
            result = block_scan(fq_list, iter_seekscan)
        except IrregularRecordLengthException:
            log.info("Failing back to line-scanning reads")
            result = block_scan(fq_list, iter_linescan)
    elif strategy == Strategy.LINE_SCAN:
        log.info("Line-scanning reads")
        result = block_scan(fq_list, iter_linescan)

    return result


def is_gzipped(filename):
    """
    Checks whether the file is gzipped
    This could probably be done in a more robust way.
    """
    return filename.endswith('.gz')


def vendorfail(fq_list: List[str], out_prefix: str, log: Logger) -> None:
    """
    Scans a list of coordinated FastQ Files to remove records marked as having
    failed vendor QC. Works on one or two (or more) FastQ files. Reads all
    files in lock step removing all concordant reads should any be marked as
    'failed'.
    """

    gzipped = is_gzipped(fq_list[0])

    # peek at the file details to pick a strategy with which to start
    log.info("Looking for required FASTQ features")
    strategy = check_format(fq_list[0], log)
    # out_filename_dict = {fq: get_outfile(fq, out_prefix) for fq in fq_list}

    # scan for bad ranges
    log.info("Scanning for failed reads.")
    bad_ranges = scan(fq_list, strategy, log)

    # merge adjacent ranges together
    log.info("Merging adjacent blocks.")
    merged_ranges = merge_ranges(bad_ranges, fq_list, log)

    # write good blocks to output fastq file
    log.info("Writing passing FASTQ reads.")
    out_file_dict = get_outfile_dict(fq_list, out_prefix)
    for fq, ranges in merged_ranges.items():
        write_pass(fq, ranges, out_file_dict[fq])
    log.info("DONE")


def main() -> None:
    """
    development CLI entrypoint for vfail
    """

    # 0.01 fractional bad data
    fq1 = 'fracfail_0.01_R1.fq'
    fq2 = 'fracfail_0.01_R2.fq'

    print("vendor fail main")
    fq_list = list(filter(lambda f: f is not None, [fq1, fq2]))
    out_prefix = 'vf.'
    vendorfail(fq_list=fq_list, out_prefix=out_prefix)


if __name__ == '__main__':
    main()