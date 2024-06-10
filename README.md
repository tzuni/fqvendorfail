# FASTQ vendor-fail filter

This tool scans FASTQ files to remove reads flagged as failing by the sequencer
software. This is a re-write of the existing tool. New features include bulk 
writes, and several attempts to avoid writing data at all. In the cases where
the FASTQ headers do not hold vendor fail flags, or all reads in a file are 
passing, the tool will avoid re-writing the data and simply create a symbolic
link. 

**Usage**

```
fqvendorfail -o <output_prefix> <fq1> [<fq2> ... <fqN>]
```

Probably this will only ever be used for single or paired end sequencing, but 
it can technically handle any number of correlated FASTQ files; filtering them
all in lock step.