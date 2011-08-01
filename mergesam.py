#!/usr/bin/env python

import atexit
import os
import os.path
import plac
import shutil
import sys
import tempfile
import logging
from subprocess import check_call, Popen, PIPE, CalledProcessError
from Bio import SeqIO
from itertools import *
from warnings import warn

logging.basicConfig(level=logging.INFO)

def fork_and_pump(input, handle):
    """Send input to handle in a child process.

This is needed to avoid deadlocks in bidirectional communication with
processes. In general, 'handle' should be int standard input of a
running command. This function forks a child process, and the child
reads lines from the iterable input and writes them to the handle.
When the input is exhausted, the child process closes the handle and
exits.

The child process ignores IOErrors, because they probably signify an
error somewhere else in the program."""
    if os.fork():
        # Parent
        handle.close()
    else:
        # Child
        try:
            handle.writelines(input)
            handle.close()
        # An IOError here means some *other* part of the program
        # crashed, so don't complain here.
        except IOError:
            pass
        exit()

class CommandPipeline(object):
    """A class that sets up a pipeline of commands.

Input is read from a filehandle or iterable object, and output can be
read from the 'stdout' property."""

    def __init__(self, command_list, input, suppress_stderr=True):
        assert len(command_list) >= 1
        logging.debug("Creating pipeline with commands: %s", command_list)
        if suppress_stderr:
            error_handle = open(os.devnull, "w")
        else:
            error_handle = sys.stderr
        commands = iter(command_list)
        first_command = next(commands)
        # The first proc needs special treatment to read from an input
        # iterable that isn't a file.
        try:
            first_proc = Popen(first_command, stdin=input, stdout=PIPE, stderr=error_handle)
        except AttributeError:
            first_proc = Popen(first_command, stdin=PIPE, stdout=PIPE, stderr=error_handle)
            input_iterator = iter(input)
            fork_and_pump(input_iterator, first_proc.stdin)
        # Append the rest of the commands to the pipeline, chaining
        # inputs to outputs
        self.procs = [first_proc]
        for cmd in commands:
            next_proc = Popen(cmd, stdin=self.procs[-1].stdout, stdout=PIPE, stderr=error_handle)
            self.procs[-1].stdout.close()
            self.procs.append(next_proc)
        # The last proc's stdout is the output of the pipeline
        self.stdout = self.procs[-1].stdout

    def close(self):
        """Close standard output and finish all processes."""
        self.stdout.close()
        for proc in self.procs:
            proc.wait()

    def poll(self):
        for proc in self.procs:
            proc.poll()

    def returncodes(self):
        """Return the list of return codes from all processes."""
        self.poll()
        return [proc.returncode for proc in self.procs]

    def failed(self):
        """Return True if any process failed."""
        return any(rc not in (0, None) for rc in self.returncodes())

    def succeeded(self):
        """Return True if every process succeeded."""
        return all(rc == 0 for rc in self.returncodes())

    def finished(self):
        """Return true if every process finished."""
        return all(rc is not None for rc in self.returncodes())

def mktempdir(*args, **kwargs):
    """Create a temporary directory with automatic cleanup.

This is identical to tempfile.mkdtemp, except that the directory and
its contents are registered for automatic deletion upon program exit."""
    path = tempfile.mkdtemp(*args, **kwargs)
    logging.debug("Created temp dir: %s", path)
    atexit.register(lambda(x): shutil.rmtree(x), path)
    return path

def is_fasta(filename):
    assert os.path.exists(filename)
    try:
        next(SeqIO.parse(filename, 'fasta'))
        logging.debug("Looks like fasta: %s (%s)", filename, "Successfuly read one fasta record")
        return True
    except:
        return False

def is_fai(filename):
    assert os.path.exists(filename)
    firstline = next(open(filename, 'r'))
    fields = firstline.split("\t")
    if len(fields) != 5:
        logging.debug("Not fai: %s (%s)", filename, "Should have 5 tab-separated fields")
        return False
    for f in fields[1:]:
        try:
            int(f)
        except ValueError:
            logging.debug("Not fai: %s (%s)", filename, "Fields 2 through 5 should be integers")
            return False
    else:
        logging.debug("Looks like fai: %s", filename)
        return True

read_devnull = open(os.devnull, "r")
write_devnull = open(os.devnull, "w")

def check_call_silent(*args, **kwargs):
    """subprocess.check_call with all handles set to /dev/null."""
    return check_call(stdin=read_devnull, stdout=write_devnull, stderr=write_devnull,
                                 *args, **kwargs)

def read_stdout_from_command(*args, **kwargs):
    """Run a command and return an iterator over its standard output.

Arguments are identical to subprocess.Popen except that stdout and
stderr are not available (stderr is discarded)."""
    return Popen(stdout=PIPE, stderr=write_devnull,
                            *args, **kwargs).stdout

def create_fasta_index(filename, samtools="samtools", *args, **kwargs):
    """Create a fasta index."""
    logging.info("Creating fasta index for %s", filename)
    check_call_silent([samtools, "faidx", filename])

def get_fai(filename):
    """Return the fai for a fasta file.

If the file is a fasta file, the corresponding fai file will be
checked for existence, and generated if it doesn't exist. It will be
generated in a temporary directory if necessary, in which case it will
be deleted upon completion.

If the file is not a fasta file, it is assumed to be a fai file and
returned directly."""
    logging.debug("Finding fai for %s", filename)
    if os.path.exists(filename):
        if is_fasta(filename):
            fai_filename = filename + ".fai"
            # Check for already-existing fai file
            if os.access(fai_filename, os.R_OK):
                logging.debug("Found existing fai file at %s", fai_filename)
                return fai_filename
            else:
                # Try to create nonexistent fai file
                try:
                    if os.path.exists(fai_filename):
                        raise Exception("Index exists but is unreadable")
                    logging.debug("Creating new fai file at %s", fai_filename)
                    create_fasta_index(filename)
                    return fai_filename
                # Try to create fai file in temporary directory where we
                # have write permission.
                except:
                    tempdir = mktempdir()
                    tempfilename = os.path.join(tempdir, os.path.basename(filename))
                    fai_filename = tempfilename + ".fai"
                    logging.debug("Creating new temp fai file at %s" + fai_filename)
                    os.symlink(os.path.abspath(filename), tempfilename)
                    create_fasta_index(tempfilename)
                    return fai_filename
        # Non-fasta file is assumed to be fai file.
        elif is_fai(filename):
            logging.debug("Found fai file at %s", filename)
            return filename
        else:
            raise ValueError("Not a fai or fasta file: %s" % filename)
    else:
        # File is nonexistent, so check if the fai file exists anyway.
        fai_filename = filename + ".fai"
        if os.path.exists(fai_filename) and is_fai(fai_filename):
            logging.debug("Found existing fai file at %s", fai_filename)
            return fai_filename
        else:
            raise IOError("Could not get fai for nonexistent file: %s" % filename)

def is_bam(filename, samtools="samtools"):
    """Returns True if filename represents a bam file."""
    # Assume standard input is not bam, since there isn't really any way to check.
    if filename == "-":
        logging.debug("Assuming standard input is sam")
        return False
    try:
        check_call_silent([samtools, "view", "-H", filename])
        return True
    except CalledProcessError:
        return False

def view_sam_header(filename, samtools="samtools", ref=None):
    """Return an iterator over the header lines of a sam or bam file."""
    command = [samtools, "view", "-H"]
    if not is_bam(filename, samtools=samtools):
        command.append("-S")
        if ref is not None:
            command.extend(["-t", ref])
    command.append(filename)
    return read_stdout_from_command(command)

def view_sam(filename, samtools="samtools", header=False, ref=None):
    """Return a handle on the output of samtools view."""
    command = [samtools, "view"]
    if header:
        command.append("-h")
    if not is_bam(filename, samtools=samtools):
        command.append("-S")
        if ref is not None:
            command.extend(["-t", ref])
    command.append(filename)
    return read_stdout_from_command(command)

def merge_sam(filenames, header=True, *args, **kwargs):
    """Return a handle on the concatenation of multiple sam files.

Unlike simply using 'cat', only the first file's header is copied to
the output. No checking is done to ensure that the files are
compatible with each other (i.e. mapped to the same reference).

Extra arguments (e.g. ref) are passed through to view_sam."""
    # Read the first header if requested, and none of the other files'
    # headers
    header_flags = chain([header], repeat(False))
    handles = imap(lambda f, h: view_sam(filename=f, header=h, *args, **kwargs), filenames, header_flags)
    return chain.from_iterable(handles)

def sam2bam_command(samtools="samtools", compress=True):
    """Take sam input and produce bam output"""
    command = [samtools, "view", "-bS", "-"]
    if not compress:
        command.append("-u")
    logging.debug("Sam2bam command: %s", " ".join(command))
    return command

def sort_bam_command(byname=False, mem=None, samtools="samtools"):
    """Take unsorted bam input and write sorted output."""
    command = [samtools, "sort", "-o"]
    if mem is not None:
        command.extend(["-m", str(int(mem))])
    if byname:
        command.append("-n")
    command.extend(["-", "-"])
    logging.debug("Sort bam command: %s", command)
    return command

def create_bam_index(filename, samtools="samtools"):
    assert os.path.exists(filename)
    check_call_silent([samtools, "index", filename])

def need_reference(files):
    """Return True if files require a reference index.

True if the first file (not including standard input) is a sam file
with no header, and False otherwise."""
    try:
        first_file = next(x for x in files if x != '-')
    except StopIteration:
        # This would mean that standard input is the only file, so we
        # can't guess whether a reference is required, and the user is
        # on her own.
        return None
    # If the file has no header, it needs a reference.
    return len(list(view_sam_header(first_file, ref=None))) == 0

@plac.annotations(
    # arg=(helptext, kind, abbrev, type, choices, metavar)
    # Flags
    nosort=("Do not sort or index bam output.", "flag", "s"),
    sortbyname=("Sort bam by name", "flag", "n"),
    noindex=("Do not index bam output file.", "flag", "i"),
    uncompressed=("Do not compress bam output.", "flag", "u"),
    sam=("Output in sam format instead of bam.", "flag", "S"),
    noheader=("Do not output sam header.", "flag", "d"),
    # Default
    ref=("Reference. This can be a fasta file or a fai file. An unindexed fasta file will be indexed automatically. This option is only required ", "option", "r"),
    samtoolspath=("Path to samtools.", "option"),
    sortmem=("Maximum memory for bam sorting (see samtools manpage).", "option", "m", int),
    outfile=("Output file name. Default is to write to standard output.", "option", "o", None, None, "FILE"),
    samfiles=("Input files to read and merge. Default is to read from standard input. Both sam and bam files can be supplied as inputs.", "positional"),
    )
def main(nosort, sortbyname, noindex, uncompressed, # BAM options
         sam, noheader,                             # SAM options
         ref=None, samtoolspath="samtools",         # External dependencies
         sortmem=500000000,                         # samtools sort option
         outfile="-", *samfiles):                   # I/O specification
    """Merge, convert, sort, index, and add headers to sam and bam files.

The primary purpose of this script is to take any number of sam or bam
files and merge them into one sorted, compressed, and indexed bam
file. Various options can produce other outputs."""
    # Read stdin by default
    if not samfiles:
        samfiles = ["-"]
    # Ensure that ref points to the fai file
    if ref:
        ref = get_fai(ref)
    # Check if a reference is needed
    elif not sam:
        needref = need_reference(samfiles)
        if needref:
            raise Exception("A reference index is required because the first input file has no header.")
        elif needref is None:
            logging.warn("A reference index may be required for proper operation.")
    # Handle all interactions between flags
    if sam:
        header = not noheader
        if nosort:
            logging.warn("--nosort option has no effect for sam output")
        if noindex:
            logging.warn("--noindex option has no effect for sam output")
        if uncompressed:
            logging.warn("--noindex option has no effect for sam output")
        if sortbyname:
            logging.warn("--sortbyname option  has no effect for sam output")
        index = False
    else:
        if noheader:
            logging.warn("--noheader option has no effect for bam output")
        if nosort and noindex:
            logging.warn("--noindex is redundant when --nosort is specified")
        if uncompressed and not nosort:
            logging.warn("--uncompressed has no effect unless --nosort is also specified")
        # Always need a header for bam
        header = True
        sort = not nosort
        index = sort and not noindex
        # Compress should be False if sort is True, because
        # "samtools sort" will do the compression.
        compress = not sort and not uncompressed
    logging.info("Input files: %s", list(samfiles))
    sam_input = merge_sam(filenames=samfiles, header=header, samtools=samtoolspath, ref=ref)
    if outfile == '-':
        output = sys.stdout
        # Cannot index standard output
        index = False
    else:
        output = open(outfile, "w")
    if sam:
        logging.info("Writing sam to %s", outfile)
        output.writelines(sam_input)
    else:
        command_list = []
        logging.info("Converting sam to bam")
        command_list.append(sam2bam_command(samtools=samtoolspath, compress=compress))
        if sort:
            logging.info("Sorting bam")
            command_list.append(sort_bam_command(samtools=samtoolspath, mem=sortmem))
        pipeline = CommandPipeline(command_list, sam_input, suppress_stderr=True)
        logging.info("Writing bam to %s", outfile)
        output.writelines(pipeline.stdout)
        pipeline.close()
        logging.debug("Pipeline return codes: %s", list(pipeline.returncodes()))
        if pipeline.failed():
            raise Exception("Pipeline failed")
    output.close()
    logging.info("Finished creating output file %s", outfile)
    if index:
        logging.info("Indexing %s", outfile)
        create_bam_index(outfile, samtools=samtoolspath)
    logging.info("Done")

# Entry point
def plac_call_main():
    return plac.call(main)

if __name__=="__main__":
    plac_call_main()
