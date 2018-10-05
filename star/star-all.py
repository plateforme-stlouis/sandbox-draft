#!/usr/bin/python
# coding: utf-8


######
#
# Variable to tune (if needed)
#
# STAR   : Where the program STAR is.
# GENOME : Directory of the option --genomeDir.
# GTF    : File of the option --sjdbGTFfile.

STAR   = "~/star/STAR-2.6.0c/bin/MacOSX_x86_64/STAR"
GENOME = "~/star/genome/"
GTF    = "~/star/Homo_sapiens.GRCh38.79.gtf"

#
######



##
#
# Python code
#

#
import sys
import os
import re

GUI = False
LAUNCHER = 'launcher_ALL-STAR.sh'

def helper():
    print("""
Usage: python {0} --in-dir=/path/to/fastq --out-dir=/path/to/bam --exclude-dir=subdir/to/exclude

\t--in-dir : the directory containing all the FASTQ (possibly in subdirs)
\t--out-dir: the output directory which will contain the BAM and relatives

\t--exclude-dir: optionnal.

\t--help   : display this.

Example:

 python {0} --in-dir=~/star --out-dir=~/star/alt_bam


Note:

    o Edit the file {1} to tune the variables:
    \tSTAR   = path/to/star/binary\n\t     default: {2}
    \tGENOME = path/to/star/genome\n\t     default: {3}
    \tGTF    = path/to/annotations\n\t     default: {4}

    You can edit the file with TextEdit or Notepad or RStudio (or Emacs).

    o This script produces the shell script {5}
                                         located in input directory (--in-dir).
      + You can go to the input directory, type in the command-line:
          cd path/to/input/directory

          # For example, if you have --in-dir=~/star
          cd ~/star

      + You can edit it with TextEdit or whatever (or Emacs).
      + You can list all the samples found, type in the command-line:
          grep Sample {5}
      + You can list all the fastq files, type in the command-line:
          grep fastq {5}
    """.format(sys.argv[0],
               os.getcwd() + '/' + sys.argv[0],
               STAR, GENOME, GTF,
               LAUNCHER))
    return sys.exit(0)


def parse_args(argv):
    config = {'in': '', 'out': '', 'exclude': '-'}
    for arg in argv:
        if arg == "--help" or arg == "-h":
            helper()
        else:
            try:
                key, val = arg.split('=')
                k, d = key[2:].split('-')
                if d != 'dir':
                    helper()
                if k == 'out' or k == 'in' or k == 'exclude':
                    config[k] = val
            except SystemExit:
                print('WRONG: {}'.format(argv))
                sys.exit(1)
    try:
        for k in config.keys():
            if config[k] == '':
                helper()
    except SystemExit:
        print('WRONG: {}'.format(argv))
        sys.exit(1)
    for k in ['in', 'out']:
        config[k] = os.path.expanduser(config[k])
    return config


def get_FASTQ(config):
    root, excluded = config['in'], config['exclude']
    files = []
    for filename in [os.path.join(dirpath, y)
                     for dirpath, _, filenames in os.walk(root)
                     if not excluded in dirpath
                     for y in filenames]:
        dirname, ext = os.path.splitext(filename)
        if ext == '.fastq' or ext == '.fastq.gz':
            files.append(filename)
    return files


def splitter(files):
    splits = []
    for f in files:
        try:
            m = re.search('\w+_([\w]+)_[0-9]+_(R[12]|[12])_[\w .]+', f)
            sample, r = m.group(1), m.group(2)
            splits.append((sample, r, f))
        except:
            raise ValueError('File name is not compliant. {}'.format(f))
    return splits


class Sample:
    def __init__(self, name):
        self.name = name
        self.R1 = ''
        self.R2 = ''
        self.reader = ''

    def __repr__(self):
        s = """
        === {0}
        R1: {1}
        R2: {2}
        """.format(self.name, self.R1, self.R2)
        if self.reader != '':
            s += "uncompressed: {}\n\t".format(self.reader)
        return s


    def _try(self, key, val):
        if key == '':
            return val
        else:
            raise ValueError("Duplicates. {0} with {1}".format(self.name, val))

    def set_R(self, tag, f):
        if '1' in tag:
            self.R1 = self._try(self.R1, f)
        elif '2' in tag:
            self.R2 = self._try(self.R2, f)
        else:
            raise ValueError("Not expected. {0} {1}".format(tag, f))


def repairer(splits):
    samples = {}
    for n, t, f in splits:
        try:
            if not n in samples.keys():
                sample = Sample(n)
                samples[n] = sample
            samples[n].set_R(t, f)
        except:
            print("Sample {0} has more than 2 files.\n{1}\n{2}\n{3}\n".format(samples[n].name,
                                                                              samples[n].R1,
                                                                              samples[n].R2,
                                                                              f))
            raise
    return [ sample for sample in samples.values() ]

def compliance_reader(samples):
    smpls = []
    for sample in samples:
        _, ext1 = os.path.splitext(sample.R1)
        _, ext2 = os.path.splitext(sample.R2)
        if ext1 != ext2:
            raise ValueError("{0}\nIncoherent extension: {1} {2}".format(sample, ext1, ext2))
        if ext1 == ".fastq":
            sample.reader = 'cat'
        elif ext1 == ".gz":
            sample.reader = 'zcat'
        else:
            raise ValueError("Format not supported. {}".format(ext1))
        smpls.append(sample)
    return smpls


def launcher(name='foo',
             outdir='.',
             R1='foo_1.fastq.gz',
             R2='foo_2.fastq.gz',
             cat='zcat'):
    cmd = """


echo '# Sample: {sample}'

{STAR} \\
\t --runThreadN 4 --genomeDir {GENOME} \\
\t --sjdbGTFfile {GTF} \\
\t --sjdbOverhang 100 --readFilesIn \\
\t\t {R1} \\
\t\t {R2} \\
\t --readFilesCommand {cat} \\
\t --outSAMtype BAM SortedByCoordinate \\
\t --quantMode GeneCounts \\
\t --outFileNamePrefix "{outdir}/{sample}"

echo {sample} > done.txt

""".format(STAR=STAR, GENOME=GENOME, GTF=GTF,
               sample=name,
               R1=R1,
               R2=R2,
               cat=cat,
           outdir=outdir)

    return cmd



def generate_script(samples, config):
    indir, outdir = config['in'], config['out']

    filename = outdir + '/' + LAUNCHER
    if os.path.exists(outdir):
        print('Warning: {} already exists.'.format(outdir))
        ans = ''
        while not ans in ['y', 'n']:
            if GUI:
                ans = 'y'
            else:
                ans = raw_input('Overwrite ? y/n  ')
            if ans == 'n':
                print('Bye.')
                sys.exit(1)
            elif ans == 'y':
                print('Ok. Let overwrite!')
            else:
                print('Just one of these letters: y n')
    else:
        os.makedirs(outdir)
    with open(filename, 'w') as f:
        f.write("#!/usr/bin/bash\n# Automatically generated.\n\n")

        f.write("\n# Create output directory (if does not already exist)\n")
        f.write("mkdir -p {}\n\n".format(outdir))

        for sample in samples:
            f.write(launcher(sample.name,
                             outdir,
                             sample.R1,
                             sample.R2,
                             sample.reader))
            f.write('\n')
    print("\n# Now type in the command-line:\n\tbash {}".format(filename))


def run(config):
    print("Config done.")
    print(config)
    samples = compliance_reader(repairer(splitter(get_FASTQ(config))))
    print(samples)
    generate_script(samples, config)
    return


if __name__ == "__main__":

    if len(sys.argv) == 1:
        try:
            print("Import Graphical User Interface...")
            # GUI tools
            from Tkinter import Tk
            from tkFileDialog import askopenfilename, askdirectory
            from tkMessageBox import showinfo, showerror

            from Tkinter import *
        except ImportError:
            raise ImportError("Tkinter module is required. Try command-line interface.")
        print("done.\nLaunch GUI...")
        GUI = True

        window = Tk()
        window.title("Shell Script generator")

        lbl_star   = Label(window, text="star binary file:").grid(column=0, row=0)
        lbl_genome = Label(window, text="directory genome:").grid(column=0, row=1)
        lbl_gtf    = Label(window, text="annotation file:").grid(column=0, row=2)
        lbl_indir  = Label(window, text="Input:").grid(column=0, row=4)
        lbl_outdir = Label(window, text="Output:").grid(column=0, row=5)

        star = Entry(window, width=len(STAR),
                     textvariable=StringVar(window, STAR))
        star.grid(column=1, row=0)

        genome = Entry(window, width=len(GENOME),
                       textvariable=StringVar(window, GENOME))
        genome.grid(column=1, row=1)

        gtf = Entry(window, width=len(GTF),
                    textvariable=StringVar(window, GTF))
        gtf.grid(column=1, row=2)

        indir = Entry(window, width=len(os.getcwd()),
                      textvariable=StringVar(window, os.getcwd()))
        indir.grid(column=1, row=4)

        outdir = Entry(window, width=len(os.getcwd()) + 6,
                       textvariable=StringVar(window, os.getcwd() + "/bam/"))
        outdir.grid(column=1, row=5)


        def go():
            config = { 'in': indir.get(), 'out': outdir.get(), 'exclude': 'subdir3'}
            run(config)
            sys.exit(0)
            return

        btn_start = Button(window, text="Ok", command=go)
        btn_start.grid(column=1, row=7)


        def askfile(entry):
            filename = askopenfilename()
            entry.delete(0, END)
            entry.insert(0, filename)
            return

        btn_star = Button(window, text="File",
                          command=lambda: askfile(star)).grid(column=2, row=0)
        btn_gtf = Button(window, text="File",
                          command=lambda: askfile(gtf)).grid(column=2, row=2)

        def askdir(entry):
            filename = askdirectory(initialdir=os.getcwd())
            entry.delete(0, END)
            entry.insert(0, filename)
            return

        btn_genome = Button(window, text="Dir",
                          command=lambda: askdir(genome)).grid(column=2, row=1)
        btn_indir = Button(window, text="Dir",
                          command=lambda: askdir(indir)).grid(column=2, row=4)
        btn_outdir = Button(window, text="Dir",
                          command=lambda: askdir(outdir)).grid(column=2, row=5)



        print("done.")
        window.mainloop()

        sys.exit(0)


    else:
        config = parse_args(sys.argv[1:])
        #print(config)

        # fastqs = get_FASTQ(config['in'], config['exclude'])
        # print(fastqs)

        # splitteds = splitter(fastqs)
        # print(splitteds)

        # repairs = repairer(splitteds)
        # print(repairs)

        # compliants = compliance_reader(repairs)
        # print(compliants)

        run(config)
