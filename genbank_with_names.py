#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""
Genbank record operations based on biopython Bio.SeqIO
https://github.com/biopython/biopython/blob/master/Bio/SeqIO/InsdcIO.py
"""

import os.path as op
import sys
import logging

from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq

from jcvi.formats.base import must_open, FileShredder, BaseFile, get_number
from jcvi.formats.gff import GffLine
from jcvi.apps.fetch import entrez
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, mkdir, glob


MT = "mol_type"
LT = "locus_tag"


class MultiGenBank (BaseFile): 
# class = blueprint for object MultiGenBank
# BaseFile = parentclass/superclass that passes on all functions and 
# variables to MultiGenBank, which inhererits from ancestor
    """
    Wrapper for parsing concatenated GenBank records.
    """

    # This is the constructor/initializer. Think of it as the factory: 
    # whenever a new object "MultiGenBank" gets made, this code will 
    # execute to set some default values for our variables.
    def __init__(self, filename, source="JCVI"):
    # The 'self' argument is required first and
    #  refers to this particular MultiGenBank object
        super(MultiGenBank, self).__init__(filename)
        # Super calls the constructor of the parent class.
        assert op.exists(filename) 
        # Assert ensures that the expression 'op.exists(filename)' is true.
        # If the expression is false (filename op does not exist), 
        # then it raises an exception error

        pf = filename.rsplit(".", 1)[0]
        # filename.rsplit separates the filename by the last period(.),
        # returns the first part of the filename before the period,
        # and assigns it to be pf 
        fastafile, gfffile = pf + ".fasta", pf + ".gff"
        # Multiple assignment defines the following code for each name
        fasta_fw = must_open(fastafile, "w")
        gff_fw = must_open(gfffile, "w")

        # Defining attributes of self (MulitGenBank object)
        self.source = source
        self.counter = defaultdict(list)

        nrecs, nfeats = 0, 0
        for rec in SeqIO.parse(filename, "gb"):
            seqid = rec.name
            rec.id = seqid
            SeqIO.write([rec], fasta_fw, "fasta")
            rf = rec.features
            for f in rf:
                self.print_gffline(gff_fw, f, seqid)
                nfeats += 1

                for sf in f.sub_features:
                    self.print_gffline(gff_fw, sf, seqid, parent=f)
                    nfeats += 1

            nrecs += 1

        logging.debug("A total of {0} records written to `{1}`.".\
                        format(nrecs, fastafile))
        fasta_fw.close()

        logging.debug("A total of {0} features written to `{1}`.".\
                        format(nfeats, gfffile))
        gff_fw.close()

    # Class method requires first argument to be "self"
    def print_gffline(self, fw, f, seqid, parent=None):

        score = phase = "."
        type = f.type
        if type == "source":
            type = "contig"

        attr = "ID=tmp"
        source = self.source

        start = get_number(f.location.start) + 1
        end = get_number(f.location.end)
        strand = '-' if f.strand < 0 else '+'
        g = "\t".join(str(x) for x in \
            (seqid, source, type, start, end, score, strand, phase, attr))
        g = GffLine(g)

        qual = f.qualifiers
        id = "tmp"
        if MT in qual:
            id = seqid
        elif LT in qual:
            id, = qual[LT]
        else:
            qual[LT] = [self.current_id]
            id, = qual[LT]

        id = id.split()[0]

        if parent:
            id, = parent.qualifiers[LT]
            id = id.split()[0]

        assert id != "tmp", f
        oid = id
        self.counter[(oid, type)].append((start, end))
        count = len(self.counter[(oid, type)])

        if type in ("mRNA", "gene"):
            if type == "gene" and count > 1:
                return
            self.start = min(a for a, b in self.counter[(id, type)])
            self.end = max(a for a, b in self.counter[(id, type)])
            self.set_attribute("gene", "Alias", qual, g)
            self.set_attribute("product", "Note", qual, g)
        else:
            suffix = ".{0}.{1}".format(type.lower(), count)
            id = id + suffix
            g.attributes["Parent"] = [oid]
            self.set_attribute("product", "Note", qual, g)

        g.attributes["ID"] = [id]
        g.update_attributes()
        print >> fw, g

        self.current_id = oid

    def set_attribute(self, gb_tag, gff_tag, qual, g):
        if gb_tag in qual:
            tag, = qual[gb_tag]
            g.attributes[gff_tag] = [tag]

# Class outlines the parameters for the new object GenBank
class GenBank(dict):
# dict is the class from which the object BenBank inherits functions and methods
    """
    Wrapper of the GenBank record object in biopython SeqIO.
    """
    # The constructor/initializer acts as a factory: 
    # whenever a new object "GenBank" gets made, this code will 
    # execute to set some default values for our variables.
    def __init__(self, filenames=None, accessions=None, idfile=None):
        self.accessions = accessions
        self.idfile = idfile
        # 'Self' now refers to the GenBank object that is created

        if filenames is not None: # If filnames exist
            self.accessions = [op.basename(f).split(".")[0] for f in filenames]
            d = dict(SeqIO.to_dict(SeqIO.parse(f, "gb")).items()[0] \
                for f in filenames)
            for (k, v) in d.iteritems():
                self[k.split(".")[0]] = v

        elif idfile is not None:
            gbdir = self._get_records()
            d = dict(SeqIO.to_dict(SeqIO.parse(f, "gb")).items()[0] \
                for f in glob(gbdir + "/*.gb"))
            for (k, v) in d.iteritems():
                self[k.split(".")[0]] = v

        else:
            sys.exit("GenBank object is initiated from either gb files or "\
                "accession IDs.")

    def __getitem__(self, accession):
        rec = self[accession]
        return rec

    def __repr__(self):
        recs = []
        for accession in self.keys():
            recs.append([accession, self.__getitem__(accession)])
        return recs

    def _get_records(self):
        gbdir = "gb"
        dirmade = mkdir(gbdir)
        if not dirmade:
            sh("rm -rf {0}_old; mv -f {0} {0}_old".format(gbdir,))
            assert mkdir(gbdir)

        entrez([self.idfile, "--format=gb", "--database=nuccore", "--outdir={0}"\
            .format(gbdir)])

        logging.debug('GenBank records written to {0}.'.format(gbdir))
        return gbdir

    @classmethod
    def write_genes_bed(cls, gbrec, outfile, emit_cds):
        seqid = gbrec.id.split(".")[0]
        if not seqid:
            seqid = gbrec.name.split(".")[0]

        genecount = 0
        consecutivecds = 0
        for feature in gbrec.features:
            if feature.type == "gene":
                # Create bedline from gene record
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = feature.qualifiers['gene'][0]
                except:
                    try: # some features only have a locus tag
                       name = feature.qualifiers['locus_tag'][0]
                    except:
                        # who knows what it's call - use the gene count to name
                        name = "{0}_{1}".format(seqid, genecount)
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = seqid + "\t{0}\t{1}\t{2}\t1000\t{3}\t{0}\t{1}\t65,105,225\n".format(start, stop, name, strand)
                outfile.write(bed_line) 
                # Return to original code               
                genecount+=1
                consecutivecds = 0
                continue # Jumps back to the for loop for next gene and skips code below
                
            if feature.type == 'CDS':
                if emit_cds == True:
                    if consecutivecds:
                        genecount+=1 # Only executed if consecutivecds is not 0
                    consecutivecds = 1
                    start = feature.location.start
                    stop = feature.location.end
                    if start > stop: start, stop = stop, start
                    if feature.strand < 0:
                        strand = "-"
                    else:
                        strand = "+"
                    # Changed from score =  "." because IGV doesn't include direction of translation without a score ("."= no score)
                    score = "1000"
                    accn = "{0}_{1}".format(seqid, genecount)

                    start = str(start).lstrip("><")
                    stop = str(stop).lstrip("><")
                    bedline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".\
                        format(seqid, start, stop, accn, score, strand)
                    outfile.write(bedline)

    @classmethod
    # @ is a decorator that labels the class method
    def write_genes_fasta(cls, gbrec, fwcds, fwpep):
        seqid = gbrec.id.split(".")[0]
        if not seqid:
            seqid = gbrec.name.split(".")[0]

        genecount = 0
        consecutivecds = 0
        for feature in gbrec.features:
            if feature.type == "gene":
                genecount+=1
                consecutivecds = 0
                continue

            if feature.type == 'CDS':
                if consecutivecds:
                    genecount+=1
                consecutivecds = 1
                accn = "{0}_{1}".format(seqid, genecount)

                if len(feature.sub_features) == 0:
                    seq = feature.extract(gbrec.seq)
                else:
                    seq = []
                    for subf in sorted(feature.sub_features, \
                        key=lambda x: x.location.start.position*x.strand):
                        seq.append(str(subf.extract(gbrec.seq)))
                    seq = "".join(seq)
                    if Seq(seq).translate().count("*")>1:
                        seq = []
                        for subf in feature.sub_features:
                            seq.append(str(subf.extract(gbrec.seq)))
                        seq = "".join(seq)
                    seq = Seq(seq)

                fwcds.write(">{0}\n{1}\n".format(accn, seq))
                fwpep.write(">{0}\n{1}\n".format(accn, seq.translate()))

    def write_genes(self, output="gbout", individual=False, pep=True, emit_cds=True):
        """ If multiple accession numbers, creates and opens 3 files in 3 formats and
        writes("w") the genomic data for each accession number in 1 file for each format.""" 
        if not individual: 
            fwbed = must_open(output+".bed", "w")
            fwcds = must_open(output+".cds", "w")
            fwpep = must_open(output+".pep", "w")
        """ For-loop that iterates the opening of a file and writing to 
        that file for each format and every accession number."""
        for recid, rec in self.iteritems():
            if individual:
                mkdir(output)
                fwbed = must_open(op.join(output, recid+".bed"), "w")
                fwcds = must_open(op.join(output, recid+".cds"), "w")
                fwpep = must_open(op.join(output, recid+".pep"), "w")

            GenBank.write_genes_bed(rec, fwbed, emit_cds)
            GenBank.write_genes_fasta(rec, fwcds, fwpep)

        if not pep:
            FileShredder([fwpep.name])

    def write_fasta(self, output="gbfasta", individual=False):
        if not individual:
            fw = must_open(output+".fasta", "w")

        for recid, rec in self.iteritems():
            if individual:
                mkdir(output)
                fw = must_open(op.join(output, recid+".fasta"), "w")

            seqid = rec.id.split(".")[0]
            if not seqid:
                seqid = rec.name.split(".")[0]
            seq = rec.seq
            fw.write(">{0}\n{1}\n".format(seqid, seq))


def main():

    actions = (
        ('tofasta', 'generate fasta file for multiple gb records'),
        ('getgenes', 'extract protein coding genes from Genbank file'),
        ('gff', 'convert Genbank file to GFF file'),
              )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def gff(args):
    """
    %prog gff seq.gbk

    Convert Genbank file to GFF and FASTA file.
    The Genbank file can contain multiple records.
    """
    p = OptionParser(gff.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gbkfile, = args
    MultiGenBank(gbkfile)


def preparegb(p, args):
    p.add_option("--gb_dir", default=None,
            help="path to dir containing GanBank files (.gb)")
    p.add_option("--id", default=None,
            help="GenBank accession IDs in a file. One ID per row, or all IDs" \
            " in one row comma separated.")
    p.add_option("--simple", default=None, type="string",
            help="GenBank accession IDs comma separated " \
            "(for lots of IDs please use --id instead).")
    p.add_option("--individual", default=False, action="store_true",
            help="parse gb accessions individually [default: %default]")
    opts, args = p.parse_args(args)
    accessions = opts.id
    filenames = opts.gb_dir

    if not (opts.gb_dir or opts.id or opts.simple):
        sys.exit(not p.print_help())

    if opts.gb_dir:
        filenames = glob(opts.gb_dir+"/*.gb")

    if opts.id:
        rows = file(opts.id).readlines()
        accessions = []
        for row in rows:
            accessions += map(str.strip, row.strip().split(","))

    if opts.simple:
        accessions = opts.simple.split(",")

    if opts.id or opts.simple:
        fw = must_open("GenBank_accession_IDs.txt", "w")
        for atom in accessions:
            print >>fw, atom
        fw.close()
        idfile = fw.name
    else:
        idfile=None

    return (filenames, accessions, idfile, opts, args)


def tofasta(args):
    """
    %prog tofasta [--options]

    Read GenBank file, or retrieve from web.
    Output fasta file with one record per file
    or all records in one file
    """
    p = OptionParser(tofasta.__doc__)
    p.add_option("--prefix", default="gbfasta",
            help="prefix of output files [default: %default]")
    filenames, accessions, idfile, opts, args = preparegb(p, args)
    prefix = opts.prefix

    GenBank(filenames=filenames, accessions=accessions, idfile=idfile).\
        write_fasta(output=prefix, individual=opts.individual)

    if opts.individual:
        logging.debug("Output written dir {0}".format(prefix))
    else:
        logging.debug("Output written to {0}.fasta".format(prefix))


def getgenes(args):
    """
    %prog getgenes [--options]

    Read GenBank file, or retrieve from web.
    Output bed, cds files, and pep file (can turn off with --nopep).
    Either --gb_dir or --id/--simple should be provided.
    """
    p = OptionParser(getgenes.__doc__)
    p.add_option("--prefix", default="gbout",
            help="prefix of output files [default: %default]")
    p.add_option("--nopep", default=False, action="store_true",
            help="Only get cds and bed, no pep [default: %default]")
    p.add_option("--emit_cds", default=False, action="store_true",
            help="only get genes, no cds or pep [default: %default]")
    filenames, accessions, idfile, opts, args = preparegb(p, args)
    prefix = opts.prefix

    GenBank(filenames=filenames, accessions=accessions, idfile=idfile).\
        write_genes(output=prefix, individual=opts.individual, \
        pep=(not opts.nopep), emit_cds=opts.emit_cds)

    if opts.individual:
        logging.debug("Output written dir {0}".format(prefix))
    elif opts.nopep:
        logging.debug("Output written to {0}.bed, {0}.cds".format(prefix,))
    else:
        logging.debug("Output written to {0}.bed, {0}.cds, {0}.pep".format(prefix,))


if __name__ == '__main__':
    main()

# $ python genbank_with_names.py getgenes --simple=NC006273