/*
 * Copyright (c) 2015. PathOS Variant Curation System. All rights reserved.
 *
 * Organisation: Peter MacCallum Cancer Centre
 * Author: Ken Doig
 */


package org.petermac.pathos.pipeline

import groovy.util.logging.Log4j
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.SAMFileWriter
import htsjdk.samtools.SAMFileWriterFactory
import htsjdk.samtools.SAMReadGroupRecord
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.samtools.fastq.FastqRecord
import htsjdk.samtools.fastq.FastqWriter
import htsjdk.samtools.fastq.FastqWriterFactory
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.petermac.annotate.MyVariant
import org.petermac.util.Fasta
import org.petermac.util.Locator
import org.petermac.util.SmithWaterman
import org.petermac.util.Vcf
import org.petermac.util.Vcf2Tsv

import java.util.regex.Pattern

/**
 * Canary is an ultra fast diagnostic report renderer which works directly
 * from FASTQ files to extract a VCF of annotated variants.
 *
 * Dependencies: Uses the fast striped Smith-Waterman C/C++ library accessed via JNI and a
 * dynamic linked library
 *
 * Author:  Kenneth Doig
 * Date:    14-Mar-2015 to 07-Apr-2015
 * Changes:
 *
 * 10   kdd 24-Nov-16       Added semantic versioning - removed IGVH mode
 */

@Log4j
class Canary
{
    private SmithWaterman           sw        = new SmithWaterman()
    private static File             debugfile = new File( "debug.log" )
    private static boolean          debug     = false
    private static boolean          complex   = false
    private static def              startTime = System.currentTimeMillis()
    private static int              flank     = 5            // size of flanking region for regexs
    private static int              min_pairs = 10           // minimum number of read pairs for variants
    private static Double           min_vaf   = 3.0          // minimum VAF
    private static int              maxMutations = 10        // maximum mutations for an alignment
    private static int              maxComplexSize = 30      // Maximum size of complex mutations
    private static int              maxMnpGap      = 15      // Maximum size of inter mutation gap for complex mutations

    //  Read pair cache
    //
    private static boolean          useCache  = true
    private static Map              cache     = [:]
    private static int              cacheQrys = 0
    private static int              cacheHits = 0

    //  Constants Todo: make arguments
    //
    private final static String     VERSION             = '1.0.0' // Canary version number
    private final static Double     COMPRESS_RATIO      = 4.8     // approximate gzip FASTQ compression ratio
    private final static int        HOMO_RUNS           = 6       // minimum size of homopolymer runs to flag
    private final static int        MIN_OVERLAP         = 10      // minimum overlap required between read pairs
    private final static int        PRIMER_BASES        = 5       // extra primer bases to add to alignment. This allows
                                                                  // for bases modified as first or second base of actual read
    //  BAM file writer (optional)
    //
    private static SAMFileWriter bamWriter = null
    private static SAMFileHeader sfh       = null

    //  FASTQ file writers (optional)
    //
    private static FastqWriter fq1Writer   = null
    private static FastqWriter fq2Writer   = null

    /**
     * main method for CLI execution
     *
     * @param args
     */
    static void main( String[] args )
    {
        //	Collect and parse command line args
        //
        def cli = new CliBuilder(   usage: "Canary [options] read1.fastq.gz read2.fastq.gz",
                                    header: '\nAvailable options (use -h for help):\n',
                                    footer: '\nExtract variants directly from FASTQ files\n')

        //	Options to command
        //
        cli.with
        {
            h(      longOpt: 'help',		        'This help message' )
            d(      longOpt: 'debug',		        'Turn on debugging (Note: will generate large file of alignments [debug.out])' )
            a(      longOpt: 'amplicon',   args: 1, required: true, 'Amplicon FASTA file [required]' )
            p(      longOpt: 'primers',    args: 1, required: true, 'Amplicon Primers file [required]' )
            v(      longOpt: 'vcf',        args: 1, 'Found variants VCF file [canary.vcf]' )
            o(      longOpt: 'output',     args: 1, 'Output report file' )
            f(      longOpt: 'flank',      args: 1, 'Size of flanking region (bp) [5]' )
            fastq(  longOpt: 'fastq',      args: 1, 'Optional FASTQ output files prefix' )
            maxmut( longOpt: 'maxmut',     args: 1, 'Maximum number of mutations allowed per read pair [10]' )
            mut(    longOpt: 'mutalyzer',  args: 1, 'Mutalyzer annotation server host [https://mutalyzer.nl]' )
            r(      longOpt: 'reads',      args: 1, 'Percent of reads to process [100]' )
            c(      longOpt: 'complex',             'Coalesce complex events aka MNPs' )
            mnpmax( longOpt: 'mnpmax',     args: 1, 'Maximum size of complex mutations [30]' )
            mnpgap( longOpt: 'mnpgap',     args: 1, 'Maximum size of inter mutation gap for complex mutations [15]' )
            n(      longOpt: 'nocache',             'Dont use read cache' )
            ver(    longOpt: 'version',             'Display Canary version and exit' )
            b(      longOpt: 'bam',        args: 1, 'Optional BAM file of alignment' )
            minpair(longOpt: 'minpair',    args: 1, 'Min read pairs for variants [10]' )
            vaf(    longOpt: 'vaf',        args: 1, 'Minimum VAF for variants [3.0%]' )
            filt(   longOpt: 'filter',     args: 1, 'List of comma separated amplicon names to use' )
            t(      longOpt: 'tsv',        args: 1, 'TSV (Tab separated variable) file of VCF output' )
            cols(   longOpt: 'columns',    args: 1, 'File of VCF field names to output to TSV (one per line with optional alias after comma)' )
            norm(   longOpt: 'normalise',  args: 1, 'Generates annotated VCF file from VCF output' )
            ts(     longOpt: 'transcript', args: 1, 'File of transcripts mapping genes -> refseq (without version)' )
            ano(    longOpt: 'annotation', args: 1, 'File of MyVariant annotation fields' )
        }

        //  Check for version option first - ignore option parsing
        //
        if ( args.size() > 0 && (args[0] =~ '-ver' || args[0] =~ '--version' ))
        {
            println version()
            return
        }

        def opt = cli.parse( args )
        if ( ! opt ) return

        if ( opt.help || opt.arguments().size() != 2 )
        {
            cli.usage()
            return
        }
        log.info( "Canary ${args}")

        //  Debug ?
        //
        if ( opt.debug )
        {
            debug = true
            Logger.getRootLogger().setLevel(Level.DEBUG)
        }
        log.debug( "Debugging turned on!" )

        //  Switch off cache mode
        //
        if ( opt.nocache ) useCache = false;

        //  Switch on complex variant mode
        //
        if ( opt.complex ) complex = true;

        //  Load variants to look for
        //
        List<List> rows = []
        if ( opt.action )
        {
            File vfile = new File(opt.action as String)
            if ( ! vfile.exists())
            {
                log.fatal( "Missing actionable variants VCF file: ${opt.action}" )
                return
            }

            //  Load in VCF file
            //
            Vcf vcf = new Vcf( vfile )
            vcf.load()
            rows = vcf.getRows()
        }

        //  Load amplicons matching FASTQ
        //
        File ampfile = new File(opt.amplicon as String)
        if ( ! ampfile.exists())
        {
            log.fatal( "Missing amplicon FASTA file: ${opt.amplicon}" )
            return
        }

        //  Load primers from file
        //
        File primfile = new File(opt.primers as String)
        if ( ! primfile.exists())
        {
            log.fatal( "Missing amplicon Primers file: ${opt.primers}" )
            return
        }

        //  Open output report file
        //
        String of  = ( opt.output ? opt.output : 'canary.out' )
        File ofile = new File(of)

        //  Set flanking region
        //
        if ( opt.flank )
        {
            try
            {
                flank = Integer.parseInt(opt.flank)
            }
            catch ( Exception e )
            {
                log.fatal( "Number formatting error: ${opt.flank} ${e}")
                return
            }
        }

        //  Set maximum mutations
        //
        if ( opt.maxmut )
        {
            try
            {
                maxMutations = Integer.parseInt(opt.maxmut)
            }
            catch ( Exception e )
            {
                log.fatal( "Number formatting error: ${opt.maxmut} ${e}")
                return
            }
        }

        //  Set minimum read pairs for a variant
        //
        if ( opt.minpair )
        {
            try
            {
                min_pairs = Integer.parseInt(opt.minpair)
            }
            catch ( Exception e )
            {
                log.fatal( "Number formatting error: ${opt.minpair} ${e}")
                return
            }
        }

        //  Setup maximum reads to process
        //
        Double pctreads = 0.0
        if ( opt.reads )
        {
            String pct = opt.reads
            if ( pct.isDouble())
                pctreads = Double.valueOf( pct )

            if ( ! (0 < pctreads && pctreads <= 100.0))
            {
                log.fatal( "Invalid percent reads to process: ${pctreads}")
                return
            }
            if ( pctreads > 99.5 ) pctreads = 0.0     // process every read - dont approximate
        }

        //  Set minimum variant allele frequency
        //
        if ( opt.vaf )
        {
            String vaf = opt.vaf
            if ( vaf.isDouble())
            {
                min_vaf = Double.valueOf( vaf )
            }
            else
            {
                log.fatal( "Invalid variant allele frequency: ${vaf}")
                return
            }
        }

        //  Set maximum size for MNP mutations
        //
        if ( opt.mnpmax )
        {
            try
            {
                maxComplexSize = Integer.parseInt(opt.mnpmax)
            }
            catch ( Exception e )
            {
                log.fatal( "Number formatting error: ${opt.mnpmax} ${e}")
                return
            }
        }

        //  Set maximum gap size for MNP mutations
        //
        if ( opt.mnpgap )
        {
            try
            {
                maxMnpGap = Integer.parseInt(opt.mnpgap)
            }
            catch ( Exception e )
            {
                log.fatal( "Number formatting error: ${opt.mnpgap} ${e}")
                return
            }
        }

        //  Extract file names
        //  Todo: support multiple pairs of read files from multiple flowcells eg NextSeq, HiSeq
        //  Todo: need to support Bzip2 files also
        //
        List<String> extra = opt.arguments()
        assert extra.size() % 2 == 0, "Input files must be paired ${extra}"

        File r1file = new File(extra[0])
        File r2file = new File(extra[1])
        if ( ! r1file.exists() || ! r2file.exists())
        {
            log.fatal( "Read files don't exist: ${r1file},${r2file}" )
            return
        }

        //  Create output VCF file
        //
        String vcfname = 'canary.vcf'
        if ( opt.vcf ) vcfname = opt.vcf as String
        File vcffile = new File( vcfname )
        if ( vcffile.exists())
            vcffile.delete()

        //  Create BAM file of alignment
        //
        if ( opt.bam )
        {
            bamWriter = openBam( opt.bam, '15K1234', 'ReadGroup' )
        }

        //  Output FASTQ files of matching reads
        //
        if ( opt.fastq )
        {
            def ff = new FastqWriterFactory()
            String prefix = opt.fastq as String
            fq1Writer = ff.newWriter( new File( prefix + "_R1.fastq"))
            fq2Writer = ff.newWriter( new File( prefix + "_R2.fastq"))
        }

        //  Check we have a transcript file if normalising
        //
        if ( opt.normalise && ! opt.transcript )
        {
            log.fatal( "--normalise option requires --transcript option")
            System.exit(1)
        }

        //  Check we have a genome if we are normlising variants
        //
        boolean haveGenome = false
        if ( opt.normalise )
        {
            Locator loc = Locator.instance

            File genome = new File( loc.genomePath )

            if ( loc.genomePath && genome.exists())
            {
                haveGenome = true
            }
            else
            {
                log.error( "--normalise option requires genome files at ${loc.genomePath}, won't perform normalisation")
            }
        }

        Map tsMap = null
        if ( opt.transcript )
        {
            File tsf = new File( opt.transcript as String )
            if ( ! tsf.exists())
            {
                log.fatal( "Transcript file doesn't exist ${opt.transcript}")
                System.exit(1)
            }

            tsMap = NormaliseVcf.loadTranscripts( tsf )
            if ( tsMap.size() < 1 )
            {
                log.fatal( "No Transcripts found in file ${opt.transcript}")
                System.exit(1)
            }

            log.info( "Loaded ${tsMap.size()} gene/transcripts from ${tsf}")
        }

        //  Annotation from MyVariants
        //
        List anoFields = []
        if ( opt.annotation )
        {
            File anof = new File( opt.annotation as String )
            if ( ! anof.exists())
            {
                log.error( "Can't read annotation file: ${opt.annotation}")
                System.exit(1)
            }
            anoFields = anof.readLines()
        }

        //  Run it all
        //
        int nlines = new Canary().runCanary( ampfile, opt.filter ?: '', primfile, rows, ofile, vcffile, r1file, r2file, pctreads, anoFields )

        //  Close BAM/FASTQ files
        //
        if ( bamWriter ) bamWriter.close()
        if ( fq1Writer ) fq1Writer.close()
        if ( fq2Writer ) fq2Writer.close()

        //  Set TSV VCF input file
        //
        File tsvVcfFile = vcffile

        //  Create optional Normalised/Annotated VCF file of output
        //
        if ( opt.normalise && nlines && haveGenome )
        {
            File normf  = new File( opt.normalise as String )
            normf.delete()

            //  Annotate with HGVS nomenclature and 3' shift all variants in the VCF
            //
            NormaliseVcf.normaliseVcfFile( vcffile, normf, tsMap, opt.mutalyzer ?: 'https://mutalyzer.nl' )

            //  use the annotated VCF if we are outputting a TSV file
            //
            tsvVcfFile = normf
        }

        //  Create optional TSV file from VCF output
        //
        if ( opt.tsv && nlines )
        {
            File tsvf  = new File( opt.tsv as String )
            tsvf.delete()

            //  Check for optional TSV column file
            //
            File colsf = null
            if ( opt.columns )
            {
                colsf = new File( opt.columns as String )
                if ( ! colsf.exists())
                {
                    log.error( "Columns file doesn't exist: ${opt.columns}")
                    colsf = null
                }
            }

            //  Convert VCF file to TSV table
            //
            Vcf2Tsv.vcf2Tsv( tsvVcfFile, tsvf, '', '', '', colsf, true )
        }

        Integer elapsed = (System.currentTimeMillis() - startTime) / 1000    // seconds elapsed
        log.info( "Done: processed ${8*nlines} lines ${nlines} read pairs into ${vcffile} in ${elapsed} seconds" )
    }

    /**
     * Run the canary in a coalmine
     *
     * @param ampfile   Amplicon FASTA file
     * @param amplist   List of Amplicon names to use
     * @param primfile  Amplicon primer file
     * @param vars      List of VCF variants
     * @param ofile     Report output file
     * @param vcffile   VCF File for output
     * @param r1        FASTQ file Read1
     * @param r2        FASTQ file Read2
     * @param pct       Number of lines processed
     * @return
     */
    int runCanary( File ampfile, String amplist, File primfile, List<List> vars, File ofile, File vcffile, File r1, File r2, Double pct, List anoFields )
    {
        //  Calculate read file stats
        //
        Map    rs = new Fasta( r1, true ).readStats()           // get read stats from FASTQ
        Double er = rs.fileSize * COMPRESS_RATIO / rs.readSize  // estimate the total no of reads in file
        Double maxReads = pct * er / 100.0                      // estimate of maximum num reads to process (0==all)
        log.info( "Found Reads file: readlen=${rs.readLen} filesize=${rs.fileSize} estimated reads=${String.format("%.0f",er)}")

        //  Load amplicons
        //
        List<Map> amplicons = loadAmplicons( ampfile, amplist, primfile, rs.readLen )

        //  Stream gzipped FASTQ files
        //
        Map stats = stream( amplicons, r1, r2, maxReads )

        //  Format stats for output
        //
        List<Map> variants = formatStats( stats, amplicons )
        log.info( "Found ${variants.size()} events" )

        //  Combine variants covered by multiple Amplicons
        //
        variants = combineVariants( variants, amplicons )
        log.info( "Found ${variants.size()} variants" )

        //  Output collected stats
        //
        outputStats( variants, stats, ofile, pct )

        //  Output a VCF file of variants
        //
        outputVcf( variants, vcffile, anoFields )

        return stats.nline
    }

    /**
     * Load in Amplicons, bases primers and genomic position
     *
     * @param   amplicons     Amplicon FASTA file
     * @param   amplist       Comma seperated list of amplicon names to use (empty==ALL)
     * @param   primers       Amplicon primer file
     * @return                List of Maps [pos: fwdlen: revlen: name: bases: ]
     */
    static List<Map> loadAmplicons( File amplicons, String amplist, File primers, int readLen )
    {
        List<Map> amplst = []
        List<Map> offamp = []
        List     ampList = amplist ? amplist.tokenize(',') : []
        Fasta fa = new Fasta(amplicons)
        List<Map> amps = fa.load()
        log.info( "Found ${amps.size()} Amplicons")

        //  Load in primers and match with FASTA amplicon bases
        //  File is in TSV format with rows = [ chr:start-end, <fwd primer len>, <rev primer len>, <name> ]
        //
        primers.splitEachLine('\t')
        {
            column ->
            if ( column.size() == 4 )
            {
                String pos     = column[0]
                String ampname = column[3]
                if ( ! amplist || ampname in ampList )
                {
                    Map a = amps.find { it.head == pos }
                    if ( a?.bases )
                    {
                        String bases = a.bases
                        int fwdlen   = column[1] as int
                        int revlen   = column[2] as int
                        Locus l      = new Locus(pos)        // Amplicon Locus

                        //  Locus for amplicon without primers
                        //
                        Locus lcap = new Locus( l.chr, l.startPos() + fwdlen as int, l.endPos() - revlen as int )

                        //  Construct an amplicon Map for processing reads
                        //
                        Map amp =   [
                                chr:        l.chr,
                                pos:        l.startPos(),
                                locus:      lcap,
                                fwdlen:     fwdlen,
                                fwdprimer:  bases[0..(fwdlen-1)],
                                revlen:     revlen,
                                revprimer:  Locus.revcom(bases[-revlen..-1]),
                                name:       ampname,
                                bases:      bases[(fwdlen-PRIMER_BASES)..(-revlen-1+PRIMER_BASES)],         // bases w/o primer
                                overlap:    2 * readLen - bases.length()        // bases overlapping bw reads
                        ]

                        //  Make sure amplicon names are unique - they are used to key variants
                        //
                        if ( amplst.name.contains(ampname))
                        {
                            log.fatal( "Duplicate amplicon name [${ampname}] exiting")
                            System.exit(1)
                        }

                        //  Ignore off target amplicons - prefixed with "Off"
                        //
                        if ( amp.name.startsWith('Off'))
                        {
                            log.debug( "Found off target amplicon ${amp.name} overlap=${amp.overlap}")
                            offamp << amp
                        }
                        else
                        {
                            amplst << amp
                            log.debug( "${amp}" )
                            if ( amp.overlap < MIN_OVERLAP ) log.warn( "Small overlap between amplicon read pair: ${amp.overlap} in ${amp.name}")
                        }
                    }
                    else
                        log.error( "Couldn't find ${pos} in Amplicons")
                }
            }
        }

        //  Find all the off target ammplicons
        //
        amplst = findOffTarget( amplst, offamp )

        log.info( "Using ${amplst.size()} Amplicons")

        return amplst
    }

    /**
     * Process the off target amplicons and add to regular amplicons
     *
     * @param amps  List of amplicons to enrich
     * @param offts List of off target amplicons
     * @return      List of enriched amplicons
     */
    static List<Map> findOffTarget( List<Map> amps, List<Map> offts )
    {
        List<Map> oamps = []

        //  Loop through the amplicons looking for matching Off targets
        //
        for ( amp in amps )
        {
            //  Look for a matching off target amplicon
            //
            List<Map> offlst = []

            for ( offt in offts )
            {
                log.debug( "Matching amp=[${amp.fwdprimer}:${amp.revprimer}] with offt=[${offt.fwdprimer}:${offt.revprimer}]")
                boolean fwd = ( offt.fwdprimer == amp.fwdprimer || offt.revprimer == amp.fwdprimer )
                boolean rev = ( offt.revprimer == amp.revprimer || offt.fwdprimer == amp.revprimer )
                if ( fwd && rev && offt.fwdprimer == amp.fwdprimer )
                {
                    log.info( "${amp.name} has ${fwd?'fwd':''} ${rev?'rev':''} matches to ${offt.name}")
                    offlst << offt
                }

                //  Have a match but with the primers on the other strand
                //
                if ( fwd && rev && offt.fwdprimer != amp.fwdprimer )
                {
                    log.info( "${amp.name} has ${fwd?'fwd':''} ${rev?'rev':''} matches to ${offt.name}")

                    //  Reverse complement bases, we don't care about other attributes for an off target
                    //
                    offlst << [ bases: Locus.revcom( offt.bases ) ]
                }
            }

            //  Add matching off target amplicons
            //
            amp.offTarget = offlst

            oamps << amp
        }

        return oamps
    }

    /**
     * Load amplicons and map variants to search for
     *
     * @param afile     Amplicon FASTA file
     * @param vars      Variants VCF file
     * @return          List of amplicon/variant Maps
     */
    static List<Map> loadVariants( List<Map> amplicons, List<List> vars, int readlen )
    {
        List amplist = []

        //  Loop through variants
        //
        log.info( "Found ${vars.size()} variants to search for")
        List<Map> varloci = []
        for ( List<String> row in vars )
        {
            //  Convert a VCF row into a Map
            //  Map of converted variant [chr:, pos:, ref:, alt:, ensvar: "chr_pos_ref/alt", hgvsg: ]
            //
            Map var = HGVS.normaliseVcfVar( row[0], row[1], row[3], row[4] )
            log.info( "Variant: ${var.hgvsg}")

            def l = new Locus( var.chr, var.pos as String, var.endpos as String )
            Map m = new HashMap(var)
            m.locus = l
            varloci << m
        }

        //  Map variants to Amplicons
        //
        for ( amp in amplicons )
        {
            //  Ignore off target amplicons
            //
            if ( amp.name.startsWith('Off'))
            {
                log.warn( "Ignoring off target amplicon ${amp.name}")
                continue
            }

            //  Locus for amplicon with primers
            //
                Locus lamp = new Locus( amp.pos )

            //  Locus for amplicon without primers
            //
            Locus lcap = new Locus( lamp.chr, lamp.startPos() + amp.fwdlen as Integer, lamp.endPos() - amp.revlen as Integer )

            //  Look for Amplicons overlapping
            //
            for ( v in varloci )
            {
                if ( lcap.contains(v.locus))
                {
                    int  fwdoff = v.locus.startPos() - lamp.startPos()
                    if ( fwdoff > readlen ) fwdoff = 0    // ignore variants not in actual read
                    int  revoff = lamp.endPos() - v.locus.startPos()
                    if ( revoff > readlen ) revoff = 0    // ignore variants not in actual read
                    if ( fwdoff==0 && revoff==0 ) continue

                    Map vamp =  [
                                var:    v.hgvsg,
                                amp:    amp.name,
                                chr:    lamp.chr,
                                pos:    lamp.startPos(),
                                bases:  amp.bases,
                                ref:    v.ref,
                                alt:    v.alt,
                                fwd:    [primer: amp.bases[0..(amp.fwdlen-1)],             offset: fwdoff, regex: regex(amp.bases, fwdoff, true,  v)],
                                rev:    [primer: Locus.revcom(amp.bases[-amp.revlen..-1]), offset: revoff, regex: regex(Locus.revcom(amp.bases), revoff, false, v)]
                                ]

                    amplist << vamp

                    log.debug( "${vamp}" )
                }
            }
        }

        return amplist
    }

    /**
     * Calculate and compile a regex for this variant and amplicon
     *
     * @param bases
     * @param dir
     * @param var
     * @return
     */
    static Pattern regex( String bases, int offset, boolean dir, Map var )
    {
        String re

        def ref = dir ? var.ref : Locus.revcom(var.ref)
        def alt = dir ? var.alt : Locus.revcom(var.alt)

        //  Insertion
        //
        if ( var.ref == '-' )
        {
            re = bases[(offset-flank)..(offset-1)] + "(${alt})?" + bases[offset..(offset+flank-1)]
        }
        //  Deletion
        //
        else if ( var.alt == '-' )
        {
            if ( dir )
                re = bases[(offset-flank)..(offset-1)] + "(${ref})?" + bases[(offset+ref.length())..(offset+ref.length()+flank-1)]
            else
                re = bases[(offset-ref.length()-flank+1)..(offset-ref.length())] + "(${ref})?" + bases[(offset+1)..(offset+flank)]
        }
        //  SNP
        //
        else
            re = bases[(offset-flank)..(offset-1)] + "(${ref}|${alt})" + bases[(offset+1)..(offset+flank)]

        return Pattern.compile( re )
    }

    /**
     * Process a streaming gzipped FASTQ file from STDIN
     *
     * @param amps      Amplicons to search
     * @param r1        Read1 gzipped file
     * @param r2        Read2 gzipped file
     * @param maxReads  Maxmimum number of reads to process
     * @return          Processed stats Map
     */
    Map stream( List<Map> amps, File r1, File r2, maxReads )
    {
        //  Initialise stuff for processing reads
        //
        int npairs = 0
        Map stats = [nline:0, mapped:0, unmapped:0, unpaired:0, offTarget:0 ]

        //  Fasta file classes for paired reads
        //
        Fasta r1fq = new Fasta( r1, true )
        Fasta r2fq = new Fasta( r2, true )

        //  Reset unmatched reads file
        //
        if ( debug ) debugfile.delete()

        //  Loop through fwd and rev reads simultaneously
        //
        Map read1, read2
        while ((read1 = r1fq.getRead()) != null && (read2 = r2fq.getRead()) != null)
        {
            ++npairs

            //  Process the reads using Smith Waterman alignment
            //
            stats = processReadSW( amps, read1, read2, stats )

            if ( maxReads && npairs > maxReads )
            {
                log.warn( "Maximum estimated read limit reached: ${npairs-1} read pairs")
                break
            }
        }

        stats.nline = npairs
        return stats
    }

    /**
     * Process a pair of amplicon reads by finding their primers and identifying all variants
     * Uses the Smith Waterman library to align the read pairs and then align pair to reference
     * Collect the stats for each mutation for reporting
     *
     * @param amps      List of Maps of amplicons to search
     * @param read1     Read 1 Map [ head:, bases:, quals: ]
     * @param read2     Read 2 Map [ head:, bases:, quals: ]
     * @param stats     Accumulated statistics before this read pair
     * @return          Accumulated statistics after  this read pair
     */
    Map processReadSW( List<Map> amps, Map read1, Map read2, Map stats )
    {
        String  readStatus  = 'unmapped'
        boolean written     = false

        //  Check each amplicon/variant against read
        //
        for ( amp in amps )
        {
            Boolean fwd = false
            Boolean rev = false
            Map fwdread = new HashMap(read1)
            Map revread = new HashMap(read2)

            //  Look for matching primers
            //
            if ( matchPrimer( fwdread.bases[0..amp.fwdlen], amp.fwdprimer ))
            {
                //  We have a match on first read
                //
                fwd = true
                rev = matchPrimer( revread.bases[0..amp.revlen], amp.revprimer )
            }
            else
            {
                //  Try other direction
                //
                if ( matchPrimer( fwdread.bases[0..amp.revlen], amp.revprimer ))
                {
                    //  We have a match on first read in reverse direction
                    //
                    fwd = true
                    if ( matchPrimer( revread.bases[0..amp.fwdlen], amp.fwdprimer ))
                    {
                        //  We have a match on both reads
                        //
                        rev = true

                        //  Swap fwd and rev reads so fwdbases matches amp.fwdprimer
                        //
                        Map tmpread = fwdread
                        fwdread     = revread
                        revread     = tmpread
                    }
                }
            }

            //  We must have at least one matching primer
            //
            if ( ! fwd && ! rev ) continue

            if ( ! fwd || ! rev )
            {
                //  Only one primer matched - unpaired reads
                //
                readStatus  = 'unpaired'
                if ( debug )
                {
                    Map res = sw.align(  revread.bases, fwdread.bases )
                    Map fmt = sw.format( revread.bases, fwdread.bases, res )
                    debugfile << "## Unpaired primers\n"
                    debugfile << "fwd: " + fwdread.bases    + '\n'
                    debugfile << "rev: " + revread.bases    + '\n'
                    debugfile << "fqc: " + fwdread.quals    + '\n'
                    debugfile << "rqc: " + revread.quals    + '\n'
                    debugfile << "${res} snp=${fmt.snps} ins=${fmt.ins} dels=${fmt.dels}\n${fmt.ref}\n${fmt.align}\n${fmt.qry}\n"
                }

                continue    // try next amplicon
            }

            //  Reverse second read to be on fwd strand for alignment
            //
            revread.bases = Locus.revcom(revread.bases)

            //  Perform alignment if we have both paired reads matching the amplicon primers
            //
            readStatus  = 'mapped'

            //  Add up how many times we've seen this amplicon
            //
            if ( ! stats[amp.name] ) stats[amp.name] = 0
            ++stats[amp.name]

            //  Perform Alignment
            //
            List<Map> muts = performAligment( amp, fwdread, revread, stats )

            for ( mut in muts )
            {
                //  Set up stats to collect data
                //
                def name = "${amp.name}-${mut.pos}-${mut.ref}-${mut.alt}"
                def pos  = amp.pos+amp.fwdlen+mut.pos-PRIMER_BASES   // genomic position of mutation
                if ( ! stats[name] )
                    stats[name] = [amp: amp.name, chr: amp.chr, pos: pos, ref: mut.ref, alt: mut.alt, cnt:0 , fsRescue:0 ]

                Map m = stats[name] as Map
                ++m.cnt
                if ( mut.complex ) m.complex = true

                if ( debug ) debugfile << "${m}\n"
            }

            //  Output reads to mapped FASTQ files
            //
            if ( fq1Writer && ! written )
            {
                fq1Writer.write( new FastqRecord( read1.head, read1.bases, '', read1.quals ))
                fq2Writer.write( new FastqRecord( read2.head, read2.bases, '', read2.quals ))
                written = true
            }

         } // amps loop

        //  Keep counts of read types
        //
        if ( readStatus == 'unpaired' ) ++stats.unpaired
        if ( readStatus == 'unmapped' ) ++stats.unmapped
        if ( readStatus == 'mapped'   ) ++stats.mapped

        return stats
    }

    /**
     * Create a BAM file from alignments
     * Todo: Set mapping quality
     *
     * @param bamFile
     * @param sample
     */
    private static SAMRecord outputBamRead( String chr, Integer pos, String head, String bases, String quals, String cigar )
    {
        /*  SAM file attributes
         
            QNAME	M01053:327:000000000-AN214:1:1101:10420:11589 2:N:0:GTGAGAGACA_Consensus
            FLAG	16  				// SEQ being reverse complemented
            RNAME	1					// chromosome
            POS		36932056			// genomic pos
            MAPQ	60					// MAPping Quality
            CIGAR	163M				// Cigar string
            RNEXT	*					// Ref. name of the mate/next read
            PNEXT	0					// Position of the mate/next read
            TLEN	0					// template length
            SEQ		GGGCTGG...AAAG		// read
            QUAL	HFGGGGH...HHGH		// Phred qual
         */

        def samRecord = new SAMRecord( sfh )
        samRecord.setReadName( head )
        samRecord.setFlags( 16 )
        samRecord.setReferenceName( chr )           // chromosome
        samRecord.setAlignmentStart( pos )          // genomic pos
        samRecord.setMappingQuality( 60 )           // constant quality for the moment
        samRecord.setCigarString( cigar )
        samRecord.setReadString( bases )
        samRecord.setBaseQualityString( quals )
        //samRecord.setAttribute( 'RG', 'ReadGroup' )
        //samRecord.setAttribute( 'XR', 'CSF3R_EX17_2' )
        
        //  Add this read to BAM
        //
        bamWriter.addAlignment( samRecord )

        return samRecord
    }

    /**
     * Criteria for matching a read with a primer
     * This is fundamental to speed as there are A * R tests where A=amplicons, R=read pairs
     *
     * @param read
     * @param primer
     * @return
     */
    private static Boolean matchPrimer( String read, String primer )
    {
        //  Try exact match
        //
        if ( read == primer ) return true

        //  Try to match the middle bases of primer - this is a tradeoff between
        //  catching every matching read and being efficient
        //
        return read.contains( primer[6..-6] )
    }

    /**
     * Perform alignment and variant calling on read pair
     *
     * @param amp       Map of amplicon to search
     * @param read1     Read 1 Map [ head:, bases:, quals: ]
     * @param read2     Read 2 Map [ head:, bases:, quals: ]
     * @param stats     Map of global stats
     * @return          List of Maps of variants
     */
    List<Map> performAligment( Map amp, Map read1, Map read2, Map stats )
    {
        //  strip off primers
        //
        String ampbases = amp.bases

        //  Trim off primer bases from reads except for a flanking region of PRIMER_BASES bases
        //
        String bases1   = read1.bases[(amp.fwdlen-PRIMER_BASES)..-1]
        String bases2   = read2.bases[0..(-amp.revlen-1+PRIMER_BASES)]
        String qual1    = read1.quals[(amp.fwdlen-PRIMER_BASES)..-1]
        String qual2    = read2.quals[0..(-amp.revlen-1+PRIMER_BASES)]

        //  Look up cache
        //
        ++cacheQrys
        def pair = new AlignPair( ampliconName: amp.name, read1: bases1, read2: bases2 )
        if ( useCache && cache[pair] )
        {
            ++cacheHits
            return cache[pair] as List<Map>
        }

        //  log the inputs to the alignment
        //
        if ( debug )
        {
            debugfile << "## Raw reads\n"
            debugfile << "amp: " + ampbases  + '\n'
            debugfile << "fwd: " + bases1    + '\n'
            debugfile << "rev: " + bases2    + '\n'
            debugfile << "fqc: " + qual1     + '\n'
            debugfile << "rqc: " + qual2     + '\n'
        }

        //  Align the read pair to each other - Todo: dont assume an overlap !
        //
        Map res = sw.align(  bases2, bases1 )
        Map fmt = sw.format( bases2, bases1, res )
        if ( debug ) debugfile << "${res} snp=${fmt.snps} ins=${fmt.ins} dels=${fmt.dels}\n${fmt.ref}\n${fmt.align}\n${fmt.qry}\n"

        //  Merge pair-aligned read pair to create a single read
        //
        String mergedRead = sw.mergePair( fmt )

        //  Align merged read to amplicon
        //
        res = sw.align(  mergedRead, ampbases )
        fmt = sw.format( mergedRead, ampbases, res )

        //  Map res = [score:324, cigar:162M, refStart:0, refEnd:161, qryStart:117, qryEnd:278] snp=0 ins=0 dels=0
        //
        if ( debug ) debugfile << "## Ref alignment\n${res} snp=${fmt.snps} ins=${fmt.ins} dels=${fmt.dels}\n${fmt.ref}\n${fmt.align}\n${fmt.qry}\n"

        //  Check for poor quality alignment or a better Off target alignment
        //
        int  changes    = fmt.snps + fmt.dels + fmt.ins                 // number of changes in alignment SNP or indels
        int  offChanges = checkOffTarget( mergedRead, amp.offTarget )   // defaults to maxMutations

        if ( changes > offChanges )
        {
            if ( debug )
            {
                if ( amp.offTarget )
                    debugfile << "Off target alignment better ${offChanges} than ${changes} for ${amp.name}\n"
                else
                    debugfile << "Ignoring crap alignment: too many changes: ${changes} for ${amp.name}\n"
            }

            if ( useCache ) cache[pair] = []
            ++stats.offTarget
            return []
        }

        //  Look for frameshifts
        //
        int dels = ( fmt.ref =~ /-/ ).count      // count deletions  in reference
        int ins  = ( fmt.qry =~ /-/ ).count      // count insertions in reference
        boolean fsRescue = ( ins && dels && ((ins-dels) % 3 == 0))   // are we still in frame ?

        //  Find any SNPs or indels
        //
        List<Map> muts = sw.variants( fmt )
        List<Map> mmap = []

        for ( mut in muts )
        {
            //  Ignore read pair discordant bases
            //
            if ( mut.alt.endsWith('N')) continue

            //  Check we are in the amplicon
            //
            if ( mut.pos >= ampbases.length()) continue

            mut.fsRescue = fsRescue

            mmap << mut
        }

        //  create complex variant if we have at least two muts per alignment
        //
        if ( complex && mmap.size() > 1 )
        {
            Map cplx = sw.complex( fmt, maxComplexSize, maxMnpGap )
            if ( cplx )
            {
                if ( debug ) debugfile << "## Complex MNP ${cplx.ref} -> ${cplx.alt} pos=${cplx.pos}\n"
                mmap << cplx
            }
        }

        //  Write to BAM file if required
        //
        if ( bamWriter )
        {
            //  Dummy quals for the moment Todo: get correct quals by merging quals like reads
            //
            String quals = (read1.quals + read2.quals)[0..(mergedRead.length()-1)]

            int pos = amp.pos + amp.fwdlen - PRIMER_BASES + res.refStart

            //  Write the merged reads as an alignment to amplicon
            //
            outputBamRead( amp.chr, pos, read1.head, mergedRead[(res.qryStart)..(res.qryEnd)], quals[(res.qryStart)..(res.qryEnd)], res.cigar )
        }

        //  Save in cache
        //
        if ( useCache ) cache[pair] = mmap

        return mmap
    }

    /**
     * Perform alignment of merged read pair with Off target amplicon
     *
     * @param mergedPair    Merged read pair
     * @param offts         List of Maps of off target amplicons to search
     * @return              best alignment of off targets or maxMutations
     */
    private int checkOffTarget( String mergedRead, List<Map> offts )
    {
        int best = maxMutations

        for ( offt in offts )
        {
            String ampbases = offt.bases

            //  Align merged read to off target amplicon
            //
            Map res = sw.align(  mergedRead, ampbases )
            Map fmt = sw.format( mergedRead, ampbases, res )
            if ( debug ) debugfile << "off: ${res} snp=${fmt.snps} ins=${fmt.ins} dels=${fmt.dels}\noff: ${fmt.ref}\noff: ${fmt.align}\noff: ${fmt.qry}\n"

            //  Check for poor quality alignment
            //
            int changes = fmt.snps + fmt.dels + fmt.ins     // number of changes in alignment SNP or indels
            if ( changes < best ) best = changes
        }

        return best
    }

    /**
     * Format the statistics for output
     *
     * @param stats   Map of processed stats
     * @return        List of variants found
     */
    static List<Map> formatStats( Map stats, List amplicons )
    {
        //  Convert the stats Map into a list of variants with annotations
        //
        List<Map> vars = []
        for ( evt in stats )
        {
            Map m = evt.value as Map
            if ( ! (m instanceof LinkedHashMap)) continue   // Only process variant Maps

            //  Get read count for this variant and amplicon
            //
            m.ampcnt = stats[m.amp] as int

            //  Ignore low VAF variants (less than 2.5%) and fewer than 10 reads
            //
            if ( m.cnt / m.ampcnt < min_vaf / 100.0 ) continue
            if ( m.cnt < min_pairs ) continue

            //  Render variant as HGVSg
            //
            String ref = m.ref
            String alt = m.alt
            String pos = String.valueOf(m.pos)

            m.hgvsg = HGVS.normaliseVcfVar( m.chr, pos, ref, alt )?.hgvsg
            m.pos   = String.valueOf(m.pos)

            //  Get total counts for this variant for all amplicons - not just amplicons with variant
            //
            m.totcnt = 0
            for ( amp in allAmps( m, amplicons ))
            {
                if ( stats[amp] )
                {
                    m.totcnt += stats[amp] as int
                }
                else
                {
                    log.error( "Missing amplicons ${allAmps(m,amplicons)} for variant ${m}")
                }
            }

            log.debug( "AmpVar: ${m.amp} amp:${m.cnt}/[${m.ampcnt}:${m.totcnt}] ${m.hgvsg} ${allAmps(m,amplicons)}" )

            vars << m
        }

        //  Sort variants on chr pos
        //
        def sorted = vars.sort{ String.format( "%5s:%12s", it.chr, it.pos)}

        return sorted
    }

    /**
     * Format the statistics for output
     *
     * @param   vars    List of variants to combine
     * @param   amps    List of amplicons
     * @return          List of combined variants
     */
    static List<Map> combineVariants( List<Map> vars, List<Map> amps )
    {
        List<Map> combined = []
        Map       last = [:]

        for ( var in vars )
        {
            //  Ignore variants with 'N' bases
            //
            if ( var.alt.contains('N')) continue

            //  Look for homopolymer runs
            //
            var = homopolymerRuns( var, amps )

            //  if this var matches previous, merge stats
            //
            if ( last && last.hgvsg == var.hgvsg )
            {
                last.vafs   << var.cnt/var.ampcnt
                last.amp     = last.amp + ',' + var.amp
                last.cnt    += var.cnt
                last.ampcnt += var.ampcnt
                continue
            }

            //  Add combined variant to final list
            //
            if ( last ) combined << varStats( last, amps )

            //  Set up next variant to combine
            //
            last = var
            last.vafs = [var.cnt/var.ampcnt]    // initialise a List of the VAFs
        }

        //  Output last element if any
        //
        if ( last ) combined << varStats( last, amps )

        //  Reject any vars below the threshold after combining
        //
        List finalVars = []
        for ( var in combined )
            if ( var.totcnt && (var.cnt / var.totcnt >= min_vaf / 100.0 ))
                finalVars << var

        return finalVars
    }

    /**
     * Look for homopolymer Runs next to variant
     *
     * @param var
     * @param amps
     * @return
     */
    static Map homopolymerRuns( Map var, List<Map> amps )
    {
        Map amp = amps.find { it.name == var.amp }

        int offset = (var.pos as int) - amp.pos - amp.fwdlen + PRIMER_BASES

        //  Upstream runs of bases
        //
        if ( offset + HOMO_RUNS < amp.bases.length())
        {
            String upstream   = amp.bases[offset+1..(offset+HOMO_RUNS)]
            if (( upstream =~ /${upstream[0]}/ ).count == HOMO_RUNS )
            {
                //println "up:  ${var.ref} ${var.alt} ${upstream}"
                var.homopolymer = upstream[0]
            }
        }

        //  Downstream runs of bases
        //
        if ( offset - HOMO_RUNS >= 0 )
        {
            String downstream = amp.bases[(offset-HOMO_RUNS)..(offset-1)]
            if (( downstream =~ /${downstream[0]}/ ).count == HOMO_RUNS )
            {
                //println "dwn: ${var.ref} ${var.alt} ${downstream}"
                var.homopolymer = downstream[0]
            }
        }

        return var
    }

    /**
     * Calculate the bias in variant across multiple Amplicons
     * This statistic is the variation about the mean VAF
     *
     * bias = sum(vaf - mean_vaf) / sqrt( N )
     *
     * @param   vafs    List of Variant allele frequencies
     * @return          Amplicon Bias statistic
     */
    static Double ampBias( List<Double> vafs )
    {
        if ( ! vafs ) return 0.0

        Double var = 0.0, avg = vafs.sum() / vafs.size()

        for ( vaf in vafs )
            var += Math.abs(avg - vaf)

        return var / Math.sqrt( vafs.size())
    }

    /**
     * Find the occurrences of this variant overlapping with the amplicons
     *
     * @param var   Map of variant
     * @param amps  List of all amplicons
     * @return      List of overlapping amplicons
     */
    private static List allAmps( Map var, List<Map> amps )
    {
        Locus vl = new Locus( var.chr, var.pos as int, (var.pos as int)+var.ref.size()-1 )

        //  Intialise the list with the variant amplicon.
        //  This fixes the bug where a variant is called but is outside the amplicon locus
        //  eg where the variant is an insertion at the beginning of an amplicon, in this case
        //  the locus is the base prior to the insertion
        //
        List allamps = [ var.amp ]
        for ( amp in amps )
            if ( amp.locus.contains(vl))
                allamps << amp.name

        return allamps.unique()
    }

    /**
     * summarise stats for this variant
     *
     * @param var   Map of variant
     * @return      Updated Map
     */
    private static Map varStats( Map var, List<Map> amps )
    {
        //  Set hit ratio for amplicons
        //
        var.numAmps = "${var.vafs.size()}/${allAmps(var,amps).size()}"

        //  Add an extra 0% VAF if not all amplicons contain the variant
        //
        if ( var.ampcnt != var.totcnt ) var.vafs << 0.0

        //  Set amplicon  bias
        //
        var.ampBias = ampBias( var.vafs )

        return var
    }


    /**
     * Output the final statistics
     *
     * @param stats   Map of processed stats
     * @param ofile   Output file
     * @return        List of variants found
     */
    static void outputStats( List<Map> vars, Map stats, File ofile, Double pct )
    {
        ofile.delete()
        if ( ! stats.nline ) return

        def now = new Date()
        Double elapsed = (System.currentTimeMillis() - startTime) / (60*1000)    // minutes elapsed

        //  Calcaulate stats
        //
        def mapPct    = outpct( stats.nline, stats.mapped )
        def unpairPct = outpct( stats.nline, stats.unpaired )
        def unmapPct  = outpct( stats.nline, stats.unmapped )
        def offtPct   = outpct( stats.nline, stats.offTarget )
        def cachePct  = outpct( cacheQrys, cacheHits )

        ofile << "##    Generated by Canary(SW) on ${now}\n"
        ofile << "##\n"
        if ( pct == 0D ) pct = 100.0
        ofile << "##    Processed   ${stats.nline}(${String.format("%.1f",pct)}%) read pairs in ${String.format("%.1f",elapsed)} minutes\n"
        ofile << "##    mapped      ${stats.mapped}/${stats.nline} ${mapPct} %\n"
        ofile << "##    unpaired    ${stats.unpaired}/${stats.nline} ${unpairPct} %\n"
        ofile << "##    unmapped    ${stats.unmapped}/${stats.nline} ${unmapPct} %\n"
        ofile << "##    off target  ${stats.offTarget}/${stats.nline} ${offtPct} %\n"
        ofile << "##    cache hits  ${cacheHits}/${cacheQrys} ${cachePct} %\n"
        ofile << "##\n"
        ofile << "#Amplicon\tHGVSg\tChr\tPosition\tCount\tRef\tAlt\tAlt%\n"

        //  Output the stats for each amplicon: one line per variant/amplicon
        //
        for ( m in vars )
        {
            //  Populate variant list for reporting
            //
            List vals = []
            vals << m.amp
            vals << m.hgvsg
            vals << m.chr
            vals << m.pos
            vals << "(${m.cnt}/${m.totcnt})"
            vals << m.ref
            vals << m.alt
            vals << outpct( m.totcnt, m.cnt )

            ofile << vals.join('\t') + '\n'
        }
    }

    /**
     * Format output by base
     *
     * @param base  Base to format
     * @param m     Stats map
     * @return      <base count> <base %>
     */
    private static String outpct( int tot, int val )
    {
        Float pct = 0
        if ( tot ) pct = val * 100 / tot

        return String.format("%.6f", pct)
    }

    /**
     * Output a VCF file of variants with attributes
     *
     * @param vars
     * @param vcffile
     */
    private static void outputVcf( List<Map> vars, File vcffile, List anoFields )
    {
        //  Use the filename without extension for sample name
        //
        vcffile << Vcf.header( 'Canary', anoFields, vcffile.name.split("\\.", 2)[0] )

        //  Annotate each variant with requested fields
        //
        annotateVars( vars, anoFields )

        //  Output each variant as a VCF line
        //
        for ( var in vars )
        {
            vcffile << vcfLine( var )
        }
    }

    /**
     * Annotate each variant with MyVariant
     *
     * @param vars      List<Map> of variants to annotate var.hgvsg is passed to annotator
     * @param fields    List of fields to find annotations for see http://myvariant.info/v1/api/
     */
    private static void annotateVars( List<Map> vars, List fields )
    {
        if ( ! fields || ! vars ) return

        //  Get slice of maps with variants
        //
        List<String> hgvsgs = vars.hgvsg

        //  Call MyVariant via REST API
        //
        List<Map> myv = MyVariant.submit( hgvsgs, fields )

        //  Get annotations as flattened Maps eg Map<String,String>
        //
        List<Map> annos = MyVariant.flatMaps( myv, '' )
        log.debug( "${fields} ${annos.size()} ${annos}" )

        for ( var in vars )
        {
            //  Look for annotations for this variant
            //
            Map ano = annos.find { Map v -> v.query == var.hgvsg }
            if ( ano && ! ano.notfound )
            {
                //  Save annotation map for VCF INFO field
                //
                var.myv = [:]
                for ( kv in ano )
                {
                    //  Only save fields of interest
                    //
                    if ( kv.key in fields )
                    {
                        var.myv << kv
                    }
                }
            }
            else
            {
                log.warn( "Couldn't find annotations for ${var.hgvsg}")
            }
        }
    }

    /**
     * Create a minimal VCF line for output
     *
     * @param row
     * @return
     */
    static String vcfLine( Map var )
    {
        String freq = outpct( var.totcnt, var.cnt )
        String info = "HGVSg=${infoClean(var.hgvsg)}"
        info += ";numAmps=${var.numAmps}"
        info += ";amps=${infoClean(var.amp)}"
        if ( var.ampBias )          info += ";ampbias=" + String.format( "%.2f", var.ampBias )
        if ( var.fsRescue > 10 )    info += ";fsRescue=${var.fsRescue}"
        if ( var.homopolymer )      info += ";homopolymer=${var.homopolymer}"

        //  Add MyVariant.info annotations
        //
        if ( var.myv )
        {
            for ( kv in var.myv )
            {
                info += ";${infoClean(kv.key)}=${infoClean(kv.value)}"
            }
        }

        //  VCF columns
        //
        List flds =     [
                        var.chr,
                        var.pos as String,
                        '.',
                        var.ref,
                        var.alt,
                        '.',
                        'PASS',
                        info,
                        'GT:DP:RD:AD:FREQ',
                        "0/1:${var.totcnt}:${var.totcnt-var.cnt}:${var.cnt}:${freq}"
                        ]

        return flds.join('\t') + '\n'
    }

    /**
     * Remove all VCF INFO field incompatible characters
     *
     * @param txt   String to sanitise
     * @return      Cleaned String
     */
    static String infoClean( String txt )
    {
        return txt.replaceAll( / |;|=/, '_' )
    }

    //  hg19 Chromosomes and their lengths
    //  Todo: allow user to overwrite for other reference datasets
    //
    private static final Map<String,Integer> chromosomes =  [
                                                                1:249250621, 2:243199373, 3:198022430, 4:191154276, 5:180915260,
                                                                6:171115067, 7:159138663, 8:146364022, 9:141213431, 10:135534747,
                                                                11:135006516, 12:133851895, 13:115169878, 14:107349540, 15:102531392,
                                                                16:90354753, 17:81195210, 18:78077248, 19:59128983, 20:63025520,
                                                                21:48129895, 22:51304566, X:155270560, Y:59373566, MT:16569
                                                            ]

    /**
     * Open a BAM file for writing
     *
     * @param bamFile
     * @param sample
     * @param readGroup
     * @return
     */
    private static SAMFileWriter openBam( String bamFile, String sample, String readGroup )
    {
        //  Create SAM header
        //
        sfh = new SAMFileHeader()
        sfh.setAttribute( 'SO', 'coordinate' )

        //  Create chromosomes
        //
        for ( chr in chromosomes )
        {
            def ssr = new SAMSequenceRecord( chr.key as String, chr.value )
            sfh.addSequence(  ssr )
        }

        //  Create Read Group
        //
        def rg  = new SAMReadGroupRecord( readGroup )
        rg.setAttribute( 'PL', 'Illumina')
        rg.setAttribute( 'SM', sample)
        sfh.addReadGroup( rg )

        return new SAMFileWriterFactory(createIndex: true).makeSAMOrBAMWriter( sfh, false, new File(bamFile))
    }

    /**
     * Version method
     *
     * @return Canary version String
     */
    public static String version()
    {
        return( "Canary ${VERSION}")
    }
}


