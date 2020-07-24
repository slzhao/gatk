package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.GermlineCNVCaller;
import org.broadinstitute.hellbender.tools.copynumber.PostprocessGermlineCNVCalls;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.IntegerCopyNumberSegment;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Helper class for {@link PostprocessGermlineCNVCalls} for single-sample postprocessing of segmented
 * {@link GermlineCNVCaller} calls.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class GermlineCNVSegmentVariantComposer extends GermlineCNVVariantComposer<IntegerCopyNumberSegment> {

    /* VCF FORMAT header keys */

    /**
     * Segment copy-number call
     */
    public static final String CN = "CN";

    /**
     * Number of points in the segment
     */
    public static final String NP = "NP";

    /**
     * Quality metric (some points called)
     */
    public static final String QS = "QS";

    /**
     * Quality metric (all points called)
     */
    public static final String QA = "QA";

    /**
     * Quality metric (segment start)
     */
    public static final String QSS = "QSS";

    /**
     * Quality metric (segment end)
     */
    public static final String QSE = "QSE";

    protected final Logger logger = LogManager.getLogger(this.getClass());

    private final IntegerCopyNumberState refAutosomalCopyNumberState;
    private final Set<String> allosomalContigSet;
    private final ReferenceSequenceFile reference;
    private final int dupeQSThreshold;
    private final int hetDelQSThreshold;
    private final int homDelQSThreshold;
    private final int siteFreqThreshold;
    private final File clusteredCohortVcf;
    private final FeatureReader<VariantContext> clusteredVCFReader;

    /**
     * Constructor.
     *
     * @param outputWriter variant context writer
     * @param sampleName sample name
     * @param refAutosomalCopyNumberState ref copy-number state on autosomal contigs
     * @param allosomalContigSet set of allosomal contigs (ref copy-number allele be chosen according to
     *                           given contig baseline copy-number states)
     * @param reference may be null
     *
     */
    public GermlineCNVSegmentVariantComposer(final VariantContextWriter outputWriter,
                                             final String sampleName,
                                             final IntegerCopyNumberState refAutosomalCopyNumberState,
                                             final Set<String> allosomalContigSet,
                                             final ReferenceSequenceFile reference,
                                             final int dupeQSThreshold,
                                             final int hetDelQSThreshold,
                                             final int homDelQSThreshold,
                                             final int siteFreqThreshold,
                                             final File clusteredCohortVcf) {
        super(outputWriter, sampleName);
        this.refAutosomalCopyNumberState = Utils.nonNull(refAutosomalCopyNumberState);
        this.allosomalContigSet = Utils.nonNull(allosomalContigSet);
        this.reference = reference;
        this.dupeQSThreshold = dupeQSThreshold;
        this.hetDelQSThreshold = hetDelQSThreshold;
        this.homDelQSThreshold = homDelQSThreshold;
        this.siteFreqThreshold = siteFreqThreshold;
        this.clusteredCohortVcf = clusteredCohortVcf;
        if (clusteredCohortVcf != null) {
            try {
                clusteredVCFReader = AbstractFeatureReader.getFeatureReader(clusteredCohortVcf.getAbsolutePath(), new VCFCodec());
            } catch ( final TribbleException ex ) {
                throw new GATKException("Error - IO problem with file " + clusteredCohortVcf, ex);
            }
        } else {
            clusteredVCFReader = null;
        }
    }

    public GermlineCNVSegmentVariantComposer(final VariantContextWriter outputWriter,
                                             final String sampleName,
                                             final IntegerCopyNumberState refAutosomalCopyNumberState,
                                             final Set<String> allosomalContigSet,
                                             final ReferenceSequenceFile reference) {
        this(outputWriter, sampleName, refAutosomalCopyNumberState, allosomalContigSet, reference, 0, 0, 0, 0, null);
    }

    @Override
    public void composeVariantContextHeader(final SAMSequenceDictionary sequenceDictionary,
                                            final Set<VCFHeaderLine> vcfDefaultToolHeaderLines) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Collections.singletonList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add default tool header lines */
        vcfDefaultToolHeaderLines.forEach(result::addMetaDataLine);

        result.setSequenceDictionary(sequenceDictionary);

        /* header lines related to genotype formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1,
                VCFHeaderLineType.Integer, "Segment genotype"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN, 1,
                VCFHeaderLineType.Integer, "Segment most-likely copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(NP, 1,
                VCFHeaderLineType.Integer, "Number of points (i.e. targets or bins) in the segment"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QS, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that at least one point " +
                "(i.e. target or bin) in the segment agrees with the segment copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QA, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that all points " +
                "(i.e. targets or bins) in the segment agree with the segment copy-number call"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QSS, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that the segment start " +
                "position is a genuine copy-number changepoint"));
        result.addMetaDataLine(new VCFFormatHeaderLine(QSE, 1,
                VCFHeaderLineType.Integer, "Complementary Phred-scaled probability that the segment end " +
                "position is a genuine copy-number changepoint"));

        /* INFO header lines */
        result.addMetaDataLine(new VCFInfoHeaderLine(VCFConstants.END_KEY, 1,
                VCFHeaderLineType.Integer, "End coordinate of the variant"));

        /*FILTER header lines */
        result.addMetaDataLine(GATKSVVCFHeaderLines.getFilterLine(GATKSVVCFConstants.LOW_QS_SCORE_FILTER_KEY));
        result.addMetaDataLine(GATKSVVCFHeaderLines.getFilterLine(GATKSVVCFConstants.FREQUENCY_FILTER_KEY));

        outputWriter.writeHeader(result);
    }

    /**
     * Compose a variant context from a given {@link IntegerCopyNumberSegment}
     *
     * @param segment an instance of {@link IntegerCopyNumberSegment}
     * @return composed variant context
     */
    @VisibleForTesting
    VariantContext composeVariantContext(final IntegerCopyNumberSegment segment) {
        final String contig = segment.getContig();
        final int start = segment.getStart();
        final int end = segment.getEnd();
        final int copyNumberCall = segment.getCallIntegerCopyNumberState().getCopyNumber();
        final Allele refAllele = reference == null ? REF_ALLELE : Allele.create(ReferenceUtils.getRefBaseAtPosition(reference, contig, start), true);

        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.chr(contig);
        variantContextBuilder.start(start);
        variantContextBuilder.stop(end);
        variantContextBuilder.id(String.format(VARIANT_PREFIX + "_%s_%d_%d", contig, start, end));

        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final IntegerCopyNumberState refCopyNumber = allosomalContigSet.contains(contig)
                ? segment.getBaselineIntegerCopyNumberState()
                : refAutosomalCopyNumberState;
        genotypeBuilder.alleles(makeGenotypeAlleles(copyNumberCall, refCopyNumber.getCopyNumber(), refAllele));
        genotypeBuilder.attribute(CN, copyNumberCall);
        genotypeBuilder.attribute(NP, segment.getNumPoints());
        genotypeBuilder.attribute(QS, FastMath.round(segment.getQualitySomeCalled()));
        genotypeBuilder.attribute(QA, FastMath.round(segment.getQualityAllCalled()));
        genotypeBuilder.attribute(QSS, FastMath.round(segment.getQualityStart()));
        genotypeBuilder.attribute(QSE, FastMath.round(segment.getQualityEnd()));
        final Genotype genotype = genotypeBuilder.make();

        final Set<Allele> uniquifiedAlleles = new HashSet<>();
        uniquifiedAlleles.add(refAllele);
        if (copyNumberCall > refCopyNumber.getCopyNumber()) {
            uniquifiedAlleles.add(DUP_ALLELE);  //dupes need additional alts since their genotypes are no-call
        } else if (copyNumberCall < refCopyNumber.getCopyNumber()) {
            uniquifiedAlleles.add(DEL_ALLELE);  //dels may be no-call
        }
        variantContextBuilder.alleles(uniquifiedAlleles);
        variantContextBuilder.attribute(VCFConstants.END_KEY, end);

        //copy over allele frequency etc.
        if (clusteredVCFReader != null) {
            try {
                final VariantContext cohortVC = clusteredVCFReader.query(
                        segment.getContig(), segment.getStart(), segment.getEnd()).stream().filter(vc -> vc.getStart() == segment.getStart() && vc.getEnd() == segment.getEnd()).collect(Collectors.toList()).get(0);
                //match segment end
                if (!(cohortVC == null)) {
                    if (cohortVC.hasAttribute(VCFConstants.ALLELE_COUNT_KEY)) {
                        final int cohortAC = cohortVC.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, -1);
                        if (cohortAC > -1) {
                            variantContextBuilder.attribute(GATKVCFConstants.ORIGINAL_AC_KEY, cohortAC);
                        }
                    }
                    if (cohortVC.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY)) {
                        final double cohortAF = cohortVC.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, -1);
                        if (cohortAF > -1) {
                            variantContextBuilder.attribute(GATKVCFConstants.ORIGINAL_AF_KEY, cohortAF);
                        }
                    }
                    if (cohortVC.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY)) {
                        final int cohortAN = cohortVC.getAttributeAsInt(VCFConstants.ALLELE_NUMBER_KEY, -1);
                        if (cohortAN > -1) {
                            variantContextBuilder.attribute(GATKVCFConstants.ORIGINAL_AN_KEY, cohortAN);
                        }
                    }
                } else {
                    logger.warn("No matching cohort VC at " + segment.getContig() + ":" + segment.getStart());
                }
            } catch (final IOException e) {
                throw new GATKException("Error querying file " + clusteredCohortVcf + " over interval " +
                        new SimpleInterval(segment.getContig(), segment.getStart(), segment.getEnd()), e);
            }
        }
        variantContextBuilder.genotypes(genotype);
        variantContextBuilder.log10PError(segment.getQualitySomeCalled()/-10.0);
        variantContextBuilder.filters(getInfoFilters(segment));
        return variantContextBuilder.make();
    }

    private List<Allele> makeGenotypeAlleles(final int copyNumberCall, final int refCopyNumber, final Allele refAllele) {
        final List<Allele> returnAlleles = new ArrayList<>();
        final Allele genotypeAllele = getAlleleForCopyNumber(copyNumberCall, refCopyNumber, refAllele);
        //some allosomes like Y can have ref copy number zero, in which case we just no-call
        if (refCopyNumber == 0) {
            return GATKVariantContextUtils.noCallAlleles(1);
        }
        //for only one haplotype we know which allele it has
        if (refCopyNumber == 1) {
           return Arrays.asList(getAlleleForCopyNumber(copyNumberCall, refCopyNumber, refAllele));
        //can't determine counts per haplotypes if there is a duplication
        } else if (genotypeAllele.equals(DUP_ALLELE)) {
            return GATKVariantContextUtils.noCallAlleles(refCopyNumber);
        //for homDels, hetDels or homRefs
        } else if (refCopyNumber == 2) {
            returnAlleles.add(genotypeAllele);
            if (copyNumberCall == 0) {
                returnAlleles.add(genotypeAllele);
            } else {
                returnAlleles.add(refAllele);
            }
            return returnAlleles;
        //multiploid dels
        } else {
            for (int i = 0; i < copyNumberCall; i++) {
                returnAlleles.add(refAllele);
            }
            for (int i = copyNumberCall; i < refCopyNumber; i++) {
                returnAlleles.add(DEL_ALLELE);
            }
            return returnAlleles;
        }
    }

    /**
     *
     * @param copyNumberCall
     * @param refCopyNumber
     * @param refAllele
     * @return variant allele if copyNumberCall != refCopyNumber, else refAllele
     */
    private Allele getAlleleForCopyNumber(final int copyNumberCall, final int refCopyNumber, final Allele refAllele) {
        final Allele allele;
        if (copyNumberCall > refCopyNumber) {
            allele = DUP_ALLELE;
        } else if (copyNumberCall < refCopyNumber) {
            allele = DEL_ALLELE;
        } else {
            allele = refAllele;
        }
        return allele;
    }

    private Set<String> getInfoFilters(final IntegerCopyNumberSegment segment) {
        final Set<String> returnFilters = new LinkedHashSet<>();
        final int qsThreshold;
        if (segment.getCallIntegerCopyNumberState().getCopyNumber() == 0) {
            qsThreshold = homDelQSThreshold;
        } else if (segment.getCallIntegerCopyNumberState().getCopyNumber() > segment.getBaselineIntegerCopyNumberState().getCopyNumber()) {
            qsThreshold = dupeQSThreshold;
        } else if (segment.getCallIntegerCopyNumberState().getCopyNumber() < segment.getBaselineIntegerCopyNumberState().getCopyNumber()) {
            qsThreshold = hetDelQSThreshold;
        } else {
            qsThreshold = 0;
        }
        if (segment.getQualitySomeCalled() < qsThreshold) {
            returnFilters.add(GATKSVVCFConstants.LOW_QS_SCORE_FILTER_KEY);
        }

        return returnFilters;
    }
}
