import argparse
import itertools
import os
import subprocess
import sys
import logging
from collections import defaultdict
from glob import glob
from shutil import copyfile
from multiprocessing import Pool

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

def reportStart(reportFile):

    reportDir = os.path.dirname(reportFile.name)+"/jquery.toc.min.js"
    fromPath=str(os.path.dirname(os.path.realpath(__file__))) + "/jquery.toc.min.js"
    copyfile(fromPath, reportDir)
    print("Copy TOC", fromPath, reportDir)

    startContent = """
    <html>
<head>
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"></script>
        <script src="jquery.toc.min.js"></script>
    <style>
        
        body {  
                background-color: whitesmoke;
                font-family: Arial, Helvetica, sans-serif;
            }
        
        #toc {
            top: 0px;
            left: 0px;
            height: 100%;
            overflow-y: scroll;
            position:fixed;
    background: #333;
    box-shadow: inset -5px 0 5px 0px #000;
    width: 300px;

    color: #fff;
}

#toc ul {
    margin: 0;
    padding: 0;
    list-style: none;
}

#toc li {
    padding: 5px 10px;
}

#toc a {
    color: #fff;
    text-decoration: none;
    display: block;
}

#toc .toc-h2 {
    padding-left: 10px;
}

#toc .toc-h3 {
    padding-left: 20px;
}

#toc .toc-active {
    background: #336699;
    box-shadow: inset -5px 0px 10px -5px #000;
}

        h1, h2, h3, h4, h5, h6 {font-stretch: bold}
        p {color: black; font-style: italic}
        tr {max-width: 90%;}
        img {max-width: 100%;}
    </style>


</head>
<body>
    
    <div style="">
        <div style="">
            <ul id="toc"></ul>
        </div>
        <div id="content" style="position:absolute; left: 400px">
    """

    reportFile.write(startContent)
    reportFile.flush()

def reportEnd(reportFile):

    endContent = """
        </div>
    </div>
        <script type="text/javascript">
            $("#toc").toc({content: "#content", headings: "h1,h2,h3,h4"});
        </script>
    
    </body>

    """

    reportFile.write(endContent)
    reportFile.flush()

def prepareDescriptions():
    plotId2Descr = {}
    """
    COUNTS

    """
    plotId2Descr[
        "FeatureCount Summary"] = "<p>For each sample/replicate, the descriptive statistics from featureCounts is shown</p>" \
                                  "<p>The number of total alignments is shown in the Total column. The Assigned column shows how many useful alignments have been contained</p>" \
                                  "<p>Ideally all Unassigned_ columns have zero counts.</p>" \
                                  "<p>First plots with absolute values are shown, followed by the same plots with relative values</p>"
    plotId2Descr[
        "Compare Replicates (msEmpiRe-normalized counts)"] = "<p>The raw counts from all samples (grouped by condition/replicates) are scattered against each other.</p>" \
                                                             "<p>In theory a diagonal is expected. The bulkier or wider the shape gets, the less consistent are the replicates.</p>" \
                                                             "<p>If particularly a few/one replicate has a wide shape compared to all others, it probably should get excluded.</p>" \
                                                             "<p>These counts are normalized by msEmpiRe (inter/intra-conditions) and thus can be directly compared.</p>"

    plotId2Descr[
        "Compare Replicates (raw counts)"] = "<p>The raw counts from all samples (grouped by condition/replicates) are scattered against each other.</p>" \
                                             "<p>In theory a diagonal is expected. The bulkier or wider the shape gets, the less consistent are the replicates.</p>" \
                                             "<p>If particularly a few/one replicate has a wide shape compared to all others, it probably should get excluded.</p>"

    plotId2Descr[
        "Compare Log Fold-Changes (raw counts)"] = "<p>These plots compare all pairwise logFCs from all replicates within each condition.</p>" \
                                                   "<p>Since these are replicates, no large logFcs are expected, nor very large ones.</p>"

    plotId2Descr[
        "Compare Inter-Log Fold-Changes (raw counts)"] = "<p>These plots compare all pairwise logFCs between the two conditions.</p>" \
                                                         "<p>The logFC distributions are expected to match closely together</p>"

    plotId2Descr[
        "Compare Inter-Log Fold-Changes (msEmpiRe-normalized counts)"] = "<p>These plots compare all pairwise logFCs between the two conditions.</p>" \
                                                                         "<p>The logFC distributions are expected to match closely together</p>"

    plotId2Descr[
        "Compare Counts Per Gene (raw counts)"] = "<p>These plots compare the count distributions of all genes with counts, count g.t. 1 and count g.t. 20</p>" \
                                                  "<p>All counts have an added pseudo-count of 1. The different samples/replicates should look similar, especially if from the same condition.</p>" \
                                                  "<p>The number of genes with a specific count threshold will drop. A rule of thumb might be, that genes with a count g.t. 20 are candidates for a differential analysis. If this count is too low, the DE analysis may not have (good) results.</p>"

    plotId2Descr["Compare Raw Counts Distribution"] = "<p>Comparison of gene counts (unnormalized).</p>" \
                                                      "<p>Here counts are taken directly from feature counts.</p>" \
                                                      "<p>The count distributions should match within conditions, and should not be too different between conditions.</p>"

    plotId2Descr["Raw Counts Per Gene"] = "<p>Comparison of gene counts (unnormalized) for the top 50 genes.</p>" \
                                          "<p>Here counts are taken directly from feature counts.</p>" \
                                          "<p>Within a condition/replicates the order and the count frequency of the listed genes should be similar.</p>" \
                                          "<p>If a single gene has a majority of the read counts (g.t. 10% for instance), this is suspicious.</p>"

    """
    DE ANALYSIS

    """


    """
    ENRICHMENT
    
    """

    plotId2Descr[
        "DAVID"] = "<p>DAVID gene set enrichment performed via <a href=\"https://david.ncifcrf.gov/\">DAVID Webclient</a></p>" \
                   "<p>The DAVID gene set enrichment performs an enrichment on all GO subsets, KEGG, REACTOME and Uniprot Keywords</p>" \
                   "<p>Significant terms have been selected by pval.adj l.t. 0.05 and logFC g.t. 1.0</p>"

    plotId2Descr["KEGG"] = "<p><a href=\"https://www.genome.jp/kegg/\">KEGG</a> gene set enrichment</p>" \
                           "<p></p>" \
                           "<p>Significant terms have been selected by pval.adj l.t. 0.05 and logFC g.t. 1.0</p>"

    plotId2Descr["REACTOME"] = "<p><a href=\"https://reactome.org/\">REACTOME</a> Pathway gene set enrichment.</p>" \
                               "<p>reactome is new in the business but aims at proving free, open-source, curated and peer-reviewed pathways.</p>" \
                               "<p>Significant terms have been selected by pval.adj l.t. 0.05 and logFC g.t. 1.0</p>"

    plotId2Descr[
        "GeneOntology (overrepresentation)"] = "<p><a href=\"http://geneontology.org/\">GeneOntology</a> based gene set overrepresentation analysis.</p>" \
                                               "<p>For this analysis, it is analysed whether significant genes are more common in a specific set than in the background.</p>" \
                                               "<p>Significant terms have been selected by pval.adj l.t. 0.05 and logFC g.t. 1.0</p>"

    plotId2Descr[
        "GeneOntology (GSEA)"] = "<p><a href=\"http://geneontology.org/\">GeneOntology</a> based gene set enrichment analysis.</p>" \
                                 "<p>The input is the ranked list of genes (here abs(logFC)).</p>" \
                                 "<p>Up/Down-regulated genes can be identified by a positive/negative NES value.</p>" \
                                 "<p>Significant terms have been selected by pval.adj l.t. 1.0 and logFC n.e. 0.0</p>"

    return plotId2Descr

def runExtProcess(x):
    print("Starting Process", x)
    subprocess.run(x, shell=True, check=True)

def readable_dir(prospective_dir):
    try:
        os.makedirs(prospective_dir)
    except:
        pass

    if not os.path.isdir(prospective_dir):
        raise Exception("readable_dir:\"{0}\" is not a valid path".format(prospective_dir))
    if os.access(prospective_dir, os.R_OK):
        return prospective_dir
    else:
        raise Exception("readable_dir:{0} is not a readable dir".format(prospective_dir))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-c', '--counts', nargs='+', type=argparse.FileType('r'), required=True, help='count files')
    parser.add_argument('-d', '--diffreg', type=readable_dir, nargs='+', required=True)

    parser.add_argument('-c1', '--cond1', nargs='+', type=str, action='append', required=True, help="columns in count files for cond1")
    parser.add_argument('-c2', '--cond2', nargs='+', type=str, action='append', required=True, help="columns in count files for cond2")

    parser.add_argument('-n', '--name', type=str, required=True)
    parser.add_argument('-o', '--organism', type=str, required=True)
    parser.add_argument('-on', '--organism-name', type=str, required=False)
    parser.add_argument('-om', '--organism-mapping', type=str, required=False)
    parser.add_argument('-s', '--save', type=readable_dir, required=True)

    parser.add_argument('-fc', '--fold_changes', action='store_true', default=False, help="run FC part")
    parser.add_argument('-ca', '--counts_analysis', action='store_true', default=False, help="run FC part")
    parser.add_argument('-enrich', '--enrichment', action='store_true', default=False, help="run FC part")
    parser.add_argument('-stats', '--stats', action='store_true', default=False, help="run FC part")

    parser.add_argument('-all', '--all', action='store_true', default=False, help="run FC part")
    parser.add_argument('-sim', '--simulate', action='store_true', default=False, help="run FC part")

    parser.add_argument('-e', '--enhance', nargs='+', type=argparse.FileType('r'), help='enhancement file', default=open("/mnt/d/dev/data/genomes/ensembl.symbol.mouse.list", "r"))
    parser.add_argument('-l', '--lengths', nargs='+', type=argparse.FileType('r'), help='lengths file', default=open("/mnt/d/dev/data/genomes/ensembl.symbol.mouse.lengths.list", "r"))

    parser.add_argument('-m', '--de_methods', nargs='+', type=str, help="differential methods for analysis. If combination, write DESeq;msEmpiRe")
    parser.add_argument('-p', '--prefixes', nargs='+', type=str, required=True, help="short names for all subruns")

    parser.add_argument('-r', '--report', type=argparse.FileType('w'), required=True, help='report files')

    parser.add_argument('--parallel', type=int, required=False, default=6)

    args = parser.parse_args()

    if args.all:
        args.fold_changes = True
        args.counts_analysis = True
        args.enrichment = True
        args.stats = True


    #args.simulate=False

    args.save = os.path.realpath(args.save)


    if args.organism_name == None:

        if args.organism == "mmu":
            args.organism_name = "mouse"
        elif args.organism == "hsa":
            args.organism_name = "human"

    if args.organism_mapping == None:

        if args.organism == "mmu":
            args.organism_mapping = "org.Mm.eg.db"


    if args.organism_mapping == None:
        raise argparse.ArgumentError("Please provide a valid organism_mapping for R")

    if args.organism_name == None:
        raise argparse.ArgumentError("Please provide a valid organism_name for Reactome")

    if not len(args.diffreg) == len(args.counts) or len(args.diffreg) > 2:
        raise argparse.ArgumentError("diffreg folders must be same length than counts and should not be more than 2")

    if not len(args.prefixes) == len(args.counts) or len(args.diffreg) > 2:
        raise argparse.ArgumentError("prefixes must be same length than counts and should not be more than 2")

    if not len(args.diffreg) == len(args.cond1) or not len(args.diffreg) == len(args.cond2):
        raise argparse.ArgumentError("cond1/2 must be same length as diffreg.")

    if any([len(x) < 2 for x in args.cond1]):
        raise argparse.ArgumentError("cond1 must have at least 3 replicates")

    if any([len(x) < 2 for x in args.cond2]):
        raise argparse.ArgumentError("cond2 must have at least 3 replicates")

    numberOfDiffSamples = len(args.prefixes)

    #CONDITIONSNAME = "macrophages"

    #STARDIFFREG = "data/macrophages/star.diffreg/"
    #HISATDIFFREG = "data/macrophages/hisat2.diffreg/"

    #STARCOUNTS = "data/macrophages/star.po.counts"
    #HISATCOUNTS = "data/macrophages/hisat2.po.counts"

    #SAVEOUT = "./save/"

    #RUNFCANALYSIS = "TRUE"
    #RUNCOUNTSANALYSIS = "TRUE"
    #RUNENRICHMENT = "TRUE"
    #RUNDESTATS = "TRUE"

    #ENHANCELIST = "--enhanced "
    #LENGTHLIST = "--lengths "

    #declare - a
    #demethods

    #demethods[0] = 'msEmpiRe'
    #demethods[1] = 'msEmpiRe;DESeq'

    #COND1PATHSTAR = "./data/macrophages/control//19083-0010_S28_R1_all_fastq.star.bam ./data/macrophages/control//19083-0011_S29_R1_all_fastq.star.bam ./data/macrophages/control//19083-0015_S33_R1_all_fastq.star.bam"
    #COND2PATHSTAR = "./data/macrophages/knockout//19083-0009_S27_R1_all_fastq.star.bam ./data/macrophages/knockout//19083-0012_S30_R1_all_fastq.star.bam ./data/macrophages/knockout//19083-0013_S31_R1_all_fastq.star.bam ./data/macrophages/knockout//19083-0014_S32_R1_all_fastq.star.bam"
    #COND1PATHHISAT = "./data/macrophages/control//19083-0010_S28_R1_all_fastq.hisat2.bam ./data/macrophages/control//19083-0011_S29_R1_all_fastq.hisat2.bam ./data/macrophages/control//19083-0015_S33_R1_all_fastq.hisat2.bam"
    #COND2PATHHISAT = "./data/macrophages/knockout//19083-0009_S27_R1_all_fastq.hisat2.bam ./data/macrophages/knockout//19083-0012_S30_R1_all_fastq.hisat2.bam ./data/macrophages/knockout//19083-0013_S31_R1_all_fastq.hisat2.bam ./data/macrophages/knockout//19083-0014_S32_R1_all_fastq.hisat2.bam"

    #COND1PATHSTARFULL = "./data/macrophages/control//18139-0002_R1_all.star.bam $COND1PATHSTAR"
    #COND2PATHSTARFULL =$COND2PATHSTAR
    #COND1PATHHISATFULL = "./data/macrophages/control//18139-0002_R1_all.hisat2.bam $COND1PATHHISAT"
    #COND2PATHHISATFULL =$COND2PATHHISAT

    cwd = os.getcwd()

    def path2rsave( inWord ):

        outw = inWord.replace("//", "/").replace("/", ".").replace("-", ".")
        outw = "X" + outw

        return outw


    cond1RPaths = []
    for cols in args.cond1:

        colRPaths = []
        for col in cols:
            ncol = col
            if ncol.startswith("./"):
                ncol = ncol[2:]

            npath = os.path.join(cwd, ncol)
            rpath = path2rsave(npath)
            #print(ncol, npath, rpath)
            colRPaths.append(rpath)
        cond1RPaths.append(colRPaths)

    cond2RPaths = []
    for cols in args.cond2:

        colRPaths = []

        for col in cols:
            ncol = col
            if ncol.startswith("./"):
                ncol = ncol[2:]

            npath = os.path.join(cwd, ncol)
            rpath = path2rsave(npath)
            #print(ncol, npath, rpath)
            colRPaths.append(rpath)

        cond2RPaths.append(colRPaths)
    logging.root.setLevel(logging.DEBUG)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s  %(name)s  %(levelname)s: %(message)s')
    consoleHandler.setFormatter(formatter)

    sysLogger = logging.getLogger("SYSTEM")
    sysLogger.addHandler(consoleHandler)
    # sysLogger.setLevel(logging.INFO)

    scriptMain = os.path.dirname(os.path.realpath(__file__)) + "/"
    sysLogger.info("Script-Main: " + str(scriptMain))

    reportStart(args.report)

    fcLogger = logging.getLogger("FOLDCHANGE")
    fcLogger.addHandler(consoleHandler)

    #args.simulate = True

    plotId2Descr = prepareDescriptions()

    if args.fold_changes:
        fcLogger.info("Starting")

        for pidx, prefix in enumerate(args.prefixes):
            fcLogger.info("Running SubSample {} ({})".format(pidx, prefix))
            sysCall = "python3 {script} foldchange_fc --output {outdir} --counts {countfile} --prefixes {prefix} --conditions {cond1} --conditions {cond2} --enhance {enhancePath} --lengths {lengthsPath} --no-rrna --fpkm --tpm".format(
                script=os.path.realpath(os.path.join(scriptMain, "../..", "scripts/poreAnalysis.py")),
                outdir=args.diffreg[pidx],
                countfile=args.counts[pidx].name,
                prefix=args.prefixes[pidx],
                cond1=" ".join(args.cond1[pidx]),
                cond2=" ".join(args.cond2[pidx]),
                enhancePath=args.enhance.name,
                lengthsPath=args.lengths.name
                )

            fcLogger.debug(sysCall)

            if not args.simulate:
                subprocess.run( sysCall, shell=True, check=True )

            fcLogger.info("Finished SubSample {} ({})".format(pidx, prefix))

    else:
        fcLogger.info("Skipping FOLDCHANGE ANALYSIS")


    def runSysCall(sysCall, descr, log, plotid, plotsPrefix, args, prefix, plotDict):
        log.info(descr)
        log.debug(sysCall)

        if not args.simulate:
            subprocess.run(sysCall, shell=True, check=True)

        if plotid != None:
            plotDict[plotid][prefix] = glob(plotsPrefix + "*.png")
            log.debug("Found Images\n" + "\n".join(plotDict[plotid][prefix]))

    caLogger = logging.getLogger("COUNTSANALYSIS")
    caLogger.addHandler(consoleHandler)

    if args.counts_analysis:
        caLogger.info("Starting")

        args.report.write("<h1>poreSTAT: COUNTS ANALYSIS</h1>\n")
        args.report.flush()

        caPlots = defaultdict(lambda: defaultdict(list))

        for pidx, prefix in enumerate(args.prefixes):
            caLogger.info("Running SubSample {} ({})".format(pidx, prefix))


            sysCall = "python3 {script} --fc {countfile} --output {outdir} --enhance {enhancePath} --lengths {lengthsPath} --no-rrna --fpkm --tpm".format(
                script=os.path.realpath(os.path.join(scriptMain, "calculateExpressionValues.py")),
                outdir=os.path.join(args.diffreg[pidx], "counts.tpm.fpkm.tsv"),
                countfile=args.counts[pidx].name,
                enhancePath=args.enhance.name,
                lengthsPath=args.lengths.name
            )

            caLogger.info("Calculate Expression Values")
            caLogger.debug(sysCall)

            if not args.simulate:
                subprocess.run(sysCall, shell=True, check=True)

            sysCall = "python3 {script} --summary {countfile} --output {outdir}".format(
                script=os.path.realpath(os.path.join(scriptMain, "visFCSummary.py")),
                outdir=os.path.join(args.diffreg[pidx], "fcsummary"),
                countfile=args.counts[pidx].name + ".summary"
            )



            runSysCall(sysCall, "Visualise featureCount Summary", caLogger, "FeatureCount Summary",
                       os.path.join(args.diffreg[pidx], "fcsummary"), args, prefix, caPlots)

            sysCall = "python3 {script} --pathname --counts {counts} --conditions {conds1} --conditions {conds2}".format(
                script=os.path.realpath(os.path.join(scriptMain, "compareReplicates.py")),
                counts=os.path.join(args.diffreg[pidx], "count_out_data_msEmpiRe.norm"),
                conds1=" ".join(cond1RPaths[pidx]),
                conds2=" ".join(cond2RPaths[pidx])
            )

            runSysCall(sysCall, "compareReplicates (normalized)", caLogger, "Compare Replicates (msEmpiRe-normalized counts)",
                       os.path.join(args.diffreg[pidx], "count_out_data_msEmpiRe.norm.replicates."), args, prefix, caPlots)

            sysCall = "python3 {script} --pathname --counts {counts} --conditions {conds1} --conditions {conds2} --output {output}".format(
                script=os.path.realpath(os.path.join(scriptMain, "compareReplicates.py")),
                counts=args.counts[pidx].name,
                conds1=" ".join(args.cond1[pidx]),
                conds2=" ".join(args.cond2[pidx]),
                output=os.path.join(args.diffreg[pidx], "orig_counts.countreplicates")
            )

            runSysCall(sysCall, "Compare Replicates (countreplicates)", caLogger, "Compare Replicates (raw counts)",
                       os.path.join(args.diffreg[pidx], "orig_counts.countreplicates"), args, prefix, caPlots)

            sysCall = "python3 {script} --counts {counts} --conditions {conds1} --conditions {conds2} --output {output}".format(
                script=os.path.realpath(os.path.join(scriptMain, "compareLogFC.py")),
                counts=args.counts[pidx].name,
                conds1=" ".join(args.cond1[pidx]),
                conds2=" ".join(args.cond2[pidx]),
                output=os.path.join(args.diffreg[pidx], "orig_counts.compLogFC")
            )


            runSysCall(sysCall, "Compare LogFC", caLogger, "Compare Log Fold-Changes (raw counts)",
                       os.path.join(args.diffreg[pidx], "orig_counts.compLogFC"), args, prefix, caPlots)

            sysCall = "python3 {script} --pathname --counts {counts} --conditions {conds1} --conditions {conds2} --output {output}".format(
                script=os.path.realpath(os.path.join(scriptMain, "compareInterLogFCs.py")),
                counts=args.counts[pidx].name,
                conds1=" ".join(args.cond1[pidx]),
                conds2=" ".join(args.cond2[pidx]),
                output=os.path.join(args.diffreg[pidx], "orig_counts.interLogFC")
            )

            runSysCall(sysCall, "Compare Inter LogFC", caLogger, "Compare Inter-Log Fold-Changes (raw counts)",
                       os.path.join(args.diffreg[pidx], "orig_counts.interLogFC"), args, prefix, caPlots)

            sysCall = "python3 {script} --counts {counts} --conditions {conds1} --conditions {conds2} --output {output}".format(
                script=os.path.realpath(os.path.join(scriptMain, "compareInterLogFCs.py")),
                counts=os.path.join(args.diffreg[pidx], "count_out_data_msEmpiRe.norm"),
                conds1=" ".join(cond1RPaths[pidx]),
                conds2=" ".join(cond2RPaths[pidx]),
                output=os.path.join(args.diffreg[pidx], "count_out_data_msEmpiRe.norm")
            )

            plotId2Descr["Compare Inter-Log Fold-Changes (msEmpiRe-normalized counts)"] = "<p>These plots compare all pairwise logFCs between the two conditions.</p>" \
                                                                          "<p>The logFC distributions are expected to match closely together</p>"


            runSysCall(sysCall, "Compare Inter LogFC (normalized)", caLogger, "Compare Inter-Log Fold-Changes (msEmpiRe-normalized counts)",
                       os.path.join(args.diffreg[pidx], "count_out_data_msEmpiRe.norm.interlogfc."), args, prefix, caPlots)

            sysCall = "python3 {script} --counts {counts} --groups {conds1} --groups {conds2} --output {output}".format(
                script=os.path.realpath(os.path.join(scriptMain, "countPlots.py")),
                counts=args.counts[pidx].name,
                conds1=" ".join(args.cond1[pidx]),
                conds2=" ".join(args.cond2[pidx]),
                output=os.path.join(args.diffreg[pidx], "counts")
            )


            runSysCall(sysCall, "Make Count Plots", caLogger, "Compare Raw Counts Distribution",
                       os.path.join(args.diffreg[pidx], "counts.countplot"), args, prefix, caPlots)

            sysCall = "python3 {script} --counts {counts} --conditions {conds1} --conditions {conds2} --output {output}".format(
                script=os.path.realpath(os.path.join(scriptMain, "compareCountsPerGene.py")),
                counts=args.counts[pidx].name,
                conds1=" ".join(args.cond1[pidx]),
                conds2=" ".join(args.cond2[pidx]),
                output=os.path.join(args.diffreg[pidx], "countspergene")
            )

            runSysCall(sysCall, "Compare Counts Per Gene", caLogger, "Raw Counts Per Gene",
                       os.path.join(args.diffreg[pidx], "countspergene.cpergenes"), args, prefix, caPlots)


            caLogger.info("Finished SubSample {} ({})".format(pidx, prefix))


        for plotId in caPlots:

            prefixCount = len(caPlots[plotId])

            tableOut = "<h3>{}</h3>".format(plotId)

            if plotId in plotId2Descr:
                tableOut += "<div>{}</div>".format(plotId2Descr[plotId])

            tableOut += "<table>"

            tableOut += "<tr>"
            for prefix in caPlots[plotId]:
                tableOut += "<td>" + str(prefix) + "</td>"
            tableOut += "</tr>"

            for filetuple in itertools.zip_longest(*[caPlots[plotId][x] for x in caPlots[plotId]]):
                tableOut += "<tr>"

                for imgFile in filetuple:
                    if imgFile != None:
                        tableOut += "<td><img src=\"" + str(
                            os.path.relpath(os.path.realpath(imgFile), os.path.dirname(args.report.name))) + "\"/></td>"
                    else:
                        tableOut += "<td>N/A</td>"
                tableOut += "</tr>"

            tableOut += "</table>"

            args.report.write(tableOut + "\n")
            args.report.flush()


    else:
        caLogger.info("Skipping FOLDCHANGE ANALYSIS")

    #args.simulate = False

    statsLogger = logging.getLogger("DE_ENRICHMENT")
    statsLogger.addHandler(consoleHandler)

    def runSysCall(sysCall, descr, log, plotid, plotsPrefix, args, prefix, method, plotDict):
        log.info(descr)
        log.debug(sysCall)

        if not args.simulate:
            subprocess.run(sysCall, shell=True, check=True)

        if plotid != None and plotsPrefix != None:
            plotDict[method][plotid][prefix] = glob(plotsPrefix + "*.png")
            log.debug("Found Images\n" + "\n".join(plotDict[method][plotid][prefix]))

    if args.stats or args.enrichment:
        statsLogger.info("Starting STATS OR ENRICHMENT")

        performMethods = [tuple(x.split(";")) for x in args.de_methods]

        statsLogger.info("Performing analysis for methods: "+" ".join([str(x) for x in performMethods]))

        deEnrichPlots = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

        if args.stats:

            args.report.write("<h1>poreSTAT: REPORTS for DE ANALYSIS</h1>\n")
            args.report.flush()

            for methods in performMethods:
                prefix2countFile = {}
                statsLogger.info("Running Methods {}".format(methods))
                methodStr = "_".join(methods)

                for pidx, prefix in enumerate(args.prefixes):
                    statsLogger.info("Running SubSample {} ({})".format(pidx, prefix))

                    countDeFile = glob("{diffreg}/count__*.tsv".format(diffreg=args.diffreg[pidx]))[0]
                    robustDeFile = os.path.join(args.save, args.name + "." + prefix + "." + methodStr + ".tsv")

                    prefix2countFile[prefix] = countDeFile
                    sysCall = "python3 {script} --de {counts} --methods {methods} --output {output}".format(
                        script=os.path.realpath(os.path.join(scriptMain, "calcRobustFCs.py")),
                        counts=countDeFile,
                        methods=" ".join(methods),
                        output=robustDeFile
                    )

                    plotId2Descr["DE Methods Overview ({})".format(" ".join(methods))] = "<p>The following plots rely on the robust values (logFC, PVAL) for the methods {}</p>" \
                                                                                         "<p>The upset plot shows how many genes have been discovered using each method.</p>"\
                                                                                         "<p>The volcano plot shows logFC and -log10(pVal) for the robustly detected genes (that are genes that are DE with all above methods).</p>".format(" ".join(methods))
                    runSysCall(sysCall, "Calculate Robust FCs", statsLogger, "DE Methods Overview ({})".format(" ".join(methods)), robustDeFile + ".rob.", args, prefix, methods, deEnrichPlots)

                    #makePCA.py --fc report_save/aortas.star.msEmpiRe_DESeq.tsv --output report_save/aortas.star.msEmpiRe_DESeq.tsv.mpca
                    sysCall = "python3 {script} --fc {counts} --output {output} --top_de ROB --num 1000".format(
                        script=os.path.realpath(os.path.join(scriptMain, "makePCA.py")),
                        counts=robustDeFile,
                        output=robustDeFile + ".mpca"
                    )

                    plotId2Descr["Cluster Raw Counts ({})".format(" ".join(
                        methods))] = "<p>The following plots cluster the raw counts for the top differential genes (methods {})</p>" \
                                     "<p>The cluster map plot shows how close the expression values (raw counts) are related..</p>" \
                                     "<p>The scatter plot has performed a UMAP transformation and displays these results.</p>".format(
                        " ".join(methods))

                    runSysCall(sysCall, "Cluster Data", statsLogger, "Cluster Raw Counts ({})".format(" ".join(
                        methods)), robustDeFile + ".mpca.", args, prefix, methods, deEnrichPlots)

                    sysCall = "python3 {script} --counts {counts} --conditions {conds1} --conditions {conds2} --last --ignoreMissing".format(
                        script=os.path.realpath(os.path.join(scriptMain, "compareCountsPerGene.py")),
                        counts=robustDeFile,
                        conds1=" ".join(args.cond1[pidx]),
                        conds2=" ".join(args.cond2[pidx])
                    )
                    plotId2Descr["Compare Counts Per Gene (DE, raw counts)"] = plotId2Descr["Compare Counts Per Gene (raw counts)"] + "" \
                                                                              "<p>In addition to the previous plot, gene names have here been replaced by gene symbols.</p>"

                    runSysCall(sysCall, "Compare Counts Per Gene (DiffReg/Robust)", statsLogger,
                               "Compare Counts Per Gene (DE, raw counts)",
                               robustDeFile + ".cpergenes.", args, prefix, methods, deEnrichPlots)

                    #python3 compareDECounts.py --de $COUNTFILESTAR --counts $STARDIFFREG/count_out_data_msEmpiRe.norm --conditions $COND1RPATHSTAR --conditions $COND2RPATHSTAR --tools ${dems[@]} > $SAVEOUT/$CONDITIONSNAME.star.$JDEMS.directcompare.tsv

                    if len(methods) >= 2:
                        sysCall = "python3 {script} --de {de} --counts {counts} --conditions {conds1} --conditions {conds2} --tools {methods} --output {output}".format(
                            de=countDeFile,
                            counts=os.path.join(args.diffreg[pidx], "count_out_data_msEmpiRe.norm"),
                            script=os.path.realpath(os.path.join(scriptMain, "compareDECounts.py")),
                            methods=" ".join(methods),
                            output=os.path.join(args.save, args.name + "." + prefix + "." + methodStr + ".directcompare.tsv"),
                            conds1=" ".join(cond1RPaths[pidx]),
                            conds2=" ".join(cond2RPaths[pidx])
                        )

                        runSysCall(sysCall, "Compare DE Counts", statsLogger, None, None, args, prefix, methods, deEnrichPlots)



                    for countType in ["TPM", "FPKM"]:

                        spConds1 = [x + "." + countType for x in args.cond1[pidx]]
                        spConds2 = [x + "." + countType for x in args.cond2[pidx]]

                        sysCall = "python3 {script} --counts {counts} --groups {conds1} --groups {conds2} --output {output}".format(
                            script=os.path.realpath(os.path.join(scriptMain, "countPlots.py")),
                            counts=countDeFile,
                            conds1=" ".join(spConds1),
                            conds2=" ".join(spConds2),
                            output=os.path.join(args.diffreg[pidx], "counts."+countType)
                        )

                        plotId2Descr["Compare {} Per Gene (raw counts)".format(countType)] = "<p>Comparison of gene expression.</p>" \
                                                                                             "<p>Here counts are taken from the {} column of the differential analysis.</p>" \
                                                                                             "<p>The count distributions should match within conditions, and should not be too different between conditions.</p>".format(countType)

                        runSysCall(sysCall, "Make Count Plots ({})".format(countType), statsLogger, "Compare {} Per Gene (raw counts)".format(countType),
                                   os.path.join(args.diffreg[pidx], "counts."+countType+".countplot"), args, prefix, methods, deEnrichPlots)

                        sysCall = "python3 {script} --counts {counts} --conditions {conds1} --conditions {conds2} --output {output}".format(
                            script=os.path.realpath(os.path.join(scriptMain, "compareCountsPerGene.py")),
                            counts=countDeFile,
                            conds1=" ".join(spConds1),
                            conds2=" ".join(spConds2),
                            output=os.path.join(args.diffreg[pidx], "countspergene."+countType)
                        )
                        plotId2Descr["Compare {} Counts Per Gene".format(countType)] = "<p>This plot shows the relative amount of all counts per gene.</p>" \
                                                                                       "<p>For all replicates the top genes should be the same. Also the order should not be changed too much, especially for high frequency genes.</p>" \
                                                                                       "<p></p>"

                        runSysCall(sysCall, "Compare Counts Per Gene", statsLogger, "Compare {} Counts Per Gene".format(countType),
                                   os.path.join(args.diffreg[pidx], "countspergene."+countType+".cpergenes"), args, prefix, methods, deEnrichPlots)



                        sysCall = "python3 {script} --fc {counts} --output {output} --top_de ROB --num 1000 --{ct}".format(
                            script=os.path.realpath(os.path.join(scriptMain, "makePCA.py")),
                            counts=robustDeFile,
                            output=robustDeFile + "." + countType + ".mpca",
                            ct=countType.lower()
                        )

                        plotId2Descr["Cluster {} Counts ({})".format(countType, " ".join(
                            methods))] = "<p>The following plots cluster the raw counts for the top differential genes (methods {})</p>" \
                                         "<p>The cluster map plot shows how close the expression values (raw counts) are related..</p>" \
                                         "<p>The scatter plot has performed a UMAP transformation and displays these results.</p>".format(
                            " ".join(methods))



                        runSysCall(sysCall, "Cluster Data", statsLogger, "Cluster {} Counts ({})".format(countType, " ".join(
                            methods)), robustDeFile + "." + countType + ".mpca", args, prefix, methods, deEnrichPlots)



                if True:

                    allPrefixes = [x for x in prefix2countFile]

                    combinedRaw = os.path.join(args.save, args.name + "." + "combined_raw" + "." + methodStr + ".tsv")
                    combinedDE = os.path.join(args.save, args.name + "." + "combined" + "." + methodStr + ".tsv")

                    sysCall = "python3 {script} --de1 {counts1} --de2 {counts2} --prefix1 {prefix1} --prefix2 {prefix2} --output {output}".format(
                        script=os.path.realpath(os.path.join(scriptMain, "mergeDiffreg.py")),
                        counts1=prefix2countFile[allPrefixes[0]],
                        counts2=prefix2countFile[allPrefixes[1]],
                        prefix1=allPrefixes[0],
                        prefix2=allPrefixes[1],
                        methods=" ".join(methods),
                        output=combinedRaw
                    )

                    runSysCall(sysCall, "Merging DE Methods", statsLogger, None, None, args, prefix, methods, deEnrichPlots)

                    prefix2countFile[prefix] = countDeFile
                    sysCall = "python3 {script} --de {counts} --methods {methods} --output {output}".format(
                        script=os.path.realpath(os.path.join(scriptMain, "calcRobustFCs.py")),
                        counts=combinedRaw,
                        methods=" ".join(methods),
                        output=combinedDE
                    )

                    runSysCall(sysCall, "Calculate Robust FCs", statsLogger, "DE Methods Overview", combinedDE + ".upset", args, prefix, methods, deEnrichPlots)


                    for countType in ["", "TPM", "FPKM"]:

                        sysCall = "python3 {script} --fc {counts} --output {output} --top_de ROB --num 1000 {ct}".format(
                            script=os.path.realpath(os.path.join(scriptMain, "makePCA.py")),
                            counts=combinedDE,
                            methods=" ".join(methods),
                            output=combinedDE + "." + countType + ".mpca",
                            ct="--"+countType.lower() if len(countType) > 0 else ""
                        )

                        plotId2Descr["Cluster Combined {} Counts ({})".format(countType, " ".join(
                            methods))] = "<p>The following plots cluster the raw counts for the top differential genes (methods {})</p>" \
                                         "<p>The cluster map plot shows how close the expression values (raw counts) are related..</p>" \
                                         "<p>The scatter plot has performed a UMAP transformation and displays these results.</p>".format(
                            " ".join(methods))


                        runSysCall(sysCall, "Cluster Data Combined", statsLogger,
                                   "Cluster Combined {} Counts ({})".format(countType, " ".join(
                                       methods)), combinedDE + "." + countType + ".mpca", args, prefix, methods,
                                   deEnrichPlots)



                #[plotid][method][prefix]
                print(deEnrichPlots)

                for method in deEnrichPlots:

                    tableOut = "<h2>{}</h2>".format("&".join(method))
                    args.report.write(tableOut + "\n")

                    for plotId in deEnrichPlots[method]:

                        prefixCount = len(deEnrichPlots[method][plotId])

                        tableOut = "<h3>{}</h3>".format(plotId)

                        if plotId in plotId2Descr:
                            tableOut += "<div>{}</div>".format(plotId2Descr[plotId])

                        tableOut += "<table>"

                        tableOut += "<tr>"
                        for prefix in deEnrichPlots[method][plotId]:
                            tableOut += "<td>" + str(prefix) + "</td>"
                        tableOut += "</tr>"

                        print(plotId)
                        for filetuple in itertools.zip_longest(*[deEnrichPlots[method][plotId][x] for x in deEnrichPlots[method][plotId]]):
                            tableOut += "<tr>"

                            print(filetuple)
                            for imgFile in filetuple:

                                if imgFile != None:
                                    tableOut += "<td><img src=\"" + str(os.path.relpath(os.path.realpath(imgFile), os.path.dirname(args.report.name))) + "\"/></td>"
                                else:
                                    tableOut += "<td>N/A</td>"

                            tableOut += "</tr>"

                        tableOut += "</table>"

                        args.report.write(tableOut + "\n")
                        args.report.flush()


        if args.enrichment:

            args.report.write("<h1>poreSTAT: ENRICHMENT ANALYSIS</h1>\n")
            args.report.flush()

            deEnrichTables = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
            deEnrichFiles = defaultdict(lambda: dict())

            for methods in performMethods:

                statsLogger.info("Running Methods {}".format(methods))
                methodStr = "_".join(methods)

                for pidx, prefix in enumerate(args.prefixes + ["combined"]):
                    statsLogger.info("Running SubSample {} ({})".format(pidx, prefix))

                    deFile = os.path.join(args.save, args.name + "." + prefix + "." + methodStr + ".tsv")

                    if prefix == "combined" and not os.path.exists(deFile):
                        print("Combined file not existing\n" + deFile + "\n", file=sys.stderr)
                        continue

                    sysCall = "python3 {script} --de {de} --output {output}".format(
                        script=os.path.realpath(os.path.join(scriptMain, "getRobustFCs.py")),
                        de=deFile,
                        output=os.path.join(args.save, args.name + "." + prefix + "." + methodStr + ".robust.tsv"),
                    )

                    runSysCall(sysCall, "Make Robust ({}, {})".format(prefix, methodStr), statsLogger, None, None, args, prefix,methods, deEnrichPlots)

                    enrichCalls = []

                    deEnrichFiles[methods][prefix] = deFile

                    for direction in ["all", "up", "down"]:

                        sysCall = "Rscript {script} {de} {org} {dir}".format(
                            script=os.path.realpath(os.path.join(scriptMain, "runDAVIDAnalysis.R")),
                            de=deFile,
                            org=args.organism,
                            dir=direction
                        )
                        enrichCalls.append(sysCall)

                        sysCall = "Rscript {script} {de} {org} {dir}".format(
                            script=os.path.realpath(os.path.join(scriptMain, "runGOAnalysis.R")),
                            de=deFile,
                            org=args.organism_mapping,
                            dir=direction)
                        enrichCalls.append(sysCall)

                        if args.organism_name != None:
                            sysCall = "Rscript {script} {de} {org} {dir}".format(
                                script=os.path.realpath(os.path.join(scriptMain, "runReactomeAnalysis.R")),
                                de=deFile,
                                org=args.organism_name,
                                dir=direction)
                            enrichCalls.append(sysCall)

                        sysCall = "Rscript {script} {de} {org} {dir}".format(
                            script=os.path.realpath(os.path.join(scriptMain, "runKeggAnalysis.R")),
                            de=deFile,
                            org=args.organism,
                            dir=direction)
                        enrichCalls.append(sysCall)

                        if direction == "all":
                            sysCall = "Rscript {script} {de} {org} {dir}".format(
                                script=os.path.realpath(os.path.join(scriptMain, "runGOGSEA.R")),
                                de=deFile,
                                org=args.organism_mapping,
                                dir=direction)
                            enrichCalls.append(sysCall)

                    for x in enrichCalls:
                        statsLogger.info(x)

                    parallel = args.parallel

                    #args.simulate = True

                    if not args.simulate:
                        with Pool(processes=parallel) as pool:
                            statsLogger.info("Started Parallel Pool with {} processes.".format(parallel))
                            pool.map(runExtProcess, enrichCalls, 1)
                        statsLogger.info("Stopped Parallel Pool with {} processes.".format(parallel))

                    #args.simulate = False

                    #[method][plotid][prefix]

            statsLogger.info("Fetching Enrichment Data")
            for methods in performMethods:
                statsLogger.info("Running Methods {}".format(methods))
                methodStr = "_".join(methods)

                for pidx, prefix in enumerate(args.prefixes):

                    deFile = deEnrichFiles[methods][prefix]

                    deEnrichTables[methods]["DAVID"][prefix] = glob(deFile + ".david*.tsv")
                    deEnrichTables[methods]["KEGG"][prefix] = glob(deFile + ".kegg*.tsv")
                    deEnrichTables[methods]["REACTOME"][prefix] = glob(deFile + ".reactome*.tsv")
                    deEnrichTables[methods]["GeneOntology (overrepresentation)"][prefix] = glob(deFile + ".GeneOntology*goenrich.tsv")
                    deEnrichTables[methods]["GeneOntology (GSEA)"][prefix] = glob(deFile + ".GeneOntology*gsea.tsv")




            for method in deEnrichTables:

                tableOut = "<h2>{}</h2>".format("&".join(method))
                args.report.write(tableOut + "\n")

                for plotId in deEnrichTables[method]:

                    prefixCount = len(deEnrichTables[method][plotId])

                    tableOut = "<h3>{}</h3>".format(plotId)

                    if plotId in plotId2Descr:
                        tableOut += "<div>{}</div>".format(plotId2Descr[plotId])

                    tableOut += "<table>"
                    tableOut += "<tr>"
                    for prefix in deEnrichTables[method][plotId]:
                        tableOut += "<td>" + str(prefix) + "</td>"
                    tableOut += "</tr>"

                    print(plotId)
                    for filetuple in itertools.zip_longest(
                            *[deEnrichTables[method][plotId][x] for x in deEnrichTables[method][plotId]]):
                        print(filetuple)
                        tableOut += "<tr>"

                        for imgFile in filetuple:
                            if imgFile != None:
                                tableOut += "<td><a href=\"" + str(os.path.relpath(os.path.realpath(imgFile),
                                                                                os.path.dirname(
                                                                                    args.report.name))) + "\">{}<a/></td>".format(os.path.basename(imgFile))
                            else:
                                tableOut += "<td>N/A</td>"
                        tableOut += "</tr>"

                    tableOut += "</table>"

                    args.report.write(tableOut + "\n")
                    args.report.flush()

    else:
        statsLogger.info("Skipping STATS and ENRICHMENT")


    reportEnd(args.report)
