import sys, os, argparse
from scipy.stats import binom


def main (argv):
    parser = argparse.ArgumentParser(description='Genotype Graphited Files', add_help=True)

    parser.add_argument('-e', '--error_rate', action="store", dest="error_rate", default=0.01, type=float, help="Error rate of the sequencing technology")
    parser.add_argument('-v', '--vcfs', action="append", dest="vcfs", default=[])
    parser.add_argument('-o', '--output_directory', action="store", dest="output_directory")

    results = parser.parse_args()
    errorRate = results.error_rate
    vcfs = results.vcfs
    outputDirectory = results.output_directory

    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)

    for vcf in vcfs:
        genotypeVCF(vcf, outputDirectory, errorRate)

def genotypeVCF(vcf, outputDirectory, errorRate):
    outFilePath = outputDirectory + '/' + vcf.split('/')[-1].replace('.vcf', '.filtered') + '.vcf'
    out = open(outFilePath, 'w')
    vcfFile = open(vcf, "r")
    sampleColumnIdxs = []
    formatInfo = {}
    formatColumnIdx = -1
    for vcfLine in vcfFile:
        vcfLine = vcfLine.replace('\n', '')
        vcfLineSplit = vcfLine.split('\t')
        if vcfLine.startswith('#'):
            if vcfLine.startswith('#CHROM'):
                headerValues = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
                for i in range(len(vcfLineSplit)):
                    if 'FORMAT' == vcfLineSplit[i]:
                        formatColumnIdx = i
                    elif vcfLineSplit[i] not in headerValues:
                        sampleColumnIdxs.append(i)
            out.write(vcfLine + '\n')
            continue
        if len(formatInfo) == 0:
            formatSplit = vcfLineSplit[formatColumnIdx].split(':')
            for i in range(len(formatSplit)):
                formatInfo[formatSplit[i]] = i

        numOfAlts = vcfLineSplit[4].count(',') + 1
        for sampleColumnIdx in sampleColumnIdxs:
            genotype = getGenotypeUpdatedSampleInfo(vcfLineSplit[sampleColumnIdx], formatInfo, errorRate, numOfAlts)
            if 'GT' not in vcfLineSplit[formatColumnIdx]:
                vcfLineSplit[formatColumnIdx] = "GT:" + vcfLineSplit[formatColumnIdx]
                vcfLineSplit[sampleColumnIdx] = genotype + ":" + vcfLineSplit[sampleColumnIdx]
            else:
                sampleSplit = vcfLineSplit[sampleColumnIdx].split(':')
                sampleSplit[formatInfo['GT']] = genotype
                vcfLineSplit[sampleColumnIdx] = ':'.join(sampleSplit)
            out.write('\t'.join(vcfLineSplit) + '\n')
    out.close()

def getGenotypeUpdatedSampleInfo(vcfSampleInfo, formatInfo, errorRate, numOfAlts):
    graphiteCounts = getCountsFromSampleInfo(vcfSampleInfo, formatInfo, numOfAlts)
    genotypeArray = [0,0]
    thresholdsCountsAndIdxs = []
    for altCount in range(1, numOfAlts + 1):
        refCounts = graphiteCounts[altCount - 1]['refCounts']
        altCounts = graphiteCounts[altCount - 1]['altCounts']
        totalCounts = refCounts + altCounts

        minAltConfirmedCounts = binom.ppf(errorRate, totalCounts, 0.5)
        minRefConfirmedCounts = binom.ppf(1-errorRate, totalCounts, 0.5)

        refThreshold = minRefConfirmedCounts
        altThreshold = max(10, minAltConfirmedCounts)
        if refCounts > refThreshold:
            if len(thresholdsCountsAndIdxs) == 0: thresholdsCountsAndIdxs.append({})
            thresholdsCountsAndIdxs[0] = {'counts': refCounts, 'index': 0}
        if altCounts > altThreshold:
            thresholdsCountsAndIdxs.append({'counts': altCounts, 'index': altCount})

    if len(thresholdsCountsAndIdxs) > 1:
        sortedList = sorted(thresholdsCountsAndIdxs, key=lambda v: v['counts'])
        genotypeArray[0] = sortedList[0]['index']
        genotypeArray[1] = sortedList[1]['index']
    elif len(thresholdsCountsAndIdxs) == 1:
        genotypeArray[0] = thresholdsCountsAndIdxs[0]['index']
        genotypeArray[1] = thresholdsCountsAndIdxs[0]['index']
    else:
        genotypeArray[0] = '.'
        genotypeArray[1] = '.'

    genotype = '/'.join(str(i) for i in genotypeArray)
    return genotype

def getCountsFromSampleInfo(vcfSampleInfo, formatInfo, numOfAlts):
    vcfSampleInfoSplit = vcfSampleInfo.split(':')
    counts = []
    for altCount in range(1, numOfAlts + 1):
        counts.append({'refCounts': 0, 'altCounts': 0})
        for field in ['DP4_NFP', 'DP4_NP', 'DP4_EP']:
            graphiteCounts = vcfSampleInfoSplit[formatInfo[field] * altCount].split(',')
            if graphiteCounts[0] == '.':
                counts[altCount - 1]['refCounts'] = -1
                counts[altCount - 1]['altCounts'] = -1
                break
            counts[altCount - 1]['refCounts'] += int(graphiteCounts[0]) + int(graphiteCounts[1])
            counts[altCount - 1]['altCounts'] += int(graphiteCounts[2]) + int(graphiteCounts[3])
    return counts

if __name__ == "__main__":
    main(sys.argv)
