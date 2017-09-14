import sys, os, argparse
from scipy.stats import binom


def main (argv):
    parser = argparse.ArgumentParser(description='Genotype Graphited Files', add_help=True)

    parser.add_argument('-e', '--error_rate', action="store", dest="error_rate", default=0.01, type=float, help="Error rate of the sequencing technology")
    parser.add_argument('-v', '--vcfs', action="append", dest="vcfs", default=[])
    parser.add_argument('-o', '--output_directory', action="store", dest="output_directory")
    parser.add_argument('-u', '--upper_bound_percent', action="store", dest="upper_bound_percent", default=0.8, type=float, help="Upper percentage bound for het calls")
    parser.add_argument('-l', '--lower_bound_percent', action="store", dest="lower_bound_percent", default=0.2, type=float, help="Lower percentage bound for het calls")
    parser.add_argument('-t', '--threshold', action="store", dest="threshold", default=10, type=float, help="Threshold for total counts")

    results = parser.parse_args()
    errorRate = results.error_rate
    vcfs = results.vcfs
    outputDirectory = results.output_directory
    upperBoundPercent = results.upper_bound_percent
    lowerBoundPercent = results.lower_bound_percent
    threshold = results.threshold

    if not os.path.exists(outputDirectory):
        os.makedirs(outputDirectory)

    for vcf in vcfs:
        genotypeVCF(vcf, outputDirectory, errorRate, lowerBoundPercent, upperBoundPercent, threshold)

def genotypeVCF(vcf, outputDirectory, errorRate, lowerBoundPercent, upperBoundPercent):
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

        insertedGeno = False
        gtIdx = -1
        if 'GT' not in vcfLineSplit[formatColumnIdx]:
            vcfLineSplit[formatColumnIdx] = "GT:" + vcfLineSplit[formatColumnIdx]
            gtIdx = 0
            insertedGeno = True
        else:
            gtIdx = formatInfo['GT']
        for sampleColumnIdx in sampleColumnIdxs:
            genotype = getGenotypeUpdatedSampleInfo(vcfLineSplit[sampleColumnIdx], formatInfo, errorRate, lowerBoundPercent, upperBoundPercent, threshold)
            sampleSplit = vcfLineSplit[sampleColumnIdx].split(':')
            if insertedGeno:
                sampleSplit.insert(0, genotype)
            else:
                sampleSplit[gtIdx] = genotype
            vcfLineSplit[sampleColumnIdx] = ':'.join(sampleSplit)
        out.write('\t'.join(vcfLineSplit) + '\n')
    out.close()

def getGenotypeUpdatedSampleInfo(vcfSampleInfo, formatInfo, errorRate, lowerBoundPercent, upperBoundPercent, threshold):
    graphiteCounts = getCountsFromSampleInfo(vcfSampleInfo, formatInfo)
    genotypeArray = ['.','.']
    thresholdsCountsAndIdxs = []
    totalCounts = sum(graphiteCounts)
    genotypeArray = [0,0]


    # minAltConfirmedCounts = binom.ppf(errorRate, totalCounts, 0.5)
    # minRefConfirmedCounts = binom.ppf(1-errorRate, totalCounts, 0.5)

    if totalCounts >= threshold:
        for idx in range(0, len(graphiteCounts)):
            alleleFreq = float(graphiteCounts[idx]) / float(totalCounts)
            if (genotypeArray[0] == '.' and alleleFreq >= lowerBoundPercent) or alleleFreq > upperBoundPercent:
                genotypeArray[0] = idx
                genotypeArray[1] = idx
            elif lowerBoundPercent <= alleleFreq and alleleFreq <= upperBoundPercent:
                genotypeArray[0] = genotypeArray[1]
                genotypeArray[1] = idx

    genotype = '/'.join(str(i) for i in genotypeArray)
    return genotype

def getCountsFromSampleInfo(vcfSampleInfo, formatInfo):
    vcfSampleInfoSplit = vcfSampleInfo.split(':')
    counts = [0] * int(len(vcfSampleInfoSplit[formatInfo['DP4_NFP']].split(',')) / 2)

    for field in ['DP4_NFP', 'DP4_NP', 'DP4_EP']:
        graphiteCounts = vcfSampleInfoSplit[formatInfo[field]].split(',')
        if graphiteCounts[0] == '.':
            return [0]
        idxCount = 0
        for alleleIdx in range(0, len(graphiteCounts), 2):
            counts[idxCount] += int(graphiteCounts[alleleIdx]) + int(graphiteCounts[alleleIdx+1])
            idxCount += 1

    return counts

if __name__ == "__main__":
    main(sys.argv)
