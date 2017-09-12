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
        if numOfAlts > 1:
            numOfAlts = 1
        gtIdx = -1
        if 'GT' not in vcfLineSplit[formatColumnIdx]:
            vcfLineSplit[formatColumnIdx] = "GT:" + vcfLineSplit[formatColumnIdx]
            gtIdx = 0
        else:
            gtIdx = formatInfo['GT']
        for sampleColumnIdx in sampleColumnIdxs:
            genotype = getGenotypeUpdatedSampleInfo(vcfLineSplit[sampleColumnIdx], formatInfo, errorRate, numOfAlts)
            sampleSplit = vcfLineSplit[sampleColumnIdx].split(':')
            sampleSplit[gtIdx] = genotype
            vcfLineSplit[sampleColumnIdx] = ':'.join(sampleSplit)
        out.write('\t'.join(vcfLineSplit) + '\n')
    out.close()

def getGenotypeUpdatedSampleInfo(vcfSampleInfo, formatInfo, errorRate, numOfAlts):
    graphiteCounts = getCountsFromSampleInfo(vcfSampleInfo, formatInfo, numOfAlts)
    genotypeArray = ['.','.']
    thresholdsCountsAndIdxs = []
    totalCounts = sum(graphiteCounts)


    # minAltConfirmedCounts = binom.ppf(errorRate, totalCounts, 0.5)
    # minRefConfirmedCounts = binom.ppf(1-errorRate, totalCounts, 0.5)

    threshold = 3
    if graphiteCounts[0] > threshold:
        genotypeArray = [0,0]
    for idx in range(1, len(graphiteCounts)):
        if graphiteCounts[idx] > threshold:
            if genotypeArray[0] == '.':
                genotypeArray[0] = idx
                genotypeArray[1] = idx
            else:
                genotypeArray[0] = genotypeArray[1]
                genotypeArray[1] = idx

    genotype = '/'.join(str(i) for i in genotypeArray)
    return genotype

def getCountsFromSampleInfo(vcfSampleInfo, formatInfo, numOfAlts):
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
