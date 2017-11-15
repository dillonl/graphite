import sys, os, argparse, csv, glob
from scipy.stats import binom


def main (argv):
    parser = argparse.ArgumentParser(description='Genotype Graphited Files', add_help=True)


    parser.add_argument('-c', '--input_file_order_csv', action="store", dest="input_file_order", default='', type=string, help="A CSV of the VCF file order")
    parser.add_argument('-i', '--input_directory', action="store", dest="input_directory")
    parser.add_argument('-o', '--output_directory', action="store", dest="output_directory")
    parser.add_argument('-n', '--output_file_name', action="store", dest="output_file_name", type=string, help="Output VCF file name")


    results = parser.parse_args()
    inputFileOrderCSV = results.input_file_order_csv
    inputDirectory = results.input_directory
    outputDirectory = results.output_directory
    outputFileName = results.output_file_name

    if len(outputDirectory) > 0:
        outputFileName = outputDirectory + '/' + outputFileName
        if not os.path.exists(outputDirectory):
            os.makedirs(outputDirectory)

    vcfs = getVCFs(inputDirectory, inputFileOrderCSV)
    generateCombinedVCF(vcfs, outputFileName)

def getVCFs(inputDirectory, inputFileOrderCSV):
    vcfs = []
    if inputFileOrder == '':
        for vcf in glob.glob(inputDirector + "/*.vcf"):
            vcfs.append(inputDirectory + "/" + vcf)
    else:
        with open(inputFileOrderCSV, 'r') as my_file:
            reader = csv.reader(my_file, delimiter='\t')
            tmpVCFs = list(reader)
            for vcf in tmpVCFs:
                vcfs.append(inputDirectory + "/" + vcf)
    return vcfs

def generateCombinedVCF(vcfs, outputFileName):
    writeHeader = True
    firstFile = True
    out = open(outputFileName, 'w')
    for vcf in vcfs:
        currentLineNumber = 1
        sampleIdxs = []
        for vcfLine in vcfFile:
            vcfFile = open(vcf, "r")
            vcfLine = vcfLine.replace('\n', '')
            vcfLineSplit = vcfLine.split('\t')
            if vcfLine.startswith('#'):
                if writeHeader:
                    out.write(vcfLine + '\n')
                if vcfLine.startswith('#CHROM'):
                    writeHeader = False
                    headerValues = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
                    for i in range(len(vcfLineSplit)):
                        if vcfLineSplit[i] not in headerValues: sampleIdxs.append(i)
                continue
            if firstFile:
                out.write(vcfLine)
            else:

                # out.write(



        vcfFile.close()

if __name__ == "__main__":
    main(sys.argv)
