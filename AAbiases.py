from Bio import SeqIO
from collections import defaultdict
from Scripts.filterPeptides import getFilteredList
import os
class AAbiasesEstimator():
    def __init__(self, fastaAAfile, groups = None):
        self.fastaAAfile = self.getWDoutfile(fastaAAfile)
        self.groupingFile = groups # two columns: seq.id \t group
        self.trnasCodeDic = {'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
    'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
    'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
    'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
    'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
    'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
    'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
    'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
    'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
    'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
    'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
    'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
    'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
    'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
    'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
    'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
   }
        self.groupPerSeq = {}
        self.whiteList = getFilteredList()
        self.outFile = self.fastaAAfile + "AAcnt"
        self.total_aa_dic = defaultdict(int) # total aa in set per group
        self.allGroups = self.__determineGroups() # all groups
        self.allAA = set([self.trnasCodeDic[i] for i in self.trnasCodeDic])
        self.aa_cnt_dic = {i:defaultdict(int) for i in self.allGroups} # group1: A: 5, G:5 ...... group2: A: 5, G:5 ......
        self.percentaGEaa_dic = {i:defaultdict() for i in self.allGroups}
        self.runCount()



    def runCount(self):
        self.__readFasta()
        self.__countPercentage()
        self.__writeOut()

    def getWDoutfile(self,filename):
        wd = r"D:\PycharmProjects\Secretome\Data"
        return wd + "\\" + filename

    def __determineGroups(self):
        allGroups = []
        if not os.path.isfile(self.groupingFile):
            print(self.groupingFile)
            self.groupPerSeq = {seq.id:self.groupingFile for seq in SeqIO.parse(self.fastaAAfile, "fasta")}
            allGroups.append(self.groupingFile)
        else:
            with open(self.groupingFile) as gr:
                for lines in gr:
                    sp = lines.rstrip().split("\t")
                    self.groupPerSeq[sp[0]] = sp[1]
                    allGroups.append(sp[1])

        print (set(allGroups))
        print(self.groupPerSeq)
        return (set(allGroups))


    def __readFasta(self):
        """
        :return: dictionary with number of each aa in total dataset and count also total aa in dataset
        """
        totalLen = 0
        for seq in SeqIO.parse(self.fastaAAfile, "fasta"):
            totalLen += len(seq.seq)
            if seq.id in self.groupPerSeq and (str(seq.seq) in self.whiteList or "random" in self.fastaAAfile):
                self.total_aa_dic[self.groupPerSeq[seq.id]] += len(seq.seq)
                for aa in self.allAA:
                    self.aa_cnt_dic[self.groupPerSeq[seq.id]][aa] += str(seq.seq).count(aa)

        print(self.aa_cnt_dic)
        print(totalLen)

    def __countPercentage(self):
        print(self.total_aa_dic)
        print(self.aa_cnt_dic)
        for groups in self.aa_cnt_dic:
            for aa in self.aa_cnt_dic[groups]:
                self.percentaGEaa_dic[groups][aa] = (self.aa_cnt_dic[groups][aa] * 100) / self.total_aa_dic[groups]

    def __writeOut(self):
        with open(self.outFile, "w") as out_tab:
            for groups in self.percentaGEaa_dic:
                for aa in self.percentaGEaa_dic[groups]:
                    out_tab.write(aa + "\t" + str(self.percentaGEaa_dic[groups][aa]) + "\t" + groups + "\n")


aamerged = open("AA_count_merged.txt", "w")
for fn in os.listdir(r'D:\PycharmProjects\Secretome\Data\fastaSets'):
    gropu = fn+ "_gr"
    aac = AAbiasesEstimator(fn, groups=gropu)
    with open(aac.outFile) as infile:
        for lines in infile:
            aamerged.write(lines)
aamerged.close()

#AAbiasesEstimator("translated_all_70095-sORF.fasta", groups="ClassifcationForAAcount.txt")
#("translated_all_70095-sORF.fasta", groups="peptide_sbs_classif.txt")




