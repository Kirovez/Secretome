from collections import defaultdict
import collections
from Bio import SeqIO
import numpy
import re
from Scripts.filterPeptides import getFilteredList
class DegaradomMS():
    def __init__(self, fastaMSpeptides, fastaProteins, portions=10):
        """
        :param fastaMSpeptides: peptides from MS in fasta format. !!!Protein name must be a pert of the name
        :param fastaProteins: full proteins
        :param portions: number of portion to which protein will be divided to estimate patterns
        """
        self.infileName = fastaMSpeptides
        self.fastaMSpeptides = fastaMSpeptides
        self.NEW_whileList = getFilteredList() ## the peptides that occuried in > 1 MS experiments for revision 02/07/2018
        self.fastaProteins = self.getWDoutfile(fastaProteins)
        self.portions = portions
        self.indProtein = SeqIO.index(self.fastaProteins, "fasta")
        self.PepperProt = defaultdict(list)
        self.buildDicPepperProt()

    def __str__(self):
        """
        Takes two files and build matrix with number of peptides un protein portion
        :return: file with matrix
        """

    def getWDoutfile(self,filename):
        wd = r"D:\PycharmProjects\Secretome\Data"
        return wd + "\\" + filename

    def getProteinId(self, pepID):
        if "Cp" in pepID or "Mp" in pepID:
            return pepID.split("_")[1]
        else:
            tmp = "P" + pepID.split("P")[1]
            return tmp.split(".p")[0] + ".p"

    def buildDicPepperProt(self):
        """
        read file with peptides and build dictionary as follow:  protein id : [pep1, pep2, ...]
        :return: dictionary peptides per protein
        """
        for pep in SeqIO.parse(self.getWDoutfile(self.fastaMSpeptides), "fasta"):
            if pep.seq in self.NEW_whileList:
                #extract protein name from pep name. CHANGE it if pep id differs from this
                #in our case pep id is like this: Secretome_Pp3c5_12820V3.1.p_1
                protein_name = self.getProteinId(pep.id)
                self.PepperProt[protein_name].append(str(pep.seq))


    def getPositionPart(self,protein, peptide):
        """
        :param protein: Seq of protein
        :param peptide: list of peptides
        :return: array number of peptides in each portion ordered by portion
        """
        pep_per_portion = {i:0 for i in range(0,self.portions)}
        row = numpy.linspace(0, len(protein.seq), self.portions+1)[1:]
        #print(row)
        for p in peptide:
            index = str(protein.seq).index(p)
            porion = numpy.searchsorted(row, index)
            pep_per_portion[porion] += 1

            #print(index, porion, pep_per_portion)


        return ([str(pep_per_portion[i]) for i in collections.OrderedDict(sorted(pep_per_portion.items()))])

    def writeMatrix(self):
        with open(self.getWDoutfile("Matrix_Degradome_" + self.fastaMSpeptides), "w") as outfile:
            for proteins in self.PepperProt:
                print(self.getPositionPart(self.indProtein[proteins], self.PepperProt[proteins]))
                outfile.write(proteins + "\t" + "\t".join(self.getPositionPart(self.indProtein[proteins], self.PepperProt[proteins])) + "\n")

    def visualizePepOnPro(self, ProteinId):
        protseq = str(self.indProtein[ProteinId].seq)
        print(protseq)

        for p in self.PepperProt[ProteinId]:
            index = protseq.index(p)
            indstr = str(index)
            print(" "*(index - len(indstr) - 1) + indstr + "-" + p)

    def run(self):
        """
        write degradation matrix
        :return:
        """
        self.writeMatrix()

    def findProteolyticSite(self, pattern):
        pat = re.compile(pattern)

        for seq in SeqIO.parse(self.fastaProteins, "fasta"):
            m = re.search(pat, str(seq.seq))
            if m:
                found = m.group(1)
                print(m)

    def countPeptides(self, protein_IDs):
        print("Number of proteins in the list", len(protein_IDs))
        cnt = 0
        total = 0
        for prot in self.PepperProt:
            total += len(self.PepperProt[prot])
            if prot in protein_IDs:
                cnt += len(self.PepperProt[prot])

        print("{} of total {} peptides covered of the protein in the list".format(cnt, total))


    def getFastaProteins(self):

        with open(self.getWDoutfile("Selected_proteins_{}.fasta".format(self.infileName.split(".")[0])), "w") as outfile:
            self.buildDicPepperProt()
            cnt = 0
            t_aa = 0
            for seq in SeqIO.parse(self.fastaProteins, "fasta"):
                if seq.id in self.PepperProt:
                    seq.description = ""
                    SeqIO.write(seq, outfile, "fasta")
                    cnt += 1
                    t_aa += len(seq.seq)

            print("total number of proteins in fasta selected:", cnt)
            print("total number of aa:", t_aa)


        """
    write file with proteins present in secretome file
    """



#DegaradomMS("Cellular_pep.txtMeJa.fasta", "Ppatens_318_v3.3.protein(1).fa").getFastaProteins()
#protein = "Pp3c22_19780V3.1.p"

#DegaradomMS("Cellular_pep.txtControl.fasta", "Ppatens_318_v3.3.protein(1).fa").visualizePepOnPro(protein)
#print("**********"*20)
#DegaradomMS("Secretome_pep.txtControl.fasta", "Ppatens_318_v3.3.protein(1).fa").visualizePepOnPro(protein)
#print("**********"*20)
# DegaradomMS("Secretome_pep.txtControl.fasta", "Ppatens_318_v3.3.protein(1).fa").run()
# DegaradomMS("Secretome_pep.txtMeJa.fasta", "Ppatens_318_v3.3.protein(1).fa").run()
#
# #DegaradomMS("Cellular_pep.txtControl.fasta", "Ppatens_318_v3.3.protein(1).fa").findProteolyticSite("")
#
# DegaradomMS("NEW_Cellular_pep.txtControl.fasta", "Ppatens_318_v3.3.protein(1).fa").run()
# DegaradomMS("NEW_Cellular_pep.txtMeJa.fasta", "Ppatens_318_v3.3.protein(1).fa").run()
