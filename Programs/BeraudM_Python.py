
A transformer en ipynb

### Lipofabrik
### Programme de traduction des coordonnées pour comparer 2 génomes annotés par 2 entités différentes. Ceci pour permettre de vérifier les mutations obtenues.

#190814-Translation of coordinates from one genome to another using BLAST

##open first genome
working_genome="C:/Users/PCPORT-LIPOFABRIK-04/Documents/Python/input/6633_curated_OneLine.fa"
genome=open(working_genome,"r")

#For ligne in genome:
 #   print (ligne)

ligne=genome.readlines()[1]


#open list of mut
coord="C:/Users/PCPORT-LIPOFABRIK-04/Documents/Python/input/input_mut.txt"
mut=open(coord,"r")
ligne_mut=mut.readlines()



#create the list
towrite=list()

##k=0
##while k<4:
##        bou=int(ligne_mut[k][0:-1])-1
##        bouend=int(ligne_mut[k][0:-1])-401
##        caption=str(bou)
##        towrite.append(caption)
##        towrite.append(ligne[bouend:bou])
##        k+=1
        
k=0
while k<4:
        bou=int(ligne_mut[k][0:-1])-1
        bouend=int(ligne_mut[k][0:-1])-401
        caption=str(bou)
        towrite.append(caption+"\n")
        towrite.append(ligne[bouend:bou]+"\n")
        k+=1

schtroumpf=list()
k=0
while k<len(towrite):
    schtroumpf.append(towrite[k])
    k+=1

#write in a file
final_fasta="C:/Users/PCPORT-LIPOFABRIK-04/Documents/Python/input/output.fa"
file_fasta=open(final_fasta,"a")


k=0
while k<len(schtroumpf):
    file_fasta.write(schtroumpf[k])
    k+=1


mut.close
genome.close
file_fasta.close

########################################################""
### Lipofabrik
#-*-coding: utf-8 -*-
"""
Created 190523 by MBE
To remove the line break in the fasta protein file
"""

# To change
chemin="C:/Users/PCPORT-LIPOFABRIK-04/Documents/Python/input/"


#Basics
from os import system
from os import chdir


fichier="contigs2019_prot.faa"

protseq=list()

chdir(chemin)
fichier1=open(fichier,"r")
ligne=fichier1.readlines()

"""
too check
len(ligne)
26270 => OK, jujuedit = 26271
print(ligne[-1])
DETKVEPDSVFRALGSHGVVRNGKAFQVIIGLSVPQMRERVEKILNQ => OK with jujuedit
bou=ligne[1].split("\n")
bou
['MDILRNQFPAQRSSKHIRERSRLLPNFSKKNAALKTSGFWAQTISGFIPKVMNSLNNLSF', '']
=> so yes, it is a problem.
bou=ligne[1].split(">")
bou
['MDILRNQFPAQRSSKHIRERSRLLPNFSKKNAALKTSGFWAQTISGFIPKVMNSLNNLSF\n']
bou
['', 'gene_00001 putative ABC transporter ATP-binding protein\n']
>>> ligne[1].split(">")
['MDILRNQFPAQRSSKHIRERSRLLPNFSKKNAALKTSGFWAQTISGFIPKVMNSLNNLSF\n']
>>> ligne[0].split(">")
['', 'gene_00001 putative ABC transporter ATP-binding protein\n']

ligne[1][-1]
'\n'
ligne[1][0:-2]
'MDILRNQFPAQRSSKHIRERSRLLPNFSKKNAALKTSGFWAQTISGFIPKVMNSLNNLS'
>>> ligne[1]
'MDILRNQFPAQRSSKHIRERSRLLPNFSKKNAALKTSGFWAQTISGFIPKVMNSLNNLSF\n'
>>> ligne[1][0:-1]
'MDILRNQFPAQRSSKHIRERSRLLPNFSKKNAALKTSGFWAQTISGFIPKVMNSLNNLSF'
>>> ligne[0].split(">")[0]==""
True
>>> ligne[1].split(">")[0]==""
False
"""



k=0
while k<len(ligne) :
	if ligne[k].split(">")[0]=="" :
		protseq.append("\n"+ligne[k])
	else:
		protseq.append(ligne[k].split("\n")[0])

	k+=1



fichier_retour="190523_fasta_retour.fa"
fichier2=open(fichier_retour,"w")
bou=str()

k=0
while k<len(protseq):
        bou=bou + protseq[k]
        k+=1

fichier2.write(bou[1:])

"""
fichier2.write(protseq[0:-1])
k=0
while k<len(protseq):
        fichier2.append(protseq[k])
        k+=1


fichier2=open(fichier_retour,"a")
bou=fichier_retour.read()
print(bou)

k=0
while k<len(protseq):
        fichier2.write(protseq[k])
        k+=1
"""

fichier2.close
fichier1.close
        






