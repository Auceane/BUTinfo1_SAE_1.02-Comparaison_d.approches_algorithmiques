import biology as b

#Partie1
#Question 1
def test_est_base():
    assert b.est_base("A")
    assert b.est_base("T")
    assert b.est_base("G")
    assert not b.est_base("a")
    assert b.est_base("C")
    assert not b.est_base("Q")
    print("Test ok")

    
#Question 2
def test_est_adn():
    assert b.est_adn("AATGC")
    assert not b.est_adn("ATBOAAT")
    assert b.est_adn("ATGTCAAA")
    assert not b.est_adn("AgTA")
    assert b.est_adn("AGATACAATT")
    print ("Test Ok")

    
    
#Question 3    
def test_arn():
    assert b.arn("AGTTTGAT")=="AGUUUGAU"
    assert b.arn("TTTTTTT")=="UUUUUUU"
    assert b.arn("AGCAACAG")=="AGCAACAG"
    assert b.arn("TTUAG")!="UUUAG"
    assert b.arn("ATTGQ")==None
    print("Test Ok")
    
    
    
#Question 4
def test_arn_to_codons():
    assert b.arn_to_codons("")==[]
    assert b.arn_to_codons("AGT")==[]
    assert b.arn_to_codons("AUGT")==[]
    assert b.arn_to_codons("AUGCC")==["AUG"]
    assert b.arn_to_codons("AGUUUUGCAACC")==["AGU","UUU","GCA","ACC"]
    assert b.arn_to_codons("U")==[]
    print("Test Ok")

    
    
#Question 5    
def test_codon_stop():
    assert b.codons_stop(b.load_dico_codons_aa('data/codons_aa.json'))==['AGA', 'AGG', 'UAA', 'UAG']
    print("Test Ok")
    
    
    
#Question 6
def test_codons_to_aa():
    assert b.codons_to_aa([])==[]
    assert b.codons_to_aa(["GUA","UAG"])==["Valine"]
    assert b.codons_to_aa(["GUA","UAC","UAG","GAC"])==["Valine","Tyrosine"]
    assert b.codons_to_aa(["GUA","CAG","UCG","GAC"])!=["Valine","Serine","Serine","Aspartic acid"]
    assert b.codons_to_aa(["ACA","CAG","UCG","GAC"])==["Threonine","Glutamine","Serine","Aspartic acid"]
    assert b.codons_to_aa(["UAA","AGA","AUA"])==[]
    print("Test Ok")
    
    
    
#Parti 2
#Question 1
def test_nextIndice():
    assert b.nextIndice(["bonjour","hello","buongiorno","ciao","bye"],0,["hello","bye"])==1
    assert b.nextIndice(["bonjour","hello","buongiorno","ciao","bye"],1,["hello","bye"])!=4
    assert b.nextIndice(["bonjour","hello","buongiorno","ciao","bye"],2,["hello","bye"])==4
    assert b.nextIndice(["bonjour","hello","buongiorno","ciao","bye"],4,["hello","bye"])==4
    assert b.nextIndice(["bonjour","buongiorno","ciao","bye"],0,["hello","bye"])==3
    assert b.nextIndice(["bonjour","hello","buongiorno","ciao"],2,["hello","bye"])==4
    print("Test Valide")
    
    
#Question 2
def test_decoupe_sequence():
    
    #Test Basic
    seq = ["val1","début","val2","val3","end","val4","fin","begin","val5","fin","val6"]
    start = ["début","begin"]
    stop = ["fin","end"]
    assert b.decoupe_sequence(seq,start,stop)==[['val2','val3'],['val5']]
    
    #Test sans "val1" (change rien)
    seq = ["début","val2","val3","end","val4","fin","begin","val5","fin","val6"]
    assert b.decoupe_sequence(seq,start,stop)==[['val2','val3'],['val5']]
    
    #Test sans le premier "end"
    seq = ["val1","début","val2","val3","val4","fin","begin","val5","fin","val6"]
    assert b.decoupe_sequence(seq,start,stop)==[['val2','val3',"val4"],['val5']]
    
    #Test sans le premier "end" et le deuxième "fin"
    seq = ["val1","début","val2","val3","val4","fin","begin","val5","val6"]
    assert b.decoupe_sequence(seq,start,stop)==[['val2','val3','val4'],['val5','val6']]
    
    #Test sans "val1" , le premier "end", et le deuxième "fin"
    seq = ["début","val2","val3","val4","fin","begin","val5","val6"]
    assert b.decoupe_sequence(seq,start,stop)==[['val2','val3','val4'],['val5','val6']]
    print("Test Valide")
    
    
#Question 3
def test_codons_to_seq_codantes():
    #Test Basic
    seq=["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC",  "CGU", "UAG", "GGG"]
    dic=b.load_dico_codons_aa("data/codons_aa.json")
    assert b.codons_to_seq_codantes(seq,dic)==[['CGU', 'AUG', 'AAU'], ['GGG', 'CCC', 'CGU']]
    
    #Test sans le premier "AUG"
    seq=["CGU", "UUU", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC",  "CGU", "UAG", "GGG"]
    assert b.codons_to_seq_codantes(seq,dic)==[['AAU'],['GGG', 'CCC', 'CGU']]
    
    #Test sans le premier "AUG" et "UAA"
    seq=["CGU", "UUU", "CGU", "AUG", "AAU", "AUG", "GGG", "CCC",  "CGU", "UAG", "GGG"]
    assert b.codons_to_seq_codantes(seq,dic)==[['AAU', 'AUG', 'GGG', 'CCC', 'CGU']]
    
    #Test sans le premier "AUG", "UAA" et "UAG"
    seq=["CGU", "UUU", "CGU", "AUG", "AAU", "AUG", "GGG", "CCC",  "CGU", "GGG"]
    assert b.codons_to_seq_codantes(seq,dic)==[['AAU', 'AUG', 'GGG', 'CCC', 'CGU','GGG']]
    
    #Test sans "UAA"
    seq=["CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "AUG", "GGG", "CCC",  "CGU", "UAG", "GGG"]
    assert b.codons_to_seq_codantes(seq,dic)==[['CGU', 'AUG', 'AAU', 'AUG', 'GGG', 'CCC', 'CGU']]
    
    #Test avec un "AUG" au début
    seq=["AUG", "CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "UAA", "AUG", "GGG", "CCC",  "CGU", "UAG", "GGG"]
    assert b.codons_to_seq_codantes(seq,dic)==[["CGU", "UUU", "AUG", 'CGU', 'AUG', 'AAU'], ['GGG', 'CCC', 'CGU']]
    
    #Test avec un "AUG" au début et sans "UAA"
    seq=["AUG","CGU", "UUU", "AUG", "CGU", "AUG", "AAU", "AUG", "GGG", "CCC",  "CGU", "UAG", "GGG"]
    assert b.codons_to_seq_codantes(seq,dic)==[["CGU","UUU","AUG",'CGU','AUG','AAU',"AUG","GGG","CCC","CGU"]]
    
    print("Test Valide")
    
    
#Question 4
def test_seq_codantes_to_seq_aas():
    #Test Basic
    dic=b.load_dico_codons_aa("data/codons_aa.json")
    seq=[['CGU', 'AUG', 'AAU'], ['GGG', 'CCC', 'CGU']]
    assert b.seq_codantes_to_seq_aas(seq,dic)== [['Arginine', 'Methionine', 'Asparagine'], ['Glycine', 'Proline', 'Arginine']]
    
    #Test avec 'UAG' en plus (UAG n'a pas de valeur correspondante donc n'apparait pas)
    seq=[['CGU', 'AUG', 'AAU','UAG'], ['GGG', 'CCC', 'CGU']]
    assert b.seq_codantes_to_seq_aas(seq,dic)== [['Arginine', 'Methionine', 'Asparagine'], ['Glycine', 'Proline', 'Arginine']]
    
    #Test en ajoutant "CAU"
    seq=[['CGU', 'AUG', 'AAU',"CAU"], ['GGG', 'CCC', 'CGU']]
    assert b.seq_codantes_to_seq_aas(seq,dic)== [['Arginine', 'Methionine', 'Asparagine', 'Histidine'], ['Glycine', 'Proline', 'Arginine']]
    
    #Test en retirant "AUG"
    seq=[['CGU', 'AAU'], ['GGG', 'CCC', 'CGU']]
    assert b.seq_codantes_to_seq_aas(seq,dic)== [['Arginine', 'Asparagine'], ['Glycine', 'Proline', 'Arginine']]
    
    print("Test Valide")
    
    
    
#Question 5
def test_adn_encode_molecule():
    molecule=["Glycine", "Proline", "Arginine"]
    dic=b.load_dico_codons_aa("data/codons_aa.json")
    #Test du brin donner
    brin="CGTTTTATGCGTATGAATTAAATGGGGCCCCGTTAGGGG"
    assert b.adn_encode_molecule(brin,dic,molecule)
    
    #A chaque fois une/des lettre/s est/sont enlevée/s
    brin="CGTTTTATGCGTATGAATTAAATGGGCCCCGTTAGGGG"
    assert not b.adn_encode_molecule(brin,dic,molecule)
    
    brin="CGTTTTATGCGTATGATTAAATGGGGCCCCGTTAGGGG"
    assert not b.adn_encode_molecule(brin,dic,molecule)
    
    brin="CGTTTTAGCGTATGAATTAAATGGGGCCCCGTTAGGGG"
    assert not b.adn_encode_molecule(brin,dic,molecule)
    
    brin="CGTTTTATGCGTATGATTAAATGGGCCCCGTTAGGG"
    assert not b.adn_encode_molecule(brin,dic,molecule)
    
    brin="CGTTTTATGCGTATGAATTAAATGGGCCCGTTAGGGG"
    assert not b.adn_encode_molecule(brin,dic,molecule)
    
    brin="CGTTTTTGCGTATGAATTAAATGGGGCCCCGTTAGGG"
    assert not b.adn_encode_molecule(brin,dic,molecule)
    
    brin="CGTTTTATGCTGAATTAAATGGGGCCCCGTTAGGGG"
    assert b.adn_encode_molecule(brin,dic,molecule)
    
    brin="CGTTTTATGCGTATGAATTAAATGCCGTTAGGGG"
    assert not b.adn_encode_molecule(brin,dic,molecule)
    
    print("Test Valide")
