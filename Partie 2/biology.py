from json import *


def est_base(c):
    return len(c) == 1 and c in "ATGC"


def est_adn(s):
    i = 0
    while i < len(s) and est_base(s[i]):
        i += 1
    return i >= len(s)


def arn(adn):
    if not est_adn(adn):
        return None
    s = ""
    i = 0
    while i < len(adn):
        if adn[i] == "T":
            s += "U"
        else:
            s += adn[i]
        i += 1
    return s


def arn_to_codons(arn):
    codons = []
    i = 0
    while i < len(arn) - 2:
        codons.append(arn[i] + arn[i + 1] + arn[i + 2])
        i += 3
    return codons


def load_dico_codons_aa(filename):
    fichier = open(filename, "r")
    strjson = fichier.read()
    fichier.close()
    return loads(strjson)


def codons_stop(dico):
    stop = []
    bases = "AUGC"
    i = 0
    while i < 4:
        j = 0
        while j < 4:
            k = 0
            while k < 4:
                if bases[i] + bases[j] + bases[k] not in dico:
                    stop.append(bases[i] + bases[j] + bases[k])
                k += 1
            j += 1
        i += 1
    return stop


def codons_to_aa(codons, dico):
    aa = []
    i = 0
    while i < len(codons) and codons[i] in dico:
        aa.append(dico[codons[i]])
        i += 1
    return aa




def nextIndice(tab, ind, elements):
    '''
        Prend en paramètre un tableau tab, un indice ind de tab, et un autre tableau d'elements. 
        La fonction recherche dans le tableau tab à partir de l'indice ind et retourne l'indice de la 
        première case du tableau tab contenant une valeur de elements
    '''
    while ind<len(tab):
        i=0
        while i<len(elements):
            if tab[ind]==elements[i]:
                return ind
            i+=1
        ind+=1
    return ind


def decoupe_sequence(seq,start,stop):
    '''
       Prend en paramètre trois tableaux seq, start et stop. 
       La fonction doit découper le tableau seq en séquences et retourner un tableau contenant 
       les différents morceaux.
    '''
    t=[]
    i=0
    while i<len(seq):
        p=[]
        if seq[i] in start:
            i+=1
            while i<len(seq) and seq[i] not in stop:
                
                p.append(seq[i])
                i+=1
            t.append(p)
        i+=1
    return t


def codons_to_seq_codantes(codons, dic):
    '''
    Prend en paramètre une séquence de codons et le dictionnaire de correspondance entre codons et 
    acides aminés.
    La fonction découpe la séquence de codons en séquences codantes. 
    Les différentes séquences sont stockées dans un tableau.
    '''
    stop=codons_stop(dic)
    start='AUG'
    fin=decoupe_sequence(codons,start,stop)
    return fin


def seq_codantes_to_seq_aas(seq, dic):
    '''
    Prend en paramètre un tableau de séquences codantes et le dictionnaire de correspondance entre codons 
    et acides aminés. 
    La fonction retourne un tableau contenant les différentes séquences d'acides aminés codées 
    par les différentes séquences codantes. 
    '''
    tabf=[]
    i=0
    while i<len(seq):
        j=0
        tab=codons_to_aa(seq[i],dic)
        i+=1
        tabf.append(tab)
    return tabf


def adn_encode_molecule(brin, dic, molecule):
    '''
    Prend en parmètre un brin d'ADN, le dictionnaire de correspondance entre codons et acides aminés et 
    une molécule.
    La fonction doit retourner True si l'ARN obtenu à partir de l'ADN puis découpé en codons contient 
    une séquence codante correspondant à la molécule, c'est-à-dire si la séquence d'acide aminée 
    correspondant à une séquence codante est la même que la molécule.
    '''
    ar=arn(brin)
    code3=arn_to_codons(ar)
    stop=codons_stop(dic)
    seq=codons_to_seq_codantes(code3,dic)
    seqfin=seq_codantes_to_seq_aas(seq,dic)
    if molecule in seqfin:
        return True
    return False

def test_fonction(adn,dic):
    '''
    Cette fonction va juste rendre plus visible le notebook et n'est pas demandé
    '''
    ar=arn(adn)
    code3=arn_to_codons(ar)
    stop=codons_stop(dic)
    seq=codons_to_seq_codantes(code3,dic)
    seqfin=seq_codantes_to_seq_aas(seq,dic)
    return seqfin