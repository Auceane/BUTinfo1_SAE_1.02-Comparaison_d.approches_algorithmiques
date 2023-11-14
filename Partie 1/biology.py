from json import *

#Question 1 
def est_base(car):
    """
    Cette fonction nous dit si le caractère traiter correspond a une base nucléique de l'adn.
        car est le caractère testé.
    """
    return car=='A' or car=='T' or car=='G' or car=='C'


#Question 2
def est_adn(ch):
    """
    Cette fonction nous dit si les caractères de la chaîne de caractère sont bien des bases nucléiques.
    vcar est la valeur de véritée de la chaîne ch qui est vérifié a chaque boucle pour voir si les caractères correspondes.
    """
    i=0
    if len(ch)==0:
        return False
    vcar=est_base(ch[i])
    while i<len(ch)-1 and vcar==True:
        i+=1
        vcar=est_base(ch[i])
    return vcar


#Question 3
def arn(ch):
    """
    La chaîne de caractère pris en compte de l'adn est transformé en arn, les thymines sont transformé en uraciles.
    ch est la chaine de caractère, t est le tableau reccupant les valeurs afin de les tranformé en une chaine de caractère qui sera retournée a la fin qui est donc f.
    """
    i=0
    t=[]
    if est_adn(ch)==False:
        return None
    f=''
    while i<len(ch):
        if ch[i]=='T':
            t.append('U')
        else:
            t.append(ch[i])
        f+=t[i]
        i+=1
    return f


#Question 4
def arn_to_codons(ch):
    """
    Retourne un tableau à partir d'une chaîne d'arn.
    Les éléments du tableau sont des chaînes de caractère a 3 caractère.
    tt est la valeur qui prend le surplus de en faisant en sorte que le resultat soit forcément de 3 caractères.
    tab est le tableau des différents codons obtenus.
    cha prend 3 caractères et les ajoute dans le tableau.
    """
    tt=len(ch)%3
    tab=[]
    i=0
    if est_arn(ch)==False:
        return tab
    while i<len(ch)-tt:
        j=0
        cha=''
        while j<3:
            cha+=ch[i+j]
            j+=1
        tab.append(cha)
        i+=3
    return tab

def est_base_arn(car):
    """
    Vérifie que le caractère est un caractère d'un arn. 
    """
    return car=='A' or car=='U' or car=='G' or car=='C'
def est_arn(ch):
    """
    Vérifie que les caractères d'une chaîne sont des caractères d'un arn. 
    """
    i=0
    if len(ch)==0:
        return False
    vcar=est_base_arn(ch[i])
    while i<len(ch)-1 and vcar==True:
        i+=1
        vcar=est_base_arn(ch[i])
    return vcar


#Question 5
def tabl(c):
    """
    Cette fonction produit le tableau des possibilités de combinaison possible en 3 caractères.
    """
    t=[]
    i=0
    
    while i<4:
        ch=c[i]
        j=0
        while j<4:
            cha=ch+c[j]
            k=0
            while k<4:
                char=cha+c[k]
                t.append(char)
                k+=1
            j+=1
        i+=1
    return t

def load_dico_codons_aa(filename):
    """
    Donne le dictionnaire du fichier codons_aa.json.
    """
    codons_aa=open(filename,'r')
    codons_aa_json=codons_aa.read()
    dic=loads(codons_aa_json)
    codons_aa.close()
    return dic

def codons_stop(dic):
    """
    Donne le tableau des codons stop existant.
    fin est le tableau retourné.
    """
    c='AUGC'
    
    fin=[]
    i=0
    t=tabl(c)
    
    while i<len(t):
        el=t[i]
        tf=False
        j=0
        while j<len(dic) and tf==False:
            tf=el==list(dic)[j]
            j+=1
        if tf==False:
            fin.append(el)
        i+=1
    return fin



#Question 6
def codons_to_aa(table):
    """
    Donne le tableau des protéines que produise les différentes combinaisons de caractères.
    Retourne un tableau t contenant les valeurs donner par le dictionnaire dic et le tableau cd. 
    """
    i=0
    t=[]
    
    dic=load_dico_codons_aa('data/codons_aa.json')
    
    tfin=True
    cd=codons_stop(dic)
    while i<len(table) and tfin==True:
        j=0
        while j<len(dic) and i==len(t):
            if list(dic)[j]==table[i]:
                t.append(list(dic.values())[j])
            elif table[i]==cd[j%len(cd)]:
                tfin=False
            j+=1
        i+=1
    return t
