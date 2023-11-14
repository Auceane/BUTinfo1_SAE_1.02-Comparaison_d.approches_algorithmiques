import biology as b


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