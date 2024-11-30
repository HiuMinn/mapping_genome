<h1 align="center">Projet d’Algorithmique pour la Génomique</h1>


$\to$ Développement d’une solution de mapping de données de séquençage à haut-débit sur un génome de référence  
Lilian Guitart Arnau, Hieu Minh Dang, Mathis Gallier Jouve et Simon Henninot

## I. L’algorithme de mapping, sa description et sa complexité

### I.1. Desciption générale de l'algorithme

Nous disposons des séquences complètes des 14 chromosomes et de l'apicoplaste de P. falciparum (parasite responsable du paludisme chez l'homme) sur lequel nous souhaitons mapper des reads obtenus par une technique de séquençage haut-débit.  
La méthode naïve consistant à parcourir le génome pour le comparer à chaque read n'est pas envisageable. En effet, la complexité temporelle d'un tel algorithme serait en  $O(N\_r \cdot T\_r \cdot T\_g)$ avec :  
$N\_r=$ le nombre de reads  
$T\_r=$ la taille d'un read  
$T\_g=$ la taille du génome complet  
Dans notre cas, à raison d'une comparaison de string toutes les 30 ns, cela représente un temps de calcul de l'ordre de $1 500 000 * 100 * 23 000 000 * 30 * 10^{-9} = 1,035 * 10^{8}$ s, c'est-à-dire plus de trois ans.
  
Pour résoudre ce problème de complexité, nous allons réaliser le mapping d'un read en le divisant en k-mers (sous-séquences de longueurs k) puis en localisant l'ensemble de ces k-mers efficacement sur la transormée de Burrows-Wheeler du génome. 
Cette stratégie repose sur la possibilité de calculer la BWT du génome via l'algorithme DC3 qui construit la table des suffixes en temps linéaire.  
Nous prendrons également en compte le fait q'un k-mer puisse être localisé sur le brin complémentaire inversé, à plusieurs endroits et/ou qu'il puisse y avoir une erreur de séquençage ou une petite mutation.  
Le reste des optimisations repose sur la mémorisation de variables globales permettant de ne pas faire les calculs à chaque itérations. On augmente raisonnablement la complexité spaciale pour diminuer grandement la complexité temporelle.

### I.2. Import des librairies :


```python
import math
from Bio import SeqIO #traitement des fichiers fastq
import pysam #traitement des fichiers bam
from collections import Counter
from bisect import bisect_left #pour la recherche dichotomique rapide
import time #permet de calculer les temps d'exécution
import cProfile #permet d'analyser finement les temps d'exécution de toutes les fonctions appellées
import json #pour stocker les résultats dans un fichier à part entière : resultats.json
```

### I.3. Import des données


```python
def read_fasta_file(fasta_path):
    fasta_seq = []
    reads = SeqIO.parse(fasta_path, "fasta")
    for read in reads:
        read_seq = str(read.seq)
        fasta_seq.append(read_seq)
    return fasta_seq
```


```python
end_of_string = '$'
```


```python
genome = read_fasta_file("GCF_000002765.5_GCA_000002765_genomic.fna")
for i in range(len(genome)):
    genome[i]=genome[i].upper()+end_of_string # on ne traitera pas ici le cas des minuscules c'est-à-dire des régions de forte simplicité. On ajoute directement le end_of_string nécessaire à la méthode de la BWT
len(genome)
```




    15



La variable _"genome"_ est la liste des 15 chaines de caractères correspondant aux 14 chromosomes de P.falciparum et à son apicoplaste (cf https://en.wikipedia.org/wiki/Plasmodium_falciparum).


```python
reads = SeqIO.parse("single_Pfal_dat.fq", "fastq") #Attention : le parcours des reads se fait en les enlevant de reads. Il faut donc bien réexécuter cette ligne à chaque utilisation de la variable. Une copie de cette ligne sera donc dans le main()
nb_reads = sum(1 for _ in reads)
print(nb_reads)
type(reads)
```

    1500000





    Bio.SeqIO.QualityIO.FastqPhredIterator



La variable _"reads"_ n'est pas convertie en liste pour des raisons de complexité spaciale.  
Les reads seront donc traités l'un après l'autre.  

Les principales variables de l'analyse (potentielement couteuses en temps de calcul) seront calculées puis stockées sous forme de variables globales. Cela ne pose pas de problème de conflit car elle ne dépendent que des données importées et ne changent pas au cours du notebook.  
Ci-dessous les premières d'entre-elles :


```python
alphabet = {'A', 'G', 'C', 'T'}
print(alphabet)
```

    {'A', 'C', 'T', 'G'}


### I.4. Création du Suffix Array

On calcule le tableau des suffixes de chaque chromosome du génome en temps linéaire (O(M) avec M la taille du génome complet) grâce à l'algorithme DC3 :


```python
def DC3(S):

    # Étape 0, Création de T
    T = []
    if len(S) == 0:
        return []
    if isinstance(S[0], str):  # On vérifie que S contient des chaînes
        for c in S:
            T.append(ord(c))
    else:
        T = S.copy()  # Si c'est déjà une liste d'entiers, on la copie
    T += [0, 0, 0]  


    # Étape 1 :

    # Création de P0, P1, P2 et P1_2
    P0, P1, P2 = [], [], []
    for i in range(len(T) - 2):
        if i % 3 == 1:
            P1.append(i)
        elif i % 3 == 2:
            P2.append(i)
        elif i % 3 == 0 and i < len(T) - 3:
            P0.append(i)
    P1_2 = P1 + P2

    # Création de R1_2
    R1_2 = []
    for p in P1_2:
        triplet = (T[p], T[p + 1], T[p + 2])
        R1_2.append(triplet)

    # Tri de R1_2 avec les positions en utilisant le nouveau critère de tri
    R1_2_with_p = list(zip(R1_2, P1_2))

    def custom_sort(triplet_p):
        triplet, p = triplet_p
        return triplet + tuple(T[p + 3:p + 10])

    R1_2_sorted = sorted(R1_2_with_p, key=custom_sort)

    # Création de Index1_2 et Order1_2
    rank = 0
    prev_key = None
    p_to_rank = {}
    order_list = []
    for triplet_p in R1_2_sorted:
        key = custom_sort(triplet_p)
        if key != prev_key:
            rank += 1
        triplet, p = triplet_p
        p_to_rank[p] = rank
        prev_key = key
        order_list.append(rank)
    Index1_2 = [p for (triplet, p) in R1_2_sorted]
    Order1_2 = order_list

    # Création de T'
    T_prime = [p_to_rank[p] for p in P1_2]

    # Vérification des duplicatas dans Order1_2 pour décider de la récursivité
    if len(set(Order1_2)) < len(Order1_2):  # Il y a des duplicatas
        Index_prime0_1_2 = DC3(T_prime)
        Index1_2 = [P1_2[i] for i in Index_prime0_1_2]
        for idx, p in enumerate(Index1_2):
            p_to_rank[p] = idx + 1


    # Étape 2 : Calcul de R0
    R0 = []
    for p in P0:
        pair = (T[p], p_to_rank.get(p + 1, 0))
        R0.append(pair)

    # Tri de R0
    R0_with_p = list(zip(R0, P0))
    R0_sorted = sorted(R0_with_p)
    Index0 = [p for (pair, p) in R0_sorted]


    # Étape 3 : Fusion de Index0 et Index1_2 pour obtenir le suffix array final

    # Fonctions auxiliaires pour la comparaison
    def leq(a1, a2, b1, b2):
        return a1 < b1 or (a1 == b1 and a2 <= b2)

    def leq3(a1, a2, a3, b1, b2, b3):
        return a1 < b1 or (a1 == b1 and leq(a2, a3, b2, b3))

    SA = []
    i = 0
    j = 0
    while i < len(Index0) and j < len(Index1_2):
        a = Index0[i]
        b = Index1_2[j]
        if b % 3 == 1:
            if leq(T[a], p_to_rank.get(a + 1, 0), T[b], p_to_rank.get(b + 1, 0)):
                SA.append(a)
                i += 1
            else:
                SA.append(b)
                j += 1
        else:
            if leq3(T[a], T[a + 1], p_to_rank.get(a + 2, 0), T[b], T[b + 1], p_to_rank.get(b + 2, 0)):
                SA.append(a)
                i += 1
            else:
                SA.append(b)
                j += 1

    # Ajout des indices restants
    while i < len(Index0):
        SA.append(Index0[i])
        i += 1
    while j < len(Index1_2):
        SA.append(Index1_2[j])
        j += 1

    # Suppression des indices des sentinelles
    SA = [idx for idx in SA if idx < len(T) - 3]
    stop = time.time()
    return SA
```


```python
a = time.time()
l_SA = [DC3(chromosome) for chromosome in genome] # ici on ne rajoute pas le end_of_string = '$'
b = time.time()
print(b-a)
```

    132.62179613113403


### I.5. Calcul de la BWT du génome

Le calcul efficace de la BWT se fait grâce au tableau des suffixes calculé par l'algorithme DC3.


```python
def BWT(T,SA=None): #on a déjà rajouté le end_of_string
    """
    Compute the BWT from the suffix table (O(n))

    Args:
        T (str): string
        end_of_string (char): end of string character to append

    Return:
        bwt (str): BWT
    """
    if SA == None:
        SA = DC3(T)
    bwt=""
    for i in SA:
        bwt += T[i-1]
    return(bwt)
```


```python
l_BWT = [BWT(genome[i],l_SA[i]) for i in range(len(genome))]
```

### I.6. Construction des kmers

Pour un read donné, on construit la liste de ses kmers simplement comme la liste des sous-mots de longueur k.


```python
def build_kmers(read,k):
    """
    construire dictionaire des indexes de kmers (O(T-k))
    :param text: genome complet
    :param k: longeur de kmers
    :return: dictionaire des sous-textes et leurs indices de kmers
    #complexite = O(n)
    """
    kmers_list  = []
    for i in range(len(read)-k+1):
        kmer = read[i:i+k]
        kmers_list.append(kmer)
    return kmers_list
```

### I.7. Mapping des kmers

Pour réaliser le mapping des kmers de manière efficace, on calcul et on mémorise deux variables globales supplémentaires :
- $count$ qui est le dictionnaire dont les clés sont les lettres de l'alphabet et les valeurs le nombre de lettre plus petite dans un chromosome
- $occurrence$ qui est le dictionnaire dont les clés sont les lettres de l'alphabet et les valeurs les listes dont l'élément i est le nombre de lettre dans BWT[:(i+1)]  

On mémorise ces deux variables pour chaque chromosome sous la forme deux de listes l_count_smaller et l_occurrences (comme pour l_SA et l_BWT)

Les fonctions de mapping efficaces des kmers qui suivent ont été inspirées par la correction de TP : Burrows-Wheeler-transform.ipynb de Sergio Peignier (https://sergiopeignier.github.io/string_search.html) et la publication github suivante : https://github.com/kemaleren/bwt



```python
def compute_count_smaller(s):
    """
    compute the number of all the smaller characters in the string (O(n))
    :param s (str): 
    :return: 
    """
    K = Counter(s)
    cumulated = 0
    count_smaller = {}
    for letter in sorted(alphabet):
        count_smaller[letter] = cumulated
        cumulated += K[letter]
    return count_smaller
```


```python
l_count_smaller = [compute_count_smaller(L) for L in l_BWT]
```


```python
def compute_occurrences(bwt):
    """
    return the occurences of each character in the BWT for each index (O(n))
    :param bwt: BWT of the genome
    :return: 
    """
    occurrences = {letter : [0] for letter in alphabet}
    occurrences[bwt[0]] = [1]
    for letter in bwt[1:]:
        for k, v in occurrences.items():
            v.append(v[-1] + (k == letter))
    for k, v in occurrences.items():
        v.extend([v[-1], 0]) #traitement des cas de dépassement d'indices dans la fonction pattern_matching_bwt
    return occurrences
```


```python
l_occurrences = [compute_occurrences(L) for L in l_BWT]
```


```python
def update_range(e, f, letter, occ, count_smaller):
    """
    update the range of the substring search based on the occurences dictionary (O(1))
    :param e: begin of the range
    :param f: end of the range
    :param letter: character to search for
    :param occ: occurences dictionary
    :param count_smaller: dictionary of number of occurences of smaller letters
    :return: new begin of the range and new end of the range
    """
    #new_e = count_smaller[letter] + occ[letter][e - 1] + 1
    #new_f = count_smaller[letter] + occ[letter][f]
    a=count_smaller[letter]
    return a + occ[letter][e - 1] + 1, a + occ[letter][f]

def pattern_matching_bwt(kmer, n, occ, count_smaller, SA):
    """
    substring search with bwt (O(k))
    :param kmer: 
    :param n: 
    :param occ: 
    :param count_smaller: 
    :param SA: 
    :return: 
    """
    results = []
    e, f = 0, n - 1
    for letter in reversed(kmer):
        e, f = update_range(e, f, letter, occ, count_smaller)
        if e > f:
            return []
    results.extend(SA[e:f+1])
    return list(dict.fromkeys(SA[e:f + 1])) #plus rapide que sorted(set(results))
```

Maintenant que l'on peut mapper un kmer efficacement, on peut construire la liste des maps de tous les kmers d'un read :


```python
def query_kmers(read, n, k, occ, count_smaller, SA):
    """
    query all the kmers of the sequences on the genome (O(kT))
    :param read: pattern for query
    :param n: 
    :param k: 
    :param occ: 
    :param count_smaller: 
    :param SA: 
    :return: 
    """
    kmers_list = build_kmers(read,k)
    res = []
    for kmer in kmers_list:
        offset = pattern_matching_bwt(kmer, n, occ, count_smaller, SA)
        res.append(offset)
    return res
```

### I.8. Analyse des mapping

On utilise une heuristique pou analyser les listes des maps des kmer d'un read.  
On se fixe une variable erreurs_max puis on parcours les différentes positions de départ possibles du read sans dépasser ce seuil (pas forcément en partant des positions du premier kmer qui map). Ensuite, on teste si les kmers suivants mappent aux positions suivantes. On prend alors la meilleure position de départ (si il en existe une sous le seuil d'erreur) et on calcul une distance de hamming entre le read et la portion du génome pour réaliser un calcul exact du pourcentage d'identité.


```python
def hamming_distance(seq1, seq2):
    return sum(1 for a, b in zip(seq1, seq2) if a != b)
```


```python
def is_in_sorted_list(sorted_list, x):
    idx = bisect_left(sorted_list, x) # cette méthode avec la fonction bisect_left est un peu plus rapide qu'une rechecher dichotomique classique
    return idx < len(sorted_list) and sorted_list[idx] == x

def cherche_succession(idx_kmer, courant, erreurs, k, erreurs_max ,res): #O((T)log(T))
    while idx_kmer < len(res) - 1: # on parcours tous les maps 
        if erreurs > erreurs_max: # cas où on a dépassé l'erreur max
            return None
        if is_in_sorted_list(res[idx_kmer + 1],courant+1):   # cas où le k_mer suivant a mappé
            idx_kmer += k
            courant += k
        else:   # cas où le k_mer suivant n'a pas mappé et donc qu'il y a un problème à la fin de celui ci
            idx_kmer += k # on saute les k kmers suivant qui pourrait comporter d'autres erreurs mais on cherche ici une approximation avant de calculer la distance.
            courant += k
            erreurs += 1
    return erreurs

def analyse(read_met,L,k,chromosome,occ,count_smaller,SA,erreurs_max=4):
    """
    O(nT^2log(T))
    :param read_met: 
    :param L: 
    :param k: 
    :param chromosome: 
    :param occ: 
    :param count_smaller: 
    :param SA: 
    :param erreurs_max: 
    :return: 
    """
    read = str(read_met.seq)
    res = query_kmers(read, len(L), k, occ, count_smaller, SA)
    n = len(res)

    erreurs = 0
    idx_res_first = 0
    # Trouver le premier k-mer qui map en comptant les erreurs
    while idx_res_first < n and not res[idx_res_first]:
        idx_res_first += k # car il se peut que l'erreur soit au k-ieme nucléotide et on ne voudrait pas la compter k fois : la variable erreurs compte le nombre d'erreur minimal de la succession en cours de construction
        erreurs += 1
        if erreurs > erreurs_max:
            return 0,None

    best_erreurs = erreurs_max
    best_start = None
    
    # Parcours des positions de départ possibles
    max_idx_res = min(idx_res_first + erreurs_max - erreurs, n)
    for i in range(idx_res_first, max_idx_res): #pire cas : T-k+1
        for courant in res[i]: # pire cas : n-T+1
            result = cherche_succession(i, courant, erreurs+i, k, erreurs_max ,res) #O((T)log(T))
            if result != None:
                current_erreurs = result
                if current_erreurs < best_erreurs:
                    best_erreurs = current_erreurs
                    best_start = courant-i
                    
    if best_erreurs >= erreurs_max:
        return 0,None
    # Une fois que l'on pense avoir trouver un candidat par les approximations ci-dessus, on réalise une vraie 
    # distance pour vérifier (il ne faut donc pas qu'elle soit réalisée trop souvent pour des problème de complexité)
    m=len(read)
    identity = int(100*(m - hamming_distance(read,chromosome[best_start:best_start+m]))/m) #O(log(T))
    return identity,best_start
```

### I.9. Fonction principale et construction du tableau des résultats

La taille des kmers a été déterminé empiriquement par tests. Un trop petit k impliquerait que les kmers mappent trop souvent. Inversement, un trop grand k impliquerait un trop petit nombre de matches et donc de ne rien détecter.

Toujours pour des raisons de complexité temporelle, lors du parcours des chromosomes dans la fonction main, on va considérer un read seulement s'il n'a pas déjà été bien mapé sur un chromosome précédent (ie à plus de 98% d'identité)

Les résultats seront stocker dans un dictionnaire similaire à celui correspondant au fichier bam (cf. partie III) mais aussi dans un fichier à part entière nommé $resultats.json$.  
(Une copie de ce fichier des résultats nommé $resultats\_copie.json$ est déjà enregistré dans le dossier de ce notebook)


```python
def main(genome=genome,k=20):
    resultats_dict = {}
    map_ok=[False for i in range(nb_reads)]
    for i in range(len(genome)): #15 chromosomes 
        start = time.time()
        reads = SeqIO.parse("single_Pfal_dat.fq", "fastq") #Lors du parcours des reads, ceux-ci sont supprimés de la variable. On la redéfinit donc pour chaque chromosome
        chromosome = genome[i]
        SA = l_SA[i]
        L = l_BWT[i]
        occ = l_occurrences[i]
        count_smaller = l_count_smaller[i]
        j=0
        for read in reads: #N reads => O(nNT^2log(T))
            name = read.name
            if not map_ok[j] :  #Si le read a déjà été correctement mappé, on ne cherche pas d'autres positions
                identity,position = analyse(read,L,k,chromosome,occ,count_smaller,SA) #O(nT^2log(T))
                if identity>=98:
                    map_ok[j]=True
                if not name in resultats_dict.keys() or identity>resultats_dict[name]["identity"] :
                    resultats_dict[name] = {"sequence": str(read.seq),      # NB : name == id donc on considère seulement le name qui est donc un identifiant
                                            "position": position,
                                            "identity": identity,
                                            "chromosome": i,
                                            "reversed_strand": False}
                identity,position = analyse(read.reverse_complement(),L,k,chromosome,occ,count_smaller,SA)
                if not name in resultats_dict.keys() or identity>resultats_dict[name]["identity"] :
                    resultats_dict[name] = {"sequence": str(read.reverse_complement().seq),
                                            "position": position,
                                            "identity": identity,
                                            "chromosome": i,
                                            "reversed_strand": True}
            j+=1
        t = time.time()-start
        if i==14:
            print("apicoplaste traité en",int(t//60),"min",int(round((time.time()-start)%60,0)),"s")
        else :
            print("chromosome",i+1,"traité en",int(t//60),"min",int(round((time.time()-start)%60,0)),"s")
    
    #enregistrement des résultats dans un fichier à part entière
    with open("resultats.json", "w", encoding="utf-8") as fichier:
        json.dump(resultats_dict, fichier, indent=4, ensure_ascii=False)
        
    return resultats_dict
```

### I.10. Calcul de la complexité pire-cas de l'algorihme

On avait noté :  
$N\_r=$ le nombre de reads  
$T\_r=$ la taille d'un read  
$T\_g=$ la taille du génome complet  
$k=$ la longeur d'un kmer  
On note également :  
$N\_c=$ le nombre de chromosome
$T\_c=$ la taille moyenne d'un chromosome


Pour des raisons de praticité, des variables globales ont été pré-calculées, en particulier les suffix-arrays (l_SA) et les BWT (l_bwt) de chaque chromosome. Ici, les complexités des fonctions qui calculent ces variables sont DC3 en $O(N\_c \cdot T\_c) = O(T\_g)$ et BWT en $O(N\_c \cdot T\_c) = O(T\_g)$.  
Aussi, les variables l_count_smaller et l_occurrences ont été calculées par les fonctions $compute_count_smaller()$ en $O(N\_c \cdot T\_c) = O(T\_g)$ et $compute_occurrences()$ en $O(N\_c \cdot T\_c) = O(T\_g)$.
Les précalculs admettent donc une complexité temporelle en $O(T\_g)$.  

La fonction $pattern\_matching\_bwt()$ admet une complexité temporelle en $O(k \cdot log(T\_c))$ donc la fonction $query_kmer()$ est en $O(T\_r \cdot k \cdot log(T\_c))$.  
De plus, la $fonction chercher\_succession()$ a une complexité de $O(T\_r \cdot log(T\_r))$.  
Ainsi, la fonction $analyse()$ admet une complexité de $O(k \cdot T\_r^3 \cdot log(T\_r) \cdot log(T\_c))$.  
La fonction $main()$ parcourt dans le pire des cas les $N\_r$ reads pour chaque chromosome et appelle la fonction $analyse()$ à chaque itération. Elle a donc une complexité temporelle pire cas en $O(N\_c \cdot N\_r \cdot k \cdot T\_r^3 \cdot log(T\_r) \cdot log(T\_c))$.
Le nombre de chromosome étant majoré, on peut réécrire cette complexité comme un $O(N\_r \cdot k \cdot T\_r^3 \cdot log(T\_r) \cdot log(T\_g))$.

## II. Exécution sur le jeu de données simulé et évaluation du temps de calcul


```python
start = time.time()
resultats = main()
stop = time.time()
```

    chromosome 1 traité en 14.97 min
    chromosome 2 traité en 18.06 min
    chromosome 3 traité en 18.62 min
    chromosome 4 traité en 17.93 min
    chromosome 5 traité en 20.48 min
    chromosome 6 traité en 18.13 min
    chromosome 7 traité en 16.31 min
    chromosome 8 traité en 16.44 min
    chromosome 9 traité en 16.65 min
    chromosome 10 traité en 15.33 min
    chromosome 11 traité en 17.19 min
    chromosome 12 traité en 16.03 min
    chromosome 13 traité en 15.83 min
    chromosome 14 traité en 14.51 min
    chromosome 15 traité en 1.51 min
    238.08559790849685 min d'exécution au totale



```python
print(int((stop-start)//3600),'h',int(((stop-start)/60)%60),"min d'exécution au totale")
```


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    Cell In[22], line 1
    ----> 1 print(int((stop-start)//3600),'h',int(((stop-start)/60)%60),"min d'exécution au totale")


    NameError: name 'stop' is not defined



```python
#cProfile.run("main()",sort="cumulative") #permet de calculer le temps d'exécution de toutes les fonctions appelées si besoin.
```

## III. Résultats obtenus 

### III.1. Explication du protocole expérimental de comparaison à la vérité terrain (fichier BAM).

Le protocole expérimental de comparaison à la vérité terrain passe par la lecture des résultats attendus notés dans le fichier BAM. Pour cela, lors du mapping de chaque read, on note dans un tableau résultat le numéro du read, sa position, sa qualité ainsi que s'il se trouve sur le brin complémentaire ou non.  
Les reads sont comparés avec la référence  en calculant la précision, la sensiblité et donc le F1-score. Ici, on considère que si le read a été bien trouvé à une position correct mais le nombre d'identité est faux, le read est un faux négatif. Il faut donc exact au niveau de position et de qualité (identity) pour être considéré comme un vrai positif.


```python
def read_bam_file(bam_path):
    """
    lecture de fichier bam (O(1))
    :param bam_path: 
    :return: dictionaire des informations de chaque read dans le fichier bam
    """
    def calculate_identity(cigar):
        matches = 0
        mismatches = 0
        total_length = 0
        import re
        operations = re.findall(r"(\d+)([=X])", cigar)

        for length, op in operations:
            length = int(length)
            if op == "=":
                matches += length
            elif op == "X":
                mismatches += length

            total_length += length

        
        if total_length > 0:
            identity = int(100* matches / total_length)
            return identity
        else:
            return 0
    bam_info = {}
    bam_name = []
    bam_file = pysam.AlignmentFile(bam_path, "rb")
    for read in bam_file:
        read_name = read.query_name
        read_seq = read.query_sequence
        read_position = read.reference_start
        read_is_reverse = read.is_reverse
        read_chromosome = read.reference_name
        cigar = read.cigarstring

        align_length = read.query_alignment_length


        bam_info[read_name] = {"sequence": read_seq,
                               "position": read_position,
                               "identity": calculate_identity(cigar),
                               "chromosome": read_chromosome,
                               "reversed_strand": read_is_reverse}  #true if seq is reversed trand
        bam_name.append(read_name)
    bam_file.close()
    return bam_info,bam_name
```


```python
def test_result(resultats_dict, bam_path):
    """
    Analyse des résultats trouvés face au fichier BAM 
    
    Args:
        reads : object iterable de fastq
        resultats_dict : dictionnaire des résultats
        bam_path (str): lien du fichier BAM
    """
    
    bam_info, bam_name = read_bam_file(bam_path)
    TP = FP = FN = 0
    TTP =quality_faux =0 #correct position, false in quality
    for read_name,read_res in resultats_dict.items() :
        if read_res["position"]: #if read is found
            
            if read_res["position"] == bam_info[read_name]["position"]:
                TP += 1
                TTP += 1
                if read_res["identity"] != bam_info[read_name]["identity"]:
                    FP += 1
                    TP -=1
                    quality_faux += 1
            else:
                FP+=1
        elif not read_res["position"] and bam_info[read_name]["position"]:
            FN +=1
    precision = TP/(TP+FP)
    recall = TP/(TP+FN)
    f1_score = 2*precision*recall/(precision+recall)
    print(f"Nombre de position correcte trouvée : {TTP} ({TTP/len(resultats_dict)*100:.2f}%), où il y a {TP} ({TP/len(resultats_dict)*100:.2f}%) vrais positives")
    print(f"Nombre de faux positive : {FP} ({FP/len(resultats_dict)*100:.2f}%), dont {quality_faux} ({quality_faux/len(resultats_dict)*100:.2f}%) sont trouvés à bonne position mais faux en l'identité")
    print(f"Nombre de read pas trouvé (faux négative) : {FN} ({FN/len(resultats_dict)*100:.2f}%), alors qu'ils existe dans le genome")
    print(f"les mappes donne un résultats avec la précision de {precision} et la sensibilité {recall} dans {len(resultats_dict)} éléments")
    print(f"F1 score : {f1_score}")
    print("Note : Les pourcentages sont par rapport au nombre total de read ")
    return precision, recall, f1_score

```

### III.2. Analyse des résultats


```python
bam_info,bam_name = read_bam_file("single_Pfal_dat.bam")
```


```python
print(bam_info["NC_037283.1-72168"]["position"])
print(bam_info["NC_037283.1-72168"]["identity"])
```

    2862965
    98



```python
resultats_path = "./resultats_copie.json"
with open(resultats_path, "r", encoding="utf-8") as fichier:
    resultats_dict = json.load(fichier)
```


```python
comparaison = test_result(resultats_dict, "single_Pfal_dat.bam")
```

    Nombre de position correcte trouvée : 1233597 (82.24%), où il y a 1233514 (82.23%) vrais positives
    Nombre de faux positive : 103549 (6.90%), dont 83 (0.01%) sont trouvés à bonne position mais faux en l'identité
    Nombre de read pas trouvé (faux négative) : 162934 (10.86%), alors qu'ils existe dans le genome
    les mappes donne un résultats avec la précision de 0.9225548833525421 et la sensibilité 0.8833225440546301 dans 1500000 éléments
    F1 score : 0.9025125561960424
    Note : Les pourcentages sont par rapport au nombre total de read 


## IV. Conclusions

Notre algorithme a significativement amélioré la complexité de la recherche de séquences dans les grands génomes par rapport aux méthodes naïves. Il est capable d'effectuer un mapping correct pour jusqu'à 80 % (précision à 0.9) des séquences. En revanche, cette méthode ainsi que l'analyse des résultats des k-mers sont approximatives, ce qui engendre une variance dans les résultats de recherche.  
Malgré une complexité temporelle en $O(N\_r \cdot k \cdot T\_r^3 \cdot log(T\_r) \cdot log(T\_g))$, la fonction $main()$ s'exécute en environ 4h. Cela peut s'expliquer par des grandes valeurs de constantes dans ce "grand O". Une piste d'amélioration serait donc de faire diminuer ces constantes via des optimisations de la structure du code ou des types d'objet utilisé.
