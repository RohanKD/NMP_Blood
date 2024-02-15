from tqdm import tqdm


def results(loc_1, loc_2, scf1,scf2):

    def score(gene, loc1, loc2):
        return loc1[gene] / scf1 - loc2[gene] / scf2

    # internal scale factor = 202.9
    # external scale factor = 5.231
    # with open("genes_.txt", "r") as f:
    #     for line in f:
    #         gene_val = line.split(" - ")
    #         gene_val[1] = float(gene_val[1].replace("\n", ""))
    #         external.append(gene_val)
    # with open("genes.txt", "r") as f:
    #     for line in f:
    #         gene_val = line.split(" - ")
    #         gene_val[1] = float(gene_val[1].replace("\n", ""))
    #         internal.append(gene_val)
    # subtract internal/scf - external/scf (scf = scale factor)
    # rank by score and create third text file with ranked score using list comp method
    # external/internal format = [[gene1, val1], [gene2, score2], ... , [gene3, score3]]
    scores = []
    internal = loc_1.keys()
    for g in tqdm(internal):
        scores.append((g, score(g, loc_1, loc_2)))

    # sort scores
    scores = sorted(scores, key=lambda tup: tup[1], reverse=True)
    # output score to text file for GUI to utilize
    result = []

    for score in scores:
        result.append(str(score[0]) + " - " + str(round(score[1], 3)) + "\n")
    return result


if __name__ == "__main__":
    results()
