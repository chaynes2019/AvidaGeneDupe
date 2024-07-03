def generateCodingSitesOfOrganism(organism):
    createAnalyzeCfg(organism)
    executeInfoAnalysis()

    #organism is a row of a pandas dataframe
    codingSites = getCodingSites(organism)

    organism["NOT Coding Sites"] = codingSites[0]
    organism["NAND Coding Sites"] = codingSites[1]
    organism["AND Coding Sites"] = codingSites[2]
    organism["ORNOT Coding Sites"] = codingSites[3]
    organism["OR Coding Sites"] = codingSites[4]
    organism["ANDNOT Coding Sites"] = codingSites[5]
    organism["NOR Coding Sites"] = codingSites[6]
    organism["XOR Coding Sites"] = codingSites[7]
    organism["EQUALS Coding Sites"] = codingSites[8]