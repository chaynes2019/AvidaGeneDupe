makeFileParameters = getMakeFileParameters("analysisMakeFile.txt")

'''
Collect the different data files over all considered replicates.
'''
passedTreatments = makeFileParameters["Treatments to Analyze"]
lineageFiles = collectDataFiles(passedTreatments)

'''
For each, extract data from each line to be put into a Pandas dataframe.
While iterating over each lineage, find slip-duplication events and task creation events; record them,
potentially for lack of a better spot, in the dataframe row of each member of that lineage.
The templates must also be recorded.
'''
analysisDataframe = pd.DataFrame()
for lineageFile in lineageFiles:
    extractDataFileToDataFrame(lineageFile, analysisDataframe)
    recordSlipDuplicationEvents(analysisDataframe)
    recordTaskCreationEvents(analysisDataframe)

'''
For each line in the Pandas dataframe, generate locus history and add that information to the dataframe.
'''
generateLocusHistory(analysisDataframe)

'''
For each line in the Pandas dataframe, generate coding sites and then replace that line with 9 lines, one for each task, detailing the coding sites.
'''
generateCodingSites(analysisDataframe)

'''
Count up the number of template sites that are coding and non-coding sites;
count up the number of de novo task coding sites that originated in
duplications.
'''
