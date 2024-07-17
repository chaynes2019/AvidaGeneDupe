with open("NOPFile.txt", 'w') as f:
    genomeStr = "acAbbbcaAAAaAAaaAbabAcbAAbcbaaabbbbaacbbbababbcAbacAcbAccAbAAAbabaAAAccacaAaaaaabAAaaAbabAAaabAAcaAA"
    for k in range(len(genomeStr)):
        if genomeStr[k] == 'a':
            f.write('nop-A      #\n')
        elif genomeStr[k] == 'b':
            f.write('nop-B     #\n')
        elif genomeStr[k] == 'c':
            f.write('nop-C     #\n')
        elif genomeStr[k] == 'A':
            f.write('nop-X     #\n')
        else:
            print("Ludicrous!")
