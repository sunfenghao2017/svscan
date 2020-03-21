|Column|Explanations
|-------|----------
|svid| sv id of this fusion, internal usage
|gene1| partner gene1 name in fusion gene
|gene2| partner gene2 name in fusion gene
|srcnt| all split read seed count supporting gene1->gene2 fusion(```measured in read```)
|dpcnt| all discordant pairs reads supporting  gene1->gene2 fusion(```measured in read```)
|mocnt| all molecules supporting gene1->gene2 fusion(```measuered in molecule```)

### Attention 
This sheet just list all number of seed read/pair supporting fusion events in COSMIC/ONCOKB, but didn't take the precise breakpoint of fusion event into account. That is to say, there may be multiple breakpoint events of gene1->gene2, and they are summarized together in this sheet, this sheet can give you some clue about gene1->gene2 fusion, but should not be used as final clinical test result.
