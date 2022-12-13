@ECHO OFF
ECHO Wilcome to this small uGene example.
python outToCSV.py file=exampleData.out -ncbi=[{'col':'ncbiID','db':'taxonomy','path':['LineageEx.[Rank=kingdom].ScientificName','LineageEx.[Rank=clade].ScientificName','LineageEx.[Rank=phylum].ScientificName','LineageEx.[Rank=order].ScientificName','LineageEx.[Rank=family].ScientificName','LineageEx.[Rank=genus].ScientificName'],'rename':['kingdom','clade','phylum','order','family','genus']}]
python main.py file=exampleData.csv -t=[{'dev_report':1, 'x_axis':'geneID', 'y_axis':'ncbiID', 'values':['FAS_F','FAS_B'],'jobs':'tax'}, {'dev_report':1, 'x_axis':'ncbiID', 'y_axis':'geneID', 'values':['FAS_F','FAS_B'],'jobs':'gene'}]
ECHO Done !
PAUSE