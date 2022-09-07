import pandas as pd

print("ciao")

names = ["analysis/genome-miRNAcounts/70459-taggedBAMcounts.txt","analysis/genome-miRNAcounts/70460-taggedBAMcounts.txt","analysis/genome-miRNAcounts/70461-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74665-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74666-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74667-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74668-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74669-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74670-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74671-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74672-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74673-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74674-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74675-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74676-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74677-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74678-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74679-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74680-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74681-taggedBAMcounts.txt","analysis/genome-miRNAcounts/74682-taggedBAMcounts.txt"]

ds_list = []

for name in names:
    tmp = pd.read_csv(name)
    print(tmp.head())
    ds_list.append(tmp)


# result = pd.concat(ds_list, axis=1, join='inner')
# print(result.head(5))
# result.to_csv("genome-counts-BR.csv")
# print(ds_list)