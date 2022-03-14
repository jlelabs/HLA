import pandas as pd
import sys
from hla_lookup_tables import flagging

output_path = sys.argv[1] + "/"
sample_name = sys.argv[2]
alternative_path = sys.argv[3] + "/"
patient_path = output_path + sample_name + "_report.d8.txt"
alternative_patient_path = alternative_path + sample_name + "_Top_2_Results.txt"

try:
    flagging(patient_path, output_path, sample_name)
except:
    alternative_output = pd.read_csv(alternative_patient_path, sep="\t", header=None).iloc[:,:1]
    gene_column = alternative_output.iloc[:,0].apply(lambda x: x.split("*")[0]).tolist()
    alternative_output.insert(0, "Gene", gene_column)
    allele1 = alternative_output.iloc[::2]
    allele2 = alternative_output.iloc[1::2]
    alternative_output = allele1.merge(allele2,on="Gene")
    alternative_output.columns = ["Gene", "Allele1", "Allele2"]
    alternative_output.to_csv(output_path + sample_name + "_top2.txt", sep="\t", index=False)
    print(alternative_output)
    flagging(alternative_output, output_path, sample_name + "_TOP2_")

    # for fillna
    filled_na = pd.read_csv(patient_path, sep="\t").fillna("")
    flagging(filled_na, output_path, sample_name + "_FILLED_NA_")
