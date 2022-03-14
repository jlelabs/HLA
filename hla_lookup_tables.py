import pandas as pd
import sys

# Database of Disease Associated HLA risk cases
path_to_reference = sys.path[0] + "/Configuration"
HLA_Database= pd.ExcelFile(path_to_reference + '/HLA_Database_August_6.xls')
hla_df = HLA_Database.parse("Sheet1")

# HLAVBSeq results

def make_allele_table(patient_path):
    """ Creates allele table
    """
    # output table
    if isinstance(patient_path, str):
        allele_table = pd.read_csv(patient_path, sep="\t")
    else:
        allele_table = patient_path

    # subset alleles by genes that are within the HLA cases
    allele_table = allele_table[allele_table["Gene"].isin(["A","B", "C", "DQA1", "DQB1", "DRB1", "DPB1"])]
    allele_table["Gene"] = "HLA-" + allele_table["Gene"]
    allele_table.set_index(["Gene"], inplace=True)
    allele_table["Risk"] = "Average Risk"
    return allele_table


def make_test_patient_table(allele_table):
    """ Create test patient allele table
    """
    # testing table
    test_patient_allele_table = pd.DataFrame()

    # create row for each allele
    test_patient_allele_table["Allele2"] = allele_table["Allele1"]
    test_patient_allele_table = test_patient_allele_table.append(allele_table[["Allele2"]])
    test_patient_allele_table = test_patient_allele_table.sort_index()
    test_patient_allele_table.columns = ["Allele"]

    # indicators for test cases
    test_patient_allele_table["Risk"] = "None"
    test_patient_allele_table["HLA ID"] = "None"

    # remove gene as index
    test_patient_allele_table = test_patient_allele_table.reset_index(drop=False)
    return test_patient_allele_table




# create dictionary with Disease Associated HLA risk cases
# each key will have a unique case with a dictionary with keys of alleles that must be present
def conditions():
    """ Creates dictionary with Disease Associated HLA Risk cases
    """
    ret = {}
    for index, row in hla_df.iterrows():
        a = row["Allele/Haplotype/Genotype_ Assignment"]
        a = a.replace("(","")
        a = a.replace(")", "")
        a = a.replace(" ","")
        a = a.replace("x1","")
        all_alleles = a.split("-")
        haplotype = []
        for i in range(len(all_alleles)):
            a = all_alleles[i]
            if "/" in a:
                a = a.split("/")
                haplotype.append(a[0])
                haplotype.append(a[1])
            else:
                haplotype.append(a)

        remove_alleles = []
        double_alleles = []
        for allele in haplotype:
            if allele[len(allele)-2:] == "x2":
                double_alleles.append(allele[0:len(allele)-2])
                remove_alleles.append(allele)
        for i in range(len(remove_alleles)):
            haplotype.remove(remove_alleles[i])
            haplotype.append(double_alleles[i])
            haplotype.append(double_alleles[i])
        new_alleles=[]
        for a in haplotype:
            new_a = a
            if "or" in a:
                new_a = a.split("or")
                new_a[1] = new_a[0][0:4] + new_a[1]
            new_alleles.append(new_a)

        lookup = {}
        for i in range(len(new_alleles)):
            if isinstance(new_alleles[i], str):
                new_alleles[i] = [new_alleles[i]]

        new_alleles.sort(key = len)
        for i in range(len(new_alleles)):
            a = new_alleles[i]
            if isinstance(a, str) != True:
                lookup_list = {}
                for j in a:
                    lookup_list[j] = ""
                lookup[i] = lookup_list
            else:
                lookup[i] = {a:""}
        ret[hla_df.iloc[index,0]] = lookup
    return ret




# flagging function
def flagging(patient_path, output_path, sample_name):
    """ Exports the three csvs for patient alleles, alleles and associated risks
    """

    risk_case_dict = conditions()
    allele_table = make_allele_table(patient_path)
    test_patient_allele_table = make_test_patient_table(allele_table)

    # returning list of HLA risk case ID's
    flag = []

    # create dictionary of all alleles present in sample
    sample_alleles = {}

    for j in range(len(test_patient_allele_table["Allele"])):
        # retrieve allele and locus to test against given HLA risk case condition dictionary
        allele = test_patient_allele_table["Allele"][j]

        if allele in sample_alleles:
            sample_alleles[allele] = 2
        else:
            sample_alleles[allele] = 1


    # iterate through each HLA risk case
    for c in risk_case_dict.keys():

        # all conditions in given HLA risk case
        conditions_dict = risk_case_dict[c]

        # length of the conditions in given HLA risk cases
        condition_length = len(conditions_dict)

        # counter to keep track of all conditions met in HLA risk case
        l = 0

        # dictionary of alleles already captured in counting for HLA risk case, prevent double counting
        hla = {}

        # dictionary of indices for test_patient_allele_table and allele_table to indicate risk
        indices= {}

        # make a copy of sample_alleles
        c_sample_alleles = sample_alleles.copy()


        # iterate through all the conditions
        for i in range(len(conditions_dict)):

            # iterate through all the HLA alleles
            for j in range(len(test_patient_allele_table["Allele"])):

                # retrieve allele and locus to test against given HLA risk case condition dictionary
                allele = test_patient_allele_table["Allele"][j]
                locus = allele.split("*")[0]

                # all HLA risk alleles are either seven or ten characters long
                seven_letter = allele[0:min(7, len(allele))]
                ten_letter = allele[0:min(10, len(allele))]

                # search seven letter abbreviation against given HLA risk case condition dictionary
                if seven_letter in conditions_dict[i]:

                    # check if the seven letter abbreviation is already accounted for in given HLA risk case
                    if c_sample_alleles[allele] > 0:
                        copies = c_sample_alleles[allele] - 1
                        c_sample_alleles[allele] = copies
                        l += 1
                        indices[j] = ""
                        break

                # search ten letter abbreviation against given HLA risk case condition dictionary
                elif ten_letter in conditions_dict[i]:
                    # check if the seven letter abbreviation is already accounted for in given HLA risk case
                    if c_sample_alleles[allele] > 0:
                        copies = c_sample_alleles[allele] - 1
                        c_sample_alleles[allele] = copies
                        l += 1
                        indices[j] = ""
                        break


        # check if all conditions are met for given risk case
        if condition_length == l:

            # if met, append case id to return
            flag.append(c)

            # for each condition, iterate over to indicate for each HLA allele risk
            for m in indices.keys():

                # generate gene name
                gene = test_patient_allele_table.iloc[m,0].split("*")[0]

                # indicate risk in returning csv
                test_patient_allele_table.iloc[m, 2] = "Risk"
                allele_table.at[gene, "Risk"] = "Increased Risk"

                # append to other case id's if one is existing for the HLA alelle
                patient_hla_id = test_patient_allele_table.iloc[m, 3]

                if patient_hla_id == "None":
                    test_patient_allele_table.iloc[m, 3] = c
                else:
                    patient_hla_id = patient_hla_id + " " + c
                    test_patient_allele_table.iloc[m, 3] = patient_hla_id

    # return list of patient risks by case ID
    patient_risks = flag

    # match patient risk cases according to their ID and create dataframe
    patient_associated_disease_table = hla_df[hla_df["HLA_ID"].isin(patient_risks)]
    patient_associated_disease_table = patient_associated_disease_table[["HLA_ID", "Associated_Disease", "Locus", "Allele/Haplotype/Genotype_ Assignment", "Risk ", "OR[95%CI]_or_RR,_p-value,_references*", "Interpretation"]]

    # change column names
    patient_associated_disease_table.columns = ["HLA ID", "Associated Disease", "Locus", "Allele/Haplotype/Genotype Assignment", "Risk", "OR[95%CI] or RR, p-value, references*", "Interpretation"]
    patient_associated_disease_table.set_index(["HLA ID", "Associated Disease", "Locus"], inplace=True)

    allele_table = allele_table.rename(columns={"Allele1":"Allele 1", "Allele2":"Allele 2"})
    test_patient_allele_table.set_index(["Gene", "Allele"], inplace=True)

    # reformat allele table
    hla_abc = allele_table[allele_table.index.isin(["HLA-A", "HLA-B", "HLA-C"])].copy()
    hla_abc["Allele 1"] = hla_abc["Allele 1"].apply(lambda x: x[:7])
    hla_abc["Allele 2"] = hla_abc["Allele 2"].apply(lambda x: x[:7])
    hla_other = allele_table[allele_table.index.isin(["HLA-DQA1", "HLA-DQB1", "HLA-DPB1", "HLA-DRB1"])].copy()
    hla_other["Allele 1"] = hla_other["Allele 1"].apply(lambda x: x[:10])
    hla_other["Allele 2"] = hla_other["Allele 2"].apply(lambda x: x[:10])
    allele_table = pd.concat([hla_abc, hla_other])


    # print the tables to the log
    print(allele_table)
    print(test_patient_allele_table)
    print(patient_associated_disease_table)

    # export dataframes to csv
    patient_associated_disease_table.to_csv(output_path + sample_name + "_patient_associated_disease.csv")
    test_patient_allele_table.to_csv(output_path + sample_name + "_patient_allele.csv")
    allele_table.to_csv(output_path + sample_name + "_allele_table.csv", index=True)
