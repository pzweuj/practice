# pzw
# 20190515

import sys
import argparse
import pandas as pd

def main(annovarFile, panelSize):

# annovarFile = "test.txt"
# panelSize = 1800000

	df = pd.read_csv(annovarFile, sep="\t", header=0)
	df.rename(columns={"Unnamed: 22":"Quality", "Unnamed: 23":"DP"}, inplace=True)

	# print df["DP"]

	df2 = df[(df["Func.refGene"] == "exonic") | (df["Func.refGene"] == "splicing")]
	del df
	df3 = df2[df2["DP"] >= 50]
	del df2

	for i in df3.index.tolist():
		if df3.loc[i, "ExAC_EAS"] == ".":
			df3.at[i, "ExAC_EAS"] = "0"

		if df3.loc[i, "1000g2015aug_eas"] == ".":
			df3.at[i, "1000g2015aug_eas"] = "0"

	df3["ExAC_EAS"] = df3["ExAC_EAS"].astype("float64", copy=False)
	df3["1000g2015aug_eas"] = df3["1000g2015aug_eas"].astype("float64", copy=False)

	df4 = df3[(df3["ExAC_EAS"] <= 0.01) & (df3["1000g2015aug_eas"] <= 0.01)]
	del df3

	for m in df4.index.tolist():
		if df4.loc[m, "cosmic70"] != ".":
			cancer_list = df4.loc[m, "cosmic70"].split("OCCURENCE=")[1].split(",")

			occur = 0
			for n in cancer_list:
				occur_tmp = int(n.split("(")[0])
				occur += occur_tmp

			# print occur
			if occur >= 4:
				df4.drop(m, inplace=True)

	tmb = float(len(df4)) / (float(panelSize) / 1000000)
	return tmb

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="calculate TMB",
		prog="calculateTMB.py",
		usage="python3 calculateTMB.py -i <annovar output> -s <panel size(bp)>")
	parser.add_argument("-v", "--version", action="version", version="Version 0.1 20190517")
	parser.add_argument("-i", "--input", type=str, help="Input the annovar file")
	parser.add_argument("-s", "--size", type=str, help="Input the panel size")
	if len(sys.argv[1:]) == 0:
		parser.print_help()
		parser.exit()
	args = parser.parse_args()
	main(annovarFile=args.input, panelSize=args.size)
