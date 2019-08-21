# pzw
# 20190819

def complementReverse(seq):
	seq_upper = seq.upper()
	seq_comp = seq_upper.replace("A", "t").replace("G", "c").replace("T", "a").replace("C", "g")
	seq_comp_degen = seq_comp.replace("R", "y").replace("Y", "r").replace("M", "k").replace("K", "m")
	seq_comp_degen2 = seq_comp_degen.replace("H", "d").replace("D", "h").replace("B", "v").replace("V", "b")

	seq_reverse = seq_comp_degen2[::-1].upper()

	return seq_reverse

