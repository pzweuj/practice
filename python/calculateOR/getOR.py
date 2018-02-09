#!/usr/bin/python
# pzw
# 20180209

name = open('desease.name', 'r')

def func(somestr,fenxingmap,ormap,freqmap):
    import gc
    gc.collect()
    newormap = {}
    newfreqmap = {}
    freqtoosmall = False
    if freqmap and freqmap.get(0) < 0.00000001:
        freqtoosmall = True
    for some in fenxingmap[somestr]:
        some_key = some
        some_or = fenxingmap[somestr][some][0]
        some_freq = fenxingmap[somestr][some][1]
        if ormap:
            for oldkey in ormap:
                oldor = ormap[oldkey]
                oldfreq = freqmap[oldkey]

                newkey = oldkey + '|' + somestr + some_key
                newor = oldor * some_or
                newfreq = oldfreq * some_freq

                if freqtoosmall:
                    #newfreq = newfreq*10
                    pass

                newormap[newkey]=newor
                newfreqmap[newkey]=newfreq
                #print newfreq
        else:
            newkey = somestr + some_key
            newor = some_or
            newfreq = some_freq
            newormap[newkey] = newor
            newfreqmap[newkey] = newfreq

    return newormap,newfreqmap

for names in name:
	if '\n' in names:
		names = names.split('\n')[0]
	desease = open('desease/'+ names + '.txt', 'r')
	results = open('results/' + names + '-results.txt', 'w')
	fenxingmap = {}
	for line in desease:
		a = line.split('\t')
		rsid = a[1]
		alt = a[2]
		altOR = float(a[3])
		altFre = float(a[4])
		ref = a[5]
		refOR = float(a[6])
		refFre = float(a[7])

		fenxingmap[rsid] = {}
		fenxingmap[rsid][alt] = [altOR, altFre]
		fenxingmap[rsid][ref] = [refOR, refFre]

	# print fenxingmap

	#
	finalormap = {}
	finalfreqmap = {}
	#

	for key in fenxingmap.keys():
		finalormap,finalfreqmap = func(key, fenxingmap, finalormap, finalfreqmap)

	# print finalfreqmap
	sortedormaplist = sorted(finalormap.items(), key = lambda item:item[1])

	totalfreq = 0.0
	newfreqlist = []
	minor = 0.0
	or1 = 0.0
	or2 = 0.0
	or3 = 0.0
	maxor = 0.0
	for item in sortedormaplist:
	    key = item[0]
	    orvalue = item[1]
	    value = finalfreqmap[key]
	    totalfreq += value

	    if minor == 0.0:
	        minor = orvalue

	    if 0.0 < totalfreq < 0.5:
	        or1 = orvalue
	    if 0.5 < totalfreq < 0.8:
	        or2 = orvalue
	    elif 0.8 <= totalfreq < 0.9:
	        or3 = orvalue
	    else:
	        maxor = orvalue

	results.write('min totalor:' + str(minor) + '\n')
	results.write('50% split totalor:' + str(or1) + '\n')
	results.write('80% split totalor:' + str(or2) + '\n')
	results.write('90% split totalor:' + str(or3) + '\n')
	results.write('max totalor:' + str(maxor) + '\n')
	results.write('totalfreq:' + str(totalfreq) + '\n')

	desease.close()
	results.close()
	print(names + ' done')

name.close()