# pzw
# 20190411
import sys
import json

############################
# QC_json = open(sys.argv[1], "r")
# outputFile = open(sys.argv[2], "w")

############################
QC_json = open(, "r")
outputFile = open(, "w")
###########################


jsonFile = json.load(QC_json)

beforeFilter = jsonFile["summary"]["before_filtering"]
afterFilter = jsonFile["summary"]["after_filtering"]

before_reads = beforeFilter["total_reads"]
before_bases = beforeFilter["total_bases"]
before_GC = "%.2f" % (beforeFilter["gc_content"] * 100) + "%"
before_Q20 = "%.2f" % (beforeFilter["q20_rate"] * 100) + "%"
before_Q30 = "%.2f" % (beforeFilter["q30_rate"] * 100) + "%"
duplicateRate = "%.2f" % (jsonFile["duplication"]["rate"] * 100) + "%"

after_reads = afterFilter["total_reads"]
after_bases = afterFilter["total_bases"]
after_GC = "%.2f" % (afterFilter["gc_content"] * 100) + "%"
after_Q20 = "%.2f" % (afterFilter["q20_rate"] * 100) + "%"
after_Q30 = "%.2f" % (afterFilter["q30_rate"] * 100) + "%"

print "raw reads: ", before_reads
print "raw bases: ", before_bases
print "raw GC content: ", before_GC
print "raw Q20: ", before_Q20
print "raw Q30: ", before_Q30
print "duplication: ", duplicateRate
print "------------------------------------------"
print "clean reads: ", after_reads
print "clean bases: ", after_bases
print "clean GC content: ", after_GC
print "clean Q20: ", after_Q20
print "clean Q30: ", after_Q30

outputFile.write("raw reads: " + str(before_reads) + "\n")
outputFile.write("raw bases: " + str(before_bases) + "\n")
outputFile.write("raw GC content: " + str(before_GC) + "\n")
outputFile.write("raw Q20: " + before_Q20 + "\n")
outputFile.write("raw Q30: " + before_Q30 + "\n")
outputFile.write("duplication: " + duplicateRate + "\n")
outputFile.write("------------------------------------------\n")
outputFile.write("clean reads: " + str(after_reads) + "\n")
outputFile.write("clean bases: " + str(after_bases) + "\n")
outputFile.write("clean GC content: " + after_GC + "\n")
outputFile.write("clean Q20: " + after_Q20 + "\n")
outputFile.write("clean Q30: " + after_Q30 + "\n")

outputFile.close()
QC_json.close()

print "task done"