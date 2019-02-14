# coding=utf-8

def isHomOrHet(gentype):
	nu = gentype.split(';')
		if nu[0] == nu[1]:
			print '纯合子'
		else:
			print '杂合子'
