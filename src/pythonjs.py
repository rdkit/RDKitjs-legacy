import json
with open('fpscores.pkl.json') as f:
    lst = json.load(f)

dico = {}
for l in lst:
	if len(l)>2:
		i=0
		for item in l:
			i=i+1
			if i>1:
				dico[item] = l[0]
	else:
		dico[l[1]] = l[0]

print len(lst)
print len(dico)

with open('fpscores.pkldico.jso', 'w') as fp:
    json.dump(dico, fp)

print "done"