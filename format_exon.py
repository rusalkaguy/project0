# format exon 
exon = { 'name': "adsf", "start":50, "stop":60}
exon2 = { 'name': "gdfd", "start":150, "stop":160}
print exon

ex_list1 = [exon]
ex_list2 = [exon, exon2]
print "name", exon["name"]

print "ex_list1 ------------------------"
for x in ex_list1:
	print x["name"]
print "====="
print map(lambda x: x["name"], ex_list1)

print "ex_list2 ------------------------"
names=map(lambda x: x["name"], ex_list2)
print names
# test functions
def format_exon(x):
	return x['name']+':'+str(x["start"])+'-'+str(x["stop"])

# "name:start-stop"
print format_exon(exon) 

# print exons lists
for xl in [ex_list1, ex_list2]:
	print xl
	print ",".join(map(format_exon,xl))